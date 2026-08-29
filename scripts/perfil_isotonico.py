# -*- coding: utf-8 -*-
"""Regrada o perfil do leito por regressao ISOTONICA suavizada, por reach.

    python scripts/perfil_isotonico.py taha_ai_novo/taha_ai.g01 --saida g32

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO QUE ISTO ENCERRA

  O leito escavado e uma ESCADARIA: patamares mortos de 0,01-0,02% de
  declividade intercalados com espelhos de 1 a 3 m, por centenas de
  quilometros. Cada espelho e uma barragenzinha; na vazao baixa o 1D
  unsteady balanca agua entre os patamares ate divergir -- 20 rodadas
  medidas, o estouro muda de endereco mas a doenca e essa.

O QUE SE FAZ

  Por reach: regressao isotonica (PAV) do talvegue contra a distancia,
  impondo leito NAO-CRESCENTE rio abaixo; em seguida media movel de
  `--janela` m re-projetada no cone isotonico. O delta (novo - velho) e
  aplicado a TODOS os pontos do canal (lb..rb), deslocando a calha sem
  deformar a forma; a planicie nao muda; o HTab acompanha. O volume do
  perfil se conserva na media (a isotonica e o ajuste de minimos
  quadrados sob a restricao de monotonia).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402


def pav_nao_crescente(y):
    """Pool Adjacent Violators para sequencia nao-crescente."""
    n = len(y)
    vals = list(y)
    pesos = [1.0] * n
    blocos = [[i] for i in range(n)]
    k = 0
    out_v, out_w, out_b = [], [], []
    for i in range(n):
        out_v.append(vals[i])
        out_w.append(pesos[i])
        out_b.append(blocos[i])
        while len(out_v) > 1 and out_v[-2] < out_v[-1]:
            v = (out_v[-2] * out_w[-2] + out_v[-1] * out_w[-1]) \
                / (out_w[-2] + out_w[-1])
            w = out_w[-2] + out_w[-1]
            b = out_b[-2] + out_b[-1]
            out_v = out_v[:-2] + [v]
            out_w = out_w[:-2] + [w]
            out_b = out_b[:-2] + [b]
    z = np.empty(n)
    for v, b in zip(out_v, out_b):
        for i in b:
            z[i] = v
    return z


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g32")
    janela = _arg(argv, "--janela", 1500.0, float)
    rios = _arg(argv, "--rios", None)
    rios = set(x.strip() for x in rios.split(",")) if rios else None
    minimo = _arg(argv, "--minimo", 0.05, float)
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    por = {}
    for i, d in enumerate(S):
        por.setdefault((d["rio"], d["reach"]), []).append((i, d))
    for k in por:
        por[k].sort(key=lambda t: -t[1]["rs"])

    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"janela : {janela:.0f} m\n")

    novos = {}
    resumo = []
    for k, secs in por.items():
        if rios is not None and k[0] not in rios:
            continue
        rs = np.array([d["rs"] for _, d in secs])
        z0 = np.array([float(np.asarray(d["z"], float).min())
                       for _, d in secs])
        # isotonica: nao-crescente rio abaixo (rs decresce no vetor)
        z1 = pav_nao_crescente(z0)
        # suavizacao por media movel em distancia + re-isotonica
        z2 = np.empty_like(z1)
        for a in range(len(rs)):
            m = np.abs(rs - rs[a]) <= janela / 2
            z2[a] = z1[m].mean()
        z2 = pav_nao_crescente(z2)
        # extremos ancorados (juncoes/contornos nao mudam de cota)
        z2[0] = z0[0]
        z2[-1] = min(z0[-1], z2[-1])
        z2 = pav_nao_crescente(z2)
        dz = z2 - z0
        n_mex = int((np.abs(dz) > minimo).sum())
        resumo.append((k, n_mex, float(np.abs(dz).max())))
        for (i, d), delta in zip(secs, dz):
            if abs(delta) <= minimo:
                continue
            st = np.asarray(d["sta"], float)
            z = np.asarray(d["z"], float).copy()
            m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
            z[m] = z[m] + delta
            novos[i] = {"sta": st, "z": z, "htab": float(z.min()) + 0.15}

    for k, n_mex, pior in resumo:
        print(f"   {k[0]:13s} {k[1]:3s}: {n_mex:4d} secoes regradadas   "
              f"|dz| max {pior:5.2f} m")
    print(f"\n   total regradado: {len(novos)} de {len(S)} secoes")

    if not novos:
        print("nada a regradar")
        return

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i_sec, saida, j = -1, [], 0
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("Type RM Length L Ch R"):
            i_sec += 1
        nv = novos.get(i_sec)
        if nv is not None:
            if l.startswith("#Sta/Elev"):
                v = []
                for a, b in zip(nv["sta"], nv["z"]):
                    v += [a, b]
                saida.append("#Sta/Elev= %d " % len(nv["sta"]))
                saida += _col(v, 8, 2)
                cnt = int(l.split("=")[1])
                j += 1
                lidos = 0
                while j < len(linhas) and lidos < 2 * cnt:
                    x = linhas[j]
                    if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                        break
                    lidos += len([1 for c in range(0, len(x), 8)
                                  if x[c:c + 8].strip()])
                    j += 1
                continue
            if l.startswith("XS HTab Starting El and Incr="):
                resto = l.split("=", 1)[1].split(",")
                saida.append("XS HTab Starting El and Incr=%.2f,%s,%s"
                             % (nv["htab"], resto[1], resto[2]))
                j += 1
                continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))

    print("\nCONFERENCIA (relendo o arquivo gravado)")
    B = ler_secoes(novo)
    print(f"   secoes: {len(S)} -> {len(B)}   (nao pode mudar)")
    por2 = {}
    for d in B:
        por2.setdefault((d["rio"], d["reach"]), []).append(d)
    sobe = 0
    for k, secs in por2.items():
        secs.sort(key=lambda d: -d["rs"])
        zz = [float(np.asarray(d["z"], float).min()) for d in secs]
        sobe += sum(1 for a, b in zip(zz, zz[1:]) if b > a + 0.05)
    print(f"   leito subindo rio abaixo: {sobe} pares   (tem de ser zero)")


if __name__ == "__main__":
    main(sys.argv[1:])
