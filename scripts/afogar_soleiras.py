# -*- coding: utf-8 -*-
"""Rebaixa soleiras ADVERSAS: leito que SOBE rio abaixo nao fica.

    python scripts/afogar_soleiras.py taha_ai_novo/taha_ai.g01 --saida g22

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  Ha secoes cujo leito e MAIS ALTO que o da vizinha de montante -- um
  contradeclive. Na foz do Taio: o leito desce a 334,43 (RS 567) e SOBE
  4,31 m ate a soleira da ultima secao. Rio nao faz isso; e artefato da
  escavacao calculada por secao. Hidraulicamente a soleira apos um poco e
  uma barragem que nao existe: na vazao baixa vira vertedouro com jusante
  seco, e o solver unsteady explode ali (medido: o estouro da rodada 4
  comecou no Taio, DZ de 94,7 m na iteracao 1).

O QUE SE FAZ

  Por reach, de montante para jusante, guarda-se o MINIMO ACUMULADO do
  talvegue. Secao com leito acima do minimo acumulado + `--tol` tem o canal
  (lb..rb) rebaixado ate o minimo acumulado. So se REBAIXA; a planicie nao
  muda; o `XS HTab Starting El` acompanha (leito novo + 0,15 m). O resultado
  e leito monotono nao-crescente rio abaixo -- pocos continuam existindo
  (viram remanso), soleiras nao.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g22")
    tol = _arg(argv, "--tol", 0.05, float)
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
    print(f"saida  : {novo}\n")

    novos = {}
    for k, secs in por.items():
        minimo = np.inf
        for i, d in secs:
            zt = float(np.asarray(d["z"], float).min())
            if zt > minimo + tol:
                delta = zt - minimo
                st = np.asarray(d["sta"], float)
                z = np.asarray(d["z"], float).copy()
                m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
                z[m] = z[m] - delta
                novos[(d["rio"], d["reach"], round(d["rs"], 2))] = {"sta": st, "z": z,
                            "htab": float(z.min()) + 0.15}
                print(f"   {k[0]:13s} {k[1]:3s} RS {d['rs']:9.1f}  soleira "
                      f"{delta:+5.2f} m rebaixada (leito {zt:.2f} -> "
                      f"{minimo:.2f})")
            else:
                minimo = min(minimo, zt)

    if not novos:
        print("nenhuma soleira adversa")
        return

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    # chave por (rio, reach, RS): com ESTRUTURA (Type 5) no arquivo o
    # indice de ler_secoes dessincroniza e o perfil ia para a secao errada
    saida, j = [], 0
    rio_c = reach_c = None
    chave = None
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio_c = p[0].strip()
            reach_c = p[1].strip() if len(p) > 1 else ""
        if l.startswith("Type RM Length L Ch R"):
            p = l.split("=", 1)[1].split(",")
            try:
                chave = (rio_c, reach_c, round(float(p[1]), 2))
            except (ValueError, IndexError):
                chave = None
        nv = novos.get(chave)
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
        mn = np.inf
        for v in zz:
            if v > mn + tol:
                sobe += 1
            mn = min(mn, v)
    print(f"   soleiras adversas restantes: {sobe}   (tem de ser zero)")


if __name__ == "__main__":
    main(sys.argv[1:])
