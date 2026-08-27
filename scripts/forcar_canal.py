# -*- coding: utf-8 -*-
"""FORCA um canal viavel num trecho doente: leito monotonico + prisma minimo.

    python scripts/forcar_canal.py taha_ai.g01 --saida g98 \
        --trecho Itajai_Norte,R2,60000,95000 [--trecho ...] \
        --profundidade 2.0

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

ORDEM DO REINALDO (noite de 26/08): "se algum rio nao der, force um canal
nessa parte" -- preferido a amputar.

O QUE SE FAZ, so DENTRO do trecho pedido:

  1. LEITO MONOTONICO: o talvegue nao pode subir rio abaixo. Regressao
     isotonica (PAV) NAO-CRESCENTE de montante para jusante sobre os
     talvegues do trecho; cada secao e rebaixada (nunca levantada) ate o
     seu valor isotonico. Fora do trecho nada muda -- a licao do global
     que piorou o modelo fica respeitada.
  2. PRISMA MINIMO: dentro do canal (lb..rb) o fundo e escavado ate ficar
     `--profundidade` m abaixo da margem MAIS BAIXA da secao, se ainda nao
     estiver -- lamina fina em rampa e onde o solver oscila.
  3. `XS HTab Starting El` acompanha o novo fundo.

  CONFERENCIA relendo o gravado: monotonicidade exata no trecho, nada
  levantado, profundidade minima atingida.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col, _arg               # noqa: E402
from ras_io import escrever                            # noqa: E402


def pav_nao_crescente(v):
    """Isotonica nao-crescente de montante (inicio) para jusante (fim)."""
    # v esta em ordem de RS DECRESCENTE (montante -> jusante):
    # exigimos v[i+1] <= v[i]; PAV sobre -v crescente
    blocos = [[-x, 1.0] for x in v]
    out = []
    for b in blocos:
        out.append(b)
        while len(out) > 1 and out[-2][0] > out[-1][0]:
            s2, n2 = out.pop()
            s1, n1 = out.pop()
            out.append([(s1 * n1 + s2 * n2) / (n1 + n2), n1 + n2])
    res = []
    for s, n in out:
        res += [-s] * int(n)
    return np.asarray(res)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "g98")
    prof = _arg(argv, "--profundidade", 2.0, float)
    trechos = []
    for k, a in enumerate(argv):
        if a == "--trecho":
            rio, reach, r0, r1 = argv[k + 1].split(",")
            trechos.append((rio.strip(), reach.strip(),
                            float(r0), float(r1)))
    if not trechos:
        raise SystemExit("faltou --trecho rio,reach,rs_min,rs_max")

    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    novos = {}
    for rio, reach, r0, r1 in trechos:
        idx = [i for i, d in enumerate(S)
               if d["rio"] == rio and d["reach"] == reach
               and r0 <= d["rs"] <= r1]
        idx.sort(key=lambda i: -S[i]["rs"])     # montante -> jusante
        if len(idx) < 3:
            print(f"   {rio} {reach} {r0:.0f}-{r1:.0f}: "
                  f"so {len(idx)} secoes -- pulado")
            continue
        tal = np.array([float(np.asarray(S[i]["z"], float).min())
                        for i in idx])
        iso = pav_nao_crescente(tal)
        alvo_fundo = np.minimum(tal, iso)       # so rebaixa
        n_mono, n_prisma = 0, 0
        for j, i in enumerate(idx):
            d = S[i]
            st = np.asarray(d["sta"], float)
            z = np.asarray(d["z"], float).copy()
            m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
            delta_iso = float(tal[j] - alvo_fundo[j])
            if delta_iso > 1e-3:
                z[m] = z[m] - delta_iso
                n_mono += 1
            # prisma minimo: fundo ate `prof` abaixo da margem mais baixa
            borda = min(float(np.interp(d["lb"], st, z)),
                        float(np.interp(d["rb"], st, z)))
            fundo = float(z[m].min())
            falta = (borda - prof) - fundo
            if falta < -1e-3:
                pass                             # ja fundo o bastante
            elif falta > 1e-3:
                z[m] = np.minimum(z[m], borda - prof)
                n_prisma += 1
            novos[i] = {"sta": st, "z": z,
                        "htab": float(z.min()) + 0.15}
        print(f"   {rio} {reach} RS {r0:.0f}-{r1:.0f}: {len(idx)} secoes, "
              f"{n_mono} rebaixadas (isotonica), {n_prisma} escavadas "
              f"(prisma {prof:.1f} m)")

    if not novos:
        print("nada a forcar")
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
    for rio, reach, r0, r1 in trechos:
        idx = [i for i, d in enumerate(B)
               if d["rio"] == rio and d["reach"] == reach
               and r0 <= d["rs"] <= r1]
        idx.sort(key=lambda i: -B[i]["rs"])
        tal = [float(np.asarray(B[i]["z"], float).min()) for i in idx]
        sobe = sum(1 for a, b in zip(tal, tal[1:]) if b > a + 0.01)
        lev = sum(1 for i in idx
                  if np.asarray(B[i]["z"], float).min()
                  > np.asarray(S[i]["z"], float).min() + 1e-3)
        print(f"   {rio} {reach}: subidas rio abaixo restantes {sobe} "
              f"(tem de ser 0); secoes levantadas {lev} (tem de ser 0)")


if __name__ == "__main__":
    main(sys.argv[1:])
