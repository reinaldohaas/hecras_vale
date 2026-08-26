# -*- coding: utf-8 -*-
"""Baixa o canal ESCAVADO NA ENCOSTA ate o talvegue da propria secao.

    python scripts/nivelar_bancadas.py taha_ai_novo/taha_ai.g01 --saida g23

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  Em 32 secoes o fundo do canal marcado (lb..rb) esta ate 14 m ACIMA do
  ponto mais baixo da propria secao: a escavacao foi feita na posicao do
  eixo, e onde o eixo passa na encosta o canal virou uma BANCADA no morro,
  enquanto o fundo do vale continua la, no meio da "planicie". O perfil
  longitudinal do rio entao sobe 14 m numa secao e desce na seguinte
  (medido no Pombas RS 17007: canal a 351,86 com talvegue a 337,90) -- e o
  solver explode nessas escadas fantasma.

O QUE SE FAZ

  Secao com fundo de canal mais de `--tol` m acima do proprio talvegue tem
  o canal inteiro (lb..rb) rebaixado pela diferenca: o fundo do canal passa
  a ser o talvegue que a secao JA tem. Nenhuma cota nova e inventada -- e o
  mesmo criterio de escavacao do modelo, aplicado na cota certa. So se
  rebaixa; planicie intacta; `XS HTab Starting El` acompanha.
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
    ext = _arg(argv, "--saida", "g23")
    tol = _arg(argv, "--tol", 1.0, float)
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    S = ler_secoes(entrada)
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}\n")

    novos = {}
    for i, d in enumerate(S):
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
        if m.sum() < 2:
            continue
        delta = float(z[m].min() - z.min())
        if delta <= tol:
            continue
        z2 = z.copy()
        z2[m] = z2[m] - delta
        novos[i] = {"sta": st, "z": z2, "htab": float(z2.min()) + 0.15}
        print(f"   {d['rio']:13s} {d['reach']:3s} RS {d['rs']:9.1f}  "
              f"bancada {delta:5.1f} m nivelada ao talvegue "
              f"{z.min():7.2f}")

    if not novos:
        print("nenhuma bancada acima da tolerancia")
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
    resta = 0
    for d in B:
        st = np.asarray(d["sta"], float)
        z = np.asarray(d["z"], float)
        m = (st >= d["lb"] - 1e-6) & (st <= d["rb"] + 1e-6)
        if m.sum() >= 2 and z[m].min() - z.min() > tol:
            resta += 1
    print(f"   bancadas restantes: {resta}   (tem de ser zero)")


if __name__ == "__main__":
    main(sys.argv[1:])
