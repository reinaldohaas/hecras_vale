# -*- coding: utf-8 -*-
"""Sanea o formato dos perfis: n na primeira estaca e estaca duplicada.

    python scripts/sanear_perfil.py taha_ai_novo/taha_ai.g01 --saida g09

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

DUAS REGRAS DE FORMATO que o preprocessador cobra em data_errors.txt:

  1. "A horizontal Manning's n value needs to be specified on first
     station."  A primeira quebra do #Mann TEM de estar na primeira estaca
     do perfil. Extensao pela esquerda desloca as quebras e deixa o trecho
     novo sem n; aqui a primeira quebra volta para a primeira estaca,
     estendendo o n que ja era o primeiro (o da planicie).

  2. "Station and elevation data contains duplicate points."  Duas estacas
     que so diferem alem da 2a casa viram duplicata quando gravadas em
     %8.2f. A segunda sai; a cota que fica e a da primeira.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from corrigir_cutlines import _col                     # noqa: E402
from ras_io import escrever                            # noqa: E402


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g09"
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")

    n_mann, n_dup = 0, 0
    saida, j = [], 0
    st0 = None          # primeira estaca do #Sta/Elev da secao corrente
    while j < len(linhas):
        l = linhas[j]
        if l.startswith("#Sta/Elev"):
            cnt = int(l.split("=")[1])
            bruto, k2 = [], j + 1
            while k2 < len(linhas) and len(bruto) < 2 * cnt:
                x = linhas[k2]
                if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                    break
                bruto += [x[c:c + 8] for c in range(0, len(x), 8)
                          if x[c:c + 8].strip()]
                k2 += 1
            v = [float(x) for x in bruto[:2 * cnt]]
            st = v[0::2]
            z = v[1::2]
            keep = [0]
            for t in range(1, len(st)):
                if round(st[t], 2) > round(st[keep[-1]], 2):
                    keep.append(t)
            if len(keep) != len(st):
                n_dup += len(st) - len(keep)
                st = [st[t] for t in keep]
                z = [z[t] for t in keep]
            st0 = st[0]
            par = []
            for a, b in zip(st, z):
                par += [a, b]
            saida.append("#Sta/Elev= %d " % len(st))
            saida += _col(par, 8, 2)
            j = k2
            continue
        if l.startswith("#Mann=") and st0 is not None:
            cnt = int(l.split("=")[1].split(",")[0])
            resto = l.split("=", 1)[1]
            bruto, k2 = [], j + 1
            while k2 < len(linhas) and len(bruto) < 3 * cnt:
                x = linhas[k2]
                if not x.strip() or x[:1].isalpha() or x[:1] == "#":
                    break
                bruto += [x[c:c + 8] for c in range(0, len(x), 8)
                          if x[c:c + 8].strip()]
                k2 += 1
            val = [float(x) for x in bruto[:3 * cnt]]
            if val and abs(val[0] - st0) > 0.005:
                val[0] = st0
                n_mann += 1
            saida.append(l)
            lin, corpo = "", []
            for t, x in enumerate(val):
                lin += ("%8.2f" % x if t % 3 == 0 else
                        "%8.3f" % x if t % 3 == 1 else "%8d" % int(x))
                if (t + 1) % 9 == 0:
                    corpo.append(lin)
                    lin = ""
            if lin:
                corpo.append(lin)
            saida += corpo
            j = k2
            continue
        saida.append(l)
        j += 1
    escrever(novo, "\n".join(saida))
    print(f"entrada: {entrada}   (intocada)")
    print(f"saida  : {novo}")
    print(f"   quebras de Manning trazidas a 1a estaca: {n_mann}")
    print(f"   estacas duplicadas removidas           : {n_dup}")

    A2 = ler_secoes(entrada)
    B2 = ler_secoes(novo)
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    print(f"   secoes: {len(A2)} -> {len(B2)}   (nao pode mudar)")
    tal = max(abs(float(a['z'].min()) - float(b['z'].min()))
              for a, b in zip(A2, B2))
    print(f"   talvegue mudou no maximo {tal:.6f}  (tem de ser zero)")


if __name__ == "__main__":
    main(sys.argv[1:])
