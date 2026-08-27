# -*- coding: utf-8 -*-
"""Prolonga as pontas do eixo de um reach, para a secao extrema CRUZAR.

    python scripts/prolongar_eixo.py modelo/mirim_t30/mirim_t30.g09 \
        --reach Canal_Retif,R1 --metros 30 --saida g10

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO QUE ISTO CORRIGE

  A primeira secao de um reach construido por script costuma cair EXATAMENTE
  sobre o primeiro vertice do eixo. Ali a cutline TOCA o eixo em vez de
  cruza-lo, e o HEC-RAS acusa "XS doesn't intersect the associated Reach" --
  medido no canal: a cutline de RS 7500 fica a 0,006 m do eixo, com o meio da
  secao na estacao 0,0 de 7.549,8 m.

  Tocar nao e cruzar: a interseccao de duas linhas que se encontram na ponta
  pode sair vazia por arredondamento, e a associacao secao-reach se perde.

  A correcao nao move secao nenhuma. Prolonga o EIXO alguns metros para fora,
  na direcao do primeiro e do ultimo segmento, e com isso a secao extrema
  passa a ter eixo dos dois lados.

  Nao altera comprimento de trecho nem estacionamento: o `Reach XY` e a
  polilinha de DESENHO e de associacao, enquanto as distancias que o solver
  usa vem de `Type RM Length L Ch R` e das juncoes.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever    # noqa: E402

L16 = 16


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    alvo = _arg(argv, "--reach")
    m = _arg(argv, "--metros", 30.0, float)
    ext = _arg(argv, "--saida", "g10")
    if not alvo or "," not in alvo:
        raise SystemExit("use --reach Rio,Reach")
    rio_a, rch_a = [x.strip() for x in alvo.split(",", 1)]
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    t = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    saida, i, feito = [], 0, False
    while i < len(t):
        l = t[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            aqui = (p[0].strip(), p[1].strip())
            if aqui == (rio_a, rch_a):
                n = int(t[i + 1].split("=")[1])
                v, j = [], i + 2
                while len(v) < 2 * n:
                    v += [float(t[j][c:c + L16])
                          for c in range(0, len(t[j]), L16)
                          if t[j][c:c + L16].strip()]
                    j += 1
                c = np.array(v).reshape(-1, 2)
                d0 = c[0] - c[1]
                d0 = d0 / max(float(np.hypot(*d0)), 1e-9)
                d1 = c[-1] - c[-2]
                d1 = d1 / max(float(np.hypot(*d1)), 1e-9)
                c2 = np.vstack([c[0] + m * d0, c, c[-1] + m * d1])
                saida.append(l)
                saida.append(f"Reach XY= {len(c2)} ")
                lin = ""
                for k, (x, y) in enumerate(c2):
                    lin += "%16.4f%16.4f" % (x, y)
                    if (k + 1) % 2 == 0:
                        saida.append(lin)
                        lin = ""
                if lin:
                    saida.append(lin)
                print(f"eixo de {rio_a},{rch_a}: {n} -> {len(c2)} vertices   "
                      f"prolongado {m:g} m em cada ponta")
                i = j
                feito = True
                continue
        saida.append(l)
        i += 1
    if not feito:
        raise SystemExit(f"nao achei o reach {rio_a},{rch_a}")
    escrever(novo, "\n".join(saida))

    # ---------------------------------------------------------- conferencia
    from qc_secoes import ler_secoes
    from qc_geometria import ler_eixos
    from corrigir_cutlines import mapa_reaches
    from shapely.geometry import LineString
    print("\nCONFERENCIA (relendo o arquivo gravado)")
    for rot, g in (("antes", entrada), ("depois", novo)):
        S = ler_secoes(g)
        eixos = ler_eixos(g)
        mapa = mapa_reaches(g)
        ruins = 0
        for k, d in enumerate(S):
            ln = LineString(np.asarray(d["cut"], float))
            x = ln.intersection(eixos[mapa[k]])
            nn = 0 if x.is_empty else (len(x.geoms)
                                       if hasattr(x, "geoms") else 1)
            if nn == 0:
                ruins += 1
        print(f"   {rot:<6}: secoes que NAO cruzam o proprio eixo: {ruins}")
    a = open(entrada, encoding="latin-1", errors="replace").read()
    b = open(novo, encoding="latin-1", errors="replace").read()
    for chave in ("#Sta/Elev=", "Bank Sta=", "XS GIS Cut Line=",
                  "Type RM Length L Ch R"):
        print(f"   {chave:<24} {a.count(chave)} -> {b.count(chave)}  "
              f"{'ok' if a.count(chave) == b.count(chave) else 'DIVERGIU'}")
    return novo


if __name__ == "__main__":
    main(sys.argv[1:])
