# -*- coding: utf-8 -*-
"""Conta cruzamentos NA EDGE LINE QUE O HEC-RAS CONSTRUIU.

    python scripts/conferir_edge_lines.py modelo/itajai_oeste/itajai_oeste.g01.hdf

POR QUE LER O HDF, E NAO REFAZER A CONTA

  A primeira versao disto ligava as pontas das cutlines numa poligonal e
  testava o cruzamento. Dava ZERO no Oeste enquanto o RAS Mapper acusava dois
  pontos de auto-interseccao -- ou seja, o conferidor concordava com o gerador
  e os dois estavam errados juntos, que e o pior arranjo possivel.

  A edge line do HEC-RAS nao e a poligonal das pontas: no Oeste ela tem 862
  vertices para 380 secoes. Ele densifica o traco, e a dobra so aparece na
  versao densificada -- a linha ia ate a ponta de duas secoes que avancavam
  cem metros para fora e voltava por cima do caminho de ida.

  Entao aqui nao se reconstroi nada: le-se `/Geometry/River Edge Lines` do
  `.gNN.hdf`, que e o traco que o RAS usa para montar a superficie de
  interpolacao, e mede-se nele. Saida 0 = limpo.
"""
import os
import sys

import h5py
from shapely.geometry import LineString
from shapely.strtree import STRtree

CAM = "/Geometry/River Edge Lines"


def cruzamentos(P):
    seg = [LineString([P[i], P[i + 1]]) for i in range(len(P) - 1)]
    if len(seg) < 3:
        return set()
    arv, par = STRtree(seg), set()
    for i, g in enumerate(seg):
        for j in arv.query(g):
            j = int(j)
            if abs(i - j) > 1 and g.intersects(seg[j]):
                par.add((min(i, j), max(i, j)))
    return par


def main(caminho):
    hdf = caminho if caminho.lower().endswith(".hdf") else caminho + ".hdf"
    if not os.path.exists(hdf):
        raise SystemExit(f"nao achei {hdf} -- rode a validacao antes, que e "
                         "quem monta o HDF pelo preprocessador geometrico")
    with h5py.File(hdf, "r") as f:
        if CAM not in f:
            raise SystemExit(f"{hdf} nao tem {CAM}")
        info = f[CAM + "/Polyline Info"][:]
        pts = f[CAM + "/Polyline Points"][:]
    print(f"geometria: {hdf}")
    total = 0
    for k, linha in enumerate(info):
        ini, n = int(linha[0]), int(linha[1])
        P = pts[ini:ini + n]
        c = cruzamentos(P)
        total += len(c)
        print(f"   edge line {k}: {n:5d} vertices   {len(c)} cruzamento(s)")
        for i, j in sorted(c)[:6]:
            print(f"      segmentos {i} e {j}, em {P[i][0]:.1f} {P[i][1]:.1f}")
    print(f"\nTOTAL: {total}")
    return total


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1]) else 0)
