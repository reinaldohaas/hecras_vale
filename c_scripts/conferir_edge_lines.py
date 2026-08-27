# -*- coding: utf-8 -*-
"""Conta defeitos NAS LINHAS QUE O HEC-RAS CONSTRUIU: edge E bank lines.

    python scripts/conferir_edge_lines.py modelo/itajai_oeste/itajai_oeste.g01.hdf

POR QUE LER O HDF, E NAO REFAZER A CONTA

  A primeira versao disto ligava as pontas das cutlines numa poligonal e
  testava o cruzamento. Dava ZERO no Oeste enquanto o RAS Mapper acusava dois
  pontos de auto-interseccao -- ou seja, o conferidor concordava com o gerador
  e os dois estavam errados juntos, que e o pior arranjo possivel.

  A edge line do HEC-RAS nao e a poligonal das pontas: no Oeste ela tem 862
  vertices para 380 secoes. Ele densifica o traco, e a dobra so aparece na
  versao densificada. Entao aqui nao se reconstroi nada: le-se do `.gNN.hdf`
  o traco que o RAS usa, e mede-se nele. Saida 0 = limpo.

POR QUE AS BANK LINES TAMBEM

  A licao valia pela metade: o pipeline conferia a edge line e deixava a BANK
  LINE sem porteiro nenhum -- e foi na bank line que o defeito do Mirim
  apareceu primeiro ("bank line atravessando o rio em todo o trecho", 82
  cruzamentos com o eixo). Aqui medem-se, por bank line:

    auto-interseccao   a linha dobra sobre si mesma (mesma conta da edge);
    cruzamento do eixo a margem passa para o outro lado do rio -- a estacao
                       de margem daquela secao esta do lado errado do
                       talvegue, e a vazao do canal fica repartida errada.

  O `TOTAL:` da saida soma tudo, e e o numero que o `construir_rio.py` le na
  etapa 9 -- bank line suja agora reprova o rio, como a edge ja reprovava.
"""
import os
import sys

import h5py
from shapely.geometry import LineString, MultiLineString
from shapely.strtree import STRtree

EDGE = "/Geometry/River Edge Lines"
BANK = "/Geometry/River Bank Lines"
EIXO = "/Geometry/River Centerlines"


def cruzamentos(P):
    """Pares de segmentos nao vizinhos que se tocam (auto-interseccao)."""
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


def linhas(f, cam):
    """[(k, Nx2)] de cada polyline do grupo `cam` do HDF (ou lista vazia)."""
    if cam not in f:
        return []
    info = f[cam + "/Polyline Info"][:]
    pts = f[cam + "/Polyline Points"][:]
    return [(k, pts[int(l[0]):int(l[0]) + int(l[1])])
            for k, l in enumerate(info)]


def n_toques(a, b):
    """Quantas vezes a linha `a` toca a linha `b` (pontos ou trechos)."""
    g = LineString(a).intersection(LineString(b))
    if g.is_empty:
        return 0
    if g.geom_type in ("Point", "LineString"):
        return 1
    return len(g.geoms)


def main(caminho):
    hdf = caminho if caminho.lower().endswith(".hdf") else caminho + ".hdf"
    if not os.path.exists(hdf):
        raise SystemExit(f"nao achei {hdf} -- rode a validacao antes, que e "
                         "quem monta o HDF pelo preprocessador geometrico")
    total = 0
    with h5py.File(hdf, "r") as f:
        print(f"geometria: {hdf}")
        eixos = linhas(f, EIXO)
        for k, P in linhas(f, EDGE):
            c = cruzamentos(P)
            total += len(c)
            print(f"   edge line {k}: {len(P):5d} vertices   "
                  f"{len(c)} cruzamento(s)")
            for i, j in sorted(c)[:6]:
                print(f"      segmentos {i} e {j}, em "
                      f"{P[i][0]:.1f} {P[i][1]:.1f}")
        for k, P in linhas(f, BANK):
            c = cruzamentos(P)
            atravessa = sum(n_toques(P, E) for _, E in eixos)
            total += len(c) + atravessa
            print(f"   bank line {k}: {len(P):5d} vertices   "
                  f"{len(c)} auto-interseccao(oes)   "
                  f"{atravessa} cruzamento(s) do eixo")
            for i, j in sorted(c)[:6]:
                print(f"      segmentos {i} e {j}, em "
                      f"{P[i][0]:.1f} {P[i][1]:.1f}")
    print(f"\nTOTAL: {total}")
    return total


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1]) else 0)
