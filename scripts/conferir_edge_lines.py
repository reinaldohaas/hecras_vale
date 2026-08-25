# -*- coding: utf-8 -*-
"""Conta cruzamentos das edge lines DIRETO NO .g01 escrito.

    python scripts/conferir_edge_lines.py modelo/itajai_acu/itajai_acu.g01

O gerador ja impoe que nao haja nenhum, mas quem confere o gerador nao pode
ser o gerador: isto le o arquivo que o HEC-RAS vai ler e refaz a conta que o
RAS Mapper faz ao montar a superficie de interpolacao -- ligar as pontas
esquerdas de todas as cutlines numa polilinha, idem as direitas, e ver se
alguma se cruza. Saida 0 = limpo.
"""
import os
import sys

from shapely.geometry import LineString
from shapely.strtree import STRtree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes   # noqa: E402


def cruzamentos(pts):
    seg = [LineString([pts[i], pts[i + 1]]) for i in range(len(pts) - 1)]
    arv, par = STRtree(seg), set()
    for i, g in enumerate(seg):
        for j in arv.query(g):
            j = int(j)
            if abs(i - j) > 1 and g.intersects(seg[j]):
                par.add((min(i, j), max(i, j)))
    return par


def main(g):
    S = ler_secoes(g)
    S.sort(key=lambda d: -d["rs"])
    print(f"geometria: {g}   {len(S)} secoes")
    total = 0
    for lado, k in (("esquerda", 0), ("direita", 1)):
        pts = [tuple(d["cut"][k]) for d in S]
        c = cruzamentos(pts)
        total += len(c)
        print(f"   edge line {lado:<9}: {len(c)} cruzamento(s)")
        for i, j in sorted(c)[:8]:
            print(f"      entre RS {S[i]['rs']:.2f} e RS {S[j]['rs']:.2f}")
    print(f"\nTOTAL: {total}")
    return total


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1]) else 0)
