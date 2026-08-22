# -*- coding: utf-8 -*-
"""Declara no .rasmap as camadas de exibicao que faltam na geometria.

    python scripts/camadas_rasmap.py modelo/mirim_t30/mirim_t30.rasmap

POR QUE OS FLOW PATHS E OS MARCADORES DE RS NAO APARECIAM

  Nao e falta de dado nem defeito de geometria: sao CAMADAS DE EXIBICAO, e o
  .rasmap deste projeto simplesmente nao as declara. O modelo de referencia
  (`itajaim_hecras`) declara as duas dentro de cada bloco `RASGeometry`:

      <Layer Type="FlowPaths" Checked="True" />
      <Layer Type="RiverStations" Checked="True" />

  Sem essas linhas o RAS Mapper nao tem o que ligar no painel, e nada se ve.

O QUE NAO E A CAUSA, medido

  `UnitsRiverStation="Feet"` num projeto `SI Units` parece errado e chama a
  atencao, mas o modelo de referencia -- que roda e desenha os marcadores --
  usa `Feet` nas suas 22 ocorrencias. E o padrao do RAS Mapper e nao explica
  os marcadores ausentes. Fica como esta.

NADA DE GEOMETRIA E TOCADO. So o .rasmap, e so acrescentando linhas.
"""
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever    # noqa: E402

FALTANTES = ("FlowPaths", "RiverStations")


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    p = argv[0]
    if not os.path.exists(p):
        raise SystemExit(f"nao achei {p}")
    t = open(p, encoding="latin-1", errors="replace").read()
    ET.parse(p)                       # se ja estiver quebrado, para aqui

    n = 0
    for tipo in FALTANTES:
        if f'Type="{tipo}"' in t:
            print(f"   '{tipo}' ja estava declarado")
            continue
        # entra logo apos a abertura de cada bloco RASGeometry, que e onde o
        # modelo de referencia as poe
        def por(m):
            nonlocal n
            n += 1
            ind = re.match(r"[ \t]*", m.group(0)).group(0)
            return (m.group(0) + "\n" + ind + "  "
                    + f'<Layer Type="{tipo}" Checked="True" />')
        t = re.sub(r'(?m)^[ \t]*<Layer[^>]*Type="RASGeometry"[^>]*>$', por, t)
        print(f"   '{tipo}' acrescentado em {n} bloco(s) RASGeometry")
        n = 0

    escrever(p, t)
    ET.parse(p)
    print(f"   XML valido: {p}")
    return p


if __name__ == "__main__":
    main(sys.argv[1:])
