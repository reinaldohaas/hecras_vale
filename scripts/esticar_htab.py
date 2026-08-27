# -*- coding: utf-8 -*-
"""Estica a tabela de propriedades (HTab) das secoes que a tem curta.

    python scripts/esticar_htab.py taha_ai.h12 --saida h15 --alcance 60

A ENTRADA NAO E TOCADA. Sai um .gXX novo.

O DEFEITO

  Secao com `XS HTab Starting El and Incr=el0,0.3, 60` so tem tabela ate
  el0+18 m. Acima disso o RAS EXTRAPOLA a conveyance -- fisica falsa, e
  instavel. As interpoladas do canion explodiam na SEGUNDA cheia (19/07)
  exatamente porque o vale cheio passava do teto da tabela (medido:
  222/222 interpoladas com alcance de 18 m).

O QUE SE FAZ

  Toda linha HTab com alcance < `--alcance` m vira
  `el0, alcance/100, 100` -- mesmo el0, tabela cobrindo o vale inteiro.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corrigir_cutlines import _arg                     # noqa: E402
from ras_io import escrever                            # noqa: E402


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = _arg(argv, "--saida", "h15")
    alcance = _arg(argv, "--alcance", 60.0, float)
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).split(".")[0]
    novo = os.path.join(raiz, f"{base}.{ext}")
    if os.path.abspath(novo) == os.path.abspath(entrada):
        raise SystemExit("saida igual a entrada -- recusado")

    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n").split("\n")
    saida, n = [], 0
    for l in linhas:
        if l.startswith("XS HTab Starting El and Incr="):
            p = l.split("=", 1)[1].split(",")
            try:
                el0 = float(p[0])
                inc = float(p[1])
                npt = int(p[2])
            except (ValueError, IndexError):
                saida.append(l)
                continue
            if inc * npt < alcance:
                # NAO engrossar o incremento: tabela grossa embaixo
                # desestabiliza a 1a cheia (testado: incr 0,6 matou em
                # 172 h). Mantem o incr e sobe o numero de pontos.
                npt2 = min(500, int(round(alcance / inc)))
                saida.append("XS HTab Starting El and Incr=%.2f,%.2f, %d "
                             % (el0, inc, npt2))
                n += 1
                continue
        saida.append(l)
    escrever(novo, "\n".join(saida))
    print(f"esticadas: {n} tabelas para {alcance:.0f} m de alcance")
    conf = open(novo, encoding="latin-1", errors="replace").read()
    print(f"CONFERENCIA: secoes {conf.count('#Sta/Elev')} (nao muda); "
          f"HTab de 100 pontos: "
          f"{conf.count(', 100 ')}")


if __name__ == "__main__":
    main(sys.argv[1:])
