# -*- coding: utf-8 -*-
"""Converte `Initial RS=` (chave FANTASMA) em `Initial Flow Loc=` no u01.

    python scripts/converter_inicial.py taha_ai.u01 [--escala 0.5]

EDITA O u01 (backup .antes_do_inicial). E o reparo que destravou o
modelo duas vezes: o gerador escreve `Initial RS=rio,reach,rs,Q`, que o
HEC-RAS IGNORA em silencio -- a chave real e `Initial Flow Loc=`. O RS
gravado e conferido contra o g01 do projeto e ajustado para a secao
EXISTENTE mais proxima (o RAS exige casamento exato).

`--escala` multiplica as vazoes (p.ex. 0.5 para iniciar no seco de
junho, quando a serie de contorno comeca antes do evento).
"""
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from ras_io import escrever                            # noqa: E402


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    u01 = argv[0]
    escala = float(argv[argv.index("--escala") + 1]) \
        if "--escala" in argv else 1.0
    g01 = u01.rsplit(".", 1)[0] + ".g01"
    shutil.copy2(u01, u01 + ".antes_do_inicial")

    S = ler_secoes(g01)
    por_reach = {}
    for d in S:
        por_reach.setdefault((d["rio"], d["reach"]), []).append(d["rs"])

    t = open(u01, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n")
    saida, n = [], 0
    for l in t.split("\n"):
        if l.startswith("Initial RS="):
            rio, reach, rs, q = [x.strip() for x in
                                 l.split("=", 1)[1].split(",")]
            rss = por_reach.get((rio, reach))
            if not rss:
                print(f"   {rio} {reach}: sem secoes no g01 -- mantido")
                saida.append(l)
                continue
            rs_ok = min(rss, key=lambda r: abs(r - float(rs)))
            q_ok = float(q) * escala
            saida.append("Initial Flow Loc=%s,%s,%s,%g"
                         % (rio, reach, ("%.2f" % rs_ok).rstrip("0")
                            .rstrip("."), q_ok))
            n += 1
        else:
            saida.append(l)
    escrever(u01, "\n".join(saida))

    print(f"convertidas: {n}")
    conf = open(u01, encoding="latin-1", errors="replace").read()
    print(f"CONFERENCIA: Initial Flow Loc={conf.count('Initial Flow Loc=')}"
          f"  Initial RS restantes={conf.count('Initial RS=')}"
          f"  (fantasma tem de ser 0)")


if __name__ == "__main__":
    main(sys.argv[1:])
