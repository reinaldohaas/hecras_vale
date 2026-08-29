# -*- coding: utf-8 -*-
"""Poe um piso nas vazoes de CABECEIRA de um .uNN, para o arranque molhado.

    python scripts/piso_de_vazao.py taha_ai_novo/taha_ai.u01 \
        --fracao 0.10 --minimo 2.0

GRAVA NO PROPRIO ARQUIVO, com backup `<arquivo>.antes_do_piso`.

POR QUE

  A hidrologia da familia parte de 2% do pico como vazao de base. Nas
  cabeceiras de montanha (Benedito com 8% de declividade, Norte, Taio) isso
  da centimetros de lamina no canal piloto -- e o solver unsteady explode na
  INICIALIZACAO, com erros de 9 a 28 m ja no primeiro passo (medido no
  taha_ai_novo, 26/08/2026: "Solution solver went unstable, iteration 9 at
  31JUL 00:00:05"). Canal quase seco e a instabilidade classica do HEC-RAS.

O QUE SE FAZ

  So nos blocos `Flow Hydrograph` (cabeceiras). Cada serie recebe
  piso = max(fracao * pico_da_serie, minimo). Laterais e mare nao mudam.
  O pico NAO muda; so as horas de estiagem sobem ate o piso. O acrescimo de
  volume sai impresso -- e pequeno (o piso vale poucas horas de subida).
"""
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ras_io import escrever            # noqa: E402


def _arg(argv, chave, padrao=None, tipo=float):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    u = argv[0]
    fracao = _arg(argv, "--fracao", 0.10)
    minimo = _arg(argv, "--minimo", 2.0)
    shutil.copy2(u, u + ".antes_do_piso")

    t = open(u, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n")
    blocos = re.split(r"(?=^Boundary Location=)", t, flags=re.M)
    print(f"arquivo: {u}   (backup .antes_do_piso)")
    print(f"piso   : max({fracao:.0%} do pico, {minimo} m3/s)\n")

    extra_total = 0.0
    for k, b in enumerate(blocos):
        m = re.search(r"^Flow Hydrograph=\s*(\d+)\s*$", b, flags=re.M)
        if not m:
            continue
        n = int(m.group(1))
        ini = m.end() + 1
        vals, pos = [], ini
        for l in b[ini:].split("\n"):
            if not l.strip() or l[:1].isalpha():
                break
            vals += [float(l[i:i + 8]) for i in range(0, len(l), 8)
                     if l[i:i + 8].strip()]
            pos += len(l) + 1
            if len(vals) >= n:
                break
        vals = vals[:n]
        pico = max(vals)
        piso = max(fracao * pico, minimo)
        novos = [max(v, piso) for v in vals]
        extra = sum(novos) - sum(vals)
        extra_total += extra
        onde = b.split("\n")[0].split("=")[1]
        print(f"   {onde.split(',')[0].strip():14s} pico {pico:7.1f}  "
              f"piso {piso:6.2f}  volume extra {extra*3600/1e6:6.3f} hm3")
        lin, corpo = "", []
        for i, x in enumerate(novos):
            lin += "%8.2f" % x
            if (i + 1) % 10 == 0:
                corpo.append(lin)
                lin = ""
        if lin:
            corpo.append(lin)
        blocos[k] = b[:ini] + "\n".join(corpo) + "\n" + b[pos:]

    escrever(u, "".join(blocos))
    print(f"\nvolume extra total: {extra_total*3600/1e6:.3f} hm3 "
          "(sobre ~109 hm3 do evento)")


if __name__ == "__main__":
    main(sys.argv[1:])
