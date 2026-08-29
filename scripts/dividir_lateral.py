# -*- coding: utf-8 -*-
"""Divide uma lateral uniforme em duas no ponto pedido (ex.: na barragem).

    python scripts/dividir_lateral.py taha_ai.u01 \
        --rio Itajai_Sul --reach R1 --em 33500 --fracao-acima 0.46

EDITA O u01 (backup .antes_da_divisao).

POR QUE

  A lateral uniforme do Sul cobria o reach inteiro (71209 -> 943): quase
  toda a agua do sistema entrava A JUSANTE da Barragem Sul, o pool so
  recebia a cabeceira (22%) e Ituporanga simulava -53%. Na realidade a
  bacia na barragem tem ~1150 km2 de 1990 (58%): o reservatorio encheu e
  verteu em julho/1983.

O QUE SE FAZ

  A lateral uniforme do reach que ATRAVESSA `--em` vira duas:

    de RS_ini ate a 1a secao ACIMA de --em, com serie x fracao-acima
    da 1a secao ABAIXO de --em ate RS_fim, com serie x (1 - fracao-acima)

  Fracao-acima = (area na barragem - area na cabeceira) / area da lateral,
  por exemplo Sul: (1150-434)/(1990-434) = 0,46.

  CONFERENCIA: relendo o gravado, soma dos picos das duas novas series ==
  pico da antiga (a agua se conserva).
"""
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                       # noqa: E402
from ras_io import escrever                            # noqa: E402


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def fmt_serie(vals):
    corpo, lin = [], ""
    for i, x in enumerate(vals):
        lin += "%8.2f" % x
        if (i + 1) % 10 == 0:
            corpo.append(lin)
            lin = ""
    if lin:
        corpo.append(lin)
    return corpo


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    u01 = argv[0]
    rio = _arg(argv, "--rio")
    reach = _arg(argv, "--reach")
    em = _arg(argv, "--em", None, float)
    fa = _arg(argv, "--fracao-acima", None, float)
    if None in (rio, reach, em, fa):
        raise SystemExit("faltou --rio/--reach/--em/--fracao-acima")

    g01 = u01.rsplit(".", 1)[0] + ".g01"
    S = [d for d in ler_secoes(g01)
         if d["rio"] == rio and d["reach"] == reach and d["tipo"] == "1"]
    acima = min((d["rs"] for d in S if d["rs"] > em),
                default=None)
    abaixo = max((d["rs"] for d in S if d["rs"] < em), default=None)
    if acima is None or abaixo is None:
        raise SystemExit(f"sem secoes dos dois lados de {em}")

    shutil.copy2(u01, u01 + ".antes_da_divisao")
    t = open(u01, encoding="latin-1", errors="replace").read() \
        .replace("\r\n", "\n").replace("\r", "\n")
    blocos = re.split(r"(?=^Boundary Location=)", t, flags=re.M)
    achou = False
    for k, b in enumerate(blocos):
        m = re.match(r"Boundary Location=([^,]+),([^,]+),([\d.]+)\s*,"
                     r"([\d.]+)\s*,", b)
        if not m or m.group(1).strip() != rio \
                or m.group(2).strip() != reach:
            continue
        if "Uniform Lateral Inflow Hydrograph" not in b:
            continue
        r_ini, r_fim = float(m.group(3)), float(m.group(4))
        if not (r_fim < em < r_ini):
            continue
        h = re.search(r"Uniform Lateral Inflow Hydrograph=\s*(\d+)", b)
        vals, resto_ini = [], None
        linhas_b = b[h.end():].split("\n")
        for li, l in enumerate(linhas_b[1:], 1):
            if not l.strip() or l[:1].isalpha():
                resto_ini = li
                break
            vals += [float(l[c:c + 8]) for c in range(0, len(l), 8)
                     if l[c:c + 8].strip()]
        resto = "\n".join(linhas_b[resto_ini:]) if resto_ini else ""
        cab = b.split("\n")[0]
        pre = b[:b.index("\n") + 1]
        meio = b[len(pre):h.start() - 0]
        # cabecalhos novos
        def bloco_novo(rs_a, rs_b, frac):
            serie = [v * frac for v in vals]
            cab_n = re.sub(r"^(Boundary Location=[^,]+,[^,]+,)[\d.]+\s*,"
                           r"[\d.]+\s*,",
                           r"\g<1>%-8s,%-8s," % (
                               ("%.2f" % rs_a).rstrip("0").rstrip("."),
                               ("%.2f" % rs_b).rstrip("0").rstrip(".")),
                           cab)
            corpo = [cab_n] + meio.strip("\n").split("\n") + \
                ["Uniform Lateral Inflow Hydrograph= %d " % len(serie)] + \
                fmt_serie(serie) + [resto.rstrip("\n")]
            return "\n".join(x for x in corpo if x != "") + "\n\n"
        blocos[k] = bloco_novo(r_ini, acima, fa) \
            + bloco_novo(abaixo, r_fim, 1.0 - fa)
        achou = True
        pico = max(vals)
        print(f"lateral {rio} {reach} {r_ini:.0f}->{r_fim:.0f} "
              f"(pico {pico:.1f}) dividida em {em:.0f}:")
        print(f"   acima : {r_ini:.0f}->{acima:.0f}  x {fa:.2f} "
              f"(pico {pico*fa:.1f})")
        print(f"   abaixo: {abaixo:.0f}->{r_fim:.0f}  x {1-fa:.2f} "
              f"(pico {pico*(1-fa):.1f})")
        break
    if not achou:
        raise SystemExit("nenhuma lateral uniforme atravessa esse ponto")
    escrever(u01, "".join(blocos))

    conf = open(u01, encoding="latin-1", errors="replace").read()
    n = len(re.findall(r"Boundary Location=%s\s*,%s\s*,[\d.]+\s*,[\d.]+"
                       % (re.escape(rio), re.escape(reach)), conf))
    print(f"CONFERENCIA: laterais de {rio} {reach} no arquivo: {n}")


if __name__ == "__main__":
    main(sys.argv[1:])
