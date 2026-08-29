# -*- coding: utf-8 -*-
"""Insere uma secao interpolada no meio de vaos graves de RS.

Uso:
    python c_scripts/preencher_vaos_graves.py modelo/_codex/itajai_acu_r6/itajai_acu_r6.g01 --saida g02

A entrada nao e tocada. A saida fica no mesmo diretorio, com a extensao dada.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes  # noqa: E402
from ras_io import escrever       # noqa: E402


def _col(v):
    out, linha = [], ""
    for i, x in enumerate(v):
        linha += f"{x:8.2f}"
        if (i + 1) % 10 == 0:
            out.append(linha)
            linha = ""
    if linha:
        out.append(linha)
    return out


def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s else "0"


def _type_line(rs, dist):
    return (f"Type RM Length L Ch R = 1 ,{rs:.2f},"
            f"{dist:8.2f},{dist:8.2f},{dist:8.2f}")


def _block(q, dist):
    A, B = np.asarray(q["cut"][0], float), np.asarray(q["cut"][1], float)
    lines = [_type_line(q["rs"], dist),
             f"Bank Sta={_fmt(q['lb'])},{_fmt(q['rb'])}",
             "XS GIS Cut Line= 2",
             "".join(f"{x:16.4f}" for x in (A[0], A[1], B[0], B[1])),
             f"#Sta/Elev= {len(q['sta'])} "]
    v = []
    for a, z in zip(q["sta"], q["z"]):
        v += [float(a), float(z)]
    lines += _col(v)
    lines += ["#Mann= 3 , 0 , 0 ",
              f"{0.0:8.2f}{0.055:8.3f}{0:8d}"
              f"{q['lb']:8.2f}{0.032:8.3f}{0:8d}"
              f"{q['rb']:8.2f}{0.055:8.3f}{0:8d}",
              f"XS HTab Starting El and Incr={min(q['z']) + 0.02:.2f},0.100, 500 ",
              "XS HTab Horizontal Distribution=-1,-1,-1",
              "XS Rating Curve= 0 ,0",
              "Exp/Cntr=0.3,0.1"]
    return lines


def _interp(a, b, rs):
    A0, A1 = np.asarray(a["cut"][0], float), np.asarray(a["cut"][1], float)
    B0, B1 = np.asarray(b["cut"][0], float), np.asarray(b["cut"][1], float)
    cut0 = 0.5 * (A0 + B0)
    cut1 = 0.5 * (A1 + B1)
    largura = float(np.hypot(*(cut1 - cut0)))
    n = max(9, int(round(largura / 4.0)) + 1)
    sta = np.linspace(0.0, largura, n)
    fa = np.linspace(0.0, 1.0, len(a["sta"]))
    fb = np.linspace(0.0, 1.0, len(b["sta"]))
    f = np.linspace(0.0, 1.0, n)
    z = 0.5 * (np.interp(f, fa, a["z"]) + np.interp(f, fb, b["z"]))
    la = float(a["sta"][-1] - a["sta"][0])
    lb = float(b["sta"][-1] - b["sta"][0])
    bank_l = largura * 0.5 * (float(a["lb"]) / la + float(b["lb"]) / lb)
    bank_r = largura * 0.5 * (float(a["rb"]) / la + float(b["rb"]) / lb)
    return {"rs": rs, "cut": (cut0, cut1), "sta": sta, "z": z,
            "lb": bank_l, "rb": bank_r}


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    ext = argv[argv.index("--saida") + 1] if "--saida" in argv else "g02"
    vao = float(argv[argv.index("--vao") + 1]) if "--vao" in argv else 1600.0
    raiz = os.path.dirname(entrada) or "."
    base = os.path.basename(entrada).rsplit(".", 1)[0]
    saida = os.path.join(raiz, f"{base}.{ext}")
    S = ler_secoes(entrada)
    linhas = open(entrada, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    starts = [i for i, l in enumerate(linhas)
              if l.startswith("Type RM Length L Ch R")]
    starts.append(len(linhas))
    rs_to_i = {round(float(s["rs"]), 2): i for i, s in enumerate(S)}
    inserts = {}
    for i in range(len(S) - 1):
        gap = float(S[i]["rs"] - S[i + 1]["rs"])
        if gap <= vao:
            continue
        rs = round(float(S[i]["rs"] - gap / 2.0), 2)
        inserts[i] = (_interp(S[i], S[i + 1], rs), gap / 2.0)

    out = linhas[:starts[0]]
    for i in range(len(S)):
        block = linhas[starts[i]:starts[i + 1]]
        if i in inserts:
            dist = inserts[i][1]
            block[0] = re.sub(r"^(Type RM Length L Ch R\s*=\s*1\s*,\s*[\d.]+,).*$",
                              lambda m: m.group(1) +
                              f"{dist:8.2f},{dist:8.2f},{dist:8.2f}",
                              block[0])
        out += block
        if i in inserts:
            q, dist = inserts[i]
            out += [""]
            out += _block(q, dist)
    escrever(saida, "\n".join(out))
    print(f"entrada: {entrada}")
    print(f"saida  : {saida}")
    print(f"vaos graves preenchidos: {len(inserts)}")
    print("RS inseridas: " + ", ".join(f"{q['rs']:.2f}" for q, _ in inserts.values()))


if __name__ == "__main__":
    main(sys.argv[1:])
