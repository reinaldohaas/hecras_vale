# -*- coding: utf-8 -*-
"""Confronta a largura de canal DO MODELO com a largura MEDIDA no SIG-SC.

    python scripts/comparar_larguras.py taha_ai.g01 \
        --medidas doc/larguras_sigsc --figura doc/figuras/larguras_modelo_vs_sigsc.png

SO DESENHA. Um painel por rio: largura entre Bank Sta do g01 (linha) contra
a lamina d'agua medida no MDT 1 m (pontos). A lamina e o PISO da largura
real; modelo muito acima da lamina = canal inventado largo demais.

A calha plena medida so entra onde e crivel (ate 3x a lamina local):
em planicie de barranco baixo (<1 m) a regra da quebra nao dispara e a
medida vira o transecto inteiro -- fora essas, e descartada.
"""
import csv
import os
import re
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)


def larguras_do_g01(g01):
    """rio -> [(rs_km, largura_bank_sta)]"""
    out = {}
    rio = None
    rs = None
    for l in open(g01, encoding="latin-1", errors="replace"):
        if l.startswith("River Reach="):
            rio = l.split("=")[1].split(",")[0].strip()
        elif l.startswith("Type RM Length L Ch R ="):
            try:
                rs = float(l.split("=")[1].split(",")[1])
            except (ValueError, IndexError):
                rs = None
        elif l.startswith("Bank Sta=") and rio and rs is not None:
            try:
                a, b = (float(x) for x in l.split("=")[1].split(","))
                out.setdefault(rio, []).append((rs / 1000.0, b - a))
            except ValueError:
                pass
    return out


def larguras_medidas(pasta, rio):
    """[(dist_km, lamina, calha_ou_nan)]"""
    arq = os.path.join(pasta, f"{rio}.csv")
    if not os.path.exists(arq):
        return []
    out = []
    for r in csv.reader(open(arq, encoding="utf-8"), delimiter=";"):
        if r[0] == "dist_foz_km" or not r[1]:
            continue
        la = float(r[1])
        ca = float(r[2])
        out.append((float(r[0]), la, ca if ca <= 3 * la else np.nan))
    return out


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g01 = argv[0]
    pasta = argv[argv.index("--medidas") + 1] if "--medidas" in argv \
        else "doc/larguras_sigsc"
    fig_arq = argv[argv.index("--figura") + 1] if "--figura" in argv \
        else "doc/figuras/larguras_modelo_vs_sigsc.png"

    modelo = larguras_do_g01(g01)
    rios = sorted(modelo)
    n = len(rios)
    ncol = 3
    nlin = -(-n // ncol)
    fig, eixos = plt.subplots(nlin, ncol, figsize=(15, 3.2 * nlin))
    eixos = np.atleast_2d(eixos)

    print(f"{'rio':16s} {'modelo med':>10s} {'lamina med':>10s} "
          f"{'modelo/lamina':>13s}")
    for k, rio in enumerate(rios):
        ax = eixos[k // ncol][k % ncol]
        mm = sorted(modelo[rio])
        xs = [p[0] for p in mm]
        ws = [p[1] for p in mm]
        ax.plot(xs, ws, "-", color="crimson", lw=1.2,
                label="modelo (Bank Sta)")
        med = larguras_medidas(pasta, rio)
        if med:
            ax.plot([p[0] for p in med], [p[1] for p in med], ".",
                    color="royalblue", ms=3, label="lamina SIG-SC")
            ca = [(p[0], p[2]) for p in med if not np.isnan(p[2])]
            if ca:
                ax.plot([p[0] for p in ca], [p[1] for p in ca], ".",
                        color="seagreen", ms=3, alpha=0.6,
                        label="calha plena (crivel)")
        arq_fbds = os.path.join("doc/larguras_fbds", f"{rio}.csv")
        if os.path.exists(arq_fbds):
            fb = [(float(r[0]), float(r[1]))
                  for r in csv.reader(open(arq_fbds, encoding="utf-8"),
                                      delimiter=";")
                  if r[0] != "dist_foz_km" and r[1]]
            if fb:
                ax.plot([p[0] for p in fb], [p[1] for p in fb], ".",
                        color="darkorange", ms=3, alpha=0.8,
                        label="poligono FBDS")
            m_mod = float(np.median(ws))
            m_lam = float(np.median([p[1] for p in med]))
            print(f"{rio:16s} {m_mod:10.0f} {m_lam:10.0f} "
                  f"{m_mod / max(m_lam, 1e-9):13.1f}x")
        ax.set_title(rio, fontsize=10)
        ax.set_xlabel("dist. da foz (km)", fontsize=8)
        ax.set_ylabel("largura (m)", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)
        if k == 0:
            ax.legend(fontsize=8)
    for k in range(n, nlin * ncol):
        eixos[k // ncol][k % ncol].axis("off")
    fig.suptitle("Largura do canal: modelo (formula de regime) vs "
                 "medido no MDT 1 m do SIG-SC", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(os.path.dirname(fig_arq) or ".", exist_ok=True)
    fig.savefig(fig_arq, dpi=130)
    print(f"\nfigura: {fig_arq}")


if __name__ == "__main__":
    main(sys.argv[1:])
