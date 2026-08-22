# -*- coding: utf-8 -*-
"""Desenha a estrutura de reaches e juncoes de uma ou mais geometrias.

    python scripts/figura_estrutura.py modelo/mirim_t30/mirim_t30.g04 \
        modelo/mirim_t30/mirim_t30.g05

Le os `River Reach=` / `Reach XY=` / `Junct` direto do .gNN e pinta um reach
por cor, sobre o MDT. Serve para conferir a topologia sem depender do
RAS Mapper.
"""
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rasterio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes   # noqa: E402

TERRENO = "modelo/Terrain/MDT_SIGSC_30m.tif"
L16 = 16
CORES = {"R1": "#1f4fd8", "R2": "#e08a00", "R3": "#7b1fa2",
         "Canal": "#d62728"}


def estrutura(g):
    t = open(g, encoding="latin-1", errors="replace").read() \
        .replace("\r", "").split("\n")
    reaches, juncs = [], []
    i = 0
    while i < len(t):
        l = t[i]
        if l.startswith("Junct Name="):
            nome = l.split("=", 1)[1].strip()
            xy = None
            j = i + 1
            up, dn = [], []
            while j < len(t) and not t[j].startswith(("Junct Name=",
                                                      "River Reach=")):
                if t[j].startswith("Junct X Y & Text X Y="):
                    v = [float(x) for x in t[j].split("=", 1)[1].split(",")]
                    xy = (v[0], v[1])
                if t[j].startswith("Up River,Reach="):
                    up.append(t[j].split(",")[-1].strip())
                if t[j].startswith("Dn River,Reach="):
                    dn.append(t[j].split(",")[-1].strip())
                j += 1
            juncs.append((nome, xy, up, dn))
            i = j
            continue
        if l.startswith("River Reach="):
            nome = l.split("=", 1)[1].split(",")[1].strip()
            n = int(t[i + 1].split("=")[1])
            v, j = [], i + 2
            while len(v) < 2 * n:
                v += [float(t[j][c:c + L16]) for c in range(0, len(t[j]), L16)
                      if t[j][c:c + L16].strip()]
                j += 1
            reaches.append((nome, np.array(v).reshape(-1, 2)))
            i = j
            continue
        i += 1
    return reaches, juncs


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    fig, axs = plt.subplots(1, len(argv), figsize=(7.6 * len(argv), 8.4))
    if len(argv) == 1:
        axs = [axs]

    # janela comum: o trecho das juncoes
    todas = [estrutura(g) for g in argv]
    P = np.vstack([xy for r, _ in todas for _, xy in r])
    J = np.array([j[1] for _, js in todas for j in js if j[1]])
    x0, x1 = J[:, 0].min() - 5000, J[:, 0].max() + 3000
    y0, y1 = J[:, 1].min() - 4000, J[:, 1].max() + 3000

    for ax, g, (reaches, juncs) in zip(axs, argv, todas):
        with rasterio.open(TERRENO) as r:
            w = rasterio.windows.from_bounds(x0, y0, x1, y1, r.transform)
            img = r.read(1, window=w).astype(float)
            img[img < -9998] = np.nan
            ext = rasterio.windows.bounds(w, r.transform)
            ax.imshow(img, extent=(ext[0], ext[2], ext[1], ext[3]),
                      origin="upper", cmap="terrain", vmin=0, vmax=45,
                      alpha=.9)
        S = ler_secoes(g)
        for d in S:
            C = np.asarray(d["cut"], float)
            if x0 < C[:, 0].mean() < x1 and y0 < C[:, 1].mean() < y1:
                ax.plot(C[:, 0], C[:, 1], "-", color="#555", lw=.45, zorder=2)
        for nome, xy in reaches:
            ax.plot(xy[:, 0], xy[:, 1], "-", lw=3.0,
                    color=CORES.get(nome, "#333"), zorder=4,
                    label=f"{nome}  ({len(xy)} vert.)")
        for nome, xy, up, dn in juncs:
            if not xy:
                continue
            ax.plot(*xy, "o", ms=15, mfc="white", mec="#111", mew=2.4,
                    zorder=6)
            ax.annotate(f"{nome.strip()}\n{len(up)} entra / {len(dn)} sai",
                        xy, xytext=(16, 12), textcoords="offset points",
                        fontsize=9, fontweight="bold", zorder=7,
                        bbox=dict(fc="white", alpha=.9, lw=.5))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ext_g = os.path.basename(g).split(".")[-1]
        ok = all(len(u) + len(d) > 2 for _, _, u, d in juncs)
        ax.set_title(f"{os.path.basename(g)}   ({len(reaches)} reaches, "
                     f"{len(juncs)} juncoes)\n"
                     + ("computa" if ok else
                        "NAO computa: juncao com 1 entra / 1 sai"),
                     fontsize=12,
                     color="#1a7f37" if ok else "#a01015")
        ax.legend(loc="lower right", fontsize=9)
        ax.ticklabel_format(style="plain", useOffset=False)
        ax.tick_params(labelsize=7)

    os.makedirs("doc", exist_ok=True)
    p = "doc/estrutura_juncoes.png"
    fig.savefig(p, dpi=130, bbox_inches="tight")
    print("figura:", p)
    return p


if __name__ == "__main__":
    main(sys.argv[1:])
