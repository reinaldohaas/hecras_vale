# -*- coding: utf-8 -*-
"""Mostra o canal preparado: o que esta pronto e o buraco do leito.

    python scripts/figura_canal_preparado.py

Le `doc/canal/canal_secoes.csv` -- nao reamostra o MDT.
"""
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rasterio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CSV = "doc/canal/canal_secoes.csv"
CUT = "doc/canal/canal_cutlines.geojson"
EIXO = "doc/canal/canal_eixo.geojson"
TERRENO = "modelo/Terrain/MDT_SIGSC_30m.tif"
Z_PONTA_MONT, Z_PONTA_JUS = -0.76, -2.68


def main():
    secs = defaultdict(list)
    for r in csv.DictReader(open(CSV, encoding="utf-8"), delimiter=";"):
        secs[float(r["rs"])].append(r)
    rss = sorted(secs, reverse=True)

    fig = plt.figure(figsize=(15, 8.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.3,
                          wspace=0.2)

    # ---------- planta
    ax = fig.add_subplot(gs[:, 0])
    cut = json.load(open(CUT))
    eixo = np.asarray(json.load(open(EIXO))["features"][0]
                      ["geometry"]["coordinates"], float)
    P = np.vstack([np.asarray(f["geometry"]["coordinates"], float)
                   for f in cut["features"]])
    bb = (P[:, 0].min() - 400, P[:, 1].min() - 400,
          P[:, 0].max() + 400, P[:, 1].max() + 400)
    with rasterio.open(TERRENO) as r:
        w = rasterio.windows.from_bounds(*bb, r.transform)
        img = r.read(1, window=w).astype(float)
        img[img < -9998] = np.nan
        ext = rasterio.windows.bounds(w, r.transform)
        ax.imshow(img, extent=(ext[0], ext[2], ext[1], ext[3]),
                  origin="upper", cmap="terrain", vmin=0, vmax=25, alpha=.92)
    for f in cut["features"]:
        c = np.asarray(f["geometry"]["coordinates"], float)
        ax.plot(c[:, 0], c[:, 1], "-", color="#2e7d32", lw=1.0, zorder=3)
    ax.plot(eixo[:, 0], eixo[:, 1], "-", color="#d62728", lw=2.6, zorder=4)
    ax.plot(*eixo[0], "o", ms=9, color="#111", zorder=5)
    ax.plot(*eixo[-1], "o", ms=9, color="#111", zorder=5)
    ax.annotate(f"montante\nleito medido {Z_PONTA_MONT:+.2f} m", eixo[0],
                xytext=(14, -34), textcoords="offset points", fontsize=8.5,
                bbox=dict(fc="white", alpha=.85, lw=0))
    ax.annotate(f"jusante\nleito medido {Z_PONTA_JUS:+.2f} m", eixo[-1],
                xytext=(-104, 16), textcoords="offset points", fontsize=8.5,
                bbox=dict(fc="white", alpha=.85, lw=0))
    ax.set_aspect("equal")
    ax.set_title(f"Canal Retificado preparado: {len(rss)} secoes a cada 150 m "
                 f"sobre 7,55 km\nplanicie cortada do MDT 1 m; leito em "
                 "aberto", fontsize=11.5)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(labelsize=7)

    # ---------- uma secao
    a = fig.add_subplot(gs[0, 1])
    r0 = rss[len(rss) // 2]
    L = secs[r0]
    st = np.array([float(x["estaca"]) for x in L])
    z = np.array([float(x["z"]) if x["z"] else np.nan for x in L])
    zl = np.array([float(x["z_mdt_lamina"]) if x["z_mdt_lamina"] else np.nan
                   for x in L])
    a.plot(st, z, "-", color="#2e7d32", lw=1.6, label="planicie -- MDT 1 m")
    a.plot(st, zl, "-", color="#4fc3f7", lw=2.4,
           label="lamina no MDT (teto do leito)")
    m = np.isfinite(zl)
    if m.any():
        x0, x1 = st[m].min(), st[m].max()
        a.axvspan(x0, x1, color="#d62728", alpha=.13)
        a.annotate("LEITO EM ABERTO\n(45 m)", (0.5 * (x0 + x1), -1.6),
                   ha="center", fontsize=8.5, color="#a01015",
                   fontweight="bold")
    a.axhline(Z_PONTA_MONT, ls="--", lw=1, color="#888")
    a.axhline(Z_PONTA_JUS, ls="--", lw=1, color="#888")
    a.annotate("faixa das duas cotas de leito medidas", (st.min() + 20, -1.9),
               fontsize=7.5, color="#666")
    a.set_xlim(st.min(), st.max())
    a.set_ylim(-3.4, 8)
    a.set_title(f"secao tipica  (RS {r0:.0f})", fontsize=10)
    a.set_xlabel("estaca (m)", fontsize=8)
    a.set_ylabel("cota (m)", fontsize=8)
    a.grid(alpha=.3)
    a.tick_params(labelsize=7)
    a.legend(fontsize=7.5, loc="upper center")

    # ---------- estado
    a = fig.add_subplot(gs[1, 1])
    a.axis("off")
    npt = len(secs[rss[0]])
    nab = sum(1 for x in secs[rss[0]] if x["origem"] == "A LEVANTAR")
    txt = (
        "  PRONTO                              \n"
        "  ----------------------------------- \n"
        f"  eixo                    7,550 km    \n"
        f"  secoes                  {len(rss):3d}         \n"
        f"  espacamento             150 m       \n"
        f"  pontos por secao        {npt:3d}         \n"
        f"  planicie do MDT 1 m     {npt-nab:3d} pontos  \n"
        "  cobertura do MDT        100%        \n"
        "  Bank Sta (canal 45 m)   posicionada \n"
        "  cutlines e orientacao   ok          \n"
        "\n"
        "  EM ABERTO                           \n"
        "  ----------------------------------- \n"
        f"  cota de leito           {nab:3d} pontos  \n"
        f"                          por secao   \n"
        f"  total a levantar        {nab*len(rss):4d} pontos\n"
        "\n"
        "  ancoras que ja existem:\n"
        f"    montante  {Z_PONTA_MONT:+.2f} m\n"
        f"    jusante   {Z_PONTA_JUS:+.2f} m\n"
        f"    queda      {Z_PONTA_MONT-Z_PONTA_JUS:.2f} m em 7,55 km\n"
        "\n"
        "  DECISAO PENDENTE\n"
        "  ----------------------------------- \n"
        "  extensao lateral: 500 m por lado,\n"
        "  por escolha. Conter a cheia de\n"
        "  6,13 m pediria ~3,3 km -- isso e\n"
        "  varzea, e pede armazenamento\n"
        "  ou 2D, nao secao mais larga."
    )
    a.text(0, 1, txt, family="monospace", fontsize=8.8, va="top")

    os.makedirs("doc", exist_ok=True)
    p = "doc/canal_preparado.png"
    fig.savefig(p, dpi=135, bbox_inches="tight")
    print("figura:", p)
    return p


if __name__ == "__main__":
    main()
