# -*- coding: utf-8 -*-
"""Mostra onde o mirim_novo termina e o que falta ate a foz.

    python scripts/figura_fim_do_canal.py

O `mirim_novo` acaba exatamente no ultimo vertice do Canal Retificado do
OpenStreetMap -- distancia ZERO --, e isso fica 1.058 m antes da confluencia
com o Itajai-Acu. No corredor entre os dois o MDT mostra um aterro de ~2,3 m
atravessando, compativel com a BR-101: o minimo dentro de +-50 m do eixo pula
de 0,05 m para 2,32 m e nao volta.

A consequencia nao e cosmetica: o contorno de jusante (hidrograma de mare do
Itajai-Acu) passa a ser imposto 1 km rio acima de onde a mare de fato entra.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rasterio
from pyproj import Transformer
from shapely.geometry import LineString
from shapely.ops import linemerge

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes                        # noqa: E402
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio     # noqa: E402

TERRENO = "modelo/Terrain/MDT_SIGSC_30m.tif"


def main():
    tr = Transformer.from_crs("EPSG:4326", "EPSG:31982", always_xy=True)
    d = json.load(open(r"C:\Users\haas\Downloads\canal_itajai_mirim.geojson",
                       encoding="utf-8"))
    segs = []
    for f in d["features"]:
        a = np.asarray(f["geometry"]["coordinates"], float)
        x, y = tr.transform(a[:, 0], a[:, 1])
        segs.append(LineString(np.c_[x, y]))
    canal = linemerge(segs)
    if hasattr(canal, "geoms"):
        canal = max(canal.geoms, key=lambda g: g.length)
    A = np.array(canal.coords[-1])

    So = ler_secoes("modelo/so_mirim.g01")
    So.sort(key=lambda x: -x["rs"])
    C = np.asarray(So[-1]["cut"], float)
    B = 0.5 * (C[0] + C[-1])

    Sn = ler_secoes("modelo/mirim_novo/mirim_novo.g01")
    Sn.sort(key=lambda x: -x["rs"])

    u = (B - A) / np.hypot(*(B - A))
    L = float(np.hypot(*(B - A)))
    n = np.array([-u[1], u[0]])
    s = np.arange(-300, L + 301, 5.0)
    off = np.arange(-300, 301, 5.0)
    pts = np.array([A + si * u + o * n for si in s for o in off])
    mdt = MosaicoSigsc(tiles=tiles_do_dominio(
        (pts[:, 0].min() - 40, pts[:, 1].min() - 40,
         pts[:, 0].max() + 40, pts[:, 1].max() + 40)))
    Z = mdt.cota(pts[:, 0], pts[:, 1]).reshape(len(s), len(off))
    fundo = np.array([np.nanmin(z[np.abs(off) <= 50])
                      if np.isfinite(z[np.abs(off) <= 50]).any() else np.nan
                      for z in Z])

    fig, (ax, a2) = plt.subplots(1, 2, figsize=(15, 7),
                                 gridspec_kw={"width_ratios": [1.25, 1]})
    bb = (min(A[0], B[0]) - 900, min(A[1], B[1]) - 700,
          max(A[0], B[0]) + 700, max(A[1], B[1]) + 700)
    with rasterio.open(TERRENO) as r:
        w = rasterio.windows.from_bounds(*bb, r.transform)
        img = r.read(1, window=w).astype(float)
        img[img < -9998] = np.nan
        e = rasterio.windows.bounds(w, r.transform)
        im = ax.imshow(img, extent=(e[0], e[2], e[1], e[3]), origin="upper",
                       cmap="terrain", vmin=0, vmax=12)
    plt.colorbar(im, ax=ax, shrink=.7, label="cota (m)")
    c = np.asarray(canal.coords)
    ax.plot(c[:, 0], c[:, 1], "-", color="#d62728", lw=3,
            label="Canal Retificado (OSM)")
    for x in Sn[-6:]:
        P = np.asarray(x["cut"], float)
        ax.plot(P[:, 0], P[:, 1], "-", color="#2e7d32", lw=1.6)
    ax.plot([], [], "-", color="#2e7d32", lw=1.6, label="secoes do mirim_novo")
    for x in So[-2:]:
        P = np.asarray(x["cut"], float)
        ax.plot(P[:, 0], P[:, 1], "-", color="#1f4fd8", lw=1.8)
    ax.plot([], [], "-", color="#1f4fd8", lw=1.8,
            label="secoes do original ate a foz")
    ax.plot(*A, "o", ms=13, mfc="#fff", mec="#111", mew=2.2, zorder=6)
    ax.annotate("fim do mirim_novo\n= ultimo vertice do OSM", A,
                xytext=(-150, -70), textcoords="offset points", fontsize=9,
                bbox=dict(fc="white", alpha=.9, lw=.4),
                arrowprops=dict(arrowstyle="->", lw=1.2))
    ax.plot(*B, "o", ms=13, mfc="#fff", mec="#111", mew=2.2, zorder=6)
    ax.annotate("foz no Itajai-Acu\n1.058 m adiante", B,
                xytext=(18, 26), textcoords="offset points", fontsize=9,
                bbox=dict(fc="white", alpha=.9, lw=.4),
                arrowprops=dict(arrowstyle="->", lw=1.2))
    k = int(np.nanargmax(np.where((s > 400) & (s < L), fundo, np.nan)))
    Pb = A + s[k] * u
    ax.plot(*Pb, "s", ms=11, color="#111", zorder=6)
    ax.annotate("aterro, %.1f m" % fundo[k], Pb, xytext=(14, -34),
                textcoords="offset points", fontsize=9, fontweight="bold",
                color="#a01015", bbox=dict(fc="white", alpha=.9, lw=.4))
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title("O modelo para onde acaba o vetor do canal,\n"
                 "e um aterro cruza o que falta", fontsize=12)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(labelsize=7)

    a2.plot(s, fundo, "-", color="#1f4fd8", lw=2,
            label="cota MINIMA a +-50 m do eixo")
    a2.axvline(0, color="#d62728", lw=2)
    a2.annotate("fim do mirim_novo", (0, 3.2), rotation=90, fontsize=9,
                color="#d62728", ha="right", va="top")
    a2.axvline(L, color="#1f4fd8", lw=2, ls="--")
    a2.annotate("foz", (L, 3.2), rotation=90, fontsize=9, color="#1f4fd8",
                ha="right", va="top")
    a2.axhspan(-0.2, 0.3, color="#4fc3f7", alpha=.35)
    a2.annotate("lamina do canal (~0,05 m)", (-260, 0.35), fontsize=8.5,
                color="#0277bd")
    a2.plot(s[k], fundo[k], "s", ms=11, color="#111")
    a2.annotate("aterro: o fundo do corredor\nsobe de 0,05 m para %.2f m" % fundo[k],
                (s[k], fundo[k]), xytext=(-30, 26),
                textcoords="offset points", fontsize=9, fontweight="bold",
                color="#a01015", bbox=dict(fc="white", alpha=.9, lw=.4))
    a2.set_xlabel("distancia a partir do fim do canal (m)", fontsize=9)
    a2.set_ylabel("cota (m)", fontsize=9)
    a2.grid(alpha=.3)
    a2.legend(fontsize=8.5, loc="lower right")
    a2.set_title("O que ha no quilometro que falta", fontsize=12)

    os.makedirs("doc", exist_ok=True)
    p = "doc/fim_do_canal.png"
    fig.savefig(p, dpi=135, bbox_inches="tight")
    print("figura:", p)
    print("fim do mirim_novo (%.0f, %.0f)   foz (%.0f, %.0f)   %.0f m"
          % (A[0], A[1], B[0], B[1], L))
    print("aterro a %.0f m do fim do canal, fundo do corredor em %.2f m"
          % (s[k], fundo[k]))
    return p


if __name__ == "__main__":
    main()
