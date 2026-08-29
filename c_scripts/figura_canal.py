# -*- coding: utf-8 -*-
"""Desenha o Canal Retificado do Itajai-Mirim contra o curso do modelo.

    python scripts/figura_canal.py modelo/mirim_t30/mirim_t30.g01

O modelo roteia os ultimos 19 km pelo curso meandrico. O Canal Retificado --
7,55 km, retilineo, encaixado 3 m no terreno -- nao esta na geometria. Esta
figura poe os dois lado a lado sobre o MDT e mostra a secao tipica de cada um.

O geojson do canal vem do OpenStreetMap (`waterway=canal`, vias 290409480 e
290766868, "Canal Retificado Rio Itajai-Mirim"), em WGS84; e reprojetado para
EPSG:31982 aqui.
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
from shapely.geometry import LineString, Point
from shapely.ops import linemerge

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_secoes import ler_secoes       # noqa: E402
from qc_geometria import ler_eixos     # noqa: E402
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio   # noqa: E402

GEOJSON = r"C:\Users\haas\Downloads\canal_itajai_mirim.geojson"
TERRENO = "modelo/Terrain/MDT_SIGSC_30m.tif"
MEIA = 200.0


def ler_canal(p=GEOJSON):
    tr = Transformer.from_crs("EPSG:4326", "EPSG:31982", always_xy=True)
    d = json.load(open(p, encoding="utf-8"))
    segs = []
    for f in d["features"]:
        a = np.asarray(f["geometry"]["coordinates"], float)
        x, y = tr.transform(a[:, 0], a[:, 1])
        segs.append(LineString(np.c_[x, y]))
    c = linemerge(segs)
    return max(c.geoms, key=lambda g: g.length) if hasattr(c, "geoms") else c


def perfis(linha, mdt, s0, s1, passo=250.0, meia=MEIA):
    s = np.arange(s0, s1, passo)
    off = np.arange(-meia, meia + 1, 4.0)
    pts = []
    for si in s:
        p0 = np.array(linha.interpolate(max(si - 20, s0)).coords[0])
        p1 = np.array(linha.interpolate(min(si + 20, s1)).coords[0])
        t = p1 - p0
        t /= max(float(np.hypot(*t)), 1e-9)
        n = np.array([-t[1], t[0]])
        c = np.array(linha.interpolate(si).coords[0])
        for o in off:
            pts.append(c + o * n)
    pts = np.array(pts)
    z = mdt.cota(pts[:, 0], pts[:, 1]).reshape(len(s), len(off))
    return off, z


def main(argv):
    geom = argv[0] if argv else "modelo/mirim_t30/mirim_t30.g01"
    saida = os.path.join("doc", "canal_retificado.png")
    os.makedirs("doc", exist_ok=True)

    canal = ler_canal()
    eixo = list(ler_eixos(geom).values())[0]
    S = ler_secoes(geom)
    S.sort(key=lambda d: -d["rs"])
    sa = eixo.project(Point(*canal.coords[0]))
    sb = eixo.project(Point(*canal.coords[-1]))
    s0, s1 = min(sa, sb), max(sa, sb)

    P = np.vstack([np.asarray(canal.coords),
                   np.asarray([eixo.interpolate(x).coords[0]
                               for x in np.arange(s0, s1, 50.0)])])
    bb = (P[:, 0].min() - 1200, P[:, 1].min() - 1200,
          P[:, 0].max() + 1200, P[:, 1].max() + 1200)
    mdt = MosaicoSigsc(tiles=tiles_do_dominio(bb))
    off, zc = perfis(canal, mdt, 0.0, canal.length)
    _, zm = perfis(eixo, mdt, s0, s1)

    fig = plt.figure(figsize=(15, 8.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.1, 1, 1], hspace=0.28,
                          wspace=0.22)

    # ---------- planta
    ax = fig.add_subplot(gs[:, 0])
    with rasterio.open(TERRENO) as r:
        w = rasterio.windows.from_bounds(*bb, r.transform)
        img = r.read(1, window=w).astype(float)
        img[img < -9998] = np.nan
        ext = rasterio.windows.bounds(w, r.transform)
        ax.imshow(img, extent=(ext[0], ext[2], ext[1], ext[3]),
                  origin="upper", cmap="terrain", vmin=0, vmax=40, alpha=.9)
    for d in S:
        C = np.asarray(d["cut"], float)
        if bb[0] < C[:, 0].mean() < bb[2] and bb[1] < C[:, 1].mean() < bb[3]:
            ax.plot(C[:, 0], C[:, 1], "-", color="#666666", lw=.5, zorder=2)
    e = np.asarray([eixo.interpolate(x).coords[0]
                    for x in np.arange(s0, s1, 25.0)])
    ax.plot(e[:, 0], e[:, 1], "-", color="#1f4fd8", lw=2.6, zorder=4,
            label=f"curso do MODELO  {(s1-s0)/1000:.1f} km")
    c = np.asarray(canal.coords)
    ax.plot(c[:, 0], c[:, 1], "-", color="#d62728", lw=3.0, zorder=5,
            label=f"Canal Retificado  {canal.length/1000:.1f} km  (AUSENTE)")
    ax.plot(*c[0], "o", ms=9, color="#111", zorder=6)
    ax.plot(*c[-1], "o", ms=9, color="#111", zorder=6)
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=10)
    ax.set_title("Os ultimos 19 km: o modelo desce pelos meandros;\n"
                 "o canal retificado corta 61% do caminho e nao esta na "
                 "geometria", fontsize=12)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(labelsize=7)

    # ---------- perfis
    for k, (Z, rot, cor) in enumerate(
            ((zc, "Canal Retificado", "#d62728"),
             (zm, "curso meandrico do modelo", "#1f4fd8"))):
        a = fig.add_subplot(gs[k, 1])
        for i in range(Z.shape[0]):
            a.plot(off, Z[i], "-", color=cor, lw=.6, alpha=.35)
        a.plot(off, np.nanmedian(Z, axis=0), "-", color="#111", lw=2.2,
               label="mediana")
        a.set_title(f"{rot}\n{Z.shape[0]} perfis a cada 250 m", fontsize=9)
        a.set_xlabel("distancia ao eixo (m)", fontsize=8)
        a.set_ylabel("cota (m)", fontsize=8)
        a.set_ylim(-1, 8)
        a.grid(alpha=.3)
        a.tick_params(labelsize=7)
        a.legend(fontsize=7)

    # ---------- tabela
    a = fig.add_subplot(gs[:, 2])
    a.axis("off")

    def stat(Z):
        cen = np.nanmin(Z[:, np.abs(off) <= 20], axis=1)
        mar = np.nanmedian(np.c_[Z[:, off <= -100], Z[:, off >= 100]], axis=1)
        p = mar - cen
        larg = []
        for i in range(Z.shape[0]):
            zz = Z[i]
            if not np.isfinite(zz).any():
                continue
            lim = np.nanmedian(np.r_[zz[off <= -100], zz[off >= 100]]) - 1.0
            dd = off[np.isfinite(zz) & (zz < lim)]
            if len(dd) > 1:
                larg.append(dd.max() - dd.min())
        return np.nanmedian(p), (np.median(larg) if larg else np.nan)

    pc, lc = stat(zc)
    pm, lm = stat(zm)
    txt = (
        "                    CANAL      MEANDROS\n"
        "                 retificado   (o modelo)\n"
        "  ------------------------------------------\n"
        f"  comprimento      {canal.length/1000:6.2f} km   {(s1-s0)/1000:6.2f} km\n"
        f"  sinuosidade        1.00        {(s1-s0)/np.hypot(*(np.array(canal.coords[-1])-np.array(canal.coords[0]))):4.2f}\n"
        f"  encaixe no MDT   {pc:6.2f} m    {pm:6.2f} m\n"
        f"  largura          {lc:6.0f} m    {lm:6.0f} m\n"
        "  ------------------------------------------\n"
        "  secoes do modelo      0          ~190\n"
        "\n"
        "  O canal encurta 11,7 km (61%).\n"
        "  Os dois estao encaixados no MDT:\n"
        "  o escoamento BIFURCA, e o modelo\n"
        "  so tem o ramo mais longo e raso.\n"
        "\n"
        "  Ressalva: o MDT nao ve dentro\n"
        "  d'agua. Os dois fundos saem em\n"
        "  ~0,05 m, que e a lamina, nao o\n"
        "  leito -- os encaixes sao PISO,\n"
        "  nao a profundidade real."
    )
    a.text(0, 1, txt, family="monospace", fontsize=9.2, va="top")

    fig.savefig(saida, dpi=135, bbox_inches="tight")
    print(f"figura: {saida}")
    return saida


if __name__ == "__main__":
    main(sys.argv[1:])
