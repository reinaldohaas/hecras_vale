# -*- coding: utf-8 -*-
"""Mede a largura REAL dos rios pelo poligono RIOS_DUPLOS da FBDS.

    python scripts/largura_de_poligono.py taha_ai.g01 \
        --fbds doc/fbds --cada 500 --saida doc/larguras_fbds

SO MEDE. Um CSV por rio, mesmo formato dos do SIG-SC.

COMO

  Em cada transecto perpendicular ao eixo (passo `--cada`), a largura e o
  comprimento da intersecao do transecto com o poligono de rio (margem
  dupla) mais proximo do eixo -- o pedaco de intersecao MAIS PERTO do
  eixo, para nao somar meandros vizinhos cortados pelo mesmo transecto.

  ARROZEIRA NAO ENTRA por construcao: o mapeamento da FBDS so poligoniza
  o rio. REPRESA e LAGO saem censurados ("massa dagua"): transecto cujo
  ponto do eixo cai numa MASSAS_DAGUA. Onde o rio nao tem margem dupla
  (mais estreito que ~10 m no RapidEye) sai "sem poligono" -- que ja e
  informacao: o rio ali e mais estreito que 10 m.
"""
import csv
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_geometria import ler_eixos                     # noqa: E402

MEIA = 250.0
PERTO = 100.0    # m: poligono de rio vale se chega a isto do eixo


def carregar(pasta, camada):
    import geopandas as gpd
    import glob
    partes = []
    for shp in sorted(glob.glob(os.path.join(pasta, "*", f"*_{camada}.shp"))):
        try:
            g = gpd.read_file(shp)
        except Exception as e:
            print(f"   pulando {shp}: {e}")
            continue
        if g.crs is None:
            g = g.set_crs(31982)
        partes.append(g.to_crs(31982)[["geometry"]])
    if not partes:
        raise SystemExit(f"nenhum shapefile *_{camada}.shp em {pasta}/*/")
    import pandas as pd
    todos = gpd.GeoDataFrame(pd.concat(partes, ignore_index=True),
                             crs=31982)
    print(f"   {camada}: {len(todos)} poligonos")
    return todos


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    from shapely.geometry import LineString, Point
    from shapely.strtree import STRtree

    g01 = argv[0]
    pasta = argv[argv.index("--fbds") + 1] if "--fbds" in argv \
        else "doc/fbds"
    cada = float(argv[argv.index("--cada") + 1]) if "--cada" in argv \
        else 500.0
    saida = argv[argv.index("--saida") + 1] if "--saida" in argv \
        else "doc/larguras_fbds"
    os.makedirs(saida, exist_ok=True)

    rios_pol = carregar(pasta, "RIOS_DUPLOS")
    massas = carregar(pasta, "MASSAS_DAGUA")
    arv_rio = STRtree(rios_pol.geometry.values)
    arv_mas = STRtree(massas.geometry.values)

    E = ler_eixos(g01)
    por_rio = {}
    for (rio, reach), ls in E.items():
        por_rio.setdefault(rio, []).append((reach, ls))

    print(f"\n{'rio':16s} {'transectos':>10s} {'com pol.':>8s} "
          f"{'largura med':>11s} {'p90':>6s}")
    for rio, partes in sorted(por_rio.items()):
        partes.sort(key=lambda t: t[0])
        coords = []
        for _, ls in partes:
            coords += list(ls.coords)
        eixo = LineString(coords)
        L = eixo.length
        linhas, larguras = [], []
        for s in np.arange(200, L - 200, cada):
            P0 = np.asarray(eixo.interpolate(s).coords[0])
            P1 = np.asarray(eixo.interpolate(min(s + 30, L)).coords[0])
            t = P1 - P0
            t = t / max(np.hypot(*t), 1e-9)
            nvec = np.array([-t[1], t[0]])
            dk = f"{(L - s) / 1000:.2f}"
            p = Point(P0)

            im = arv_mas.query(p, predicate="intersects")
            if len(im):
                linhas.append([dk, "", "massa dagua"])
                continue

            ir = arv_rio.query(p.buffer(PERTO), predicate="intersects")
            if not len(ir):
                linhas.append([dk, "", "sem poligono"])
                continue
            geoms = [rios_pol.geometry.values[k] for k in ir]
            pol = min(geoms, key=lambda g: g.distance(p))
            tran = LineString([P0 - MEIA * nvec, P0 + MEIA * nvec])
            corte = tran.intersection(pol)
            if corte.is_empty:
                linhas.append([dk, "", "sem poligono"])
                continue
            pedacos = getattr(corte, "geoms", [corte])
            perto = min(pedacos, key=lambda g: g.distance(p))
            w = float(perto.length)
            linhas.append([dk, f"{w:.1f}", ""])
            larguras.append(w)

        arq = os.path.join(saida, f"{rio}.csv")
        with open(arq, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(["dist_foz_km", "largura_fbds_m", "obs"])
            w.writerows(linhas)
        if larguras:
            print(f"{rio:16s} {len(linhas):10d} {len(larguras):8d} "
                  f"{np.median(larguras):10.0f}m {np.percentile(larguras, 90):5.0f}m")
        else:
            print(f"{rio:16s} {len(linhas):10d} {0:8d} {'-':>11s} {'-':>6s}")
    print(f"\nCSVs em {saida}/")


if __name__ == "__main__":
    main(sys.argv[1:])
