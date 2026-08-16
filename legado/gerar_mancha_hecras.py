# -*- coding: utf-8 -*-
"""
Gera a MANCHA DE INUNDACAO a partir dos resultados do HEC-RAS.

Substitui o HAND sintetico do app por hidraulica de verdade:

    .p01.hdf  ->  cota d'agua por secao  ->  superficie d'agua interpolada
              ->  subtracao do DEM  ->  profundidade  ->  poligonos

Diferencas para o `manchas_inundacao_hand_dinamico.geojson` atual do app:
  - a cota vem do solver do HEC-RAS, nao de uma cota sintetica;
  - `area_km2` e a area REAL do poligono (no arquivo atual a soma dava
    15.111 km2, mais que a bacia inteira);
  - a profundidade tem gradiente: a franja rasa passa a ser a MAIOR classe,
    como numa cheia real (no arquivo atual 98% da area era "> 2,5 m");
  - so entra area conectada hidraulicamente ao canal (remove pocas isoladas
    em encostas).

Saida (mesmo schema do app, para leitura direta pelo index.html):
    app/manchas_inundacao_hecras.geojson

Uso:  python gerar_mancha_hecras.py
"""
import json
import os
import numpy as np
import rasterio
from rasterio import features, windows
from scipy.spatial import cKDTree
from scipy import ndimage
from shapely.geometry import LineString, shape, mapping
import h5py

import sys
# Projeto de origem. Padrao: a cheia de 1983, que e o evento de referencia
# do vale (regua de Blumenau em 15,34 m) e o que valida melhor:
#   Blumenau +5,4% | foz -6,8% contra a vazao observada.
PROJECT   = sys.argv[1] if len(sys.argv) > 1 else "Itajai_Rede_1983"
DEM       = "dem_itajai.tif"
UTM_EPSG  = 31982
SAIDA     = os.path.join("app", "manchas_inundacao_hecras.geojson")

BUFFER    = 2500.0     # corredor em torno do rio (m). DEVE casar com o
                       # HALFWIDTH de gerar_rede_hecras: um modelo 1D so tem
                       # validade dentro da secao que ele roteou -- extrapolar
                       # alem disso inventa area, e cortar aquem dela descarta
                       # area que o modelo de fato calculou.
DMIN      = 0.10       # profundidade minima considerada inundada (m)
PASSOS    = None      # None = escolhe automaticamente 8 instantes
                      # distribuidos ate o pico e um pouco alem
CLASSES = [                              # (limite inf, limite sup, nome, cor)
    (0.1, 1.0, "Lâmina Baixa (0.1 - 1.0m)",  "#38bdf8"),
    (1.0, 2.5, "Lâmina Média (1.0 - 2.5m)",  "#0284c7"),
    (2.5, 1e9, "Lâmina Severa (> 2.5m)",     "#1e3a8a"),
]


def ler_resultados():
    """Cota d'agua por secao e por instante + coordenadas do talvegue."""
    with h5py.File(f"{PROJECT}.p01.hdf", "r") as f:
        g = f["Results/Unsteady/Output/Output Blocks/Base Output/"
              "Unsteady Time Series/Cross Sections"]
        ws = g["Water Surface"][:]
        at = f["Geometry/Cross Sections/Attributes"][:]
        riv = [r["River"].decode().strip() for r in at]
        rch = [r["Reach"].decode().strip() for r in at]
        rs = np.array([float(r["RS"].decode()) for r in at])
    return ws, riv, rch, rs


def eixos_do_g01():
    """Eixo (Reach XY) de cada trecho do .g01, em UTM."""
    txt = open(f"{PROJECT}.g01", encoding="utf-8", errors="ignore").read().splitlines()
    eixos, i = {}, 0
    while i < len(txt):
        if txt[i].startswith("River Reach="):
            p = txt[i].split("=", 1)[1].split(",")
            rio, rea = p[0].strip(), p[1].strip()
            n = int(txt[i + 1].split("=")[1]); v = []; j = i + 2
            while len(v) < 2 * n and j < len(txt) and not txt[j].startswith("Type RM"):
                v += [float(x) for x in txt[j].split()]
                j += 1
            eixos[(rio, rea)] = LineString(list(zip(v[0::2], v[1::2])))
            i = j
            continue
        i += 1
    return eixos


def pontos_com_cota(ws, riv, rch, rs, eixos, t):
    """Densifica cada eixo e carrega a cota d'agua do instante t interpolada
    pela estaca (RS)."""
    px, py, pz = [], [], []
    for (rio, rea), ln in eixos.items():
        idx = [k for k in range(len(rs)) if riv[k] == rio and rch[k] == rea]
        if len(idx) < 2:
            continue
        idx.sort(key=lambda k: rs[k])
        rr = rs[idx]                      # crescente (jusante -> montante)
        zz = ws[t, idx]
        # o RS e medido a partir da foz -> distancia ao longo do eixo
        rs_max = rr.max()
        for s in np.arange(0, ln.length, 60.0):
            p = ln.interpolate(s)
            rs_p = rs_max - s if ln.length >= rs_max - rr.min() else rs_max - s
            px.append(p.x); py.append(p.y)
            pz.append(float(np.interp(np.clip(rs_p, rr.min(), rr.max()), rr, zz)))
    return np.array(px), np.array(py), np.array(pz)


def main():
    ws, riv, rch, rs = ler_resultados()
    eixos = eixos_do_g01()
    print(f"{ws.shape[0]} instantes, {ws.shape[1]} secoes, {len(eixos)} trechos")

    # --- recorte do DEM ao corredor dos rios
    todos = [p for ln in eixos.values() for p in ln.coords]
    xs = np.array([p[0] for p in todos]); ys = np.array([p[1] for p in todos])
    ds = rasterio.open(DEM)
    from pyproj import Transformer
    tr = Transformer.from_crs(UTM_EPSG, ds.crs.to_epsg(), always_xy=True)
    lo, la = tr.transform(xs, ys)
    # NAO usar rasterio.windows: from_bounds() derruba o processo nesta versao.
    # Le a banda inteira (isso funciona) e recorta com numpy pela transform.
    full = ds.read(1)
    T = ds.transform
    inv = ~T
    c0, r0 = inv * (lo.min() - 0.05, la.max() + 0.05)
    c1, r1 = inv * (lo.max() + 0.05, la.min() - 0.05)
    c0 = max(int(np.floor(c0)), 0); r0 = max(int(np.floor(r0)), 0)
    c1 = min(int(np.ceil(c1)), full.shape[1]); r1 = min(int(np.ceil(r1)), full.shape[0])
    dem = full[r0:r1, c0:c1].astype(float)
    if ds.nodata is not None:
        dem[dem == ds.nodata] = np.nan
    dem[dem < -500] = np.nan
    tw = T * __import__("affine").Affine.translation(c0, r0)
    print(f"recorte do DEM: {dem.shape}")

    # grade de coordenadas UTM dos pixels
    rows, cols = np.indices(dem.shape)
    lon = tw.c + (cols + 0.5) * tw.a
    lat = tw.f + (rows + 0.5) * tw.e
    tr2 = Transformer.from_crs(ds.crs.to_epsg(), UTM_EPSG, always_xy=True)
    gx, gy = tr2.transform(lon, lat)

    # instantes exportados: do inicio da subida ate depois do pico
    passos = PASSOS
    if passos is None:
        qmax = ws.max(axis=1)
        hp = int(np.argmax(qmax))
        passos = sorted(set(int(round(x)) for x in
                            np.linspace(max(hp - 36, 6), min(hp + 12, ws.shape[0] - 1), 8)))
        print(f"pico da lamina em h={hp}; instantes exportados: {passos}")

    feats = []
    for h in passos:
        t = min(h, ws.shape[0] - 1)
        px, py, pz = pontos_com_cota(ws, riv, rch, rs, eixos, t)
        tree = cKDTree(np.c_[px, py])
        d, k = tree.query(np.c_[gx.ravel(), gy.ravel()], workers=-1)
        wse = pz[k].reshape(dem.shape)
        dist = d.reshape(dem.shape)

        prof = wse - dem
        prof[dist > BUFFER] = -9999          # fora do corredor
        prof[~np.isfinite(dem)] = -9999
        molhado = prof > DMIN

        # -- conectividade: so mantem area ligada ao canal
        lab, n = ndimage.label(molhado)
        canal = np.unique(lab[(dist < 120) & molhado])
        canal = canal[canal > 0]
        conect = np.isin(lab, canal)
        prof = np.where(conect, prof, -9999)

        area_px = abs(tw.a * tw.e) * (111320.0 ** 2) * np.cos(np.radians(lat.mean()))
        for lo_, hi_, nome, cor in CLASSES:
            m = (prof >= lo_) & (prof < hi_)
            if not m.any():
                continue
            area = float(m.sum() * area_px / 1e6)
            vol = float(np.where(m, prof, 0).sum() * area_px / 1e6)
            for geom, val in features.shapes(m.astype(np.uint8), mask=m,
                                             transform=tw):
                if val != 1:
                    continue
                g = shape(geom)
                if g.area * (111320.0 ** 2) < 20000:     # descarta ruido
                    continue
                feats.append({
                    "type": "Feature",
                    "properties": {"event": "sim", "time_hour": h,
                                   "class_name": nome, "fill_color": cor,
                                   "fill_opacity": 0.5,
                                   "area_km2": round(area, 2),
                                   "volume_hm3": round(vol, 2),
                                   "fonte": "HEC-RAS " + PROJECT},
                    "geometry": mapping(g)})
            print(f"  t={h:>2}h  {nome:<30} {area:8.2f} km2   vol {vol:8.2f} hm3")

    with open(SAIDA, "w", encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": feats}, f)
    print(f"\n[OK] {SAIDA}  ({len(feats)} poligonos)")


if __name__ == "__main__":
    main()
