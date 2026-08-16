# -*- coding: utf-8 -*-
"""
Limites da planicie de inundacao, a partir das secoes transversais.

Para cada secao, procura ao longo da cutline onde a lamina d'agua encontra o
terreno -- a borda esquerda e a direita da agua. Ligando as bordas de secoes
consecutivas sai o POLIGONO da mancha, que e o limite da planicie inundada.

Isto nao depende do RAS Mapper nem de terreno importado: usa a mesma geometria
do .g01 (as cutlines ja estao georreferenciadas) e a cota d'agua de qualquer
uma das duas fontes:

    --fonte motor    <PROJETO>_motor.npz    (192 passos, roda sempre)
    --fonte hecras   <PROJETO>.p01.hdf      (ate onde o solver chegou)

Saidas:
    <PROJETO>_planicie.shp        para abrir no RAS Mapper / QGIS (EPSG:31982)
    <PROJETO>_planicie.geojson    em lat/lon, para a interface web

Uso:
    python gerar_planicie.py Itajai_Rede_1983 --fonte motor
    python gerar_planicie.py Itajai_Rede_1983 --fonte motor --instante 120
"""
import argparse
import os

import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union

UTM = 31982
CLASSES = [(0.0, 0.5, "Rasa"), (0.5, 2.0, "Media"),
           (2.0, 5.0, "Profunda"), (5.0, 1e9, "Severa")]


def ler_geometria(projeto):
    """Secoes com estacas, cotas e as PONTAS georreferenciadas da cutline."""
    txt = open(f"{projeto}.g01", encoding="utf-8", errors="ignore").read().splitlines()
    secoes = []
    rio = reach = rs = cut = None
    i = 0
    while i < len(txt):
        l = txt[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio, reach = p[0].strip(), p[1].strip()
        elif l.startswith("Type RM"):
            try:
                rs = float(l.split(",")[1])
            except ValueError:
                rs = None
        elif l.startswith("XS GIS Cut Line"):
            n = int(l.split("=")[1])
            v, j = [], i + 1
            while len(v) < 2 * n and j < len(txt):
                v += [float(x) for x in txt[j].split()]
                j += 1
            cut = np.array(v[:2 * n]).reshape(-1, 2)
            i = j
            continue
        elif l.startswith("#Sta/Elev="):
            n = int(l.split("=")[1])
            v = []
            i += 1
            while i < len(txt) and len(v) < 2 * n:
                s = txt[i]
                v += [float(s[c:c + 8]) for c in range(0, len(s.rstrip()), 8)
                      if s[c:c + 8].strip()]
                i += 1
            if cut is not None and rs is not None:
                secoes.append({"rio": rio, "reach": reach, "rs": rs,
                               "sta": np.array(v[0::2]), "z": np.array(v[1::2]),
                               "cut": cut})
            continue
        i += 1
    return secoes


def ponto_na_cutline(cut, sta, frac):
    """Coordenada UTM na fracao 'frac' (0..1) do comprimento da cutline."""
    d = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(cut[:, 0]),
                                                  np.diff(cut[:, 1])))])
    alvo = frac * d[-1]
    return (float(np.interp(alvo, d, cut[:, 0])),
            float(np.interp(alvo, d, cut[:, 1])))


def bordas(sec, zw):
    """(x,y) das bordas esquerda e direita da agua, e a profundidade maxima.

    A borda e a estaca mais externa, de cada lado do talvegue, onde o terreno
    ainda esta abaixo da lamina. Percorrer do talvegue PARA FORA e nao pegar o
    minimo global evita saltar para outro canal que a secao atravesse.
    """
    sta, z = sec["sta"], sec["z"]
    i0 = int(np.nanargmin(z))
    if z[i0] >= zw:
        return None
    e = i0
    while e > 0 and z[e - 1] <= zw:
        e -= 1
    d = i0
    while d < len(z) - 1 and z[d + 1] <= zw:
        d += 1
    if d <= e:
        return None
    L = sta[-1] - sta[0]
    if L <= 0:
        return None
    pe = ponto_na_cutline(sec["cut"], sta, (sta[e] - sta[0]) / L)
    pd = ponto_na_cutline(sec["cut"], sta, (sta[d] - sta[0]) / L)
    return pe, pd, float(zw - z[i0]), float(sta[d] - sta[e])


def poligonos(secoes, ws, idx_por_secao, instante):
    """Um poligono por trecho, ligando as bordas de secoes consecutivas."""
    por_trecho = {}
    for k, sec in enumerate(secoes):
        j = idx_por_secao.get(k)
        if j is None:
            continue
        b = bordas(sec, float(ws[instante, j]))
        if b is None:
            continue
        por_trecho.setdefault((sec["rio"], sec["reach"]), []).append(
            {"rs": sec["rs"], "e": b[0], "d": b[1], "prof": b[2], "larg": b[3]})
    feats = []
    for (rio, reach), lst in por_trecho.items():
        lst.sort(key=lambda x: -x["rs"])
        if len(lst) < 2:
            continue
        # quebra em blocos contiguos: onde a agua some, o poligono termina
        bloco = [lst[0]]
        blocos = []
        for a, b in zip(lst, lst[1:]):
            if a["rs"] - b["rs"] > 6000.0:      # vao grande = trecho seco
                blocos.append(bloco); bloco = [b]
            else:
                bloco.append(b)
        blocos.append(bloco)
        for bl in blocos:
            if len(bl) < 2:
                continue
            anel = [x["e"] for x in bl] + [x["d"] for x in reversed(bl)]
            p = Polygon(anel)
            if not p.is_valid:
                p = p.buffer(0)
            if p.is_empty or p.area <= 0:
                continue
            prof = float(np.mean([x["prof"] for x in bl]))
            cls = next(c for lo, hi, c in CLASSES if lo <= prof < hi)
            feats.append({"rio": rio, "reach": reach, "classe": cls,
                          "prof_med": round(prof, 2),
                          "prof_max": round(max(x["prof"] for x in bl), 2),
                          "larg_max": round(max(x["larg"] for x in bl), 1),
                          "rs_ini": round(bl[0]["rs"] / 1000, 2),
                          "rs_fim": round(bl[-1]["rs"] / 1000, 2),
                          "geometry": p})
    return feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("projeto", nargs="?", default="Itajai_Rede_1983")
    ap.add_argument("--fonte", choices=["motor", "hecras"], default="motor")
    ap.add_argument("--instante", type=int, default=-1,
                    help="indice do instante; -1 = envelope maximo da cheia")
    a = ap.parse_args()

    secoes = ler_geometria(a.projeto)
    print(f"{a.projeto}: {len(secoes)} secoes lidas do .g01")

    if a.fonte == "motor":
        d = np.load(f"{a.projeto}_motor.npz", allow_pickle=True)
        ws = d["ws"]
        riv = [str(x) for x in d["river"]]
        rch = [str(x) for x in d["reach"]]
        rs_h = np.array(d["rs"], dtype=float)
    else:
        import h5py
        with h5py.File(f"{a.projeto}.p01.hdf", "r") as f:
            g = f["Results/Unsteady/Output/Output Blocks/Base Output/"
                  "Unsteady Time Series/Cross Sections"]
            ws = g["Water Surface"][:]
            at = f["Geometry/Cross Sections/Attributes"][:]
            riv = [x["River"].decode().strip() for x in at]
            rch = [x["Reach"].decode().strip() for x in at]
            rs_h = np.array([float(x["RS"].decode()) for x in at])
    print(f"  fonte {a.fonte}: {ws.shape[0]} instantes x {ws.shape[1]} secoes")

    # casa cada secao do .g01 com sua coluna na serie
    idx = {}
    for k, s in enumerate(secoes):
        cand = [j for j in range(len(riv))
                if riv[j] == s["rio"] and rch[j] == s["reach"]]
        if cand:
            idx[k] = min(cand, key=lambda j: abs(rs_h[j] - s["rs"]))

    if a.instante < 0:
        # ENVELOPE: a maior lamina que cada secao viu -- e o limite da planicie
        # inundada pelo evento inteiro, nao a foto de um instante.
        ws = ws.max(axis=0)[None, :]
        inst, rotulo = 0, "envelope maximo"
    else:
        inst = min(a.instante, ws.shape[0] - 1)
        rotulo = f"instante {inst}"

    feats = poligonos(secoes, ws, idx, inst)
    if not feats:
        raise SystemExit("nenhum poligono gerado (secoes sem agua?)")
    gdf = gpd.GeoDataFrame(feats, crs=f"EPSG:{UTM}")
    gdf["area_km2"] = (gdf.geometry.area / 1e6).round(3)

    shp = f"{a.projeto}_planicie.shp"
    gdf.to_file(shp)
    gdf.to_crs(4326).to_file(f"{a.projeto}_planicie.geojson", driver="GeoJSON")

    area = float(unary_union(gdf.geometry).area / 1e6)
    print(f"\n[OK] {shp}  ({len(gdf)} poligonos, {rotulo})")
    print(f"[OK] {a.projeto}_planicie.geojson")
    print(f"     area inundada: {area:.1f} km2")
    print(f"     largura maxima: {gdf['larg_max'].max():.0f} m"
          f"   profundidade maxima: {gdf['prof_max'].max():.2f} m")
    for c in gdf["classe"].unique():
        s = gdf[gdf["classe"] == c]
        print(f"     {c:<10} {len(s):3d} poligonos  {s['area_km2'].sum():7.1f} km2")


if __name__ == "__main__":
    main()
