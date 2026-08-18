# -*- coding: utf-8 -*-
"""
Ponte com a geometria do HEC-RAS: le um .g01 e devolve eixos e cutlines.

Sem isto o programa so serviria para quem ja tem as secoes em shapefile. O
caso de uso real -- auditar a geometria de um modelo existente -- comeca no
.g01, que e onde as secoes de fato estao.

Le-se:
  River Reach=<rio>,<trecho>   seguido de Reach XY= n  e as coordenadas
  Type RM Length L Ch R = 1 ,<RS>,...
  XS GIS Cut Line=<n>          seguido das coordenadas da linha de corte

As coordenadas vem em colunas de 16 caracteres, 4 por linha (2 pares).
"""
import geopandas as gpd
from shapely.geometry import LineString


def _coords(linhas, i, n_pares):
    """Le n_pares coordenadas a partir da linha i, em colunas de 16."""
    vals = []
    while i < len(linhas) and len(vals) < 2 * n_pares:
        s = linhas[i].rstrip()
        if not s or "=" in s:
            break
        vals += [float(s[c:c + 16]) for c in range(0, len(s), 16)
                 if s[c:c + 16].strip()]
        i += 1
    pares = list(zip(vals[0::2], vals[1::2]))
    return pares, i


def ler_g01(caminho, com_perfil=False):
    """Devolve (eixos, secoes) como listas de dicionarios com geometria."""
    txt = open(caminho, encoding="utf-8", errors="ignore").read().splitlines()
    eixos, secoes = [], []
    rio = reach = None
    rs = None
    i = 0
    while i < len(txt):
        l = txt[i]
        if l.startswith("River Reach="):
            p = l.split("=", 1)[1].split(",")
            rio, reach = p[0].strip(), p[1].strip()
        elif l.startswith("Reach XY="):
            n = int(l.split("=")[1])
            pares, i = _coords(txt, i + 1, n)
            if len(pares) >= 2:
                eixos.append({"river": rio, "reach": reach,
                              "geometry": LineString(pares)})
            continue
        elif l.startswith("Type RM"):
            try:
                rs = float(l.split(",")[1])
            except (IndexError, ValueError):
                rs = None
        elif l.startswith("XS GIS Cut Line"):
            n = int(l.split("=")[1])
            pares, i = _coords(txt, i + 1, n)
            if len(pares) >= 2 and rs is not None:
                secoes.append({"river": rio, "reach": reach, "rs": rs,
                               "geometry": LineString(pares)})
            continue
        elif com_perfil and l.startswith("#Sta/Elev="):
            # station/elevation em colunas de 8, 10 por linha -- e o que o
            # modelo.py compara contra o terreno
            n = int(l.split("=")[1])
            v = []
            i += 1
            while i < len(txt) and len(v) < 2 * n:
                t = txt[i].rstrip()
                if not t or "=" in t:
                    break
                v += [float(t[c:c + 8]) for c in range(0, len(t), 8)
                      if t[c:c + 8].strip()]
                i += 1
            if secoes and len(v) >= 6:
                secoes[-1]["sta"] = v[0::2]
                secoes[-1]["z"] = v[1::2]
            continue
        i += 1
    return eixos, secoes


def exportar(caminho_g01, crs, prefixo=None, rios=None):
    """Grava <prefixo>_eixo.geojson e <prefixo>_secoes.geojson."""
    import os
    eixos, secoes = ler_g01(caminho_g01)
    if rios:
        alvo = set(rios)
        eixos = [e for e in eixos if e["river"] in alvo]
        secoes = [s for s in secoes if s["river"] in alvo]
    prefixo = prefixo or os.path.splitext(caminho_g01)[0]
    pe = f"{prefixo}_eixo.geojson"
    ps = f"{prefixo}_secoes.geojson"
    gpd.GeoDataFrame(eixos, crs=crs).to_file(pe, driver="GeoJSON")
    gpd.GeoDataFrame(secoes, crs=crs).to_file(ps, driver="GeoJSON")
    return pe, ps


if __name__ == "__main__":
    import sys
    g01 = sys.argv[1]
    crs = sys.argv[2] if len(sys.argv) > 2 else "EPSG:31982"
    rios = sys.argv[3].split(",") if len(sys.argv) > 3 else None
    print(*exportar(g01, crs, rios=rios), sep="\n")
