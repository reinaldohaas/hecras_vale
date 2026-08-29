# -*- coding: utf-8 -*-
"""Mapa de cobertura do SIG-SC: que folhas existem e quais FALTAM trazer.

    python scripts/mapa_cobertura_sigsc.py

Dominio necessario = corredor de 2 km em torno de TODOS os rios da bacia
(vale_itajai_full_network.geojson -- inclui Alfredo Wagner no Sul e as
cabeceiras do Oeste rumo a Santa Cecilia). Grade = molde das folhas
existentes em C:/Users/haas/Downloads/sigsc. Sai:

  doc/figuras/cobertura_sigsc.png       o mapa (verde=tem, vermelho=falta)
  doc/qgis/folhas_faltantes.geojson     molduras faltantes p/ QGIS
"""
import glob
import json
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
sys.path.insert(0, DIR)
os.chdir(RAIZ)
PASTA = r'C:\Users\haas\Downloads\sigsc'


def codigo_carta(lon, lat):
    """Nomenclatura sistematica brasileira ate 1:10.000 (ex.
    SG-22-Z-B-V-1-SO-F) a partir de um ponto lon/lat (graus)."""
    faixa = 'ABCDEFGHIJ'[int(-lat // 4)]        # hemisferio sul
    fuso = int((lon + 180) // 6) + 1
    x0 = (fuso - 1) * 6 - 180
    y0 = -4 * (int(-lat // 4) + 1)
    partes = [f'S{faixa}-{fuso}']
    # (ncol, nlin, rotulos em ordem de leitura NW->SE)
    niveis = [(2, 2, ['V', 'X', 'Y', 'Z']),
              (2, 2, ['A', 'B', 'C', 'D']),
              (3, 2, ['I', 'II', 'III', 'IV', 'V', 'VI']),
              (2, 2, ['1', '2', '3', '4']),
              (2, 2, ['NO', 'NE', 'SO', 'SE']),
              (2, 3, ['A', 'B', 'C', 'D', 'E', 'F'])]
    dx, dy = 6.0, 4.0
    for ncol, nlin, rot in niveis:
        dx, dy = dx / ncol, dy / nlin
        i = min(int((lon - x0) / dx), ncol - 1)
        j = min(int((y0 + nlin * dy - lat) / dy), nlin - 1)
        partes.append(rot[j * ncol + i])
        x0 += i * dx
        y0 += (nlin - 1 - j) * dy
    return '-'.join(partes)


def main():
    import rasterio
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from shapely.geometry import LineString, box
    from shapely.ops import unary_union
    from pyproj import Transformer

    # folhas existentes
    folhas = []
    for p in sorted(glob.glob(os.path.join(PASTA, '*.tif'))):
        with rasterio.open(p) as s:
            b = s.bounds
        folhas.append((b.left, b.bottom, b.right, b.top,
                       os.path.basename(p)))
    if not folhas:
        raise SystemExit(f'nenhuma folha em {PASTA}')
    W = float(np.median([f[2] - f[0] for f in folhas]))
    H = float(np.median([f[3] - f[1] for f in folhas]))
    x0 = min(f[0] for f in folhas)
    y0 = min(f[1] for f in folhas)
    print(f'{len(folhas)} folhas existentes; molde {W:.0f} x {H:.0f} m')

    # dominio: TODOS os fios da FBDS dos municipios da bacia (a rede
    # geojson e recortada -- nem chega a Alfredo Wagner/Santa Cecilia)
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    import geopandas as gpd
    from scipy.spatial import cKDTree
    pontos_rio = []
    linhas = []
    for shp in (glob.glob('doc/fbds/*/*_RIOS_SIMPLES.shp')
                + glob.glob('doc/fbds/*/*_RIOS_DUPLOS.shp')):
        try:
            g = gpd.read_file(shp)
            for geom in g.geometry:
                if geom is None:
                    continue
                gs = [geom] if geom.geom_type in ('LineString',
                                                  'Polygon') \
                    else list(geom.geoms)
                for gg in gs:
                    ls = LineString(gg.exterior.coords) \
                        if gg.geom_type == 'Polygon' else gg
                    linhas.append(ls)
                    for s in np.arange(0, ls.length, 300):
                        p = ls.interpolate(s)
                        pontos_rio.append((p.x, p.y))
        except Exception:
            pass
    print(f'fios FBDS (desenho): {len(linhas)}')
    # dominio OFICIAL: Micro Regiao Hidrografica Itajai (ANA/IBGE),
    # baixada de snirh.gov.br -> doc/qgis/bacia_itajai_ana.geojson
    from shapely.geometry import shape
    bac = json.load(open('doc/qgis/bacia_itajai_ana.geojson',
                         encoding='utf-8'))
    corredor = shape(bac['features'][0]['geometry'])
    print(f'dominio: bacia oficial ANA ({corredor.area/1e6:.0f} km2)')

    # grade ancorada no molde das folhas existentes
    bx = corredor.bounds
    i0 = int(np.floor((bx[0] - x0) / W))
    i1 = int(np.ceil((bx[2] - x0) / W))
    j0 = int(np.floor((bx[1] - y0) / H))
    j1 = int(np.ceil((bx[3] - y0) / H))
    tem = {(round((f[0] - x0) / W), round((f[1] - y0) / H))
           for f in folhas}
    presentes, faltantes = [], []
    for i in range(i0, i1 + 1):
        for j in range(j0, j1 + 1):
            cel = box(x0 + i * W, y0 + j * H,
                      x0 + (i + 1) * W, y0 + (j + 1) * H)
            if not cel.intersects(corredor):
                continue
            (presentes if (i, j) in tem else faltantes).append(cel)
    print(f'celulas no corredor: {len(presentes)} presentes, '
          f'{len(faltantes)} FALTANTES')

    # valida codigo_carta contra os nomes das folhas existentes
    inv = Transformer.from_crs(31982, 4326, always_xy=True)
    ok = err = 0
    for f in folhas:
        lon, lat = inv.transform((f[0] + f[2]) / 2, (f[1] + f[3]) / 2)
        esperado = f[4].replace('MDT_', '').replace('.tif', '')
        if codigo_carta(lon, lat) == esperado:
            ok += 1
        else:
            err += 1
            if err <= 3:
                print(f'  DIVERGE: {esperado} != '
                      f'{codigo_carta(lon, lat)}')
    print(f'validacao da nomenclatura: {ok} ok, {err} divergentes')

    # geojson das faltantes, com o codigo da carta a baixar
    os.makedirs('doc/qgis', exist_ok=True)
    cartas = []
    for c in faltantes:
        lon, lat = inv.transform(c.centroid.x, c.centroid.y)
        cartas.append(codigo_carta(lon, lat))
    geo = {'type': 'FeatureCollection',
           'crs': {'type': 'name',
                   'properties': {'name': 'EPSG:31982'}},
           'features': [{'type': 'Feature',
                         'properties': {'n': k + 1, 'carta': cartas[k]},
                         'geometry': json.loads(json.dumps(
                             c.__geo_interface__))}
                        for k, c in enumerate(faltantes)]}
    json.dump(geo, open('doc/qgis/folhas_faltantes.geojson', 'w'))
    with open('doc/qgis/folhas_faltantes.txt', 'w') as fh:
        fh.write('\n'.join(sorted(cartas)) + '\n')
    print('lista p/ portal: doc/qgis/folhas_faltantes.txt')

    # cidades de referencia
    cidades = {'Alfredo Wagner': (-49.3344, -27.7000),
               'Santa Cecília': (-50.4269, -26.9608),
               'Taió': (-50.1156, -27.1225),
               'Ituporanga': (-49.5983, -27.4142),
               'Rio do Sul': (-49.6430, -27.2144),
               'Blumenau': (-49.0661, -26.9186),
               'Itajaí': (-48.6650, -26.9078),
               'Brusque': (-48.9142, -27.0989)}

    fig, ax = plt.subplots(figsize=(16, 13))
    for geoms, cor, alfa in [(presentes, '#2d6a4f', 0.35),
                             (faltantes, '#d62828', 0.55)]:
        for c in geoms:
            xs, ys = c.exterior.xy
            ax.fill(xs, ys, color=cor, alpha=alfa, lw=0)
            ax.plot(xs, ys, color=cor, lw=0.4)
    for l in linhas[::4]:
        x, y = l.xy
        ax.plot(x, y, '-', color='navy', lw=0.15, alpha=0.35)
    bx_, by_ = corredor.exterior.xy
    ax.plot(bx_, by_, 'k-', lw=1.6, alpha=0.8,
            label='bacia do Itajaí (ANA/IBGE)')
    for c, carta in zip(faltantes, cartas):
        ax.annotate(carta.replace('SG-22-', ''),
                    (c.centroid.x, c.centroid.y), ha='center',
                    va='center', fontsize=5.5, color='#7a0000')
    ax.legend(loc='lower left', fontsize=10)
    for nome, (lon, lat) in cidades.items():
        x, y = tr.transform(lon, lat)
        ax.plot(x, y, 'k*', ms=12)
        ax.annotate(nome, (x, y), textcoords='offset points',
                    xytext=(6, 6), fontsize=10, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(alpha=0.2)
    ax.set_title(f'Cobertura do SIG-SC 1 m no corredor da bacia: '
                 f'{len(presentes)} folhas presentes (verde), '
                 f'{len(faltantes)} FALTANTES (vermelho)', fontsize=13)
    ax.set_xlabel('E (m, SIRGAS2000 UTM 22S)')
    ax.set_ylabel('N (m)')
    fig.tight_layout()
    fig.savefig('doc/figuras/cobertura_sigsc.png', dpi=120)
    print('figura: doc/figuras/cobertura_sigsc.png')
    print('molduras: doc/qgis/folhas_faltantes.geojson')


if __name__ == '__main__':
    main()
