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

    # dominio: corredor da rede completa
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    rede = json.load(open('vale_itajai_full_network.geojson',
                          encoding='utf-8'))
    linhas = []
    for f in rede['features']:
        g = f['geometry']
        cs = [g['coordinates']] if g['type'] == 'LineString' \
            else g['coordinates']
        for c in cs:
            if len(c) > 1:
                linhas.append(LineString([tr.transform(x, y)
                                          for x, y in c]))
    corredor = unary_union([l.buffer(2000, resolution=4)
                            for l in linhas])
    print('corredor da bacia pronto')

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

    # geojson das faltantes
    os.makedirs('doc/qgis', exist_ok=True)
    geo = {'type': 'FeatureCollection',
           'crs': {'type': 'name',
                   'properties': {'name': 'EPSG:31982'}},
           'features': [{'type': 'Feature',
                         'properties': {'n': k + 1},
                         'geometry': json.loads(json.dumps(
                             c.__geo_interface__))}
                        for k, c in enumerate(faltantes)]}
    json.dump(geo, open('doc/qgis/folhas_faltantes.geojson', 'w'))

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
    for l in linhas:
        x, y = l.xy
        ax.plot(x, y, '-', color='navy', lw=0.3, alpha=0.5)
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
