# -*- coding: utf-8 -*-
"""Retrato do corredor v2 com os rios da FBDS por cima.

    python scripts/figura_mosaico_rios.py

Sai doc/figuras/mosaico_v2_com_rios.png: terreno decimado (~30 m) com
sombreamento, rios simples (linha) e duplos/massas (contorno) da FBDS
em azul, cidades de referencia.
"""
import glob
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
os.chdir(RAIZ)

MOSAICO = os.path.join('taha_ai_novo', 'Terrain',
                       'taha_ai_corredor_1m_v2.tif')


def main():
    import rasterio
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    import pyogrio
    from pyproj import Transformer

    with rasterio.open(MOSAICO) as s:
        fac = 32
        h, w = s.height // fac, s.width // fac
        a = s.read(1, out_shape=(h, w)).astype(np.float32)
        if s.nodata is not None:
            a[np.isclose(a, s.nodata)] = np.nan
        b = s.bounds
    print(f'mosaico lido: {w}x{h} (decimado {fac}x)')

    ls = LightSource(azdeg=315, altdeg=45)
    rgb = ls.shade(np.where(np.isfinite(a), a, np.nanmedian(a)),
                   cmap=plt.cm.terrain, blend_mode='overlay',
                   vert_exag=2, dx=fac, dy=fac,
                   vmin=np.nanpercentile(a, 2),
                   vmax=np.nanpercentile(a, 98))
    rgb[~np.isfinite(a)] = (1, 1, 1, 1)

    fig, ax = plt.subplots(figsize=(19, 15))
    ax.imshow(rgb, extent=(b.left, b.right, b.bottom, b.top),
              interpolation='nearest')

    # rios da FBDS (ja em EPSG:31982)
    n = 0
    for shp in (glob.glob('doc/fbds/*/*_RIOS_SIMPLES.shp')
                + glob.glob('doc/fbds/*/*_RIOS_DUPLOS.shp')
                + glob.glob('doc/fbds/*/*_MASSAS_DAGUA.shp')):
        try:
            g = pyogrio.read_dataframe(shp)
        except Exception:
            continue
        fino = 'SIMPLES' in shp
        for geom in g.geometry:
            if geom is None:
                continue
            gs = list(geom.geoms) if hasattr(geom, 'geoms') else [geom]
            for gg in gs:
                xy = (gg.exterior.xy if gg.geom_type == 'Polygon'
                      else gg.xy)
                ax.plot(*xy, '-', color='#1050a0',
                        lw=0.25 if fino else 0.6,
                        alpha=0.5 if fino else 0.8)
                n += 1
    print(f'{n} tracos de rio')

    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    cidades = {'Taió': (-50.1156, -27.1225),
               'Rio do Sul': (-49.6430, -27.2144),
               'Ituporanga': (-49.5983, -27.4142),
               'Ibirama': (-49.5217, -27.0561),
               'Blumenau': (-49.0661, -26.9186),
               'Brusque': (-48.9142, -27.0989),
               'Itajaí': (-48.6650, -26.9078),
               'Timbó': (-49.2719, -26.8236)}
    for nome, (lon, lat) in cidades.items():
        x, y = tr.transform(lon, lat)
        ax.plot(x, y, '*', color='k', ms=13, mec='w', mew=0.8)
        ax.annotate(nome, (x, y), textcoords='offset points',
                    xytext=(7, 7), fontsize=11, fontweight='bold',
                    color='k',
                    path_effects=None)
    ax.set_xlim(b.left, b.right)
    ax.set_ylim(b.bottom, b.top)
    ax.set_aspect('equal')
    ax.set_title('Corredor 1 m v2 (QA/QC: zero=vazio + Copérnico '
                 'esponjado) com a hidrografia FBDS', fontsize=15)
    ax.set_xlabel('E (m, SIRGAS2000 UTM 22S)')
    ax.set_ylabel('N (m)')
    fig.tight_layout()
    fig.savefig('doc/figuras/mosaico_v2_com_rios.png', dpi=110)
    print('figura: doc/figuras/mosaico_v2_com_rios.png')


if __name__ == '__main__':
    main()
