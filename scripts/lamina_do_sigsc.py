# -*- coding: utf-8 -*-
"""Vetoriza a lamina d'agua do SIG-SC (pixels 0,0 do MDT = agua do voo
2010-12) para comparar com FBDS/OSM/modelo no painel.

    python scripts/lamina_do_sigsc.py

Le cada folha decimada 8x (~8 m), poligoniza o zero, descarta area
< 5000 m2 (ruido e valas), simplifica 10 m e grava
doc/painel/lamina_sigsc.geojson em EPSG:4326. O colar de borda das
folhas tambem e zero: um anel de 30 px e descartado da mascara.
"""
import csv
import glob
import json
import os

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RAIZ)
PASTA = r'C:\Users\haas\Downloads\sigsc'
SAIDA = os.path.join('doc', 'painel', 'lamina_sigsc.geojson')
REL = os.path.join('doc', 'qaqc_sigsc', 'relatorio.csv')
DECIM = 8
MIN_AREA = 5000.0
COLAR_MIN = 2         # px decimados sempre descartados na borda


def colares():
    """Colar medido por folha (m) no relatorio do qaqc, se existir."""
    if not os.path.exists(REL):
        return {}
    out = {}
    for r in csv.DictReader(open(REL)):
        out[r['folha']] = float(r.get('colar_px', 0) or 0)
    return out


def main():
    import rasterio
    from rasterio.features import shapes as rshapes
    from rasterio.transform import from_bounds
    from shapely.geometry import shape as sshape
    from shapely.ops import transform as stransform
    from pyproj import Transformer
    tr = Transformer.from_crs(31982, 4326, always_xy=True)

    def para4326(x, y):
        return tr.transform(x, y)

    col = colares()
    feats = []
    folhas = sorted(glob.glob(os.path.join(PASTA, 'MDT_*.tif')))
    for k, p in enumerate(folhas):
        with rasterio.open(p) as s:
            h, w = s.height // DECIM, s.width // DECIM
            a = s.read(1, out_shape=(h, w)).astype(np.float32)
            b = s.bounds
        agua = (np.isfinite(a) & (np.abs(a) <= 0.05)).astype(np.uint8)
        # descarta so o colar MEDIDO desta folha (relatorio qaqc),
        # nunca menos que COLAR_MIN px
        cpx = int(col.get(os.path.basename(p), 0)) // DECIM + COLAR_MIN
        agua[:cpx, :] = 0
        agua[-cpx:, :] = 0
        agua[:, :cpx] = 0
        agua[:, -cpx:] = 0
        if not agua.any():
            continue
        t = from_bounds(b.left, b.bottom, b.right, b.top, w, h)
        for geom, val in rshapes(agua, mask=agua.astype(bool),
                                 transform=t):
            g = sshape(geom)
            if g.area < MIN_AREA:
                continue
            g = g.simplify(10.0)
            if g.is_empty:
                continue
            feats.append({'type': 'Feature', 'properties': {},
                          'geometry': stransform(
                              para4326, g).__geo_interface__})
        if (k + 1) % 100 == 0:
            print(f'  {k + 1}/{len(folhas)}: {len(feats)} poligonos',
                  flush=True)
    os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
    gj = {'type': 'FeatureCollection', 'features': feats}
    with open(SAIDA, 'w') as fh:
        json.dump(gj, fh, separators=(',', ':'))
    print(f'{len(feats)} poligonos de lamina -> {SAIDA} '
          f'({os.path.getsize(SAIDA) // 1024} kB)')


if __name__ == '__main__':
    main()
