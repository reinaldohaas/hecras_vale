# -*- coding: utf-8 -*-
"""Exporta o relevo do projeto em camadas prontas para o QGIS.

    python scripts/exportar_relevo_qgis.py

  doc/qgis/relevo_copernicus_30m_utm.tif   bacia inteira, EPSG:31982
  doc/qgis/hand_45m.tif                    altura acima do rio mais
                                           proximo (a pintura das
                                           planicies), EPSG:31982
  pirâmides .ovr no corredor 1 m do SIG-SC (arquivo original intocado)

O corredor 1 m ja e GeoTIFF EPSG:31982 com nodata: abre direto no QGIS
(taha_ai_novo/Terrain/taha_ai_corredor_1m_completo.tif).
"""
import os
import subprocess
import sys

import numpy as np
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.windows import from_bounds

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RAIZ)
BB = (565000, 6933000, 740000, 7070000)   # oeste ate cobrir o corredor
                                          # inteiro (v2 comeca em 572489)

# 1. Copernicus em UTM
src = rasterio.open('dem_bacia_itajai.tif')
vrt = WarpedVRT(src, crs='EPSG:31982',
                resampling=rasterio.enums.Resampling.bilinear)
w = from_bounds(*BB, transform=vrt.transform)
NY, NX = int((BB[3]-BB[1])/30), int((BB[2]-BB[0])/30)
Z = vrt.read(1, window=w, out_shape=(NY, NX)).astype('float32')
Z = np.where(Z < -100, -9999, Z)
from rasterio.transform import from_bounds as tfb
perfil = {'driver': 'GTiff', 'height': NY, 'width': NX, 'count': 1,
          'dtype': 'float32', 'crs': 'EPSG:31982', 'nodata': -9999.0,
          'transform': tfb(*BB, NX, NY),
          'compress': 'deflate', 'tiled': True}
with rasterio.open('doc/qgis/relevo_copernicus_30m_utm.tif', 'w',
                   **perfil) as dst:
    dst.write(Z, 1)
print('copernicus UTM ok')

# 2. HAND 45 m (mesma receita do editor de linhas)
sys.path.insert(0, 'scripts')
from qc_secoes import ler_secoes
from qc_geometria import ler_eixos
from shapely.geometry import LineString
from scipy.spatial import cKDTree
G = 'taha_ai.g01.estado_mes_relevo1m'
S = ler_secoes(G)
E = ler_eixos(G)
eixos = {}
for (rio, reach), lsg in sorted(E.items()):
    eixos.setdefault(rio, []).extend(list(lsg.coords))
pontos, cotas = [], []
for rio, coords in eixos.items():
    ls = LineString(coords)
    secs = sorted([(d['rs'], float(np.asarray(d['z'], float).min()))
                   for d in S if d['rio'] == rio and d['tipo'] == '1'])
    if not secs:
        continue
    rss = np.array([s[0] for s in secs])
    tals = np.array([s[1] for s in secs])
    for s in np.arange(0, ls.length, 100):
        p = ls.interpolate(s)
        pontos.append((p.x, p.y))
        cotas.append(float(np.interp(ls.length - s, rss, tals)))
arv = cKDTree(np.array(pontos))
cotas = np.array(cotas)
RES = 45.0
NY2, NX2 = int((BB[3]-BB[1])/RES), int((BB[2]-BB[0])/RES)
w2 = from_bounds(*BB, transform=vrt.transform)
Z2 = vrt.read(1, window=w2, out_shape=(NY2, NX2)).astype('float32')
gx, gy = np.meshgrid(np.linspace(BB[0], BB[2], NX2),
                     np.linspace(BB[3], BB[1], NY2))
dist, idx = arv.query(np.column_stack([gx.ravel(), gy.ravel()]),
                      workers=-1)
hand = (Z2.ravel() - cotas[idx]).astype('float32')
hand[dist > 12000] = -9999
hand = hand.reshape(Z2.shape)
perfil.update({'height': NY2, 'width': NX2,
               'transform': tfb(*BB, NX2, NY2)})
with rasterio.open('doc/qgis/hand_45m.tif', 'w', **perfil) as dst:
    dst.write(hand, 1)
print('hand ok')

# 3. piramides externas no corredor 1 m
alvo = 'taha_ai_novo/Terrain/taha_ai_corredor_1m_completo.tif'
gdaladdo = os.path.join(os.path.dirname(sys.executable),
                        'Library', 'bin', 'gdaladdo.exe')
if os.path.exists(gdaladdo) and not os.path.exists(alvo + '.ovr'):
    subprocess.run([gdaladdo, '-ro', '-r', 'average', alvo,
                    '4', '16', '64', '256'], check=True)
    print('piramides .ovr do corredor 1 m ok')
else:
    print('piramides: gdaladdo ausente ou .ovr ja existe')
print('doc/qgis/ pronto')
