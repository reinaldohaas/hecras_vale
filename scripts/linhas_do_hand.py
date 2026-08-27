# -*- coding: utf-8 -*-
"""Gera as linhas INICIAIS de edge e bank a partir do terreno, por HAND.

    python scripts/linhas_do_hand.py taha_ai.g01.estado_mes_relevo1m \
        --hand-edge 15 --saida doc/linhas_hand

SO GERA LINHAS (GeoJSON por grupo + JSON do editor). Nada de g01.

COMO (ideia do Reinaldo, 27/08)

  EDGES  a pintura HAND vira o limite: altura acima do talvegue do rio
         mais proximo (Copernicus 45 m + talvegues do modelo). Cada
         pixel pertence ao rio mais proximo, entao onde duas planicies
         se tocam (Acu-Mirim) as edges dos dois rios se ENCONTRAM na
         linha media -- o divisor sai de graca. O contorno da mancha
         HAND <= --hand-edge de cada rio, separado em margem N e S pelo
         lado do eixo e ordenado pela quilometragem, e a edge inicial.

  BANKS  comecam no rio DESENHADO no SIG-SC/FBDS: o contorno dos
         poligonos RIOS_DUPLOS perto do eixo, tambem separado N/S e
         ordenado. E a margem d'agua real -- o editor ajusta dali.

  As linhas vao para o editor interativo (artifact) para ajuste a mouse
  e aprovacao; so depois viram secoes.
"""
import glob
import json
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_secoes import ler_secoes                       # noqa: E402
from qc_geometria import ler_eixos                     # noqa: E402

BB = (580000, 6933000, 740000, 7070000)
RES = 45.0
CORES = {'Itajai_Acu': '#d62828', 'Itajai_Mirim': '#f77f00',
         'Itajai_Norte': '#5f0f40', 'Itajai_Oeste': '#9d4edd',
         'Itajai_Sul': '#0077b6', 'Rio_Benedito': '#2d6a4f',
         'Rio_dos_Cedros': '#74a57f', 'Rio_do_Testo': '#b5838d',
         'Rio_Iraputa': '#6c757d', 'Rio_Taio': '#c9184a',
         'Rio_Trombudo': '#7f5539', 'Rio_das_Pombas': '#3a5a40'}
ROTULOS = {'Itajai_Acu': 'Açu', 'Itajai_Mirim': 'Mirim',
           'Itajai_Norte': 'Norte', 'Itajai_Oeste': 'Oeste',
           'Itajai_Sul': 'Sul', 'Rio_Benedito': 'Benedito',
           'Rio_dos_Cedros': 'Cedros', 'Rio_do_Testo': 'Testo',
           'Rio_Iraputa': 'Iraputá', 'Rio_Taio': 'Taió',
           'Rio_Trombudo': 'Trombudo', 'Rio_das_Pombas': 'Pombas'}


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def eixos_e_talvegues(g01):
    from shapely.geometry import LineString
    S = ler_secoes(g01)
    E = ler_eixos(g01)
    eixos = {}
    for (rio, reach), lsg in sorted(E.items()):
        eixos.setdefault(rio, []).extend(list(lsg.coords))
    saida = {}
    for rio, coords in eixos.items():
        ls = LineString(coords)
        secs = sorted([(d['rs'], float(np.asarray(d['z'], float).min()))
                       for d in S if d['rio'] == rio and d['tipo'] == '1'])
        if not secs:
            continue
        rss = np.array([s[0] for s in secs])
        tals = np.array([s[1] for s in secs])
        pontos, cotas = [], []
        for s in np.arange(0, ls.length, 100):
            p = ls.interpolate(s)
            pontos.append((p.x, p.y))
            cotas.append(float(np.interp(ls.length - s, rss, tals)))
        saida[rio] = (ls, np.array(pontos), np.array(cotas))
    return saida


def lado_e_km(ls, P):
    """para cada ponto: lado (+1 N/esq, -1 S/dir) e quilometragem."""
    from shapely.geometry import Point
    lados, kms = [], []
    for x, y in P:
        s = ls.project(Point(x, y))
        a = ls.interpolate(max(s - 30, 0))
        b = ls.interpolate(min(s + 30, ls.length))
        cz = (b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x)
        lados.append(1 if cz > 0 else -1)
        kms.append(ls.length - s)
    return np.array(lados), np.array(kms)


def em_linhas(P, lados, kms, lado, simpl):
    """pontos de um lado -> polilinha ordenada pela km e simplificada."""
    from shapely.geometry import LineString
    m = lados == lado
    if m.sum() < 4:
        return None
    o = np.argsort(-kms[m])
    pts = P[m][o]
    # corta saltos > 3 km (bracos soltos do contorno)
    partes, atual = [], [pts[0]]
    for p in pts[1:]:
        if np.hypot(p[0] - atual[-1][0], p[1] - atual[-1][1]) > 3000:
            partes.append(atual)
            atual = [p]
        else:
            atual.append(p)
    partes.append(atual)
    pts = max(partes, key=len)
    if len(pts) < 4:
        return None
    ls = LineString(pts).simplify(simpl)
    # teto de vertices por linha: editor a mouse nao quer milhares
    while len(ls.coords) > 700:
        simpl *= 1.6
        ls = ls.simplify(simpl)
    return [[round(x, 1), round(y, 1)] for x, y in ls.coords]


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    g01 = argv[0]
    hand_edge = _arg(argv, '--hand-edge', 15.0, float)
    pasta = _arg(argv, '--saida', 'doc/linhas_hand')
    os.makedirs(pasta, exist_ok=True)
    import rasterio
    from rasterio.vrt import WarpedVRT
    from rasterio.windows import from_bounds
    from scipy.spatial import cKDTree
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import LineString

    rios = eixos_e_talvegues(g01)
    todos_p, todos_c, todos_r = [], [], []
    for k, (rio, (ls, P, C)) in enumerate(sorted(rios.items())):
        todos_p.append(P)
        todos_c.append(C)
        todos_r += [rio] * len(P)
    arv = cKDTree(np.vstack(todos_p))
    cotas = np.concatenate(todos_c)
    rotulo = np.array(todos_r)

    src = rasterio.open('dem_bacia_itajai.tif')
    vrt = WarpedVRT(src, crs='EPSG:31982',
                    resampling=rasterio.enums.Resampling.bilinear)
    w = from_bounds(*BB, transform=vrt.transform)
    NY, NX = int((BB[3] - BB[1]) / RES), int((BB[2] - BB[0]) / RES)
    Z = vrt.read(1, window=w, out_shape=(NY, NX)).astype('float32')
    Z = np.where(Z < -10, np.nan, Z)
    gx = np.linspace(BB[0], BB[2], NX)
    gy = np.linspace(BB[3], BB[1], NY)
    GX, GY = np.meshgrid(gx, gy)
    dist, idx = arv.query(np.column_stack([GX.ravel(), GY.ravel()]),
                          workers=-1)
    hand = (Z.ravel() - cotas[idx]).reshape(Z.shape)
    perto = (dist < 12000).reshape(Z.shape)
    dono = rotulo[idx].reshape(Z.shape)
    print(f'HAND {NX}x{NY} pronto')

    # rio desenhado (FBDS) para os banks
    partes = []
    for shp in glob.glob('doc/fbds/*/SC_*_RIOS_DUPLOS.shp'):
        try:
            g = gpd.read_file(shp)
            if g.crs is None:
                g = g.set_crs(31982)
            partes.append(g.to_crs(31982)[['geometry']])
        except Exception:
            pass
    fbds = gpd.GeoDataFrame(pd.concat(partes, ignore_index=True),
                            crs=31982)

    linhas = []
    geo = {'type': 'FeatureCollection',
           'crs': {'type': 'name', 'properties': {'name': 'EPSG:31982'}},
           'features': []}
    for rio, (ls, P, C) in sorted(rios.items()):
        # ---- EDGES pelo HAND (contorno via matplotlib, sem skimage)
        mancha = (hand <= hand_edge) & perto & (dono == rio)
        figc, axc = plt.subplots()
        cs = axc.contour(gx, gy, mancha.astype(float), levels=[0.5])
        segs = [s for s in cs.allsegs[0] if len(s) > 3]
        plt.close(figc)
        if segs:
            PC = max(segs, key=len)
            lados, kms = lado_e_km(ls, PC)
            for lado, tag in [(1, 'N'), (-1, 'S')]:
                lin = em_linhas(PC, lados, kms, lado, 80.0)
                if lin:
                    nome = f'{ROTULOS.get(rio, rio)} — edge {tag} (HAND {hand_edge:.0f} m)'
                    linhas.append({'nome': nome, 'cor': CORES.get(rio, '#333'),
                                   'grupo': ROTULOS.get(rio, rio),
                                   'pontos': lin})
                    geo['features'].append({'type': 'Feature',
                                            'properties': {'nome': nome},
                                            'geometry': {'type': 'LineString',
                                                         'coordinates': lin}})
        # ---- BANKS pelo rio desenhado
        try:
            faixa = ls.buffer(600)
            sel = fbds[fbds.intersects(faixa)]
            if len(sel):
                uniao = sel.union_all() if hasattr(sel, 'union_all') \
                    else sel.unary_union
                uniao = uniao.intersection(faixa)
                borda = []
                geoms = getattr(uniao, 'geoms', [uniao])
                for gpoly in geoms:
                    if gpoly.geom_type != 'Polygon':
                        continue
                    borda += list(gpoly.exterior.coords)
                if len(borda) > 10:
                    PB = np.array(borda)
                    lados, kms = lado_e_km(ls, PB)
                    for lado, tag in [(1, 'N'), (-1, 'S')]:
                        lin = em_linhas(PB, lados, kms, lado, 25.0)
                        if lin:
                            nome = f'{ROTULOS.get(rio, rio)} — bank {tag} (rio SIG-SC)'
                            linhas.append({'nome': nome,
                                           'cor': CORES.get(rio, '#333'),
                                           'grupo': ROTULOS.get(rio, rio),
                                           'banco': True, 'pontos': lin})
                            geo['features'].append(
                                {'type': 'Feature',
                                 'properties': {'nome': nome},
                                 'geometry': {'type': 'LineString',
                                              'coordinates': lin}})
        except Exception as e:
            print(f'   {rio}: banks falharam ({e})')
        print(f'   {rio}: ok')

    dados = {'bbox': list(BB), 'linhas': linhas}
    json.dump(dados, open(os.path.join(pasta, 'linhas_editor.json'), 'w'))
    json.dump(geo, open(os.path.join(pasta, 'linhas_hand.geojson'), 'w'))
    n_v = sum(len(l['pontos']) for l in linhas)
    print(f'\n{len(linhas)} linhas, {n_v} vertices -> {pasta}/')


if __name__ == '__main__':
    main(sys.argv[1:])
