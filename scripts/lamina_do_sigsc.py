# -*- coding: utf-8 -*-
"""Vetoriza a lamina d'agua do SIG-SC (pixels 0,0 do MDT = agua do voo
2010-12) para comparar com FBDS/OSM/modelo no painel.

    python scripts/lamina_do_sigsc.py

Le cada folha decimada 8x (~8 m), poligoniza o zero, descarta area
< 5000 m2 (ruido e valas), simplifica 10 m e grava
doc/painel/lamina_sigsc.geojson em EPSG:4326.

TRES LICOES PAGARAM ESTE ALGORITMO:
  1. zero-de-colar e zero-de-rio se tocam onde o rio cruza a folha --
     teste de forma por componente apaga os rios (rodada 2);
  2. o voto entre folhas sobrepostas mata o colar, mas so o baixo vale
     tem lamina = 0,0: rio acima o aplainamento guarda a COTA REAL da
     agua (340 m em Taio) e o criterio zero nada ve (rodada 3);
  3. dai a FONTE EXTRA (ordem do professor): candidatos = zero-votado
     OU celula PLANA (desnivel 3x3 < 0,15 m, lamina em qualquer cota);
     sementes = candidatos sobre agua da FBDS (massas + rios duplos +
     simples); agua final = propagacao das sementes pelos candidatos.
     A propagacao atravessa as cidades onde a FBDS para (o rio plano e
     continuo); pasto plano sem rio por perto fica de fora.
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

    folhas = sorted(glob.glob(os.path.join(PASTA, 'MDT_*.tif')))
    # grade global de votos a DECIM m
    RES = float(DECIM)
    bs = []
    for p in folhas:
        with rasterio.open(p) as s:
            bs.append(s.bounds)
    x0 = min(b.left for b in bs)
    y0 = min(b.bottom for b in bs)
    x1 = max(b.right for b in bs)
    y1 = max(b.top for b in bs)
    W = int(np.ceil((x1 - x0) / RES))
    H = int(np.ceil((y1 - y0) / RES))
    print(f'grade de votos: {W} x {H} ({W * H / 1e6:.0f} M celulas)')
    zero = np.zeros((H, W), np.uint8)     # alguem disse agua (0,0)
    terra = np.zeros((H, W), np.uint8)    # alguem disse terreno valido
    from scipy.ndimage import maximum_filter, minimum_filter, \
        binary_propagation, binary_closing
    plano = np.zeros((H, W), np.uint8)    # celula plana (lamina alta)
    for k, (p, b) in enumerate(zip(folhas, bs)):
        with rasterio.open(p) as s:
            h, w = s.height // DECIM, s.width // DECIM
            a = s.read(1, out_shape=(h, w)).astype(np.float32)
        j0 = int(round((b.left - x0) / RES))
        i0 = int(round((y1 - b.top) / RES))
        h = min(h, H - i0)
        w = min(w, W - j0)
        a = a[:h, :w]
        z = np.isfinite(a) & (np.abs(a) <= 0.05)
        v = np.isfinite(a) & (a > 0.05)
        desnivel = maximum_filter(a, 3) - minimum_filter(a, 3)
        f = v & (desnivel < 0.15)
        zero[i0:i0 + h, j0:j0 + w] |= z.astype(np.uint8)
        terra[i0:i0 + h, j0:j0 + w] |= v.astype(np.uint8)
        plano[i0:i0 + h, j0:j0 + w] |= f.astype(np.uint8)
        if (k + 1) % 200 == 0:
            print(f'  votos {k + 1}/{len(folhas)}', flush=True)
    cand = ((zero == 1) & (terra == 0)) | (plano == 1)
    del zero, terra, plano
    print(f'candidatos: {int(cand.sum())} celulas')

    # fonte extra = FBDS DENTRO DO DIVISOR da bacia (ANA): fora dele
    # e rio de outra bacia e nao semeia nada
    from rasterio.features import rasterize
    from shapely.geometry import shape as sshape2
    from shapely.prepared import prep
    import pyogrio
    bac = json.load(open('doc/qgis/bacia_itajai_ana.geojson',
                         encoding='utf-8'))
    bacia = sshape2(bac['features'][0]['geometry']).buffer(300.0)
    pbacia = prep(bacia)
    t = from_bounds(x0, y0, x1, y1, W, H)
    fontes, duplos = [], []
    for shp in (glob.glob('doc/fbds/*/*_MASSAS_DAGUA.shp')
                + glob.glob('doc/fbds/*/*_RIOS_DUPLOS.shp')
                + glob.glob('doc/fbds/*/*_RIOS_SIMPLES.shp')):
        eh_duplo = 'RIOS_DUPLOS' in shp or 'MASSAS' in shp
        try:
            g = pyogrio.read_dataframe(shp)
        except Exception:
            continue
        for geom in g.geometry:
            if geom is None or geom.is_empty:
                continue
            if not pbacia.intersects(geom):
                continue
            if not pbacia.contains(geom):
                geom = geom.intersection(bacia)
                if geom.is_empty:
                    continue
            if geom.geom_type in ('LineString', 'MultiLineString'):
                geom = geom.buffer(12.0)
            fontes.append(geom)
            if eh_duplo:
                duplos.append(geom)
    print(f'fonte extra (FBDS na bacia): {len(fontes)} feicoes, '
          f'{len(duplos)} primarias (duplos+massas)')
    fbds = rasterize(((g, 1) for g in fontes), out_shape=(H, W),
                     transform=t, fill=0, dtype='uint8')
    semente = cand & (fbds == 1)
    del fbds
    print(f'sementes: {int(semente.sum())} celulas')
    agua = binary_propagation(semente, mask=cand)
    del semente
    # duplos+massas sao FONTE PRIMARIA, mas com revisao: duplo sobre
    # SOLO tambem existe (rio que mudou, mapeamento velho). O duplo so
    # entra onde se conecta a apoio do MDT: semente = duplo sobre
    # candidato dilatado, propagada DENTRO do proprio duplo. Corredeira
    # ligada ao rio entra; mancha de duplo isolada no seco cai.
    from scipy.ndimage import binary_dilation
    prim = rasterize(((g, 1) for g in duplos), out_shape=(H, W),
                     transform=t, fill=0, dtype='uint8') == 1
    apoio = binary_dilation(cand, np.ones((3, 3), bool), iterations=2)
    del cand
    prim_ok = binary_propagation(prim & apoio, mask=prim)
    print(f'duplos: {int(prim.sum())} celulas mapeadas, '
          f'{int(prim_ok.sum())} com apoio do MDT '
          f'({int(prim.sum() - prim_ok.sum())} sobre solo, fora)')
    del prim, apoio
    agua |= prim_ok
    del prim_ok
    # divisor tambem recorta a SAIDA (mar e bacias vizinhas fora)
    dentro = rasterize([(bacia, 1)], out_shape=(H, W), transform=t,
                       fill=0, dtype='uint8')
    agua &= dentro == 1
    del dentro
    agua = binary_closing(agua, np.ones((3, 3), bool))
    print(f'agua final: {int(agua.sum())} celulas')
    t = from_bounds(x0, y0, x1, y1, W, H)
    feats = []
    for geom, val in rshapes(agua.astype(np.uint8), mask=agua,
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
    os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
    gj = {'type': 'FeatureCollection', 'features': feats}
    with open(SAIDA, 'w') as fh:
        json.dump(gj, fh, separators=(',', ':'))
    print(f'{len(feats)} poligonos de lamina -> {SAIDA} '
          f'({os.path.getsize(SAIDA) // 1024} kB)')


if __name__ == '__main__':
    main()
