# -*- coding: utf-8 -*-
"""Cursos d'agua LINEARES da ANA, nomeados por ancora, p/ comparacao.

    python scripts/gerar_ana_rios.py

O shapefile da ANA nao traz nome, mas COCURSODAG e constante ao longo
de cada curso. Ancoramos cada rio num ponto conhecido (cidade/foz),
pegamos o codigo do segmento mais proximo e reconstruimos o curso
inteiro. Sai doc/painel/ana_rios.geojson (4326, nome em properties).

E o desempate de identidade: o que o legado chama de alto Itajai do
Oeste pode ser o Mirim Doce -- aqui cada linha tem dono oficial.
"""
import json
import os

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RAIZ)
SHP = r'C:\Users\haas\Downloads\ANA_curso_dagua\curso_dagua.shp'
SAIDA = os.path.join('doc', 'painel', 'ana_rios.geojson')

# ancoras (lon, lat) SO onde a identificacao conferiu; os demais
# cursos ficam com o codigo Otto (melhor sem nome que nome errado)
ANCORAS = {
    'Rio Mirim Doce':    (-50.3600, -27.1980),   # perto de Mirim Doce
    'Rio Trombudo':      (-49.7930, -27.3010),   # Trombudo Central
    'Itajaí do Sul':     (-49.6050, -27.4090),   # Ituporanga
    'Rio Taió':          (-50.1130, -27.1650),   # sul de Taio
}


def main():
    import pyogrio
    from pyproj import Transformer
    from shapely.ops import transform as stransform
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    inv = Transformer.from_crs(31982, 4326, always_xy=True)

    bac = json.load(open('doc/qgis/bacia_itajai_ana.geojson',
                         encoding='utf-8'))
    from shapely.geometry import shape
    bacia = shape(bac['features'][0]['geometry']).buffer(2000)
    bx = bacia.bounds
    g = pyogrio.read_dataframe(SHP, bbox=bx)
    print(f'{len(g)} segmentos na caixa da bacia')

    # todos os cursos com >= 15 km dentro da bacia, por codigo
    comp = g.assign(L=g.geometry.length).groupby('COCURSODAG')['L'] \
            .sum().sort_values(ascending=False)
    cursos = comp[comp >= 15000]
    print(f'{len(cursos)} cursos >= 15 km')

    # nomes por ancora, quando a ancora acerta o curso
    from shapely.geometry import Point
    nomes = {}
    for nome, (lon, lat) in ANCORAS.items():
        x, y = tr.transform(lon, lat)
        d = g.geometry.distance(Point(x, y))
        k = int(np.argmin(d.values))
        if d.values[k] < 800:
            cod = g.iloc[k]['COCURSODAG']
            if cod in cursos.index and cod not in nomes:
                nomes[cod] = nome
    feats = []
    for cod, L in cursos.items():
        sel = g[g['COCURSODAG'] == cod]
        rotulo = nomes.get(cod, '')
        etiqueta = (f'{rotulo} ' if rotulo else '') + \
            f'[{cod}] {L / 1000:.0f} km'
        if rotulo:
            print(f'  {etiqueta}')
        for geom in sel.geometry:
            if geom is None or geom.is_empty:
                continue
            s = geom.simplify(20.0)
            s = stransform(lambda a, b: inv.transform(a, b), s)
            feats.append({'type': 'Feature',
                          'properties': {'nome': etiqueta,
                                         'curso': str(cod)},
                          'geometry': s.__geo_interface__})
    os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
    json.dump({'type': 'FeatureCollection', 'features': feats},
              open(SAIDA, 'w'), separators=(',', ':'))
    print(f'{len(feats)} tracos -> {SAIDA} '
          f'({os.path.getsize(SAIDA) // 1024} kB)')


if __name__ == '__main__':
    main()
