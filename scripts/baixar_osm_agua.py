# -*- coding: utf-8 -*-
"""Baixa do OpenStreetMap (Overpass) a agua da bacia do Itajai:
rios (linhas), espelhos (poligonos), PONTES sobre rios e REPRESAS.

    python scripts/baixar_osm_agua.py

Sai em doc/osm/: osm_rios.geojson, osm_agua.geojson,
osm_pontes.geojson, osm_represas.geojson -- tudo recortado pelo
divisor oficial ANA (doc/qgis/bacia_itajai_ana.geojson) e em 4326.
"""
import json
import os
import time
import urllib.request

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
os.chdir(RAIZ)
APIS = ['https://overpass.kumi.systems/api/interpreter',
        'https://lz4.overpass-api.de/api/interpreter',
        'https://overpass-api.de/api/interpreter']
# bacia em 4 sub-caixas (sul, oeste, norte, leste) p/ nao estourar
LAT0, LON0, LAT1, LON1 = -27.85, -50.65, -26.30, -48.55
LATM, LONM = (LAT0 + LAT1) / 2, (LON0 + LON1) / 2
CAIXAS = [f'{LAT0},{LON0},{LATM},{LONM}',
          f'{LAT0},{LONM},{LATM},{LON1}',
          f'{LATM},{LON0},{LAT1},{LONM}',
          f'{LATM},{LONM},{LAT1},{LON1}']

CONSULTAS = {
    'osm_rios': 'way["waterway"~"^(river|canal)$"]({bb});',
    'osm_agua': ('way["natural"="water"]({bb});'
                 'way["waterway"="riverbank"]({bb});'),
    'osm_pontes': ('way["bridge"]["highway"]({bb});'
                   'way["bridge"]["railway"]({bb});'),
    'osm_represas': ('way["waterway"~"^(dam|weir)$"]({bb});'
                     'node["waterway"~"^(dam|weir)$"]({bb});'),
}


def overpass(q):
    corpo = f'[out:json][timeout:120];({q});out geom;'
    ultimo = None
    for api in APIS:
        for tentativa in range(2):
            try:
                req = urllib.request.Request(
                    api, data=corpo.encode(),
                    headers={'User-Agent': 'hecras-vale-itajai/1.0'})
                with urllib.request.urlopen(req, timeout=200) as r:
                    return json.load(r)
            except Exception as e:
                ultimo = e
                print(f'    ({api.split("/")[2]}: {e}; tentando de '
                      f'novo)', flush=True)
                time.sleep(10)
    raise ultimo


def para_geojson(dado, bacia):
    from shapely.geometry import LineString, Polygon, Point, shape
    from shapely.prepared import prep
    pb = prep(bacia)
    feats = []
    for el in dado.get('elements', []):
        tags = el.get('tags', {})
        if el['type'] == 'way' and 'geometry' in el:
            pts = [(p['lon'], p['lat']) for p in el['geometry']]
            if len(pts) < 2:
                continue
            fechado = pts[0] == pts[-1] and len(pts) >= 4
            g = Polygon(pts) if fechado else LineString(pts)
        elif el['type'] == 'node':
            g = Point(el['lon'], el['lat'])
        else:
            continue
        if not pb.intersects(g):
            continue
        feats.append({'type': 'Feature',
                      'properties': {'nome': tags.get('name', ''),
                                     'tipo': (tags.get('waterway')
                                              or tags.get('natural')
                                              or tags.get('bridge'))},
                      'geometry': g.__geo_interface__})
    return {'type': 'FeatureCollection', 'features': feats}


def main():
    from shapely.geometry import shape
    from shapely.ops import transform as stransform
    from pyproj import Transformer
    inv = Transformer.from_crs(31982, 4326, always_xy=True)
    bac = json.load(open('doc/qgis/bacia_itajai_ana.geojson',
                         encoding='utf-8'))
    bacia = stransform(lambda x, y: inv.transform(x, y),
                       shape(bac['features'][0]['geometry'])
                       ).buffer(0.005)
    os.makedirs('doc/osm', exist_ok=True)
    for nome, q in CONSULTAS.items():
        destino = f'doc/osm/{nome}.geojson'
        print(f'{nome}...', flush=True)
        vistos, elementos = set(), []
        for bb in CAIXAS:
            dado = overpass(q.format(bb=bb))
            for el in dado.get('elements', []):
                chave = (el['type'], el['id'])
                if chave in vistos:
                    continue
                vistos.add(chave)
                elementos.append(el)
            time.sleep(3)      # cortesia com o servidor publico
        gj = para_geojson({'elements': elementos}, bacia)
        json.dump(gj, open(destino, 'w'), separators=(',', ':'))
        print(f'  {len(gj["features"])} feicoes -> {destino} '
              f'({os.path.getsize(destino) // 1024} kB)')


if __name__ == '__main__':
    main()
