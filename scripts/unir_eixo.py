# -*- coding: utf-8 -*-
"""Une um eixo de rio partido, costurando o meio pela hidrografia FBDS.

    python scripts/unir_eixo.py --rio Rio_Luis_Alves \
        --rede vale_itajai_full_network.geojson --fbds doc/fbds

Pega TODOS os componentes nomeados do rio na rede completa (lat/lon),
ordena da cabeceira para a foz e costura os vaos caminhando pela
hidrografia da FBDS (RIOS_SIMPLES + contornos dos RIOS_DUPLOS): a cada
passo entra o segmento cuja ponta esta a < TOL m da ponta atual e cuja
OUTRA ponta mais se aproxima do proximo componente. Grava o eixo unido
em eixos_do_relevo.geojson (substitui a entrada do rio).
"""
import json
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)

TOL = 400.0


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    import glob
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import LineString, Point
    from shapely.ops import linemerge
    from pyproj import Transformer

    nome_mod = _arg(argv, '--rio')
    rede_arq = _arg(argv, '--rede', 'vale_itajai_full_network.geojson')
    fbds = _arg(argv, '--fbds', 'doc/fbds')
    nome_fonte = {'Rio_Luis_Alves': 'Rio Luís Alves',
                  'Rio_Krauel': 'Rio Krauel'}.get(nome_mod, nome_mod)

    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    rede = json.load(open(rede_arq, encoding='utf-8'))
    geoms = []
    for f in rede['features']:
        if (f['properties'].get('nome') or
                f['properties'].get('name')) != nome_fonte:
            continue
        g = f['geometry']
        cs = [g['coordinates']] if g['type'] == 'LineString' \
            else g['coordinates']
        for c in cs:
            geoms.append(LineString([tr.transform(x, y) for x, y in c]))
    m = linemerge(geoms)
    comps = [m] if m.geom_type == 'LineString' else \
        sorted(m.geoms, key=lambda g: -g.length)
    print(f'{nome_fonte}: {len(comps)} componentes nomeados')

    # ordena da cabeceira (mais alto? usa o mais longo como ancora e o
    # resto pela proximidade da foz -- aqui: por Y nao serve; ordena
    # pela distancia ao 1o componente) -- pratico: encadeia guloso
    atual = comps[0]
    resto = comps[1:]
    cadeia = [atual]
    while resto:
        fim = Point(cadeia[-1].coords[-1])
        prox = min(resto, key=lambda c: min(fim.distance(Point(c.coords[0])),
                                            fim.distance(Point(c.coords[-1]))))
        resto.remove(prox)
        if fim.distance(Point(prox.coords[-1])) < \
                fim.distance(Point(prox.coords[0])):
            prox = LineString(list(prox.coords)[::-1])
        cadeia.append(prox)

    def costurar(p_de, p_ate):
        """carrega a hidrografia SO do corredor do vao (bbox), funde e
        acha o caminho entre as pontas por BFS de componentes."""
        from shapely.ops import substring
        corredor = LineString([p_de, p_ate]).buffer(1500)
        dentro = []
        for camada in ('RIOS_SIMPLES', 'RIOS_DUPLOS'):
            for shp in glob.glob(os.path.join(fbds, '*',
                                              f'*_{camada}.shp')):
                try:
                    g = gpd.read_file(shp, bbox=corredor.bounds)
                    for geom in g.geometry:
                        if geom is None:
                            continue
                        if geom.geom_type == 'LineString':
                            dentro.append(geom)
                        elif geom.geom_type == 'MultiLineString':
                            dentro += list(geom.geoms)
                        elif geom.geom_type == 'Polygon':
                            dentro.append(
                                LineString(geom.exterior.coords))
                except Exception:
                    pass
        dentro = [f for f in dentro if f.intersects(corredor)]
        if not dentro:
            return None
        fundido = linemerge(dentro)
        comps2 = [fundido] if fundido.geom_type == 'LineString' \
            else list(fundido.geoms)
        # BFS entre componentes (conexos se distancia < 120 m): a FBDS
        # quebra o tronco na transicao SIMPLES <-> DUPLOS
        n = len(comps2)
        ini = min(range(n), key=lambda i: comps2[i].distance(p_de))
        fim = min(range(n), key=lambda i: comps2[i].distance(p_ate))
        if comps2[ini].distance(p_de) > TOL \
                or comps2[fim].distance(p_ate) > TOL:
            return None
        de_onde = {ini: None}
        fila = [ini]
        while fila and fim not in de_onde:
            i = fila.pop(0)
            for j in range(n):
                if j in de_onde:
                    continue
                if comps2[i].distance(comps2[j]) < 300:
                    de_onde[j] = i
                    fila.append(j)
        if fim not in de_onde:
            return None
        rota = []
        k = fim
        while k is not None:
            rota.append(k)
            k = de_onde[k]
        rota.reverse()
        trechos = []
        pos = p_de
        for idx, k in enumerate(rota):
            c = comps2[k]
            alvo2 = p_ate if idx == len(rota) - 1 \
                else Point(np.asarray(comps2[rota[idx + 1]]
                                      .interpolate(0.5,
                                                   normalized=True).coords[0]))
            s0, s1 = c.project(pos), c.project(alvo2)
            if abs(s1 - s0) < 1:
                pos = c.interpolate(s1)
                continue
            if s0 > s1:
                t2 = substring(c, s1, s0)
                t2 = LineString(list(t2.coords)[::-1])
            else:
                t2 = substring(c, s0, s1)
            trechos.append(t2)
            pos = Point(t2.coords[-1])
        return trechos if trechos else None

    coords = list(cadeia[0].coords)
    for prox in cadeia[1:]:
        p_de = Point(coords[-1])
        p_ate = Point(prox.coords[0])
        vao = p_de.distance(p_ate)
        if vao > TOL:
            fio = costurar(p_de, p_ate)
            if fio is None:
                raise SystemExit(f'nao consegui costurar vao de '
                                 f'{vao:.0f} m -- hidrografia nao liga')
            for f in fio:
                coords += list(f.coords)
            print(f'   vao de {vao:.0f} m costurado com '
                  f'{len(fio)} fios da FBDS')
        coords += list(prox.coords)
    eixo = LineString(coords).simplify(15)
    print(f'eixo unido: {eixo.length/1000:.1f} km, {len(eixo.coords)} pts')

    eixos = json.load(open('eixos_do_relevo.geojson', encoding='utf-8'))
    eixos['features'] = [f for f in eixos['features']
                         if f['properties'].get('nome') != nome_mod]
    eixos['features'].append({'type': 'Feature',
        'properties': {'rio': nome_mod.lower(), 'nome': nome_mod,
                       'km': round(eixo.length / 1000, 1)},
        'geometry': {'type': 'LineString',
                     'coordinates': [[round(x, 1), round(y, 1)]
                                     for x, y in eixo.coords]}})
    json.dump(eixos, open('eixos_do_relevo.geojson', 'w'),
              ensure_ascii=False)
    print('eixos_do_relevo.geojson atualizado')


if __name__ == '__main__':
    main(sys.argv[1:])
