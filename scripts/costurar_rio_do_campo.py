# -*- coding: utf-8 -*-
"""Traca o eixo real de Taio ate a Barragem Oeste por CONECTIVIDADE
ESPACIAL na malha bruta da ANA (curso_dagua.shp) -- o agrupamento por
codigo Otto (775499) pegou um ramo errado, a 3,1 km da ancora real.

Mesma tecnica do unir_eixo.py (ja provada no projeto): funde segmentos
que se tocam, depois costura por BFS entre componentes a < 300 m.

Retorna um shapely LineString ancora->barragem(+folga), em 31982.
"""
import os

import numpy as np

ANCORA = (599894.10, 7000513.33)
TOL_HOP = 300.0


def eixo_por_conectividade(dam_xy, alem_da_barragem=3000.0):
    import pyogrio
    from shapely.geometry import Point, LineString
    from shapely.ops import linemerge

    ax, ay = ANCORA
    dx, dy = dam_xy
    pad = 4000
    bbox = (min(ax, dx) - pad, min(ay, dy) - pad,
            max(ax, dx) + pad, max(ay, dy) + pad)
    g = pyogrio.read_dataframe(
        r'C:\Users\haas\Downloads\ANA_curso_dagua\curso_dagua.shp',
        bbox=bbox)
    print(f'{len(g)} segmentos brutos na caixa')
    geoms = [geo for geo in g.geometry if geo is not None
             and not geo.is_empty]
    fundidos = linemerge(geoms)
    comps = (list(fundidos.geoms) if fundidos.geom_type
             == 'MultiLineString' else [fundidos])
    print(f'{len(comps)} componentes apos fundir toques diretos')

    p_de = Point(ax, ay)
    p_ate = Point(dx, dy)
    ini = min(range(len(comps)), key=lambda i: comps[i].distance(p_de))
    fim = min(range(len(comps)), key=lambda i: comps[i].distance(p_ate))
    print(f'componente da ancora: dist {comps[ini].distance(p_de):.0f} m '
          f'({comps[ini].length:.0f} m de comprimento)')
    print(f'componente da barragem: dist {comps[fim].distance(p_ate):.0f} m '
          f'({comps[fim].length:.0f} m de comprimento)')

    n = len(comps)
    adj = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            d = comps[i].distance(comps[j])
            if d < TOL_HOP:
                adj[i].append((j, d))
                adj[j].append((i, d))

    # BFS (menor numero de saltos; empate por distancia total)
    import heapq
    dist0 = {ini: 0.0}
    prev = {ini: None}
    fila = [(0.0, ini)]
    visto = set()
    while fila:
        dcur, u = heapq.heappop(fila)
        if u in visto:
            continue
        visto.add(u)
        if u == fim:
            break
        for v, w in adj[u]:
            nd = dcur + w
            if v not in dist0 or nd < dist0[v]:
                dist0[v] = nd
                prev[v] = u
                heapq.heappush(fila, (nd, v))
    if fim not in prev and fim != ini:
        raise SystemExit(
            f'sem caminho ancora->barragem em {TOL_HOP} m de tolerancia')
    caminho = [fim]
    while prev[caminho[-1]] is not None:
        caminho.append(prev[caminho[-1]])
    caminho = caminho[::-1]
    print(f'caminho: {len(caminho)} componentes, saltos: '
          f'{[round(comps[caminho[k]].distance(comps[caminho[k+1]]), 0) for k in range(len(caminho)-1)]}')

    # concatena SO O TRECHO USADO de cada componente (substring entre
    # a posicao de entrada e o alvo desta etapa) -- um componente pode
    # ter dezenas de km dos quais so um pedaco pertence ao caminho
    from shapely.ops import substring
    trechos = []
    pos = p_de
    for k, idx in enumerate(caminho):
        c = comps[idx]
        alvo = (p_ate if k == len(caminho) - 1
                else comps[caminho[k + 1]].interpolate(
                    0.5, normalized=True))
        s0, s1 = c.project(pos), c.project(alvo)
        if abs(s1 - s0) < 1.0:
            pos = c.interpolate(s1)
            continue
        if k == len(caminho) - 1 and alem_da_barragem > 0:
            sinal = 1.0 if s1 >= s0 else -1.0
            s1 = float(np.clip(s1 + sinal * alem_da_barragem,
                               0.0, c.length))
        t = (LineString(list(substring(c, s1, s0).coords)[::-1])
             if s0 > s1 else substring(c, s0, s1))
        trechos.append(t)
        pos = Point(t.coords[-1])
    pontos = list(trechos[0].coords)
    for t in trechos[1:]:
        pontos.extend(list(t.coords)[1:])
    linha = LineString(pontos)
    print(f'eixo costurado: {linha.length:.0f} m total')
    return linha


if __name__ == '__main__':
    from pyproj import Transformer
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    dam_xy = tr.transform(-50.03813, -27.09795)
    linha = eixo_por_conectividade(dam_xy)
    print('primeiro ponto (ancora):', linha.coords[0])
    print('ultimo ponto:', linha.coords[-1])
