# -*- coding: utf-8 -*-
"""Converte as linhas do editor (edges/banks/divisor) em .geojson.

    python scripts/exportar_linhas_geojson.py doc/linhas_hand/linhas_editor.json
    python scripts/exportar_linhas_geojson.py editor_salvo.html --tipo edge
    python scripts/exportar_linhas_geojson.py ... --epsg 4326 --saida edges.geojson

ACEITA como entrada:
  - o JSON do editor (doc/linhas_hand/linhas_editor.json), ou
  - o PROPRIO HTML salvo do artifact (extrai o bloco <script id="dados">)
    -- e assim que as linhas EDITADAS no navegador voltam para o SIG.

FILTROS
  --tipo edge|bank|modelo|divisor|todas   (padrao: edge -- o pedido)
  --rio Acu|Mirim|...                     (fragmento do nome do grupo)
  --epsg 31982|4326                       (padrao 31982; 4326 reprojeta
                                           para lat/lon, ex.: QGIS/web)
"""
import json
import os
import re
import sys


def _arg(argv, chave, padrao=None):
    return argv[argv.index(chave) + 1] if chave in argv else padrao


def carregar(caminho):
    t = open(caminho, encoding='utf-8', errors='replace').read()
    if caminho.lower().endswith('.html') or '<script' in t[:2000]:
        m = re.search(r'<script id="dados" type="application/json">(.*?)'
                      r'</script>', t, flags=re.S)
        if not m:
            raise SystemExit('nao achei o bloco de dados no HTML')
        return json.loads(m.group(1))
    return json.loads(t)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    entrada = argv[0]
    tipo = _arg(argv, '--tipo', 'edge').lower()
    rio = _arg(argv, '--rio')
    epsg = int(_arg(argv, '--epsg', '31982'))
    saida = _arg(argv, '--saida')
    if saida is None:
        base = os.path.splitext(os.path.basename(entrada))[0]
        saida = f'{base}_{tipo}.geojson'

    dados = carregar(entrada)

    def quer(l):
        nome = l['nome'].lower()
        grupo = (l.get('grupo') or '').lower()
        if tipo != 'todas':
            if tipo == 'divisor' and not nome.startswith('divisor'):
                return False
            if tipo == 'modelo' and '(modelo)' not in nome:
                return False
            if tipo in ('edge', 'bank'):
                if f'{tipo} ' not in nome and not nome.startswith(tipo):
                    return False
                if '(modelo)' in nome:
                    return False
        if rio and rio.lower() not in grupo and rio.lower() not in nome:
            return False
        return True

    escolhidas = [l for l in dados['linhas'] if quer(l)]
    if not escolhidas:
        raise SystemExit(f'nenhuma linha casa com tipo={tipo} rio={rio}')

    reproj = None
    if epsg != 31982:
        from pyproj import Transformer
        reproj = Transformer.from_crs(31982, epsg, always_xy=True)

    feats = []
    for l in escolhidas:
        pts = l['pontos']
        if reproj:
            pts = [list(reproj.transform(x, y)) for x, y in pts]
            pts = [[round(x, 6), round(y, 6)] for x, y in pts]
        feats.append({'type': 'Feature',
                      'properties': {'nome': l['nome'],
                                     'grupo': l.get('grupo', ''),
                                     'cor': l.get('cor', '')},
                      'geometry': {'type': 'LineString',
                                   'coordinates': pts}})
    geo = {'type': 'FeatureCollection',
           'crs': {'type': 'name',
                   'properties': {'name': f'EPSG:{epsg}'}},
           'features': feats}
    with open(saida, 'w', encoding='utf-8') as f:
        json.dump(geo, f, ensure_ascii=False, indent=1)
    print(f'{len(feats)} linhas ({tipo}'
          f"{', rio ' + rio if rio else ''}) -> {saida}  [EPSG:{epsg}]")


if __name__ == '__main__':
    main(sys.argv[1:])
