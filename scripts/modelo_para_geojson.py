# -*- coding: utf-8 -*-
"""Converte as linhas DO MODELO taha_ai (g01.hdf) em .geojson.

    python scripts/modelo_para_geojson.py taha_ai.g01.hdf --tipo edge
    python scripts/modelo_para_geojson.py taha_ai.g01.hdf --tipo bank --rio Mirim
    python scripts/modelo_para_geojson.py taha_ai.g01.hdf --tipo cutline --epsg 4326

FONTE AUTORITATIVA: o proprio HDF da geometria que o HEC-RAS computa
(`Geometry/Cross Sections`): cutlines georreferenciadas (Polyline
Points), `Left Bank`/`Right Bank` e as estacoes (Station Elevation).
A posicao do banco na cutline segue a regra do RAS: o vao de estacoes
da secao e mapeado LINEARMENTE no comprimento da cutline.

  --tipo edge     pontas das secoes (limite do modelo), N e S por rio
  --tipo bank     Bank Sta reais, N e S por rio
  --tipo cutline  cada secao como sua propria LineString
  --rio <frag>    filtra pelo nome do rio
  --epsg 4326     reprojeta (padrao: 31982, UTM 22S SIRGAS2000)
"""
import json
import os
import sys

import numpy as np


def _arg(argv, chave, padrao=None):
    return argv[argv.index(chave) + 1] if chave in argv else padrao


def partes(P, kms, salto=2500.0):
    o = np.argsort(-np.asarray(kms))
    pts = np.asarray(P)[o]
    blocos, atual = [], [pts[0]]
    for p in pts[1:]:
        if np.hypot(p[0] - atual[-1][0], p[1] - atual[-1][1]) > salto:
            blocos.append(atual)
            atual = [p]
        else:
            atual.append(p)
    blocos.append(atual)
    return [b for b in blocos if len(b) >= 2]


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    import h5py
    hdf = argv[0]
    tipo = _arg(argv, '--tipo', 'edge').lower()
    rio_f = (_arg(argv, '--rio') or '').lower()
    epsg = int(_arg(argv, '--epsg', '31982'))
    saida = _arg(argv, '--saida')
    if saida is None:
        saida = os.path.splitext(os.path.basename(hdf))[0].split('.')[0] \
            + f'_{tipo}.geojson'

    f = h5py.File(hdf, 'r')
    g = f['Geometry/Cross Sections']
    at = g['Attributes'][()]
    pinfo = g['Polyline Info'][()]
    ppts = g['Polyline Points'][()]
    sinfo = g['Station Elevation Info'][()]
    svals = g['Station Elevation Values'][()]

    reproj = None
    if epsg != 31982:
        from pyproj import Transformer
        reproj = Transformer.from_crs(31982, epsg, always_xy=True)

    def repro(pts):
        if not reproj:
            return [[round(float(x), 2), round(float(y), 2)]
                    for x, y in pts]
        return [[round(v, 6) for v in reproj.transform(x, y)]
                for x, y in pts]

    por_rio = {}
    cutlines = []
    for i in range(len(at)):
        rio = at['River'][i].decode().strip()
        reach = at['Reach'][i].decode().strip()
        try:
            rs = float(at['RS'][i].decode().strip().rstrip('*'))
        except ValueError:
            continue
        if rio_f and rio_f not in rio.lower():
            continue
        i0, n = int(pinfo[i][0]), int(pinfo[i][1])
        if n < 2:
            continue
        C = ppts[i0:i0 + n].astype(float)
        seg = np.hypot(*np.diff(C, axis=0).T)
        L = float(seg.sum())
        s_acum = np.concatenate([[0.0], np.cumsum(seg)])
        j0, m = int(sinfo[i][0]), int(sinfo[i][1])
        stas = svals[j0:j0 + m, 0].astype(float)
        s_ini, s_fim = float(stas[0]), float(stas[-1])

        def xy(sta):
            # regra do RAS: vao de estacoes mapeado linearmente na cutline
            fr = (sta - s_ini) / max(s_fim - s_ini, 1e-9)
            s = min(max(fr, 0.0), 1.0) * L
            x = float(np.interp(s, s_acum, C[:, 0]))
            y = float(np.interp(s, s_acum, C[:, 1]))
            return (x, y)
        r = por_rio.setdefault(rio, {'eS': [], 'eN': [], 'bS': [],
                                     'bN': [], 'rs': []})
        r['eS'].append(tuple(C[0]))
        r['eN'].append(tuple(C[-1]))
        r['bS'].append(xy(float(at['Left Bank'][i])))
        r['bN'].append(xy(float(at['Right Bank'][i])))
        r['rs'].append(rs)
        cutlines.append((rio, reach, rs, C))

    feats = []
    if tipo == 'cutline':
        for rio, reach, rs, C in cutlines:
            feats.append({'type': 'Feature',
                          'properties': {'rio': rio, 'reach': reach,
                                         'rs': rs},
                          'geometry': {'type': 'LineString',
                                       'coordinates': repro(C)}})
    elif tipo in ('edge', 'bank'):
        chaves = {'edge': [('eN', 'N'), ('eS', 'S')],
                  'bank': [('bN', 'N'), ('bS', 'S')]}[tipo]
        for rio, r in sorted(por_rio.items()):
            for chave, tag in chaves:
                for k, bloco in enumerate(partes(r[chave], r['rs'])):
                    feats.append({'type': 'Feature',
                                  'properties': {'rio': rio, 'lado': tag,
                                                 'parte': k + 1,
                                                 'tipo': tipo},
                                  'geometry': {'type': 'LineString',
                                               'coordinates': repro(bloco)}})
    else:
        raise SystemExit(f'tipo desconhecido: {tipo}')

    geo = {'type': 'FeatureCollection',
           'crs': {'type': 'name',
                   'properties': {'name': f'EPSG:{epsg}'}},
           'features': feats}
    with open(saida, 'w', encoding='utf-8') as fo:
        json.dump(geo, fo, ensure_ascii=False, indent=1)
    print(f'{len(feats)} feicoes ({tipo}'
          f"{', rio ' + rio_f if rio_f else ''}) de {len(at)} secoes "
          f'-> {saida}  [EPSG:{epsg}]')


if __name__ == '__main__':
    main(sys.argv[1:])
