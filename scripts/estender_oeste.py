# -*- coding: utf-8 -*-
"""Variante g02/u02/p02: Oeste estendido pelo Rio do Campo (curso ANA
775499) ate alem da Barragem Oeste, com a barragem inline.

    python scripts/estender_oeste.py            constroi g98 (so extensao)
    python scripts/estender_oeste.py --tudo     g98 -> barragem -> g02,
                                                u02, p02, registro no prj

NUNCA TOCA g01/u01/p01. A barragem entra pelo construir_barragem.py
(padrao BaldEagle Type 5), com crista MEDIDA do terreno v2 (o voo
2010-12 fotografou o aterro) e fenda de 1,1 m ~ 163 m3/s dos 7
condutos (JICA 2011, Anexo A).

CONTORNO: o hidrograma observado de Taio (que embute a operacao real
de 1983) e MOVIDO para o novo topo, serie inalterada -- a massa que
entrava em Taio agora atravessa o reservatorio. Dupla atenuacao
assumida e documentada: esta variante serve a CENARIOS de operacao,
nao recalibra julho/83.
"""
import json
import os
import re
import shutil
import subprocess
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
sys.path.insert(0, DIR)
os.chdir(RAIZ)

MOSAICO = os.path.join('taha_ai_novo', 'Terrain',
                       'taha_ai_corredor_1m_v2.tif')
ANA_RIOS = os.path.join('doc', 'painel', 'ana_rios.geojson')
RS_TOPO = 56723.30
ANCORA = (599894.10, 7000513.33)        # 1o ponto do Reach XY do R1
DAM_XY = None                            # projetado do SNISB adiante
SNISB_DAM = (-50.03813, -27.09795)
PASSO = 250.0
ALCANCE = 15500.0
MEIA_LARG = 300.0
DX = 2.0
CRLF = '\r\n'


def eixo_extensao():
    """LineString 31982 do curso 775499, da ancora rumo a montante."""
    from shapely.geometry import shape, LineString, Point
    from shapely.ops import linemerge, transform as stransform, \
        substring
    from pyproj import Transformer
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    gj = json.load(open(ANA_RIOS, encoding='utf-8'))
    linhas = [shape(f['geometry']) for f in gj['features']
              if f['properties'].get('curso') == '775499']
    assert linhas, 'curso 775499 ausente em ana_rios.geojson'
    m = linemerge(linhas) if len(linhas) > 1 else linhas[0]
    if m.geom_type == 'MultiLineString':
        m = max(m.geoms, key=lambda g: g.length)
    eixo = stransform(lambda x, y: tr.transform(x, y), m)
    pa = Point(*ANCORA)
    dx, dy = tr.transform(*SNISB_DAM)
    pd = Point(dx, dy)
    sa = eixo.project(pa)
    sd = eixo.project(pd)
    print(f'ancora s={sa:.0f} m, barragem s={sd:.0f} m '
          f'(dist eixo-barragem: {eixo.distance(pd):.0f} m)')
    if sd > sa:
        ext = substring(eixo, sa, min(sa + ALCANCE, eixo.length))
        s_dam = sd - sa
    else:
        ext = substring(eixo, max(0.0, sa - ALCANCE), sa)
        ext = LineString(list(ext.coords)[::-1])
        s_dam = sa - sd
    print(f'extensao: {ext.length:.0f} m; barragem a s={s_dam:.0f} m '
          f'do topo velho (RS {RS_TOPO + s_dam:.1f})')
    return ext, s_dam


def cortar_secoes(ext, s_dam):
    """Secao a cada PASSO m; nada a menos de 120 m da barragem."""
    import rasterio
    src = rasterio.open(MOSAICO)
    nod = src.nodata or -9999.0
    secoes = []
    s = PASSO
    while s <= min(ALCANCE - PASSO, ext.length - PASSO):
        if abs(s - s_dam) < 120.0:
            s += PASSO
            continue
        p = ext.interpolate(s)
        p2 = ext.interpolate(min(s + 10.0, ext.length))
        tx, ty = p2.x - p.x, p2.y - p.y
        n = np.hypot(tx, ty)
        nx, ny = -ty / n, tx / n
        estacoes = np.arange(0.0, 2 * MEIA_LARG + DX / 2, DX)
        xs = p.x + (estacoes - MEIA_LARG) * nx
        ys = p.y + (estacoes - MEIA_LARG) * ny
        z = np.fromiter((v[0] for v in src.sample(zip(xs, ys))),
                        np.float32, count=len(xs))
        z[np.isclose(z, nod)] = np.nan
        ok = np.isfinite(z)
        if ok.sum() < 50:
            s += PASSO
            continue
        z = np.interp(estacoes, estacoes[ok], z[ok])
        secoes.append({'s': s, 'rs': RS_TOPO + s,
                       'x0': xs[0], 'y0': ys[0],
                       'x1': xs[-1], 'y1': ys[-1],
                       'sta': estacoes, 'z': z})
        s += PASSO
    print(f'{len(secoes)} secoes cortadas do v2')
    return secoes


def bloco_secao(sec, comprimento):
    """Texto CRLF de um bloco de secao no dialeto do g01."""
    sta, z = sec['sta'], np.round(sec['z'], 2)
    # bancas: minimo central +- 40 m, na grade
    c0, c1 = len(sta) // 4, 3 * len(sta) // 4
    imin = c0 + int(np.argmin(z[c0:c1]))
    ib0 = max(0, imin - int(40 / DX))
    ib1 = min(len(sta) - 1, imin + int(40 / DX))
    linhas = [f'Type RM Length L Ch R = 1 ,{sec["rs"]:8.2f},'
              f'{comprimento:8.2f},{comprimento:8.2f},'
              f'{comprimento:8.2f}',
              f'Bank Sta={sta[ib0]:.2f},{sta[ib1]:.2f}',
              'XS GIS Cut Line= 2',
              f'{sec["x0"]:16.2f}{sec["y0"]:16.2f}'
              f'{sec["x1"]:16.2f}{sec["y1"]:16.2f}'.rstrip()]
    # rasterio da cutline com 2 pontos: uma linha com os 2 pares
    linhas[3] = (f'{sec["x0"]:16.2f}{sec["y0"]:16.2f}')
    linhas.append(f'{sec["x1"]:16.2f}{sec["y1"]:16.2f}')
    linhas.append(f'#Sta/Elev= {len(sta)} ')
    corpo, lin = [], ''
    for i in range(len(sta)):
        lin += f'{sta[i]:8.2f}{z[i]:8.2f}'
        if (i + 1) % 5 == 0:
            corpo.append(lin)
            lin = ''
    if lin:
        corpo.append(lin)
    linhas.extend(corpo)
    linhas.append('#Mann= 3 , 0 , 0 ')
    linhas.append(f'{0.0:8.3f}{0.06:8.3f}{0.0:8.3f}'
                  f'{sta[ib0]:8.3f}{0.035:8.3f}{0.0:8.3f}'
                  f'{sta[ib1]:8.3f}{0.06:8.3f}{0.0:8.3f}')
    linhas.append('XS Rating Curve= 0 ,0')
    linhas.append('Exp/Cntr=0.3,0.1')
    linhas.append('')
    return CRLF.join(linhas) + CRLF


def montar_g98(ext, secoes):
    txt = open('taha_ai.g01', encoding='latin-1', newline='').read()
    # 1) Reach XY do R1: prepende o eixo da extensao (montante 1o)
    pts_ext = list(ext.coords)[::-1]          # montante -> jusante??
    # ext esta ancora->montante; Reach XY comeca no MONTANTE:
    pts_ext = list(ext.coords)[1:][::-1]      # montante ... ate ancora
    m = re.search(r'(River Reach=Itajai_Oeste {4},R1 +\r?\n'
                  r'Reach XY= *)(\d+)( *\r?\n)((?:.+\r?\n)+?)'
                  r'(Rch Text X Y=)', txt)
    assert m, 'Reach XY do Oeste R1 nao achado'
    velhos = m.group(4)
    novos_pts = pts_ext
    corpo, lin = [], ''
    for i, (x, y) in enumerate(novos_pts):
        lin += f'{x:16.4f}{y:16.4f}'
        if (i + 1) % 2 == 0:
            corpo.append(lin)
            lin = ''
    if lin:
        corpo.append(lin)
    n_total = int(m.group(2)) + len(novos_pts)
    novo_xy = (m.group(1) + str(n_total) + m.group(3)
               + CRLF.join(corpo) + CRLF + velhos + m.group(5))
    txt = txt[:m.start()] + novo_xy + txt[m.end():]
    # 2) blocos das secoes novas antes da secao RS_TOPO
    alvo = 'Type RM Length L Ch R = 1 ,56723.30'
    k = txt.find(alvo)
    assert k > 0
    blocos = ''
    ordenadas = sorted(secoes, key=lambda s: -s['rs'])
    for i, sec in enumerate(ordenadas):
        comprimento = (ordenadas[i]['rs']
                       - (ordenadas[i + 1]['rs']
                          if i + 1 < len(ordenadas) else RS_TOPO))
        blocos += bloco_secao(sec, comprimento)
    txt = txt[:k] + blocos + txt[k:]
    with open('taha_ai.g98', 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(txt)
    print(f'g98 escrito ({len(secoes)} secoes novas)')


def medir_barragem(ext, s_dam):
    """Crista do ATERRO = p99 do quadrado de 400 m no ponto SNISB
    (a travessia pelo eixo pega o canal de desvio, que e baixo --
    media rodada 1: topo 346 m = barragem de 5 m, absurdo)."""
    import rasterio
    from rasterio.windows import from_bounds
    from pyproj import Transformer
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    bx, by = tr.transform(*SNISB_DAM)
    src = rasterio.open(MOSAICO)
    w = from_bounds(bx - 200, by - 200, bx + 200, by + 200,
                    src.transform)
    a = src.read(1, window=w).astype(np.float32)
    a[np.isclose(a, src.nodata or -9999.0)] = np.nan
    crista = float(np.nanpercentile(a, 99))
    print(f'crista do aterro (p99 do quadrado SNISB): {crista:.1f} m; '
          f'leito local {float(np.nanmin(a)):.1f} m')
    return crista


def main(argv):
    ext, s_dam = eixo_extensao()
    secoes = cortar_secoes(ext, s_dam)
    montar_g98(ext, secoes)
    if '--tudo' not in argv:
        print('parou no g98 (rode com --tudo para barragem+u02+p02)')
        return
    crista = medir_barragem(ext, s_dam)
    rs_dam = RS_TOPO + s_dam
    py = sys.executable
    r = subprocess.run(
        [py, os.path.join('scripts', 'construir_barragem.py'),
         'taha_ai.g98', '--saida', 'g02', '--rio', 'Itajai_Oeste',
         '--reach', 'R1', '--rs', f'{rs_dam:.2f}',
         '--crista', f'{crista - 3.0:.1f}', '--topo', f'{crista:.1f}',
         '--larg-vertedouro', '120', '--fenda', '1.1',
         '--nome', 'Barragem Oeste (Taio)'],
        capture_output=True, text=True)
    print(r.stdout[-600:])
    if r.returncode != 0:
        print(r.stderr[-800:])
        raise SystemExit('construir_barragem falhou')
    os.remove('taha_ai.g98')
    # u02: contorno do Oeste movido para o novo topo
    topo_novo = max(s['rs'] for s in secoes)
    u = open('taha_ai.u01', encoding='latin-1', newline='').read()
    u = u.replace('Initial Flow Loc=Itajai_Oeste,R1,56723.3,',
                  f'Initial Flow Loc=Itajai_Oeste,R1,'
                  f'{topo_novo:.2f},')
    u = u.replace('Boundary Location=Itajai_Oeste    ,'
                  'R1              ,56723.30',
                  f'Boundary Location=Itajai_Oeste    ,'
                  f'R1              ,{topo_novo:8.2f}')
    with open('taha_ai.u02', 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(u)
    # p02
    p = open('taha_ai.p01', encoding='latin-1', newline='').read()
    p = p.replace('Geom File=g01', 'Geom File=g02')
    p = p.replace('Flow File=u01', 'Flow File=u02')
    p = re.sub(r'Plan Title=.*', 'Plan Title=1983 + Barragem Oeste',
               p, count=1)
    p = re.sub(r'Short Identifier=.*',
               'Short Identifier=1983_bOeste', p, count=1)
    with open('taha_ai.p02', 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(p)
    # registro no prj (se ainda nao houver)
    prj = open('taha_ai.prj', encoding='latin-1', newline='').read()
    for chave, arq in [('Geom File=g02', 'Geom File=g01'),
                       ('Unsteady File=u02', 'Unsteady File=u01'),
                       ('Plan File=p02', 'Plan File=p01')]:
        if chave not in prj:
            prj = prj.replace(arq + CRLF, arq + CRLF + chave + CRLF, 1)
    with open('taha_ai.prj', 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(prj)
    print(f'g02/u02/p02 prontos; topo novo RS {topo_novo:.2f}, '
          f'barragem RS {rs_dam:.2f}, crista {crista - 3:.1f}/'
          f'{crista:.1f} m')


if __name__ == '__main__':
    main(sys.argv[1:])
