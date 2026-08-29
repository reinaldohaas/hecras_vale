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
    """LineString 31982 real (por conectividade espacial na malha
    bruta da ANA) da ancora ate alem da Barragem Oeste.

    O metodo antigo agrupava por codigo Otto (775499) e projetava a
    ancora na linha mais longa com esse codigo -- a ancora ficou a
    3131 m dessa linha (codigo ERRADO; o de verdade, 775495882, passa
    a 20 m). Root cause do "Overflow"/"Error plotting cross section
    lines": um teleporte de 3 km disfarcado de reach continuo.
    """
    from pyproj import Transformer
    from shapely.geometry import Point
    from costurar_rio_do_campo import eixo_por_conectividade
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    dam_xy = tr.transform(*SNISB_DAM)
    ext = eixo_por_conectividade(dam_xy)   # coords[0]=ancora (Taio)
    s_dam = ext.project(Point(*dam_xy))     # distancia ancora->barragem
    print(f'extensao real: {ext.length:.0f} m; barragem a s={s_dam:.0f} m '
          f'da ancora (RS {RS_TOPO + s_dam:.1f})')
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
    # bancas: minimo central +- largura MEDIDA no SIG-SC (sec['largura_
    # alvo']/2); sem medida, cai no chute antigo de 40 m
    meia = sec.get('largura_alvo', 80.0) / 2.0
    c0, c1 = len(sta) // 4, 3 * len(sta) // 4
    imin = c0 + int(np.argmin(z[c0:c1]))
    ib0 = max(0, imin - int(meia / DX))
    ib1 = min(len(sta) - 1, imin + int(meia / DX))
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
    # HTab explicito -- padrao do gerador comprovado (gerar_mirim_do_
    # zero.py); secoes novas sem isto foram a causa provavel do
    # "Overflow" (junto com a ancora errada, ja corrigida acima)
    z_min = float(z.min())
    linhas.append(f'XS HTab Starting El and Incr={z_min + 0.02:.2f},'
                  '0.100, 500 ')
    linhas.append('XS HTab Horizontal Distribution=-1,-1,-1')
    linhas.append('XS Rating Curve= 0 ,0')
    linhas.append('Exp/Cntr=0.3,0.1')
    linhas.append('')
    return CRLF.join(linhas) + CRLF


GAP_BARRAGEM = 300.0
MIN_DECL = 0.0001


def largura_alvo_serie(csv_path, janela_km=1.0, cada_km=0.25,
                       fator=1.0):
    """(d_km ordenado, largura suavizada) do CSV de largura_do_sigsc,
    ou None. Mediana movel de `janela_km`, igual ao encolher_canal.py."""
    import csv
    if not os.path.exists(csv_path):
        return None
    pares = []
    for r in csv.reader(open(csv_path, encoding='utf-8'),
                        delimiter=';'):
        if r[0] == 'dist_foz_km' or len(r) < 3 or not r[2]:
            continue
        pares.append((float(r[0]), float(r[2]) * fator))
    if len(pares) < 3:
        return None
    pares.sort()
    d = np.array([p[0] for p in pares])
    w = np.array([p[1] for p in pares])
    meia = max(1, int(round(janela_km / cada_km / 2)))
    suave = np.array([np.median(w[max(0, i - meia):i + meia + 1])
                      for i in range(len(w))])
    return d, suave


def aplicar_largura_medida(ordenadas, serie, minimo=20.0):
    """Preenche sec['largura_alvo'] por interpolacao na quilometragem
    (sec['s'] = distancia da FOZ/ancora, m)."""
    if serie is None:
        return ordenadas
    d_km, w = serie
    for sec in ordenadas:
        alvo = float(np.interp(sec['s'] / 1000.0, d_km, w))
        sec['largura_alvo'] = max(alvo, minimo)
    return ordenadas


def suavizar_talvegue(ordenadas):
    """Desloca cada secao (montante->jusante) p/ nunca subir o
    talvegue, exceto no vao da barragem (degrau real, nao mexe)."""
    for i in range(1, len(ordenadas)):
        prev, cur = ordenadas[i - 1], ordenadas[i]
        dx = prev['rs'] - cur['rs']
        if dx > GAP_BARRAGEM:
            continue
        teto = float(np.min(prev['z'])) - MIN_DECL * dx
        z_min_cur = float(np.min(cur['z']))
        if z_min_cur > teto:
            cur['z'] = cur['z'] + (teto - z_min_cur)
    return ordenadas


def montar_g98(ext, secoes):
    txt = open('taha_ai.g01', encoding='latin-1', newline='').read()
    # 1) Reach XY do R1: prepende o eixo da extensao. ext vai
    # ancora(coords[0])->cabeceira(coords[-1]); Reach XY comeca no
    # MONTANTE (cabeceira), entao inverte e tira a ancora (ja e o 1o
    # ponto do bloco velho, nao duplicar)
    pts_ext = list(ext.coords)[1:][::-1]
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
    ordenadas = suavizar_talvegue(ordenadas)
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
