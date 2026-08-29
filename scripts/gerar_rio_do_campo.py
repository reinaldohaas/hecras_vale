# -*- coding: utf-8 -*-
"""Gera Rio_do_Campo (Taio -> Barragem Oeste) como RIO AVULSO, no
padrao provado de gerar_mirim_do_zero.py: HTab explicito por secao e
talvegue forcado monotonico (sem contra-declividade, exceto no vao da
barragem, onde o degrau e real). Depois junta com incorporar_rio.py.

    python scripts/gerar_rio_do_campo.py

Sai projeto_sigsc/Rio_do_Campo.g01 (rio avulso, RS = eixo_extensao()).

POR QUE (dois defeitos achados no metodo anterior, que so prependia
seCoes no MESMO reach do Oeste via regex):
  1. Nenhuma secao nova tinha "XS HTab Starting El and Incr=" nem
     "XS HTab Horizontal Distribution="; o gerador COMPROVADO do Mirim
     sempre escreve as duas em toda secao nova.
  2. O talvegue bruto (extraido direto do MDT) tinha 16 de 58 trechos
     em contra-declividade (a maioria ruido <0.5 m, uma de +3.1 m no
     vao da barragem). O gerador do Mirim IMPOE um talvegue PCHIP
     estritamente decrescente; aqui a correcao e mais simples (sem
     pontos de controle conhecidos): um piso de declividade minima
     (0,0001 m/m) aplicado por deslocamento vertical da secao inteira,
     PULANDO o vao da barragem (gap > 300 m), onde o degrau e real.
"""
import os
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
sys.path.insert(0, DIR)

import estender_oeste as base                          # noqa: E402
from ras_io import escrever                             # noqa: E402
from costurar_rio_do_campo import eixo_por_conectividade  # noqa: E402

base.MOSAICO = os.path.join(RAIZ, 'taha_ai_novo', 'Terrain',
                            'taha_ai_corredor_10m_v2.tif')
SAIDA = os.path.join(RAIZ, 'projeto_sigsc', 'Rio_do_Campo.g01')
MIN_DECL = 0.0001       # m/m -- mesmo piso do gerador do Mirim
GAP_BARRAGEM = 300.0    # m: acima disso e o vao da barragem (pula)


def suavizar(ordenadas):
    """Desloca cada secao (montante->jusante) para nunca subir o
    talvegue, exceto no vao da barragem."""
    for i in range(1, len(ordenadas)):
        prev, cur = ordenadas[i - 1], ordenadas[i]
        dx = prev['rs'] - cur['rs']
        if dx > GAP_BARRAGEM:
            continue          # vao da barragem: degrau real, nao mexe
        z_min_prev = float(np.min(prev['z']))
        teto = z_min_prev - MIN_DECL * dx
        z_min_cur = float(np.min(cur['z']))
        if z_min_cur > teto:
            delta = teto - z_min_cur
            cur['z'] = cur['z'] + delta
    return ordenadas


def bloco_secao_v2(sec, comprimento):
    """Como bloco_secao(), + HTab explicito (padrao do gerador do
    Mirim, unica diferenca estrutural comprovada contra o metodo que
    deu Overflow)."""
    sta, z = sec['sta'], np.round(sec['z'], 2)
    c0, c1 = len(sta) // 4, 3 * len(sta) // 4
    imin = c0 + int(np.argmin(z[c0:c1]))
    ib0 = max(0, imin - int(40 / base.DX))
    ib1 = min(len(sta) - 1, imin + int(40 / base.DX))
    z_min = float(z.min())
    linhas = [f'Type RM Length L Ch R = 1 ,{sec["rs"]:8.2f},'
              f'{comprimento:8.2f},{comprimento:8.2f},'
              f'{comprimento:8.2f}',
              f'Bank Sta={sta[ib0]:.2f},{sta[ib1]:.2f}',
              'XS GIS Cut Line= 2',
              f'{sec["x0"]:16.2f}{sec["y0"]:16.2f}',
              f'{sec["x1"]:16.2f}{sec["y1"]:16.2f}',
              f'#Sta/Elev= {len(sta)} ']
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
    linhas.append(f'XS HTab Starting El and Incr={z_min + 0.02:.2f},'
                  '0.100, 500 ')
    linhas.append('XS HTab Horizontal Distribution=-1,-1,-1')
    linhas.append('XS Rating Curve= 0 ,0')
    linhas.append('Exp/Cntr=0.3,0.1')
    linhas.append('')
    return base.CRLF.join(linhas) + base.CRLF


def main():
    from pyproj import Transformer
    from shapely.geometry import Point
    tr = Transformer.from_crs(4326, 31982, always_xy=True)
    dam_xy = tr.transform(*base.SNISB_DAM)
    ext = eixo_por_conectividade(dam_xy)
    s_dam = ext.project(Point(*dam_xy))
    print(f'eixo real (conectividade): {ext.length:.0f} m; '
          f'barragem a s={s_dam:.0f} m da ancora')
    secoes = base.cortar_secoes(ext, s_dam)
    ordenadas = sorted(secoes, key=lambda s: -s['rs'])
    antes = [float(np.min(s['z'])) for s in ordenadas]
    ordenadas = suavizar(ordenadas)
    depois = [float(np.min(s['z'])) for s in ordenadas]
    ajustes = sum(1 for a, b in zip(antes, depois) if abs(a - b) > 0.01)
    print(f'talvegue suavizado: {ajustes} de {len(ordenadas)} secoes '
          f'ajustadas (piso {MIN_DECL} m/m)')

    linhas = ['Geom Title=Rio_do_Campo', 'Program Version=7.01', '']
    n_pts = len(ext.coords)
    linhas.append('River Reach=Rio_do_Campo   ,R1              ')
    linhas.append(f'Reach XY= {n_pts} ')
    # incorporar_rio.py espera a FOZ no ULTIMO ponto do Reach XY; ext
    # vai da ancora (foz, Taio) ate a cabeceira -- precisa inverter
    pts = list(ext.coords)[::-1]
    corpo, lin = [], ''
    for i, (x, y) in enumerate(pts):
        lin += f'{x:16.4f}{y:16.4f}'
        if (i + 1) % 2 == 0:
            corpo.append(lin)
            lin = ''
    if lin:
        corpo.append(lin)
    linhas.extend(corpo)
    linhas.append('Rch Text X Y=0,0,0,0')
    linhas.append('')

    corpo2 = ''
    for i, sec in enumerate(ordenadas):
        comprimento = (ordenadas[i]['rs']
                       - (ordenadas[i + 1]['rs']
                          if i + 1 < len(ordenadas) else base.RS_TOPO))
        corpo2 += bloco_secao_v2(sec, comprimento)
    texto = base.CRLF.join(linhas) + base.CRLF + corpo2
    os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
    escrever(SAIDA, texto)
    print(f'{len(ordenadas)} secoes -> {SAIDA}')
    return ext, s_dam, ordenadas


if __name__ == '__main__':
    main()
