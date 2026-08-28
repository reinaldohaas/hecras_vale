# -*- coding: utf-8 -*-
"""Incorpora um rio avulso (g01 proprio) na rede, com juncao e contorno.

    python scripts/incorporar_rio.py taha_ai.g02 \
        --rio-novo modelos/krauel/krauel.g01 --nome Rio_Krauel \
        --alvo Itajai_Norte,R2 --saida g04 \
        --u01 taha_ai.u01 --u01-saida taha_ai.u04 \
        --serie doc/ana_1983/83440000_Ibirama_vazao.csv --escala 0.11

A REDE DE ENTRADA NAO E TOCADA (sai gNN novo; u01 novo em --u01-saida).

O QUE SE FAZ

  1. A foz do rio novo e projetada no eixo do reach alvo; o reach e
     PARTIDO na fronteira de secoes mais proxima: a parte de montante
     mantem o nome, a de jusante vira o proximo RN livre (Reach XY
     partida no vertice mais proximo; secoes e comprimentos intactos).
  2. Referencias ao reach de jusante mudam de nome onde preciso:
     juncao existente a jusante e laterais do u01 (lateral que atravessa
     o corte e dividida por comprimento, como dividir_lateral).
  3. O rio novo entra como bloco proprio (rio,R1) e a juncao nova e
     escrita no padrao das irmas (Junct Name/Desc/XY/Up/Dn/Junc L&A).
  4. u01: cabeceira do rio novo ganha Flow Hydrograph diario da
     `--serie` x `--escala` (janela 01JUL-05AGO/1983), Initial Flow Loc
     do dia 1, e a parte de jusante do alvo ganha Initial Flow Loc.

  CONFERENCIA relendo o gravado: secoes da rede + do novo = total;
  juncao presente; nenhum 'Initial RS='; contagem de reaches.
"""
import csv
import datetime
import os
import re
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
from qc_secoes import ler_secoes                       # noqa: E402
from qc_geometria import ler_eixos                     # noqa: E402
from ras_io import escrever                            # noqa: E402

INICIO = datetime.date(1983, 7, 1)
FIM = datetime.date(1983, 8, 5)


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def fmt8(vals):
    corpo, lin = [], ""
    for i, x in enumerate(vals):
        lin += "%8.2f" % x
        if (i + 1) % 10 == 0:
            corpo.append(lin)
            lin = ""
    if lin:
        corpo.append(lin)
    return corpo


def serie_diaria(arq):
    vals = {}
    for r in csv.reader(open(arq, encoding='utf-8'), delimiter=';'):
        if r[0] == 'data' or len(r) < 2 or not r[1]:
            continue
        d = datetime.date.fromisoformat(r[0])
        if INICIO <= d <= FIM:
            vals[d] = float(r[1])
    dias = (FIM - INICIO).days + 1
    s = np.full(dias, np.nan)
    for d, v in vals.items():
        s[(d - INICIO).days] = v
    idx = np.arange(dias)
    ok = ~np.isnan(s)
    return np.interp(idx, idx[ok], s[ok])


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    rede = argv[0]
    g_novo = _arg(argv, '--rio-novo')
    nome_rio = _arg(argv, '--nome')
    alvo_rio, alvo_reach = [x.strip() for x in
                            _arg(argv, '--alvo').split(',')]
    ext = _arg(argv, '--saida', 'g04')
    u01_in = _arg(argv, '--u01', 'taha_ai.u01')
    u01_out = _arg(argv, '--u01-saida', 'taha_ai.u04')
    serie_arq = _arg(argv, '--serie')
    escala = _arg(argv, '--escala', 1.0, float)
    raiz = os.path.dirname(rede) or '.'
    base = os.path.basename(rede).split('.')[0]
    saida = os.path.join(raiz, f'{base}.{ext}')

    # ---------------- foz e ponto de corte
    E = ler_eixos(rede)
    En = ler_eixos(g_novo)
    from shapely.geometry import Point
    ls_novo = list(En.values())[0]
    ls_alvo = E[(alvo_rio, alvo_reach)]
    foz = Point(ls_novo.coords[-1])
    s_foz = ls_alvo.project(foz)
    dist_foz = foz.distance(ls_alvo.interpolate(s_foz))
    rs_foz = ls_alvo.length - s_foz
    S = ler_secoes(rede)
    secs_alvo = sorted([d['rs'] for d in S if d['rio'] == alvo_rio
                        and d['reach'] == alvo_reach and d['tipo'] == '1'],
                       reverse=True)
    # corte entre a secao imediatamente ACIMA da foz e a de baixo
    acima = [r for r in secs_alvo if r >= rs_foz]
    abaixo = [r for r in secs_alvo if r < rs_foz]
    if not acima or not abaixo:
        raise SystemExit(f'foz (RS {rs_foz:.0f}) fora do miolo do reach')
    rs_corte_up, rs_corte_dn = min(acima), max(abaixo)
    reaches_rio = sorted({d['reach'] for d in S if d['rio'] == alvo_rio})
    novo_reach = f'R{max(int(r[1:]) for r in reaches_rio) + 1}'
    print(f'{nome_rio}: foz a {dist_foz:.0f} m do eixo de '
          f'{alvo_rio} {alvo_reach}, RS {rs_foz:.0f}')
    print(f'   corte entre RS {rs_corte_up:.2f} e {rs_corte_dn:.2f}; '
          f'jusante vira {novo_reach}')

    # ---------------- texto da rede
    t = open(rede, encoding='latin-1', errors='replace').read() \
        .replace('\r\n', '\n').replace('\r', '\n')
    linhas = t.split('\n')

    # acha o bloco do reach alvo e o indice da linha do 1o Type de jusante
    i_reach = i_fim = i_split = None
    rio_c = reach_c = None
    for i, l in enumerate(linhas):
        if l.startswith('River Reach='):
            p = l.split('=', 1)[1].split(',')
            r_novo, rc = p[0].strip(), p[1].strip()
            if i_reach is not None and i_fim is None:
                i_fim = i
            if r_novo == alvo_rio and rc == alvo_reach:
                i_reach = i
            rio_c, reach_c = r_novo, rc
        if (i_reach is not None and i_fim is None
                and l.startswith('Type RM Length L Ch R')):
            try:
                rs_l = float(l.split('=', 1)[1].split(',')[1])
            except (ValueError, IndexError):
                continue
            if abs(rs_l - rs_corte_dn) < 0.01 and i_split is None:
                i_split = i
    if i_fim is None:
        i_fim = len(linhas)
    if i_split is None:
        raise SystemExit('nao achei a secao de corte no texto')

    # Reach XY do alvo: partir no vertice mais proximo do corte
    p_corte = ls_alvo.interpolate(ls_alvo.length - (rs_corte_up
                                                    + rs_corte_dn) / 2)
    verts = np.array(ls_alvo.coords)
    k_corte = int(np.argmin(np.hypot(verts[:, 0] - p_corte.x,
                                     verts[:, 1] - p_corte.y)))
    up_xy = verts[:k_corte + 1]
    dn_xy = verts[k_corte:]

    def bloco_xy(V):
        cab = ['Reach XY= %d ' % len(V)]
        vals = []
        for x, y in V:
            vals += [x, y]
        lin = ''
        corpo = []
        for i, v in enumerate(vals):
            lin += '%16.4f' % v
            if (i + 1) % 4 == 0:
                corpo.append(lin)
                lin = ''
        if lin:
            corpo.append(lin)
        return cab + corpo

    # reescreve: cabecalho do alvo com XY de montante; no corte, novo
    # cabecalho de reach com XY de jusante
    novo_txt = []
    i = 0
    while i < len(linhas):
        l = linhas[i]
        if i == i_reach:
            novo_txt.append(l)
            i += 1
            # pula o Reach XY antigo
            if linhas[i].startswith('Reach XY='):
                n = int(linhas[i].split('=')[1])
                i += 1
                pulados = 0
                while i < len(linhas) and pulados < 2 * n:
                    pulados += len([1 for c in range(0, len(linhas[i]), 16)
                                    if linhas[i][c:c + 16].strip()])
                    i += 1
            novo_txt += bloco_xy(up_xy)
            continue
        if i == i_split:
            novo_txt.append(f'River Reach={alvo_rio:<16s},{novo_reach:<16s}')
            novo_txt += bloco_xy(dn_xy)
            novo_txt.append('')
        novo_txt.append(l)
        i += 1
    linhas = novo_txt

    # juncao existente a jusante que referencia o alvo_reach: atualiza
    t2 = '\n'.join(linhas)
    padrao_up = f'Up River,Reach={alvo_rio:<16s},{alvo_reach:<16s}'
    novo_up = f'Up River,Reach={alvo_rio:<16s},{novo_reach:<16s}'
    n_jup = t2.count(padrao_up)
    t2 = t2.replace(padrao_up, novo_up)

    # rio novo: bloco inteiro do g01 avulso (do River Reach ao fim)
    tn = open(g_novo, encoding='latin-1', errors='replace').read() \
        .replace('\r\n', '\n').replace('\r', '\n')
    m = re.search(r'^River Reach=.*$', tn, flags=re.M)
    bloco_novo = tn[m.start():].strip('\n')
    # renomeia para o nome pedido
    bloco_novo = re.sub(r'^River Reach=[^,]+,',
                        f'River Reach={nome_rio:<16s},', bloco_novo,
                        count=1, flags=re.M)
    # '#Mann= N ,-1,0' liga o modo "n horizontal variado" e o RAS exige
    # outra estrutura ("needs n on first station"); as irmas usam 0
    bloco_novo = re.sub(r'^#Mann= *(\d+) *, *-1 *, *0',
                        r'#Mann= \1 , 0 , 0 ', bloco_novo, flags=re.M)
    # e o gerador grava os valores em campos de 12 colunas -- o RAS le
    # em 8 e acha a tabela desalinhada ("n value not set"): reescreve
    # cada bloco de valores no formato canonico de 8 colunas
    ls_n2 = bloco_novo.split('\n')
    saida_n = []
    i2 = 0
    grade = np.array([0.0])
    while i2 < len(ls_n2):
        l2 = ls_n2[i2]
        if l2.startswith('#Sta/Elev='):
            # guarda a grade de estacoes da secao corrente, para ancorar
            # Bank Sta e quebras de Manning (exigencia do RAS)
            cnt3 = int(l2.split('=')[1])
            vals3 = []
            j3 = i2 + 1
            while j3 < len(ls_n2) and len(vals3) < 2 * cnt3:
                x3 = ls_n2[j3]
                if not x3.strip() or x3.lstrip()[0].isalpha():
                    break
                vals3 += [float(x3[c:c + 8])
                          for c in range(0, len(x3), 8)
                          if x3[c:c + 8].strip()]
                j3 += 1
            grade = np.array(vals3[0::2]) if vals3 else grade

        def ancora(v):
            return float(grade[int(np.argmin(np.abs(grade - v)))])
        if l2.startswith('Bank Sta='):
            a3, b3 = (float(x) for x in l2.split('=')[1].split(','))
            l2 = 'Bank Sta=%.2f,%.2f' % (ancora(a3), ancora(b3))
        saida_n.append(l2)
        if l2.startswith('#Mann='):
            cnt2 = int(l2.split('=')[1].split(',')[0])
            vals2 = []
            i2 += 1
            while i2 < len(ls_n2) and len(vals2) < 3 * cnt2:
                x2 = ls_n2[i2]
                if not x2.strip() or x2.lstrip()[0].isalpha():
                    break
                vals2 += [float(v) for v in x2.split()]
                i2 += 1
            # quebras do Manning ancoradas na grade (posicoes 0, 3, 6...)
            for j2 in range(0, len(vals2), 3):
                vals2[j2] = ancora(vals2[j2])
            lin2 = ''
            for j2, v2 in enumerate(vals2):
                lin2 += '%8.3f' % v2
                if (j2 + 1) % 10 == 0:
                    saida_n.append(lin2)
                    lin2 = ''
            if lin2:
                saida_n.append(lin2)
            continue
        i2 += 1
    bloco_novo = '\n'.join(saida_n)
    # anexa ao fim
    t2 = t2.rstrip('\n') + '\n\n' + bloco_novo + '\n'

    # juncao nova, no padrao das irmas, apos a ultima juncao existente
    nome_j = f'Foz_{nome_rio}'[:15]
    jx, jy = p_corte.x, p_corte.y
    juncao = '\n'.join([
        f'Junct Name={nome_j:<16s}',
        'Junct Desc=Confluencia, 0 , 0 , 0 ,0',
        f'Junct X Y & Text X Y={jx:.2f},{jy:.2f},{jx+800:.2f},{jy+800:.2f}',
        f'Up River,Reach={alvo_rio:<16s},{alvo_reach:<16s}',
        f'Up River,Reach={nome_rio:<16s},R1              ',
        f'Dn River,Reach={alvo_rio:<16s},{novo_reach:<16s}',
        'Junc L&A=150.00,0',
        'Junc L&A=150.00,0', '', ''])
    ultima_j = [mm.end() for mm in
                re.finditer(r'^Junc L&A=.*$\n\n', t2, flags=re.M)]
    if ultima_j:
        pos = ultima_j[len([1 for _ in
                            re.finditer(r'^Junct Name=', t2[:max(ultima_j)],
                                        flags=re.M)]) - 1] \
            if False else max(ultima_j)
        t2 = t2[:pos] + juncao + t2[pos:]
    else:
        raise SystemExit('nao achei juncoes existentes para ancorar')
    escrever(saida, t2)

    # ---------------- u01
    u = open(u01_in, encoding='latin-1', errors='replace').read() \
        .replace('\r\n', '\n').replace('\r', '\n')
    # laterais do alvo_reach que atravessam o corte: dividir
    blocos = re.split(r'(?=^Boundary Location=)', u, flags=re.M)
    for k, b in enumerate(blocos):
        mm = re.match(r'Boundary Location=([^,]+),([^,]+),([\d.]+)\s*,'
                      r'([\d.]+)\s*,', b)
        if not mm or mm.group(1).strip() != alvo_rio \
                or mm.group(2).strip() != alvo_reach:
            continue
        r_ini, r_fim2 = float(mm.group(3)), float(mm.group(4))
        if 'Uniform Lateral Inflow' in b and r_fim2 < rs_corte_dn < r_ini:
            h = re.search(r'Uniform Lateral Inflow Hydrograph=\s*(\d+)', b)
            vals, corte_i = [], None
            ls_b = b[h.end():].split('\n')
            for li, l in enumerate(ls_b[1:], 1):
                if not l.strip() or l[:1].isalpha():
                    corte_i = li
                    break
                vals += [float(l[c:c + 8]) for c in range(0, len(l), 8)
                         if l[c:c + 8].strip()]
            resto = '\n'.join(ls_b[corte_i:]) if corte_i else ''
            fr_up = (r_ini - rs_corte_up) / max(r_ini - r_fim2, 1e-9)
            # o RAS proibe lateral terminando na ULTIMA secao do reach
            # ou comecando na PRIMEIRA: encosta uma secao para dentro
            fim_up = sorted(acima)[1] if len(acima) > 1 else rs_corte_up
            ini_dn = sorted(abaixo, reverse=True)[1] \
                if len(abaixo) > 1 else rs_corte_dn
            cab = b.split('\n')[0]
            meio = b[len(cab) + 1:h.start()]

            def novo_bloco(rioX, reachX, ra, rb, frac):
                sN = [v * frac for v in vals]
                cabN = (f'Boundary Location={rioX:<16s},{reachX:<16s},'
                        f'{ra:<8s},{rb:<8s},                ,'
                        '                ')
                return '\n'.join([cabN] + meio.strip('\n').split('\n')
                                 + ['Uniform Lateral Inflow Hydrograph= '
                                    f'{len(sN)} '] + fmt8(sN)
                                 + [resto.rstrip('\n')]) + '\n\n'
            ra_up = ('%.2f' % r_ini).rstrip('0').rstrip('.')
            rb_up = ('%.2f' % fim_up).rstrip('0').rstrip('.')
            ra_dn = ('%.2f' % ini_dn).rstrip('0').rstrip('.')
            rb_dn = ('%.2f' % r_fim2).rstrip('0').rstrip('.')
            blocos[k] = (novo_bloco(alvo_rio, alvo_reach, ra_up, rb_up,
                                    fr_up)
                         + novo_bloco(alvo_rio, novo_reach, ra_dn, rb_dn,
                                      1 - fr_up))
            print(f'   lateral {r_ini:.0f}->{r_fim2:.0f} dividida no corte '
                  f'({fr_up:.2f}/{1-fr_up:.2f})')
    u = ''.join(blocos)

    # cabeceira do rio novo
    Sn = ler_secoes(g_novo)
    rs_topo = max(d['rs'] for d in Sn if d['tipo'] == '1')
    rs_topo_txt = ('%.2f' % rs_topo).rstrip('0').rstrip('.')
    if serie_arq:
        s = serie_diaria(serie_arq) * escala
        bloco_q = '\n'.join(
            [f'Boundary Location={nome_rio:<16s},R1              ,'
             f'{rs_topo_txt:<8s},        ,                ,'
             '                ',
             'Interval=1DAY',
             f'Flow Hydrograph= {len(s)} '] + fmt8(s)
            + ['DSS Path=', 'Use DSS=False', 'Use Fixed Start Time=False',
               'Fixed Start Date/Time=,', 'Flow Hydrograph Slope= 0.001 ',
               '', ''])
        u = u.rstrip('\n') + '\n\n' + bloco_q
        q0 = float(s[0])
    else:
        q0 = 10.0
    # iniciais: rio novo e reach de jusante
    ult_ini = list(re.finditer(r'^Initial Flow Loc=.*$', u, flags=re.M))
    ins = (f'\nInitial Flow Loc={nome_rio:<16s},R1              ,'
           f'{rs_topo_txt:<8s},{q0:g}')
    # jusante herda o Q inicial do alvo
    m_alvo = re.search(rf'^Initial Flow Loc={re.escape(alvo_rio)}\s*,'
                       rf'{re.escape(alvo_reach)}\s*,[^,]+,([\d.]+)',
                       u, flags=re.M)
    q_alvo = float(m_alvo.group(1)) if m_alvo else 50.0
    rs_dn_txt = ('%.2f' % rs_corte_dn).rstrip('0').rstrip('.')
    ins += (f'\nInitial Flow Loc={alvo_rio:<16s},{novo_reach:<16s},'
            f'{rs_dn_txt:<8s},{q_alvo + q0:g}')
    if ult_ini:
        fim_ini = ult_ini[-1].end()
        u = u[:fim_ini] + ins + u[fim_ini:]
    escrever(u01_out, u)

    # ---------------- conferencia
    print('\nCONFERENCIA (relendo o gravado)')
    B = ler_secoes(saida)
    n_rede = len([d for d in S if d['tipo'] == '1'])
    n_novo = len([d for d in Sn if d['tipo'] == '1'])
    n_fim = len([d for d in B if d['tipo'] == '1'])
    print(f'   secoes: {n_rede} + {n_novo} = {n_fim}  '
          f'({"OK" if n_fim == n_rede + n_novo else "ERRO"})')
    tt = open(saida, encoding='latin-1', errors='replace').read()
    print(f'   juncao {nome_j}: {"presente" if nome_j in tt else "FALTA"}')
    print(f'   reaches de {alvo_rio}: '
          f'{sorted({d["reach"] for d in B if d["rio"] == alvo_rio})}')
    uu = open(u01_out, encoding='latin-1', errors='replace').read()
    print(f'   u01: Initial RS fantasma={uu.count("Initial RS=")} '
          f'(tem de ser 0); Flow Hydrograph do {nome_rio}: '
          f'{"presente" if nome_rio in uu else "FALTA"}')
    print(f'   juncoes a jusante atualizadas: {n_jup}')


if __name__ == '__main__':
    main(sys.argv[1:])
