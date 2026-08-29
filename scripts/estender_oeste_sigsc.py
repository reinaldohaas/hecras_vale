# -*- coding: utf-8 -*-
"""Como estender_oeste.py, mas no projeto_sigsc (rede 100% SIG-SC,
Mirim ja completo) e terreno 10 m -- sem a costura de resolucao contra
o Copernico que provavelmente causou o "Overflow" de hoje no raiz.

    python scripts/estender_oeste_sigsc.py            g96 (so extensao)
    python scripts/estender_oeste_sigsc.py --tudo     g03/u03/p03 finais

NUNCA TOCA g01/g02/u01/u02/p01/p02 do projeto_sigsc.
"""
import os
import re
import subprocess
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(DIR)
sys.path.insert(0, DIR)

import estender_oeste as base                          # noqa: E402

PASTA = os.path.join(RAIZ, 'projeto_sigsc')
base.RAIZ = PASTA
base.MOSAICO = os.path.join(RAIZ, 'taha_ai_novo', 'Terrain',
                            'taha_ai_corredor_10m_v2.tif')
base.G01_BASE = 'taha_ai_sigsc.g02'      # a rede com Mirim completo


def caminho(nome):
    return os.path.join(PASTA, nome)


def montar_g96():
    ext, s_dam = base.eixo_extensao()      # cwd continua em RAIZ
    secoes = base.cortar_secoes(ext, s_dam)
    txt = open(caminho('taha_ai_sigsc.g02'), encoding='latin-1',
              newline='').read()
    pts_ext = list(ext.coords)[1:][::-1]
    m = re.search(r'(River Reach=Itajai_Oeste {4},R1 +\r?\n'
                  r'Reach XY= *)(\d+)( *\r?\n)((?:.+\r?\n)+?)'
                  r'(Rch Text X Y=)', txt)
    assert m, 'Reach XY do Oeste R1 nao achado no sigsc'
    velhos = m.group(4)
    corpo, lin = [], ''
    for i, (x, y) in enumerate(pts_ext):
        lin += f'{x:16.4f}{y:16.4f}'
        if (i + 1) % 2 == 0:
            corpo.append(lin)
            lin = ''
    if lin:
        corpo.append(lin)
    n_total = int(m.group(2)) + len(pts_ext)
    novo_xy = (m.group(1) + str(n_total) + m.group(3)
               + base.CRLF.join(corpo) + base.CRLF + velhos
               + m.group(5))
    txt = txt[:m.start()] + novo_xy + txt[m.end():]
    alvo = 'Type RM Length L Ch R = 1 ,56723.30'
    k = txt.find(alvo)
    assert k > 0
    blocos = ''
    ordenadas = sorted(secoes, key=lambda s: -s['rs'])
    serie = base.largura_alvo_serie(
        os.path.join(RAIZ, 'doc', 'larguras_sigsc', 'Rio_do_Campo.csv'))
    ordenadas = base.aplicar_largura_medida(ordenadas, serie)
    if serie:
        alvos = [s['largura_alvo'] for s in ordenadas]
        print(f'largura medida aplicada: {min(alvos):.0f}-'
              f'{max(alvos):.0f} m (mediana {np.median(alvos):.0f} m)')
    ordenadas = base.suavizar_talvegue(ordenadas)
    for i, sec in enumerate(ordenadas):
        comprimento = (ordenadas[i]['rs']
                       - (ordenadas[i + 1]['rs']
                          if i + 1 < len(ordenadas) else base.RS_TOPO))
        blocos += base.bloco_secao(sec, comprimento)
    txt = txt[:k] + blocos + txt[k:]
    with open(caminho('taha_ai_sigsc.g96'), 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(txt)
    print(f'g96 escrito ({len(secoes)} secoes novas, terreno 10m)')
    return ext, s_dam, secoes


def main(argv):
    ext, s_dam, secoes = montar_g96()
    if '--tudo' not in argv:
        print('parou no g96 (rode com --tudo p/ barragem+u03+p03)')
        return
    crista = base.medir_barragem(ext, s_dam)
    rs_dam = base.RS_TOPO + s_dam
    py = sys.executable
    r = subprocess.run(
        [py, os.path.join(RAIZ, 'scripts', 'construir_barragem.py'),
         'taha_ai_sigsc.g96', '--saida', 'g03', '--rio',
         'Itajai_Oeste', '--reach', 'R1', '--rs', f'{rs_dam:.2f}',
         '--crista', f'{crista - 3.0:.1f}', '--topo', f'{crista:.1f}',
         '--larg-vertedouro', '120', '--fenda', '1.1',
         '--nome', 'Barragem Oeste (Taio)'],
        capture_output=True, text=True, cwd=PASTA)
    print(r.stdout[-600:])
    if r.returncode != 0:
        print(r.stderr[-800:])
        raise SystemExit('construir_barragem falhou')
    os.remove(caminho('taha_ai_sigsc.g96'))
    topo_novo = max(s['rs'] for s in secoes)
    u = open(caminho('taha_ai_sigsc.u02'), encoding='latin-1',
            newline='').read()
    u = u.replace('Boundary Location=Itajai_Oeste    ,'
                  'R1              ,56723.30',
                  f'Boundary Location=Itajai_Oeste    ,'
                  f'R1              ,{topo_novo:8.2f}')
    u = u.replace('Initial Flow Loc=Itajai_Oeste,R1,56723.3,',
                  f'Initial Flow Loc=Itajai_Oeste,R1,'
                  f'{topo_novo:.2f},')
    # escala pela area (854 km2 barragem / 1570 km2 Taio)
    linhas = u.split('\r\n')
    i = next(k for k, l in enumerate(linhas)
             if l.startswith('Boundary Location=Itajai_Oeste')
             and f'{topo_novo:8.2f}' in l)
    j = i + 2
    n = int(linhas[j].split('=')[1].strip())
    nlin = (n + 9) // 10
    vals = []
    for k in range(j + 1, j + 1 + nlin):
        for c in range(0, len(linhas[k]), 8):
            campo = linhas[k][c:c + 8]
            if campo.strip():
                vals.append(float(campo))
    fracao = 854.0 / 1570.0
    novos = [v * fracao for v in vals]
    linhas[j + 1:j + 1 + nlin] = [
        ''.join(f'{v:8.2f}' for v in novos[c:c + 10])
        for c in range(0, len(novos), 10)]
    u2 = '\r\n'.join(linhas)
    with open(caminho('taha_ai_sigsc.u03'), 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(u2)
    p = open(caminho('taha_ai_sigsc.p02'), encoding='latin-1',
            newline='').read()
    p = p.replace('Geom File=g02', 'Geom File=g03')
    p = p.replace('Flow File=u02', 'Flow File=u03')
    p = re.sub(r'Plan Title=[^\r\n]*',
               'Plan Title=1983 + Barragem Oeste (SIG-SC 10m)',
               p, count=1)
    p = re.sub(r'Short Identifier=[^\r\n]*',
               'Short Identifier=1983_bOeste_sc', p, count=1)
    with open(caminho('taha_ai_sigsc.p03'), 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(p)
    prj = open(caminho('taha_ai_sigsc.prj'), encoding='latin-1',
              newline='').read()
    for chave, arq in [('Geom File=g03', 'Geom File=g02'),
                       ('Unsteady File=u03', 'Unsteady File=u02'),
                       ('Plan File=p03', 'Plan File=p02')]:
        if chave not in prj:
            prj = prj.replace(arq + '\r\n', arq + '\r\n' + chave
                              + '\r\n', 1)
    with open(caminho('taha_ai_sigsc.prj'), 'w', encoding='latin-1',
              newline='') as fh:
        fh.write(prj)
    print(f'g03/u03/p03 prontos no projeto_sigsc; topo RS {topo_novo:.2f}, '
          f'barragem RS {rs_dam:.2f}, crista {crista - 3:.1f}/'
          f'{crista:.1f} m, pico escalado {max(novos):.1f} m3/s')


if __name__ == '__main__':
    main(sys.argv[1:])
