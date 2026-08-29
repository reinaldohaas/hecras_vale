# -*- coding: utf-8 -*-
"""Separa as variantes do taha_ai em PROJETOS HEC-RAS proprios.

    python scripts/separar_projetos.py

  raiz (taha_ai.prj)          so o consagrado: Copernicus g01/u01/p01
  projeto_sigsc/              taha_ai_sigsc.prj = SIG-SC calibrado
                              (g02/p02 de hoje viram g01/p01 de la)
  projeto_experimentos/       taha_ai_lab.prj = rede plena (p01),
                              Krauel (p02), Luis Alves (p03)

Cada arquivo copiado tem as referencias reescritas (Geom/Flow File) em
CRLF; resultados .pNN.hdf acompanham. A raiz e PODADA (g02-g05/u04-u05/
p02-p05 saem do .prj e do disco -- ficam nos projetos novos). Nada dos
backups .estado_* / .antes_* muda. [[nunca-sobrescrever-o-projeto]]
"""
import os
import re
import shutil
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def ler(a):
    return open(a, encoding='latin-1', errors='replace').read() \
        .replace('\r\n', '\n').replace('\r', '\n')


def gravar(a, t):
    open(a, 'w', encoding='latin-1', newline='\r\n').write(t)


def prj_novo(titulo):
    return (f'Proj Title={titulo}\nCurrent Plan=p01\n'
            'Default Exp/Contr=0.3,0.1\nSI Units\n')


def plano(p_orig, geom, flow, titulo, ident):
    t = ler(p_orig)
    t = re.sub(r'Geom File=g\d+', f'Geom File={geom}', t, count=1)
    t = re.sub(r'Flow File=u\d+', f'Flow File={flow}', t, count=1)
    t = re.sub(r'Plan Title=.*', f'Plan Title={titulo}', t, count=1)
    t = re.sub(r'Short Identifier=.*', f'Short Identifier={ident}', t,
               count=1)
    return t


def main(argv):
    os.chdir(RAIZ)

    # ------------------------------------------------ projeto SIG-SC
    d = 'projeto_sigsc'
    os.makedirs(d, exist_ok=True)
    b = os.path.join(d, 'taha_ai_sigsc')
    shutil.copy2('taha_ai.g02', b + '.g01')
    shutil.copy2('taha_ai.u01', b + '.u01')
    gravar(b + '.p01', plano('taha_ai.p02', 'g01', 'u01',
                             '1983 relevo SIG-SC 1 m', '1983_sigsc'))
    gravar(b + '.prj', prj_novo('taha_ai_sigsc')
           + 'Geom File=g01\nUnsteady File=u01\nPlan File=p01\n'
           + 'Y Axis Title=Elevation\n'
           + 'X Axis Title(PF)=Main Channel Distance\n'
           + 'X Axis Title(XS)=Station\n'
           + 'BEGIN DESCRIPTION:\n1983 observado, relevo SIG-SC 1 m, '
             'Barragem Sul\nEND DESCRIPTION:\n')
    for orig, dest in [('taha_ai.p02.hdf', b + '.p01.hdf'),
                       ('taha_ai.g02.hdf', b + '.g01.hdf')]:
        if os.path.exists(orig):
            shutil.copy2(orig, dest)
    print(f'{d}/: g01/u01/p01 (+resultados)')

    # ------------------------------------------- projeto experimentos
    d = 'projeto_experimentos'
    os.makedirs(d, exist_ok=True)
    b = os.path.join(d, 'taha_ai_lab')
    mapa = [
        # (geom origem, u origem, p origem, geom, u, p, titulo, ident)
        ('taha_ai.g03', 'taha_ai.u01', 'taha_ai.p03', 'g01', 'u01',
         'p01', 'rede plena SIG-SC (lab)', 'lab_rede_plena'),
        ('taha_ai.g04', 'taha_ai.u04', 'taha_ai.p04', 'g02', 'u02',
         'p02', 'Krauel (lab)', 'lab_krauel'),
        ('taha_ai.g05', 'taha_ai.u05', 'taha_ai.p05', 'g03', 'u03',
         'p03', 'Luis Alves (lab)', 'lab_luis_alves'),
    ]
    linhas_prj = [prj_novo('taha_ai_lab').rstrip('\n')]
    for g_o, u_o, p_o, g_n, u_n, p_n, tit, idt in mapa:
        if not os.path.exists(g_o):
            print(f'   {g_o} nao existe -- pulado')
            continue
        shutil.copy2(g_o, f'{b}.{g_n}')
        shutil.copy2(u_o, f'{b}.{u_n}')
        gravar(f'{b}.{p_n}', plano(p_o, g_n, u_n, tit, idt))
        hdf = p_o + '.hdf'
        if os.path.exists(hdf):
            shutil.copy2(hdf, f'{b}.{p_n}.hdf')
        linhas_prj += [f'Geom File={g_n}', f'Unsteady File={u_n}',
                       f'Plan File={p_n}']
    linhas_prj += ['Y Axis Title=Elevation',
                   'X Axis Title(PF)=Main Channel Distance',
                   'X Axis Title(XS)=Station', 'BEGIN DESCRIPTION:',
                   'experimentos: rede plena, Krauel, Luis Alves',
                   'END DESCRIPTION:']
    gravar(b + '.prj', '\n'.join(linhas_prj) + '\n')
    print(f'{d}/: 3 planos de laboratorio')

    # ------------------------------------------------- poda da raiz
    t = ler('taha_ai.prj')
    for n in range(2, 6):
        t = t.replace(f'Geom File=g0{n}\n', '')
        t = t.replace(f'Plan File=p0{n}\n', '')
        t = t.replace(f'Unsteady File=u0{n}\n', '')
    gravar('taha_ai.prj', t)
    removidos = []
    for arq in ['taha_ai.g02', 'taha_ai.g03', 'taha_ai.g04',
                'taha_ai.g05', 'taha_ai.u04', 'taha_ai.u05',
                'taha_ai.p02', 'taha_ai.p03', 'taha_ai.p04',
                'taha_ai.p05', 'taha_ai.g02.hdf', 'taha_ai.g04.hdf',
                'taha_ai.g05.hdf', 'taha_ai.p02.hdf', 'taha_ai.p03.hdf',
                'taha_ai.p04.hdf', 'taha_ai.p05.hdf']:
        if os.path.exists(arq):
            os.remove(arq)
            removidos.append(arq)
    print(f'raiz podada: {len(removidos)} arquivos movidos aos projetos')

    print('\nCONFERENCIA')
    for prj in ['taha_ai.prj', 'projeto_sigsc/taha_ai_sigsc.prj',
                'projeto_experimentos/taha_ai_lab.prj']:
        t = ler(prj)
        gs = re.findall(r'Geom File=(g\d+)', t)
        ps = re.findall(r'Plan File=(p\d+)', t)
        pasta = os.path.dirname(prj) or '.'
        base = os.path.basename(prj)[:-4]
        falta = [x for x in gs + ps
                 if not os.path.exists(os.path.join(pasta,
                                                    f'{base}.{x}'))]
        print(f'   {prj}: geoms {gs} planos {ps}'
              f'{"  FALTAM: " + str(falta) if falta else "  arquivos OK"}')


if __name__ == '__main__':
    main(sys.argv[1:])
