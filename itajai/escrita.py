# -*- coding: utf-8 -*-
"""
Escrita dos arquivos do HEC-RAS.

Formatos validados contra os projetos-exemplo oficiais e, agora, conferidos de
volta pelo ras-commander -- ele le o .g01 que escrevemos e devolve as secoes,
station-elevation, bank stations e as tres zonas de Manning corretas.

O que so se descobre apanhando, e que esta aqui:

  - series e #Sta/Elev em COLUNAS DE 8 CARACTERES, 10 por linha. Escrever numa
    linha unica separada por espacos faz o RAS ler "dados faltando";
  - Boundary Location tem SEIS campos com padding (16/16/8 + 3 vazios);
  - Bank Sta precisa casar com um sta da tabela na MESMA precisao (.2f);
  - Junc L&A e o caminho ATRAVES da juncao, da ultima secao de montante a
    primeira de jusante. Gravar 500 m fixo onde a geometria da 75 desequilibra
    a continuidade exatamente na secao que o solver reporta;
  - Viewing Rectangle e a extensao real. Com o placeholder "0,1,1,0" o HDF sai
    com Extents=[0,1,0,1] e o RAS Mapper abre VAZIO, procurando a bacia inteira
    dentro de um metro quadrado na origem;
  - o terreno no .rasmap so pode ser o .hdf. Declarar o GeoTIFF faz o HEC-RAS
    tentar abri-lo como HDF5 e despejar HDF5-DIAG a cada execucao.
"""
import datetime
import os

import numpy as np

from .config import EPSG, WKT
from .perfil import cota_talvegue

MESES = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
         "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
EDITADO = "Node Last Edited Time= Aug\\15\\2026 00:00:00"


def p16(s):
    return f"{str(s)[:16]:<16}"


def f8(v):
    """Um valor em coluna de 8 caracteres, sem estourar."""
    for fmt in ("{:8.2f}", "{:8.1f}", "{:8.0f}"):
        t = fmt.format(v)
        if len(t) <= 8:
            return t
    return f"{v:8.3g}"[:8]


def serie8(vals):
    """Serie em colunas de 8, 10 por linha."""
    v = list(vals)
    return [" ".join([]) or "".join(f8(x) for x in v[i:i + 10])
            for i in range(0, len(v), 10)]


def data_ras(dt):
    return f"{dt.day:02d}{MESES[dt.month-1]}{dt.year},{dt.hour:02d}00"


def contorno(rio, reach, rs):
    return (f"Boundary Location={p16(rio)},{p16(reach)},{str(rs)[:8]:<8}"
            f",        ,                ,                ")


# ------------------------------------------------------------------ GEOMETRIA
def _secao_manual(rio, reach, d, dx):
    """Formatacao propria da secao. Reserva, se o ras-commander faltar."""
    c = d["cut"]
    linhas = [f"Type RM Length L Ch R = 1 ,{d['rs']:<8.2f},{dx},{dx},{dx}",
              "XS GIS Cut Line=2",
              f"{c[0]:16.4f}{c[1]:16.4f}{c[2]:16.4f}{c[3]:16.4f}",
              EDITADO,
              f"#Sta/Elev= {len(d['sta'])} "]
    linhas += serie8([v for p in zip(d["sta"], d["z"]) for v in p])
    linhas += ["#Mann= 3 ,-1,0",
               "".join(f"{v:>8}" for v in
                       [f"{d['sta'][0]:.2f}", f"{d['n_planicie']:.3f}", "0",
                        f"{d['lb']:.2f}", f"{d['n']:.3f}", "0",
                        f"{d['rb']:.2f}", f"{d['n_planicie']:.3f}", "0"]),
               f"Bank Sta={d['lb']:.2f},{d['rb']:.2f}",
               "XS Rating Curve= 0 ,0",
               "Exp/Cntr=0.3,0.1", ""]
    return linhas


def secao_texto(rio, reach, d, dx):
    """Bloco da secao, montado pelo ras-commander quando disponivel.

    Ganho concreto sobre a formatacao propria: o builder INSERE pontos nas
    estacas das margens, garantindo que Bank Sta coincida com um sta da tabela.
    Fazer isso a mao exige casar a precisao (.2f) dos dois lados, e foi fonte
    de erro. Ele tambem resolve reach lengths e o bloco de Manning no formato
    que o HEC-RAS espera.

    A formatacao propria fica como reserva: sem o ras-commander instalado o
    modelo continua sendo gerado.
    """
    try:
        import pandas as pd
        from ras_commander.geom import GeomCrossSection as G
    except ImportError:
        return _secao_manual(rio, reach, d, dx)
    c = d["cut"]
    try:
        r = G.build_cross_section(
            river=rio, reach=reach, rs=f"{d['rs']:.2f}",
            station_elevation=pd.DataFrame({"Station": d["sta"],
                                            "Elevation": d["z"]}),
            cut_line=[(c[0], c[1]), (c[2], c[3])],
            bank_left=float(d["lb"]), bank_right=float(d["rb"]),
            n_lob=float(d["n_planicie"]), n_channel=float(d["n"]),
            n_rob=float(d["n_planicie"]),
            length_left=float(dx), length_channel=float(dx),
            length_right=float(dx))
        linhas = r.text.splitlines()
        # O builder grava Manning com DUAS casas: 0,035 vira 0,04 (14% a mais
        # de rugosidade em todo o modelo) e os valores de Jarrett -- 0,052,
        # 0,066, 0,070 -- colapsam em tres niveis. Substitui-se so essa linha,
        # preservando tudo o mais que ele resolve.
        for k, l in enumerate(linhas):
            if l.startswith("#Mann=") and k + 1 < len(linhas):
                linhas[k + 1] = "".join(
                    f"{v:>8}" for v in
                    [f"{d['sta'][0]:.2f}", f"{d['n_planicie']:.3f}", "0",
                     f"{d['lb']:.2f}", f"{d['n']:.3f}", "0",
                     f"{d['rb']:.2f}", f"{d['n_planicie']:.3f}", "0"])
                break
        return linhas + ["XS Rating Curve= 0 ,0", "Exp/Cntr=0.3,0.1", ""]
    except Exception:
        return _secao_manual(rio, reach, d, dx)


def geometria(projeto, trechos, juncoes, titulo=None):
    xs_all, ys_all = [], []
    for t in trechos:
        for x, y in t["linha"].coords:
            xs_all.append(x); ys_all.append(y)
        for d in t["xs"]:
            c = d["cut"]
            xs_all += [c[0], c[2]]
            ys_all += [c[1], c[3]]
    folga = 0.02 * max(max(xs_all) - min(xs_all), max(ys_all) - min(ys_all), 1.0)
    x0, x1 = min(xs_all) - folga, max(xs_all) + folga
    y0, y1 = min(ys_all) - folga, max(ys_all) + folga

    g = [f"Geom Title={titulo or projeto}",
         "Program Version=7.01",
         f"Viewing Rectangle= {x0:.6f} , {x1:.6f} , {y1:.6f} , {y0:.6f} ",
         f"Spatial Reference System={WKT}", ""]

    for j in juncoes:
        g += [f"Junct Name={p16(j['nome'])}",
              "Junct Desc=Confluencia, 0 , 0 , 0 ,0",
              f"Junct X Y & Text X Y={j['x']:.2f},{j['y']:.2f},"
              f"{j['x']+800:.2f},{j['y']+800:.2f}"]
        for rio, rea in j["up"]:
            g.append(f"Up River,Reach={p16(rio)},{p16(rea)}")
        g.append(f"Dn River,Reach={p16(j['dn'][0])},{p16(j['dn'][1])}")
        for d in j["dists"]:
            g.append(f"Junc L&A={d:.2f},0")
        g.append("")

    for t in trechos:
        pts = list(t["linha"].coords)
        g += [f"River Reach={p16(t['rio'])},{p16(t['reach'])}",
              f"Reach XY= {len(pts)} "]
        par = [v for p in pts for v in p]
        g += ["".join(f"{v:16.4f}" for v in par[i:i + 4])
              for i in range(0, len(par), 4)]
        g += ["Rch Text X Y=0,0,0,0", ""]

        for i, d in enumerate(t["xs"]):
            prox = t["xs"][i + 1] if i + 1 < len(t["xs"]) else None
            dx = round(d["rs"] - prox["rs"], 2) if prox else 0.0
            g += secao_texto(t["rio"], t["reach"], d, dx)

    with open(f"{projeto}.g01", "w", encoding="ascii", errors="replace") as f:
        f.write("\n".join(g) + "\n")
    return f"{projeto}.g01"


# --------------------------------------------------------------------- FLUXO
def fluxo(projeto, trechos, cabeceiras, saida, mare, n_horas):
    """Contornos e condicao inicial.

    'Initial RS' vai para TODOS os trechos, nao so para as cabeceiras. Um
    trecho sem vazao inicial parte indefinido, e o HEC-RAS acusa isso como
    erro de volume gigantesco ja no primeiro passo -- com 12 rios eram dez
    trechos nessa situacao (Acu R1-R5, Oeste R2-R4, Norte R2, Benedito R2) e o
    erro de balanco subiu de 5.372% para 73.310% conforme a rede cresceu.

    A vazao de cada trecho e a ACUMULADA: o que ja entrou pelas cabeceiras a
    montante dele. Comecar um trecho de jusante com a vazao de uma cabeceira
    isolada cria um degrau na juncao, e o solver nao converge.
    """
    u = [f"Flow Title=Cheia_{projeto}", "Program Version=7.01", "Use Restart= 0 "]
    for t in trechos:
        u += [f"Initial RS={p16(t['rio'])},{p16(t['reach'])},"
              f"{t['xs'][0]['rs']:<8.0f},{t.get('q_base', 20.0):.0f}"]
    u.append("")
    for t in cabeceiras:
        u += [contorno(t["rio"], t["reach"], f"{t['xs'][0]['rs']:.2f}"),
              "Interval=1HOUR",
              f"Flow Hydrograph= {len(t['serie'])} "]
        u += serie8(t["serie"])
        u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=False",
              "Fixed Start Date/Time=,", "Is Critical Boundary=False",
              "Critical Boundary Flow=", ""]
    u += [contorno(saida["rio"], saida["reach"], f"{saida['xs'][-1]['rs']:.2f}"),
          "Interval=1HOUR",
          f"Stage Hydrograph= {len(mare)} "]
    u += serie8(mare)
    u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=False",
          "Fixed Start Date/Time=,", ""]
    with open(f"{projeto}.u01", "w", encoding="ascii", errors="replace") as f:
        f.write("\n".join(u) + "\n")
    return f"{projeto}.u01"


# ---------------------------------------------------------- PLANO E PROJETO
def plano(projeto, inicio, n_horas):
    fim = inicio + datetime.timedelta(hours=n_horas - 1)
    p = ["Plan Title=Rede_do_Relevo", "Program Version=7.01",
         "Short Identifier=TAJAI", "Geom File=g01", "Flow File=u01",
         f"Simulation Date={data_ras(inicio)},{data_ras(fim)}",
         "Mixed Flow Regime",
         # 'Mixed Flow Regime' e o LPI sao coisas distintas: o primeiro so
         # permite regime misto; o amortecimento dos termos de inercia perto do
         # critico e o Froude Reduction, que vem DESLIGADO de fabrica.
         "UNET Froude Reduction=True",
         "UNET Froude Limit= 0.8 ",
         "UNET Froude Power= 4 ",
         "UNET MxIter= 40 ",
         "UNET Max Iter WO Improvement= 20 ",
         "UNET Theta= 1 ", "UNET Theta Warmup= 1 ",
         "UNET ZTol= 0.01 ", "UNET ZSATol= 0.01 ",
         "UNET DZMax Abort= 30 ",
         "UNET MaxInSteps= 200 ", "UNET DtIC= 0 ",
         "Flow Smoothing Iterations=10",
         "Unsteady Friction Slope Method= 2 ",
         "UNET 1D Methodology=Finite Difference",
         "Computation Interval=1MIN", "Output Interval=1HOUR",
         "Instantaneous Interval=1HOUR", "Mapping Interval=1HOUR",
         "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1",
         "Run RASMapper=-1"]
    open(f"{projeto}.p01", "w", encoding="ascii").write("\n".join(p) + "\n")
    open(f"{projeto}.prj", "w", encoding="ascii").write("\n".join([
        f"Proj Title={projeto}", "Current Plan=p01",
        "Default Exp/Contr=0.3,0.1", "SI Units", "Geom File=g01",
        "Unsteady File=u01", "Plan File=p01", "Y Axis Title=Elevation",
        "X Axis Title(PR)=Distance", "X Axis Title(CS)=Station",
        f"RASMap Filename={projeto}.rasmap"]) + "\n")
    open(f"{projeto}.projection", "w", encoding="utf-8").write(WKT)
    return f"{projeto}.prj"


def rasmap(projeto, terreno_hdf=None):
    """Configuracao do RAS Mapper.

    A GEOMETRIA nao e declarada: o RAS Mapper ja a descobre pelo projeto, com a
    arvore completa e o nome vindo do Geom Title. Declarar tambem, com outro
    nome, fazia aparecerem DUAS entradas para o mesmo arquivo.
    """
    terr = ('  <Terrains Checked="True" Expanded="True">\n'
            f'    <Layer Name="{os.path.splitext(os.path.basename(terreno_hdf))[0]}" '
            f'Type="TerrainLayer" Checked="True" '
            f'Filename=".\\{terreno_hdf}" />\n  </Terrains>\n'
            ) if terreno_hdf else '  <Terrains Checked="True" />\n'
    x = ('<?xml version="1.0" encoding="utf-8"?>\n<RASMapper>\n'
         '  <Version>2.00</Version>\n'
         f'  <RASProjectionFilename Filename=".\\{projeto}.projection" />\n'
         '  <Geometries Checked="True" Expanded="True" />\n'
         '  <Results Checked="True" Expanded="True">\n'
         f'    <Layer Name="TAJAI" Type="RASResults" Checked="True" '
         f'Expanded="True" Filename=".\\{projeto}.p01.hdf">\n'
         '      <Layer Name="Depth" Type="RASResultsMap" Checked="True">\n'
         '        <MapParameters MapType="depth" ProfileIndex="2147483647" '
         'ProfileName="Max" />\n      </Layer>\n'
         '      <Layer Name="WSE" Type="RASResultsMap" Checked="True">\n'
         '        <MapParameters MapType="elevation" ProfileIndex="2147483647" '
         'ProfileName="Max" />\n      </Layer>\n'
         '    </Layer>\n  </Results>\n' + terr +
         '  <MapLayers Checked="True" />\n</RASMapper>\n')
    open(f"{projeto}.rasmap", "w", encoding="utf-8").write(x)
    return f"{projeto}.rasmap"


def validar(projeto):
    """Le de volta com o ras-commander. Erro de formato aparece AQUI."""
    try:
        from ras_commander.geom import GeomCrossSection as G
    except ImportError:
        try:
            from ras_commander import RasGeometry as G
        except ImportError:
            return None
    g = f"{projeto}.g01"
    xs = G.get_cross_sections(g)
    r = xs.iloc[0]
    se = G.get_station_elevation(g, r["River"], r["Reach"], r["RS"])
    bk = G.get_bank_stations(g, r["River"], r["Reach"], r["RS"])
    mn = G.get_mannings_n(g, r["River"], r["Reach"], r["RS"])
    return {"secoes": len(xs), "pontos": len(se),
            "margens": bk, "zonas_manning": len(mn)}
