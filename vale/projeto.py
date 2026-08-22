# -*- coding: utf-8 -*-
"""
Escrita dos arquivos do HEC-RAS.

Formatos conferidos contra os projetos-exemplo oficiais e, quando o
ras-commander esta disponivel, montados por ele. O que so se descobre
apanhando, e que esta aqui:

  - series e #Sta/Elev em COLUNAS DE 8 CARACTERES, 10 por linha. Numa linha
    unica separada por espacos o RAS le "dados faltando";
  - Boundary Location tem SEIS campos com padding (16/16/8 + 3 vazios);
  - Bank Sta precisa casar com um sta da tabela na MESMA precisao (.2f);
  - Junc L&A e o caminho ATRAVES da juncao, da ultima secao de montante a
    primeira de jusante. Gravar 500 m fixo onde a geometria da 75 desequilibra
    a continuidade exatamente na secao que o solver reporta;
  - Reach XY e o trecho, nao o rio inteiro. Gravar o rio todo em cada trecho
    faz a juncao virar teia de aranha e as bank lines saem erradas;
  - Viewing Rectangle e a extensao REAL. Com o placeholder "0,1,1,0" o HDF sai
    com Extents=[0,1,0,1] e o RAS Mapper abre vazio, procurando a bacia dentro
    de um metro quadrado na origem;
  - no .rasmap o terreno so pode ser o .hdf. Declarar o GeoTIFF faz o HEC-RAS
    tentar abri-lo como HDF5 e despejar HDF5-DIAG a cada execucao.
"""
import datetime
import os

import numpy as np

from .config import EPSG, NOME_SRS, WKT

MESES = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
         "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
EDITADO = "Node Last Edited Time= Aug\\15\\2026 00:00:00"


def p16(s):
    return f"{str(s)[:16]:<16}"


def f8(v):
    for fmt in ("{:8.2f}", "{:8.1f}", "{:8.0f}"):
        t = fmt.format(v)
        if len(t) <= 8:
            return t
    return f"{v:8.3g}"[:8]


def serie8(vals):
    v = list(vals)
    return ["".join(f8(x) for x in v[i:i + 10]) for i in range(0, len(v), 10)]


def data_ras(dt):
    return f"{dt.day:02d}{MESES[dt.month-1]}{dt.year},{dt.hour:02d}00"


def contorno(rio, reach, rs):
    return (f"Boundary Location={p16(rio)},{p16(reach)},{str(rs)[:8]:<8}"
            f",        ,                ,                ")


def contorno_faixa(rio, reach, rs_hi, rs_lo):
    return (f"Boundary Location={p16(rio)},{p16(reach)},"
            f"{str(rs_hi)[:8]:<8},{str(rs_lo)[:8]:<8},                ,"
            f"                ")


# ------------------------------------------------------------------ SECAO
def sem_duplicatas(d):
    """Tira estacas repetidas da secao. REDE DE PROTECAO, nao a correcao.

    O HEC-RAS recusa a geometria inteira por causa de um par de pontos com a
    mesma estaca -- "Station and elevation data contains duplicate points" --
    e nao computa nada. A causa foi resolvida no corte (espacamento de pontos
    proporcional a largura), mas o construtor tambem INSERE pontos nas estacas
    das margens, e uma insercao sobre uma amostra existente reintroduz o
    problema. Uma checagem de duas linhas aqui vale mais que um dia perdido
    procurando por que o solver nao roda.
    """
    sta = np.asarray(d["sta"], float)
    if len(sta) < 2:
        return d
    manter = np.concatenate([[True], np.diff(sta) > 1e-4])
    if manter.all():
        return d
    d = dict(d)
    d["sta"] = sta[manter]
    d["z"] = np.asarray(d["z"], float)[manter]
    d["i_thal"] = int(np.argmin(d["z"]))
    return d


def sem_duplicatas_no_texto(linhas):
    """Tira estacas repetidas do TEXTO ja gerado da secao.

    O `sem_duplicatas` acima opera sobre os floats e nao alcanca este caso por
    duas razoes, e o solver recusou o modelo inteiro pelas duas juntas:

    PRECISAO. Ele compara com tolerancia de 1e-4, mas a estaca e GRAVADA com
    duas casas. Dois pontos a 1 mm de distancia sao distintos na memoria e
    identicos no arquivo -- e o HEC-RAS le o arquivo.

    ORDEM. Ele roda ANTES do construtor, e e o construtor que insere as
    estacas das margens. Uma insercao sobre uma amostra existente cria a
    duplicata depois que a limpeza ja passou.

    Aqui a checagem e no que o solver vai ler, campo a campo, entao nao ha
    caminho que escape: construtor, secao manual ou insercao de margem. Custou
    uma rodada dos 12 rios -- Benedito RS 34604.89, Cedros RS 27106.84,
    Itajai_Sul RS 76808.17 e Pombas RS 513.28, quatro secoes de 11.251.
    """
    for i, l in enumerate(linhas):
        if not l.startswith("#Sta/Elev="):
            continue
        try:
            n = int(l.split("=", 1)[1].strip())
        except ValueError:
            return linhas
        campos, j = [], i + 1
        while j < len(linhas) and len(campos) < 2 * n:
            campos += [linhas[j][k:k + 8] for k in range(0, len(linhas[j]), 8)]
            j += 1
        campos = campos[:2 * n]
        if len(campos) < 2 * n:
            return linhas                       # bloco truncado: nao mexe
        pares, ultima = [], None
        for k in range(0, len(campos), 2):
            if campos[k] == ultima:             # estacas sao monotonicas:
                continue                        # comparar com a anterior basta
            ultima = campos[k]
            pares.append((campos[k], campos[k + 1]))
        if len(pares) == n:
            return linhas
        planos = [v for p in pares for v in p]
        novas = ["".join(planos[k:k + 10]) for k in range(0, len(planos), 10)]
        return (linhas[:i] + [f"#Sta/Elev= {len(pares)} "] + novas
                + linhas[j:])
    return linhas


def _secao_manual(d, dx):
    linhas = [f"Type RM Length L Ch R = 1 ,{d['rs']:<8.2f},{dx},{dx},{dx}",
              "XS GIS Cut Line=2",
              "{:16.4f}{:16.4f}{:16.4f}{:16.4f}".format(*d["cut"]),
              EDITADO,
              f"#Sta/Elev= {len(d['sta'])} "]
    linhas += serie8([v for p in zip(d["sta"], d["z"]) for v in p])
    linhas += ["#Mann= 3 ,-1,0",
               "".join(f"{v:>8}" for v in
                       [f"{d['sta'][0]:.2f}", f"{d['n_planicie']:.3f}", "0",
                        f"{d['lb']:.2f}", f"{d['n']:.3f}", "0",
                        f"{d['rb']:.2f}", f"{d['n_planicie']:.3f}", "0"]),
               f"Bank Sta={d['lb']:.2f},{d['rb']:.2f}",
               "XS Rating Curve= 0 ,0", "Exp/Cntr=0.3,0.1", ""]
    return linhas


def htab(d, usar=True):
    """Tabela hidraulica desta secao, pelo GeomHtabUtils.

    Sem isto toda secao fica no padrao do HEC-RAS, e o log enche de
    "Extrapolated above Cross Section Table" -- que e a agua passando do topo
    da TABELA, nao do topo da secao. As duas coisas se confundem facilmente: a
    secao pode ter 100 m de altura util e a tabela cobrir 5.

    A funcao recebe ESCALARES (invert, max_wse), nao um caminho de geometria.
    Chamada com o .g01 ela falha com "missing 1 required positional argument:
    'max_wse'", e a otimizacao simplesmente nao acontece -- foi o que ocorreu
    na primeira execucao, sem que nada mais acusasse.
    """
    if not usar:
        return []
    try:
        from ras_commander.geom import GeomHtabUtils
    except ImportError:
        return []
    z = np.asarray(d["z"], float)
    # SEM try/except aqui, de proposito. Engolir a excecao devolvendo lista
    # vazia faz a otimizacao nao acontecer sem que nada acuse -- e o modelo
    # roda com as tabelas no padrao, extrapolando. Se a chamada quebrar, e
    # melhor quebrar na cara de quem roda.
    r = GeomHtabUtils.calculate_optimal_xs_htab(
        invert=float(z.min()), max_wse=float(z.max()))
    d["htab"] = r
    # PISO CONTRA O FUNDO COMO ELE E GRAVADO. O otimizador devolve fundo+0,01,
    # e 1 cm e da ordem do arredondamento das duas casas com que a estaca e
    # escrita: qualquer diferenca entre o fundo que ele recebeu e o fundo que o
    # HEC-RAS le no arquivo cabe dentro dessa folga. O sintoma e a advertencia
    # de carregamento "XS Htab Starting Elevations ... below the XS invert",
    # que apareceu em 2 secoes de 819 no Benedito -- o RAS repoe o padrao e
    # segue, mas a tabela daquelas secoes deixa de ser a otimizada.
    piso = round(float(z.min()), 2) + 0.02
    inicial = max(float(r["starting_el"]), piso)
    return [f"XS HTab Starting El and Incr={inicial:.2f},"
            f"{r['increment']:.3f}, {int(r['num_points'])} ",
            "XS HTab Horizontal Distribution=-1,-1,-1"]


def secao(rio, reach, d, dx, usar_builder=True, usar_htab=True):
    """Bloco da secao, ja sem estacas repetidas.

    A limpeza fica AQUI, no unico ponto por onde os tres caminhos passam
    (construtor, secao manual e o recuo por excecao). Pendurada em cada
    `return` de dentro, bastaria um caminho novo para reabrir o buraco.
    """
    return sem_duplicatas_no_texto(
        _secao(rio, reach, d, dx, usar_builder, usar_htab))


def _secao(rio, reach, d, dx, usar_builder=True, usar_htab=True):
    """Bloco da secao. Pelo ras-commander quando disponivel.

    Ganho concreto do builder: ele INSERE pontos nas estacas das margens,
    garantindo que Bank Sta coincida com um sta da tabela. A mao isso exige
    casar a precisao dos dois lados, e foi fonte de erro.

    A linha de valores do Manning e reescrita depois: o builder grava com DUAS
    casas, entao 0,035 vira 0,04 (14% a mais de rugosidade no modelo inteiro) e
    os valores de Jarrett -- 0,052, 0,066, 0,070 -- colapsam em tres niveis.
    """
    d = sem_duplicatas(d)
    if not usar_builder:
        return _secao_manual(d, dx)
    try:
        import pandas as pd
        from ras_commander.geom import GeomCrossSection as G
    except ImportError:
        return _secao_manual(d, dx)
    c = d["cut"]
    try:
        r = G.build_cross_section(
            river=rio, reach=reach, rs=f"{d['rs']:.2f}",
            station_elevation=pd.DataFrame({"Station": d["sta"],
                                            "Elevation": d["z"]}),
            cut_line=[(c[0], c[1]), (c[2], c[3])],
            bank_left=float(d["lb"]), bank_right=float(d["rb"]),
            # A COTA das margens, explicita. Sem ela o builder procura um
            # terreno que nao passamos e registra, POR SECAO, um ERROR:
            # "terrain unavailable for bank elevations; interpolated LOB/ROB
            # elevations from station/elevation profile". Sao ~2.000 linhas de
            # ERROR num modelo em que nada esta errado -- ele acaba
            # interpolando do proprio perfil, que e o que queremos. O problema
            # nao e o valor: e que 2.000 falsos ERROR escondem os verdadeiros,
            # e este projeto ja perdeu horas por erro mascarado por outro.
            bank_left_elevation=float(np.interp(d["lb"], d["sta"], d["z"])),
            bank_right_elevation=float(np.interp(d["rb"], d["sta"], d["z"])),
            n_lob=float(d["n_planicie"]), n_channel=float(d["n"]),
            n_rob=float(d["n_planicie"]),
            length_left=float(dx), length_channel=float(dx),
            length_right=float(dx))
        linhas = r.text.splitlines()
        for k, l in enumerate(linhas):
            if l.startswith("#Mann=") and k + 1 < len(linhas):
                linhas[k + 1] = "".join(
                    f"{v:>8}" for v in
                    [f"{d['sta'][0]:.2f}", f"{d['n_planicie']:.3f}", "0",
                     f"{d['lb']:.2f}", f"{d['n']:.3f}", "0",
                     f"{d['rb']:.2f}", f"{d['n_planicie']:.3f}", "0"])
                break
        return (linhas + htab(d, usar_htab)
                + ["XS Rating Curve= 0 ,0", "Exp/Cntr=0.3,0.1", ""])
    except Exception:                                        # noqa: BLE001
        return _secao_manual(d, dx) + htab(d, usar_htab)


# -------------------------------------------------------------- GEOMETRIA
def geometria(op, trechos, juncoes, titulo=None):
    xs_all, ys_all = [], []
    for t in trechos:
        for x, y in t["linha"].coords:
            xs_all.append(x)
            ys_all.append(y)
        for d in t["xs"]:
            xs_all += [d["cut"][0], d["cut"][2]]
            ys_all += [d["cut"][1], d["cut"][3]]
    folga = 0.02 * max(max(xs_all) - min(xs_all), max(ys_all) - min(ys_all), 1.0)
    x0, x1 = min(xs_all) - folga, max(xs_all) + folga
    y0, y1 = min(ys_all) - folga, max(ys_all) + folga

    g = [f"Geom Title={titulo or op.projeto}", "Program Version=7.01",
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
        for dist in j["dists"]:
            g.append(f"Junc L&A={dist:.2f},0")
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
            g += secao(t["rio"], t["reach"], d, dx, op.usar_build_xs,
                       op.usar_htab)

    caminho = op.caminho(f"{op.projeto}.g01")
    with open(caminho, "w", encoding="ascii", errors="replace") as f:
        f.write("\n".join(g) + "\n")
    return caminho


# ------------------------------------------------------------------ FLUXO
def fluxo(op, trechos, cabeceiras, saida, mare, laterais, inicio,
          decl_saida=None):
    u = [f"Flow Title={op.projeto}", "Program Version=7.01",
         f"Use Restart= 0 ", ""]

    # Initial RS para TODO trecho. Um trecho sem vazao inicial parte indefinido
    # e o HEC-RAS acusa. E o valor tem de ser a vazao ACUMULADA pela rede em
    # t=0 -- ver vale/hidrologia.py: com valores chutados, o remanso inicial e
    # calculado para uma vazao e o passo 1 recebe outra, e o sistema inteiro
    # leva um choque.
    for t in trechos:
        # DUAS casas na vazao: com "%.0f" uma cabeceira de 0,21 m3/s virava
        # "0", e trecho com vazao inicial zero parte seco -- exatamente o que
        # se esta tentando evitar aqui.
        u.append(f"Initial RS={p16(t['rio'])},{p16(t['reach'])},"
                 f"{t['xs'][0]['rs']:.0f},{t['q_base']:.2f}")
    u.append("")

    for t in cabeceiras:
        u += [contorno(t["rio"], t["reach"], f"{t['xs'][0]['rs']:.2f}"),
              "Interval=1HOUR",
              f"Flow Hydrograph= {len(t['serie'])} "]
        u += serie8(t["serie"])
        u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=True",
              f"Fixed Start Date/Time={data_ras(inicio)}",
              "Is Critical Boundary=False", "Critical Boundary Flow=", ""]

    for l in laterais:
        u += [contorno_faixa(l["rio"], l["reach"],
                             f"{l['rs_hi']:.2f}", f"{l['rs_lo']:.2f}"),
              "Interval=1HOUR",
              f"Uniform Lateral Inflow Hydrograph= {len(l['serie'])} "]
        u += serie8(l["serie"])
        u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=True",
              f"Fixed Start Date/Time={data_ras(inicio)}", ""]

    # Contorno de jusante: mare SO quando a foz esta no mar. Rodando um rio
    # isolado a foz dele vira a saida do modelo, e a mare de 0,3 m foi imposta
    # a uma secao com fundo em 50 m -- o HEC-RAS recusou os dados antes de
    # computar ("Stage(s) in time series data are below the cross section
    # minimum") e a rodada terminou anunciando exito. Foz terra adentro usa
    # profundidade normal, que e o contorno correto para trecho que continua.
    u.append(contorno(saida["rio"], saida["reach"],
                      f"{saida['xs'][-1]['rs']:.2f}"))
    if mare is not None:
        u += ["Interval=1HOUR", f"Stage Hydrograph= {len(mare)} "]
        u += serie8(mare)
        u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=True",
              f"Fixed Start Date/Time={data_ras(inicio)}", ""]
    else:
        u += [f"Friction Slope={max(float(decl_saida or 0.0), 1e-4):.6f},0", ""]

    caminho = op.caminho(f"{op.projeto}.u01")
    open(caminho, "w", encoding="ascii", errors="replace").write("\n".join(u) + "\n")
    return caminho


# ------------------------------------------------------------------ PLANO
def plano(op, inicio):
    fim = inicio + datetime.timedelta(hours=op.horas - 1)
    p = [f"Plan Title={op.projeto}", "Program Version=7.01",
         f"Short Identifier={op.projeto[:12]}", "Geom File=g01", "Flow File=u01",
         f"Simulation Date={data_ras(inicio)},{data_ras(fim)}"]
    if op.lpi:
        # O solver nao permanente resolve Saint-Venant completo, estavel apenas
        # em regime FLUVIAL. Esta rede tem trechos de serra com 6 a 10% de
        # declividade, onde o escoamento e torrencial por fisica. Sem LPI o
        # solver oscila desde o aquecimento, com erros de 10 a 18 m de nivel
        # batendo o teto de iteracoes aos 20 minutos de simulacao.
        p += ["Mixed Flow Regime", "UNET Froude Reduction=True",
              "UNET Froude Limit= 0.8 ", "UNET Froude Power= 4 "]
    # ZTol: o padrao sao 6 mm. Exigir isso num modelo construido sobre DEM e
    # pedir precisao muito abaixo da que o dado tem.
    p += [f"UNET ZTol= {op.ztol} ", f"UNET ZSATol= {op.ztol} ",
          f"UNET MxIter= {op.max_iter} ",
          f"Computation Interval={op.intervalo}",
          "Output Interval=1HOUR", "Instantaneous Interval=1HOUR",
          "Mapping Interval=1HOUR",
          "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1", "Run RASMapper=-1"]
    caminho = op.caminho(f"{op.projeto}.p01")
    open(caminho, "w", encoding="ascii").write("\n".join(p) + "\n")

    prj = op.caminho(f"{op.projeto}.prj")
    open(prj, "w", encoding="ascii").write("\n".join([
        f"Proj Title={op.projeto}", "Current Plan=p01", "Default Exp/Contr=0.3,0.1",
        "Geom File=g01", "Unsteady File=u01", "Plan File=p01",
        "Y Axis Title=Elevation", "X Axis Title(s)=Main Channel Distance",
        "BEGIN DESCRIPTION:", "Vale do Itajai -- MDT 1 m do SIG-SC",
        "END DESCRIPTION:", "DSS Start Date=", "DSS Start Time=",
        "DSS End Date=", "DSS End Time=", "DSS Export Filename=",
        "DSS Export Rating Curves= 0 ", "DSS Export Rating Curve Sorted= 0 ",
        "DSS Export Volume Flow Curves= 0 ", "DXF Filename=",
        "DXF OffsetX= 0 ", "DXF OffsetY= 0 ", "DXF ScaleX= 1 ",
        "DXF ScaleY= 10 ", "GIS Export Profiles= 0 "]) + "\n")
    return caminho, prj


# ---------------------------------------------------------------- RASMAP
def rasmap(op, terreno_hdf=None):
    caminho = op.caminho(f"{op.projeto}.rasmap")
    t = ""
    if terreno_hdf:
        rel = os.path.relpath(terreno_hdf, op.saida).replace("\\", "/")
        t = (f'  <Terrains>\n'
             f'    <Layer Name="Terreno" Type="TerrainLayer" '
             f'Filename="{rel}" />\n  </Terrains>\n')
    open(caminho, "w", encoding="utf-8").write(
        '<?xml version="1.0" encoding="utf-8"?>\n<RASMapper>\n'
        f'  <Version>2.0.0</Version>\n'
        # O ARQUIVO DE PROJECAO, e nao o de PROJETO. Isto apontava para
        # "<projeto>.prj", que e o arquivo de PROJETO do HEC-RAS -- comeca com
        # "Proj Title=..." e nao tem WKT nenhum. O RAS Mapper tenta ler aquilo
        # como projecao e abre com "Corrupt Projection: GDAL issued a warning
        # that this projection file is corrupted", e dali em diante nada que
        # dependa de georreferencia funciona -- inclusive as modificacoes de
        # terreno, que somem da arvore de camadas. O WKT esta no
        # NOME_SRS, escrito logo abaixo nesta mesma funcao.
        #
        # Os dois arquivos terem extensao .prj no HEC-RAS (projeto e projecao)
        # e a armadilha; o nome parecido escondeu isto por toda a
        # reconstrucao, com o aviso aparecendo nos logs como ruido:
        # "Could not parse WKT from HDF file: Invalid WKT string: Proj Title=".
        #
        # A saida NAO e trocar a extensao para ".projection": isso derruba o
        # RAS Mapper em DisplaySRSFileText com "Referencia de objeto nao
        # definida" assim que a pagina de projecao abre. O que muda e o NOME.
        f'  <RASProjectionFilename Filename=".{os.sep}{NOME_SRS}" />\n'
        f'{t}'
        f'  <Geometries>\n    <Layer Name="{op.projeto}" Type="RASGeometry" '
        f'Filename=".{os.sep}{op.projeto}.g01.hdf" />\n  </Geometries>\n'
        f'  <Results>\n    <Layer Name="{op.projeto}" Type="RASResults" '
        f'Filename=".{os.sep}{op.projeto}.p01.hdf" />\n  </Results>\n'
        '</RASMapper>\n')
    open(op.caminho(NOME_SRS), "w", encoding="ascii").write(WKT)
    # O ".projection" de rodadas anteriores fica para tras e o RAS Mapper
    # ainda o encontra pela arvore de camadas. Sair sem ele evita reabrir o
    # projeto velho e ver o mesmo estouro.
    velho = op.caminho(f"{op.projeto}.projection")
    if os.path.exists(velho):
        os.remove(velho)
    return caminho
