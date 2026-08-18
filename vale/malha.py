# -*- coding: utf-8 -*-
"""
Modelo bidimensional: area de escoamento, contornos e os arquivos do HEC-RAS.

POR QUE 2D. Em 1D toda falha desta reconstrucao esteve em ESTABELECER E MANTER
a condicao inicial de vazao baixa -- o entalhe piloto, a vazao inicial por
trecho, a profundidade normal na foz, o remanso de partida. O modelo 2D nao
precisa de nada disso: comeca SECO e a agua chega pelos contornos. A assinatura
que ficou sem explicacao no Benedito (um terco das leituras de vazao NEGATIVAS,
identica sobre Copernicus e sobre SIG-SC) e de metodo, nao de dado.

O QUE NAO MUDA. 2D sobre modelo de SUPERFICIE nao resolve nada: a calha tem de
estar no terreno, e o Copernicus so enxerga a lamina d'agua. Este modulo so faz
sentido com fonte='sigsc', e `conferir()` recusa a outra.

A ORDEM E OBRIGATORIA, e nao e a que parece natural:

    1. .prj + .p01 + .u01 + .rasmap    o projeto tem de existir antes
    2. perimetro no .g01               GeomStorage.set_2d_flow_area_perimeter
    3. regioes de refino               GeomMesh.add_flowline_refinement_regions
    4. linhas de contorno              GeomBcLines.add_bc_lines
    5. compile_geometry                Ras.exe compila o .g01 -> .g01.hdf
    6. GeomMesh.generate               os centros de celula VOLTAM para o texto
    7. compute_property_tables         tabelas de sub-grade sobre o terreno

O texto e a fonte da verdade; o .g01.hdf e area de trabalho. Gerar a malha
antes de existir o HDF compilado falha, e compilar antes de haver perimetro
compila o vazio -- as duas coisas em silencio.

ONDE A AGUA ENTRA. Um rio isolado em 2D tem tres tipos de contorno, e os tres
sao arcos DO PROPRIO PERIMETRO -- o HEC-RAS nao aceita contorno por dentro da
area:

    montante   a tampa de cima do buffer, com o hidrograma da cabeceira
    foz        a tampa de baixo, com a mare
    laterais   arcos ao longo dos dois lados, com o afluxo difuso

As laterais nao sao um detalhe. No Itajai-Mirim a cabeceira e os tres afluentes
nomeados da BHO somam 338 dos 1.676 km2: sem elas o modelo chegaria em Brusque
com um quinto da agua. Elas sao a traducao 2D do "Uniform Lateral Inflow" do
modelo 1D, e ficam melhor do que la -- a agua desce a encosta ate o rio em vez
de aparecer dentro da calha.
"""
import datetime
import os

import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import substring

from .config import EPSG, WKT
from .projeto import data_ras, p16, serie8


def conferir(op, log=print):
    """Recusa a combinacao que nao faz sentido, antes de gastar a rodada."""
    if op.superficie():
        raise SystemExit(
            "malha 2D sobre modelo de SUPERFICIE nao vale: o Copernicus ve a "
            "lamina d'agua, nao o leito, e a calha simplesmente nao esta no "
            "terreno. Rode com fonte=sigsc.")
    if op.corredor_m < op.buffer_2d:
        # Nao e aviso: a borda da area 2D cairia fora do terreno de 1 m e a
        # malha assentaria sobre o fundo grosseiro JUSTAMENTE onde a agua
        # espalha -- a planicie de inundacao, que e o resultado do modelo.
        raise SystemExit(
            f"corredor_m={op.corredor_m:.0f} m e menor que "
            f"buffer_2d={op.buffer_2d:.0f} m: a borda da area 2D ficaria sem "
            f"terreno fino. Rode com corredor_m={op.buffer_2d + 100:.0f} "
            f"(ou maior).")


# ================================================================ PERIMETRO
def area_2d(eixo, buffer_m, max_pontos=1200, log=print):
    """Poligono da area de escoamento: o eixo do rio dilatado.

    Tres cuidados, os tres pagos com geometria que o RAS recusa:

    ILHAS FORA. Um meandro que se fecha deixa um buraco no buffer. Area 2D e
    poligono SIMPLES; o buraco viraria uma parede dentro da malha.

    UM SO POLIGONO. Eixo com um salto minusculo produz MultiPolygon, e o
    HEC-RAS le so o primeiro -- meio modelo desapareceria calado.

    VERTICES DE MENOS. Buffer de 1.500 m sobre 162 km de rio meandrico sai com
    dezenas de milhares de vertices. Simplifica-se ate caber, conferindo a
    validade a cada passo: simplify() pode cruzar a propria linha.
    """
    pol = eixo.buffer(float(buffer_m), resolution=8, cap_style=1, join_style=1)
    pol = pol.buffer(0)
    if pol.is_empty:
        raise ValueError("o buffer do eixo saiu vazio")
    if pol.geom_type == "MultiPolygon":
        partes = sorted(pol.geoms, key=lambda p: p.area, reverse=True)
        log(f"      buffer saiu em {len(partes)} pedacos; fica o maior "
            f"({partes[0].area / 1e6:.1f} km2)")
        pol = partes[0]
    if list(pol.interiors):
        log(f"      {len(pol.interiors)} ilha(s) de meandro descartada(s)")
        pol = Polygon(pol.exterior)

    tol = 0.0
    while len(pol.exterior.coords) > max_pontos:
        tol = (tol * 1.8) or (float(buffer_m) / 40.0)
        cand = pol.simplify(tol, preserve_topology=True).buffer(0)
        if cand.geom_type == "MultiPolygon":
            cand = max(cand.geoms, key=lambda p: p.area)
        if cand.is_empty or not cand.is_valid or not cand.exterior.is_simple:
            log(f"      simplificacao parou em tol={tol:.0f} m "
                f"({len(pol.exterior.coords)} vertices)")
            break
        pol = Polygon(cand.exterior)
    log(f"      area 2D: {pol.area / 1e6:.1f} km2, "
        f"{len(pol.exterior.coords)} vertices no perimetro")
    return pol


def _anel(pol):
    """Vertices do perimetro sem o ponto repetido do fechamento."""
    c = [tuple(p) for p in pol.exterior.coords]
    if c and c[0] == c[-1]:
        c = c[:-1]
    return c


def _corridas(marca):
    """Todas as sequencias CIRCULARES de True, da maior para a menor.

    Circulares porque o anel nao tem comeco: a tampa de montante quase sempre
    cai em cima do vertice 0, e uma varredura linear a partiria em duas.
    """
    n = len(marca)
    if n == 0 or not any(marca):
        return []
    if all(marca):
        return [list(range(n))]
    i0 = marca.index(False)
    runs, atual = [], []
    for k in range(n):
        i = (i0 + k) % n
        if marca[i]:
            atual.append(i)
        elif atual:
            runs.append(atual)
            atual = []
    if atual:
        runs.append(atual)
    return sorted(runs, key=len, reverse=True)


def _tampa(anel, eixo, buffer_m, montante, frac=0.08):
    """Arco do anel que fecha a ponta de montante (ou de jusante) do buffer.

    Duas condicoes, e as duas sao necessarias. Perto da PONTA do eixo, porque
    a tampa e um semicirculo centrado nela; e projetando no INICIO (ou no fim)
    do eixo, porque um meandro pode trazer o perimetro de volta para perto da
    nascente 40 km rio abaixo, e so a distancia o aceitaria como tampa.
    """
    p = Point(eixo.coords[0] if montante else eixo.coords[-1])
    L = eixo.length
    marca = []
    for x, y in anel:
        q = Point(x, y)
        s = eixo.project(q)
        na_ponta = (s <= frac * L) if montante else (s >= (1.0 - frac) * L)
        marca.append(bool(q.distance(p) <= 1.15 * buffer_m and na_ponta))
    runs = _corridas(marca)
    return runs[0] if runs else []


def _linha_dos(anel, idx):
    return LineString([anel[i] for i in idx]) if len(idx) >= 2 else None


def _centrar(ls, alvo):
    """Recorta `alvo` metros de linha, pelo meio.

    A tampa inteira de um buffer de 1.500 m e um semicirculo de 4,7 km. Lancar
    a vazao de entrada nos 4,7 km poe agua na meia encosta dos dois lados; o
    meio do arco e o que olha rio acima.
    """
    if ls is None or not alvo or ls.length <= alvo:
        return ls
    m = ls.length / 2.0
    return substring(ls, m - alvo / 2.0, m + alvo / 2.0)


# ============================================================= AREA DRENADA
def perfil_area(bho, chave, eixo, tol=600.0):
    """Area drenada ao longo do eixo, dos trechos da BHO. Devolve (s, km2).

    E daqui que sai a distribuicao do afluxo lateral -- e nao de uma constante
    por quilometro. Os degraus do perfil sao os afluentes: onde a area salta,
    entra agua de verdade naquele ponto.
    """
    import geopandas as gpd

    from .rios import normalizar

    g = gpd.read_file(bho).to_crs(EPSG)
    g["chave"] = g["NORIOCOMP"].map(normalizar)
    pares = []
    for r in g[g["chave"] == chave].itertuples():
        geom = r.geometry
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "LineString":
            pontas = [Point(geom.coords[0]), Point(geom.coords[-1])]
        else:
            pontas = [Point(geom.geoms[0].coords[0]),
                      Point(geom.geoms[-1].coords[-1])]
        # so os trechos da CADEIA PRINCIPAL: os bracos com o mesmo nome ficam
        # longe do eixo e trariam area que nao passa por aqui.
        if min(q.distance(eixo) for q in pontas) > tol:
            continue
        pares.append((max(float(eixo.project(q)) for q in pontas),
                      float(r.NUAREAMONT)))
    if not pares:
        return None
    pares.sort()
    s = np.array([p[0] for p in pares], float)
    a = np.maximum.accumulate(np.array([p[1] for p in pares], float))
    # ancora as duas pontas: o eixo comeca na nascente e termina na foz
    return (np.concatenate(([0.0], s, [eixo.length])),
            np.concatenate(([a[0]], a, [a[-1]])))


def area_entre(perfil, s0, s1):
    """Area drenada incremental entre duas chainages, em km2."""
    s, a = perfil
    return float(max(np.interp(s1, s, a) - np.interp(s0, s, a), 0.0))


# ================================================================ CONTORNOS
def contornos(pol, eixo, op, log=print):
    """Os arcos do perimetro que recebem agua. Devolve lista de dicts.

    Cada um traz `nome`, `coords` e `tipo` ('montante'|'foz'|'lateral'); as
    laterais trazem tambem a faixa de chainage (s0, s1) que representam.
    """
    anel = _anel(pol)
    n = len(anel)
    i_mont = _tampa(anel, eixo, op.buffer_2d, True)
    i_foz = _tampa(anel, eixo, op.buffer_2d, False)
    if len(i_mont) < 2 or len(i_foz) < 2:
        raise ValueError(
            f"nao achei as tampas do buffer (montante {len(i_mont)} vertices, "
            f"jusante {len(i_foz)}); o eixo pode estar dobrado sobre si")

    saida = [{"nome": "BC_Montante", "tipo": "montante",
              "coords": list(_centrar(_linha_dos(anel, i_mont),
                                      op.bc_largura).coords)},
             {"nome": "BC_Foz", "tipo": "foz",
              "coords": list(_centrar(_linha_dos(anel, i_foz),
                                      op.bc_largura).coords)}]
    if op.n_laterais <= 0:
        return saida

    # os dois lados sao o que sobra do anel entre as tampas
    usados = set(i_mont) | set(i_foz)
    lados = _corridas([i not in usados for i in range(n)])[:2]
    if len(lados) < 2:
        log("      AVISO: nao separei os dois lados do perimetro; "
            "o modelo fica sem contorno lateral")
        return saida

    L = eixo.length
    nlat = int(op.n_laterais)
    for k, lado in enumerate(lados):
        rot = "E" if k == 0 else "D"
        s_lado = np.array([eixo.project(Point(anel[i])) for i in lado], float)
        for j in range(nlat):
            s0, s1 = j * L / nlat, (j + 1) * L / nlat
            i_meio = int(np.argmin(np.abs(s_lado - 0.5 * (s0 + s1))))
            ini = fim = i_meio
            comp = 0.0
            while comp < op.bc_lateral_largura and (ini > 0 or
                                                    fim < len(lado) - 1):
                ini = max(ini - 1, 0)
                fim = min(fim + 1, len(lado) - 1)
                comp = LineString([anel[i] for i in lado[ini:fim + 1]]).length
            ls = _centrar(_linha_dos(anel, lado[ini:fim + 1]),
                          op.bc_lateral_largura)
            if ls is None or ls.length < op.bc_minima:
                log(f"      AVISO: contorno lateral {rot}{j+1} saiu com "
                    f"{0.0 if ls is None else ls.length:.0f} m "
                    f"(minimo {op.bc_minima:.0f} m); descartado")
                continue
            saida.append({"nome": f"BC_Lat_{rot}{j+1:02d}", "tipo": "lateral",
                          "coords": list(ls.coords), "s0": s0, "s1": s1,
                          "lado": rot})
    return saida


def sem_sobreposicao(bcs, folga=1.0, log=print):
    """Nenhum contorno pode encostar no outro.

    Dois arcos que compartilham face fazem o HEC-RAS somar as duas vazoes na
    mesma celula. Como todos saem do MESMO anel, basta conferir distancia.
    """
    fora = []
    for i, a in enumerate(bcs):
        la = LineString(a["coords"])
        for b in bcs[i + 1:]:
            if la.distance(LineString(b["coords"])) < folga:
                fora.append((a["nome"], b["nome"]))
                log(f"      AVISO: contornos {a['nome']} e {b['nome']} "
                    f"se tocam")
    return fora


# =============================================================== HIDROLOGIA
def vazoes(op, d, bcs, perfil, log=print):
    """Um hidrograma para cada contorno de entrada, pela area que ele drena.

    A vazao especifica e a MESMA do modelo 1D (Q_REF_FOZ / AREA_REF_FOZ), para
    que os dois sejam comparaveis. A area de cada lateral e a metade do
    incremento da faixa dela -- metade por lado.

    CONFERE A SOMA. Um contorno lateral descartado por ser curto demais leva a
    area dele junto, calado, e o modelo passa a receber menos agua do que a
    bacia produz. Aqui a soma e comparada com a area do rio e a diferenca vira
    aviso: e o mesmo tipo de sumico que fez a cabeceira do Benedito receber
    sete vezes menos agua sem nada acusar.
    """
    from .hidrologia import AREA_REF_FOZ, Q_REF_FOZ, hidrograma

    q_esp = Q_REF_FOZ / AREA_REF_FOZ
    bf = float(getattr(op, "base_frac", 0.02))
    a_cab = (float(np.interp(0.0, perfil[0], perfil[1])) if perfil
             else d["area"] * op.fracao_cabeceira)
    laterais = [b for b in bcs if b["tipo"] == "lateral"]

    conferir_km2 = 0.0
    for b in bcs:
        if b["tipo"] == "montante":
            b["area_km2"] = a_cab
        elif b["tipo"] == "lateral":
            inc = (area_entre(perfil, b["s0"], b["s1"]) if perfil
                   else (d["area"] - a_cab) / max(len(laterais), 1))
            b["area_km2"] = inc / 2.0          # metade por lado
        else:
            continue                            # a foz nao recebe vazao
        conferir_km2 += b["area_km2"]
        b["serie"] = hidrograma(q_esp * b["area_km2"], op.horas, base_frac=bf)

    log(f"      cabeceira {a_cab:.0f} km2 + {len(laterais)} laterais = "
        f"{conferir_km2:.0f} km2 de {d['area']:.0f} km2 do rio")
    if conferir_km2 < 0.80 * d["area"]:
        log(f"      AVISO: {d['area'] - conferir_km2:.0f} km2 ficaram de "
            f"fora -- o modelo recebe menos agua do que a bacia produz")
    return bcs


# ============================================================== ARQUIVOS 2D
def geometria(op, pol, bcs, nome_area, eixo=None, titulo=None, log=print):
    """Escreve o .g01 2D: cabecalho, perimetro, refino e contornos."""
    from ras_commander.geom import GeomBcLines, GeomMesh, GeomStorage

    x0, y0, x1, y1 = pol.bounds
    folga = 0.02 * max(x1 - x0, y1 - y0, 1.0)
    cab = [f"Geom Title={titulo or op.projeto}", "Program Version=7.01",
           f"Viewing Rectangle= {x0-folga:.6f} , {x1+folga:.6f} , "
           f"{y1+folga:.6f} , {y0-folga:.6f} ",
           f"Spatial Reference System={WKT}", ""]
    caminho = op.caminho(f"{op.projeto}.g01")
    with open(caminho, "w", encoding="ascii", errors="replace") as f:
        f.write("\n".join(cab) + "\n")

    GeomStorage.set_2d_flow_area_perimeter(
        caminho, nome_area, coordinates=list(pol.exterior.coords),
        point_generation_data=f",,{op.celula:.1f},{op.celula:.1f}",
        create_backup=False)
    GeomStorage.set_2d_flow_area_settings(
        caminho, nome_area, mannings_n=op.n_2d, create_backup=False)
    log(f"      area 2D '{nome_area}': celula {op.celula:g} m, n={op.n_2d:g}")

    if eixo is not None and op.refino_2d > 0:
        import geopandas as gpd
        GeomMesh.add_flowline_refinement_regions(
            caminho, gpd.GeoDataFrame(geometry=[eixo], crs=EPSG),
            buffer_width=op.refino_largura, spacing_dx=op.refino_2d,
            spacing_dy=op.refino_2d, name_prefix="Calha", project_crs=EPSG)
        log(f"      refino ao longo do eixo: {op.refino_2d:g} m numa faixa "
            f"de {2*op.refino_largura:g} m")

    GeomBcLines.add_bc_lines(
        caminho, [{"name": b["nome"], "storage_area": nome_area,
                   "coordinates": b["coords"]} for b in bcs],
        replace_existing=True)
    log(f"      {len(bcs)} linhas de contorno escritas")
    return caminho


def fluxo(op, nome_area, bcs, mare, inicio):
    """Escreve o .u01 2D.

    O `Boundary Location` de contorno 2D tem OITO campos, e nao os seis do 1D:
    o nome da area 2D vai no indice 5 e o nome da linha no indice 7. Com seis
    campos o HEC-RAS le a linha como contorno de RIO, procura um rio chamado
    "" e recusa os dados -- o mesmo desfecho da mare na foz terra adentro.

    Nao ha `Initial RS` e nao ha condicao inicial: a area comeca SECA. Era
    justamente estabelecer e manter essa condicao que derrubava o 1D.
    """
    u = [f"Flow Title={op.projeto}", "Program Version=7.01", "Use Restart= 0 ",
         ""]
    for b in bcs:
        u.append(f"Boundary Location={p16('')},{p16('')},{'':<8},{'':<8},"
                 f"{p16('')},{p16(nome_area)},{p16('')},{p16(b['nome'])}")
        if b["tipo"] == "foz":
            u += ["Interval=1HOUR", f"Stage Hydrograph= {len(mare)} "]
            u += serie8(mare)
        else:
            u += ["Interval=1HOUR", f"Flow Hydrograph= {len(b['serie'])} "]
            u += serie8(b["serie"])
        u += ["DSS Path=", "Use DSS=False", "Use Fixed Start Time=True",
              f"Fixed Start Date/Time={data_ras(inicio)}"]
        if b["tipo"] != "foz":
            # A declividade da linha de energia distribui a vazao pelas celulas
            # do contorno. Sem ela o HEC-RAS pergunta na tela -- e numa rodada
            # sem tela, simplesmente nao roda.
            u += ["Is Critical Boundary=False", "Critical Boundary Flow=",
                  f"Flow Hydrograph Slope= {op.decl_bc:.5f} "]
        u.append("")
    caminho = op.caminho(f"{op.projeto}.u01")
    open(caminho, "w", encoding="ascii", errors="replace").write(
        "\n".join(u) + "\n")
    return caminho


def plano(op, inicio, nome_area, log=print):
    """Escreve o .p01, o .prj e as opcoes 2D do plano."""
    fim = inicio + datetime.timedelta(hours=op.horas - 1)
    p = [f"Plan Title={op.projeto}", "Program Version=7.01",
         f"Short Identifier={op.projeto[:12]}", "Geom File=g01",
         "Flow File=u01",
         f"Simulation Date={data_ras(inicio)},{data_ras(fim)}",
         f"Computation Interval={op.intervalo_2d}",
         "Output Interval=1HOUR", "Instantaneous Interval=1HOUR",
         "Mapping Interval=1HOUR",
         "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1",
         "Run RASMapper=-1"]
    caminho = op.caminho(f"{op.projeto}.p01")
    open(caminho, "w", encoding="ascii").write("\n".join(p) + "\n")

    prj = op.caminho(f"{op.projeto}.prj")
    open(prj, "w", encoding="ascii").write("\n".join([
        f"Proj Title={op.projeto}", "Current Plan=p01",
        "Default Exp/Contr=0.3,0.1", "Geom File=g01", "Unsteady File=u01",
        "Plan File=p01", "Y Axis Title=Elevation",
        "X Axis Title(s)=Main Channel Distance", "BEGIN DESCRIPTION:",
        "Vale do Itajai -- modelo 2D sobre MDT 1 m do SIG-SC",
        "END DESCRIPTION:", "DSS Start Date=", "DSS Start Time=",
        "DSS End Date=", "DSS End Time=", "DSS Export Filename=",
        "DSS Export Rating Curves= 0 ", "DSS Export Rating Curve Sorted= 0 ",
        "DSS Export Volume Flow Curves= 0 ", "DXF Filename=",
        "DXF OffsetX= 0 ", "DXF OffsetY= 0 ", "DXF ScaleX= 1 ",
        "DXF ScaleY= 10 ", "GIS Export Profiles= 0 "]) + "\n")

    from ras_commander import RasPlan
    # DWE (onda difusiva) de proposito na primeira rodada. O conjunto completo
    # (SWE-ELM) e mais fisico e MUITO menos tolerante: se o modelo nao fecha
    # volume em DWE, nao ha por que atribuir a diferenca a fisica.
    RasPlan.set_2d_flow_options(
        caminho, mesh_name=nome_area, equation_set=op.equacao_2d,
        water_surface_tolerance=op.ztol_2d, volume_tolerance=op.ztol_2d,
        max_iterations=op.max_iter_2d, theta=1.0, theta_warmup=1.0,
        initial_conditions_time_hours=op.ic_horas,
        time_step_use_courant=bool(op.courant_2d),
        time_step_max_courant=op.courant_2d or None,
        time_step_min_courant=(op.courant_2d / 3.0) if op.courant_2d else None,
        computation_interval=op.intervalo_2d, include_default=True)
    log(f"      plano 2D: {op.equacao_2d}, dt={op.intervalo_2d}, "
        f"Courant<={op.courant_2d or 'fixo'}, aquecimento {op.ic_horas} h")
    return caminho, prj


def malhar(op, nome_area, log=print):
    """Compila o .g01, gera as celulas e as tabelas de sub-grade.

    `generate` extrai os centros de celula do HDF compilado e os grava DE VOLTA
    no texto. `compute_property_tables` e o passo que faz o MDT de 1 m valer
    sob celula de 100 m: cada celula guarda a curva cota-volume e cada face a
    curva cota-area do terreno FINO que ha dentro dela. Sem essas tabelas a
    celula vira um bloco plano e a calha do rio desaparece -- que e exatamente
    o defeito que o 2D veio corrigir. Por isso falhar aqui e ERRO, nao aviso.
    """
    from ras_commander import init_ras_project
    from ras_commander.geom import GeomMesh

    ras = init_ras_project(op.caminho(f"{op.projeto}.prj"), op.ras_exe)
    hdf = GeomMesh.compile_geometry("01", ras_object=ras)
    log(f"      geometria compilada -> {os.path.basename(str(hdf))}")

    r = GeomMesh.generate("01", mesh_name=nome_area, cell_size=op.celula,
                          min_face_length_ratio=op.face_minima, ras_object=ras)
    n = getattr(r, "cell_count", None) or getattr(r, "n_cells", None)
    log(f"      malha gerada: {n if n is not None else '?'} celulas")

    conflitos = GeomMesh.detect_bc_conflicts(hdf, op.celula)
    for c in conflitos:
        log(f"      AVISO: conflito de contorno -- {c}")

    if not GeomMesh.compute_property_tables("01", mesh_name=nome_area,
                                            force=True, ras_object=ras):
        raise RuntimeError(
            "compute_property_tables falhou: sem as tabelas de sub-grade a "
            "celula fica plana e a calha do rio some do modelo")
    log("      tabelas de sub-grade: ok")
    return hdf, r, conflitos


# =================================================================== MONTAR
def montar(op, d, mare, inicio, log=print):
    """Constroi o modelo 2D de um rio, do eixo aos arquivos prontos.

    `d` e um eixo de vale.eixos.montar(): precisa de 'linha', 'area', 'ras' e
    'chave'. Devolve o dicionario de estado para as figuras e a auditoria.
    """
    conferir(op)
    eixo = d["linha"]
    nome_area = str(d["ras"])[:16]

    pol = area_2d(eixo, op.buffer_2d, op.perimetro_max_pontos, log)
    bcs = contornos(pol, eixo, op, log)
    tocam = sem_sobreposicao(bcs, log=log)
    perfil = perfil_area(op.bho, d["chave"], eixo)
    if perfil is None:
        log("      AVISO: sem perfil de area da BHO; afluxo lateral uniforme")
    bcs = vazoes(op, d, bcs, perfil, log)

    plano(op, inicio, nome_area, log)
    fluxo(op, nome_area, bcs, mare, inicio)
    geometria(op, pol, bcs, nome_area, eixo, log=log)
    hdf, res, conflitos = malhar(op, nome_area, log)

    return {"area": pol, "eixo": eixo, "bcs": bcs, "perfil_area": perfil,
            "nome_area": nome_area, "hdf": str(hdf), "malha": res,
            "conflitos": conflitos, "tocam": tocam}
