# -*- coding: utf-8 -*-
"""Monta o modelo 2D da REDE: o vale inteiro em malha, cabeceiras 1D nos arcos.

    python scripts/vale_2d_rede.py taha_ai_novo --saida taha_ai_2d

E a extensao do vale/malha.py (que e por rio) para a rede dos 12 rios:

    perimetro   uniao dos buffers dos 12 eixos RETIFICADOS do 1D
    cabeceiras  11 arcos no perimetro, um onde cada rio entra (o Acu nao
                tem: ele NASCE dentro da area, da juncao Sul+Oeste)
    laterais    12 arcos, um por rio, no ponto do anel mais proximo do
                MEIO do eixo do rio -- o "Uniform Lateral Inflow" do 1D
                entra ali, inteiro (v0: um arco por rio; a distribuicao
                fina vem depois se o modelo fechar volume)
    foz         a mare do Atlantico, no arco mais proximo do fim do Acu

AS SERIES SAO AS MESMAS DO 1D: lidas do taha_ai.u01 da rede (com os pisos
de estiagem e o toco de cabeceira do Acu ja fundido na lateral). Nada e
re-derivado; um modelo e a traducao do outro, comparavel numero a numero.

POR QUE ISTO E O "1D/2D" DESTA BACIA. As 33 rodadas do 1D morreram todas em
ESTABELECER E MANTER a vazao baixa na rede densa; o 2D comeca seco e nao tem
o problema. O que o 1D de montanha fazia -- rotear minutos de defasagem ate
o vale -- os arcos de contorno fazem por hidrograma. A calha levantada de
1983 continua no 1D para o que o 1D faz bem; a MANCHA, que e o resultado,
sai da malha sobre o terreno de 1 m via tabelas de sub-grade.
"""
import datetime
import json
import os
import re
import sys

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DIR)
sys.path.insert(0, os.path.dirname(DIR))

from shapely.geometry import LineString, MultiLineString, Point   # noqa: E402
from shapely.ops import linemerge, substring                      # noqa: E402

from qc_geometria import ler_eixos                                # noqa: E402
from vale.config import Opcoes                                    # noqa: E402
from vale import malha, projeto as vprojeto                       # noqa: E402


def _arg(argv, chave, padrao=None, tipo=str):
    return tipo(argv[argv.index(chave) + 1]) if chave in argv else padrao


def eixos_da_rede(g01):
    """Um LineString por RIO (reaches emendados), na ordem montante->foz."""
    E = ler_eixos(g01)
    por_rio = {}
    for (rio, reach), ls in E.items():
        por_rio.setdefault(rio, []).append((reach, ls))
    saida = {}
    for rio, partes in por_rio.items():
        partes.sort(key=lambda t: t[0])          # R1, R2, ...
        m = linemerge(MultiLineString([ls for _, ls in partes]))
        if m.geom_type == "MultiLineString":
            # emenda com folga: liga na ordem dos reaches mesmo com hiato
            coords = []
            for _, ls in partes:
                coords += list(ls.coords)
            m = LineString(coords)
        saida[rio] = m
    return saida


def series_do_u01(u01):
    """{(rio, tipo): serie} com tipo em 'cab'|'lat'|'mare', do u01 da rede."""
    t = open(u01, encoding="latin-1", errors="replace").read() \
        .replace("\r", "")
    blocos = re.split(r"(?=^Boundary Location=)", t, flags=re.M)
    out = {}
    for b in blocos:
        m = re.match(r"Boundary Location=([^,]+),", b)
        if not m:
            continue
        rio = m.group(1).strip()
        for chave, tipo in (("Flow Hydrograph", "cab"),
                            ("Uniform Lateral Inflow Hydrograph", "lat"),
                            ("Lateral Inflow Hydrograph", "lat"),
                            ("Stage Hydrograph", "mare")):
            mm = re.search(r"^%s=\s*(\d+)\s*$" % chave, b, flags=re.M)
            if not mm:
                continue
            n = int(mm.group(1))
            vals = []
            for l in b[mm.end() + 1:].split("\n"):
                if not l.strip() or l[:1].isalpha():
                    break
                vals += [float(l[i:i + 8]) for i in range(0, len(l), 8)
                         if l[i:i + 8].strip()]
                if len(vals) >= n:
                    break
            out[(rio, tipo)] = np.array(vals[:n])
            break
    return out


def arco_no_anel(anel_ls, ponto, largura):
    """Arco de `largura` m do anel, centrado no vertice mais proximo."""
    s = anel_ls.project(Point(ponto))
    L = anel_ls.length
    a = max(s - largura / 2.0, 0.0)
    b = min(s + largura / 2.0, L)
    ls = substring(anel_ls, a, b)
    return list(ls.coords)


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    pasta_1d = argv[0].rstrip("/\\")
    saida = _arg(argv, "--saida", "taha_ai_2d")
    celula = _arg(argv, "--celula", None, float)

    # ------------------------------------------------------------- opcoes
    bruto = json.load(open(os.path.join(pasta_1d, "opcoes.json"),
                           encoding="utf-8"))
    import dataclasses
    campos = {f.name for f in dataclasses.fields(Opcoes)}
    bruto = {k: v for k, v in bruto.items() if k in campos}
    op = Opcoes(**bruto)
    op.projeto = os.path.basename(saida)
    op.saida = saida
    op.fonte = "sigsc"            # o corredor de 1 m E SIG-SC; o campo
    op.buffer_2d = min(op.buffer_2d, op.corredor_m - 100.0)
    # refino exige o HDF compilado ANTES (ovo-galinha na biblioteca);
    # v0 roda sem ele -- a calha vem da tabela de sub-grade de 1 m
    op.refino_2d = 0.0
    if celula:
        op.celula = celula
    os.makedirs(saida, exist_ok=True)
    malha.conferir(op)
    print(f"1D de origem : {pasta_1d}")
    print(f"projeto 2D   : {saida}   celula {op.celula:g} m   "
          f"buffer {op.buffer_2d:g} m")

    # -------------------------------------------------------------- eixos
    g01 = os.path.join(pasta_1d, "taha_ai.g01")
    eixos = eixos_da_rede(g01)
    print(f"eixos        : {len(eixos)} rios, "
          f"{sum(e.length for e in eixos.values())/1000:.0f} km")
    linhas_rede = list(eixos.values())
    # canais reais que o eixo 1D ignora (ex.: canal retificado do baixo
    # Mirim, visivel no relevo e fora do buffer do eixo) entram no
    # perimetro por doc/canais_extras.geojson -- em 2D a agua os acha
    # sozinha pelo terreno, desde que a area os cubra
    extras = os.path.join(os.path.dirname(DIR), "doc",
                          "canais_extras.geojson")
    if os.path.exists(extras):
        gj = json.load(open(extras, encoding="utf-8"))
        for f in gj.get("features", []):
            if f["geometry"]["type"] == "LineString":
                linhas_rede.append(
                    LineString(f["geometry"]["coordinates"]))
                print(f"canal extra  : "
                      f"{f['properties'].get('nome', '?')} no perimetro")
    rede = MultiLineString(linhas_rede)
    pol = malha.area_2d(rede, op.buffer_2d, op.perimetro_max_pontos)
    # perimetro em passo UNIFORME >= celula/2: vertices apertados geram
    # "face errors on the perimeter of the mesh" no escritor do RAS
    # (medido: 8 faces com o anel simplificado cru)
    from shapely.geometry import Polygon as _Pol
    # fechamento morfologico: entalhes e estrangulamentos concavos do
    # poligono-uniao viram faces impossiveis no mesher ("5 face errors")
    pol = pol.buffer(2 * op.celula).buffer(-2 * op.celula).buffer(0)
    if pol.geom_type == "MultiPolygon":
        pol = max(pol.geoms, key=lambda p: p.area)
    pol = _Pol(pol.exterior)
    anel = LineString(list(pol.exterior.coords))
    passo = max(1.5 * op.celula, 150.0)
    n_pts = max(int(anel.length / passo), 32)
    pts = [anel.interpolate(k * anel.length / n_pts).coords[0]
           for k in range(n_pts)]
    pol = _Pol(pts).buffer(0)
    if pol.geom_type == "MultiPolygon":
        pol = max(pol.geoms, key=lambda p: p.area)
    pol = _Pol(pol.exterior)

    # o escritor recusa perimetro fora do TERRENO; puxar vertice a vertice
    # serrilha e cria face errors (medido: 37). O caminho limpo: recortar
    # pelo FOOTPRINT de validade do terreno (erodido), e alisar por
    # ABERTURA morfologica, que nao re-expande alem do dado
    import rasterio
    from rasterio import features as rfeatures
    from shapely.geometry import shape as _shape
    vrt = rasterio.open(os.path.join(pasta_1d, "Terrain",
                                     "taha_ai_Terreno_v2.vrt"))
    dec = 20                       # footprint a ~20 m
    alt = vrt.height // dec
    larg = vrt.width // dec
    banda = vrt.read(1, out_shape=(alt, larg))
    valido = (banda != vrt.nodata).astype(np.uint8)
    Tdec = vrt.transform * vrt.transform.scale(vrt.width / larg,
                                               vrt.height / alt)
    pes = [_shape(g) for g, v in
           rfeatures.shapes(valido, mask=valido.astype(bool),
                            transform=Tdec) if v == 1]
    import shapely.ops as _ops
    foot = _ops.unary_union(pes).buffer(0)
    foot = foot.buffer(-150.0)     # margem de seguranca
    pol = pol.intersection(foot).buffer(0)
    if pol.geom_type != "Polygon":
        pol = max(pol.geoms, key=lambda p: p.area)
    pol = pol.buffer(-op.celula).buffer(op.celula).buffer(0)  # abertura
    if pol.geom_type != "Polygon":
        pol = max(pol.geoms, key=lambda p: p.area)
    pol = _Pol(pol.exterior)
    # reamostra de novo apos o recorte
    anel2 = LineString(list(pol.exterior.coords))
    n2 = max(int(anel2.length / passo), 32)
    pol = _Pol([anel2.interpolate(k * anel2.length / n2).coords[0]
                for k in range(n2)]).buffer(0)
    if pol.geom_type != "Polygon":
        pol = max(pol.geoms, key=lambda p: p.area)
    pol = _Pol(pol.exterior)
    print(f"perimetro    : recortado ao terreno e reamostrado a "
          f"{passo:.0f} m ({len(pol.exterior.coords)} vertices, "
          f"{pol.area/1e6:.0f} km2)")

    # ----------------------------------------------------------- contornos
    series = series_do_u01(os.path.join(pasta_1d, "taha_ai.u01"))
    anel_ls = LineString(list(pol.exterior.coords))
    bcs = []
    for rio, eixo in eixos.items():
        if (rio, "cab") in series and rio != "Itajai_Acu":
            p0 = np.asarray(eixo.coords[0], float)
            bcs.append({"nome": f"BC_{rio[:11]}"[:16], "tipo": "montante",
                        "coords": arco_no_anel(anel_ls, p0, op.bc_largura),
                        "serie": series[(rio, "cab")]})
        if (rio, "lat") in series:
            meio = np.asarray(
                eixo.interpolate(0.5 * eixo.length).coords[0], float)
            bcs.append({"nome": f"BC_L_{rio[:10]}"[:16], "tipo": "lateral",
                        "coords": arco_no_anel(anel_ls, meio,
                                               op.bc_lateral_largura),
                        "serie": series[(rio, "lat")]})
    mare = series.get(("Itajai_Acu", "mare"))
    if mare is None:
        for (rio, tp), v in series.items():
            if tp == "mare":
                mare = v
    fim_acu = np.asarray(eixos["Itajai_Acu"].coords[-1], float)
    bcs.append({"nome": "BC_Foz", "tipo": "foz",
                "coords": arco_no_anel(anel_ls, fim_acu, op.bc_largura)})
    print(f"contornos    : {len(bcs)} arcos "
          f"({sum(1 for b in bcs if b['tipo']=='montante')} cabeceiras, "
          f"{sum(1 for b in bcs if b['tipo']=='lateral')} laterais, 1 foz)")
    malha.sem_sobreposicao(bcs)

    q_total = sum(b["serie"].max() for b in bcs if b["tipo"] != "foz")
    print(f"pico somado  : {q_total:.0f} m3/s entrando na area")

    # ------------------------------------------------------------ arquivos
    # a ordem do malha.py: prj+p01+u01+rasmap+g01 existem ANTES de
    # inicializar o projeto; so entao as opcoes 2D do plano e a malha
    inicio = datetime.datetime(2026, 8, 1, 0, 0)
    nome_area = "Vale_Itajai"
    fim = inicio + datetime.timedelta(hours=op.horas - 1)
    d_ras = vprojeto.data_ras
    p = [f"Plan Title={op.projeto}", "Program Version=7.01",
         f"Short Identifier={op.projeto[:12]}", "Geom File=g01",
         "Flow File=u01",
         f"Simulation Date={d_ras(inicio)},{d_ras(fim)}",
         f"Computation Interval={op.intervalo_2d}",
         "Output Interval=1HOUR", "Instantaneous Interval=1HOUR",
         "Mapping Interval=1HOUR",
         "Run HTab=-1", "Run UNet=-1", "Run PostProcess=-1",
         "Run RASMapper=-1"]
    open(op.caminho(f"{op.projeto}.p01"), "w", encoding="ascii") \
        .write("\n".join(p) + "\n")
    open(op.caminho(f"{op.projeto}.prj"), "w", encoding="ascii").write(
        "\n".join([f"Proj Title={op.projeto}", "Current Plan=p01",
                   "Default Exp/Contr=0.3,0.1", "Geom File=g01",
                   "Unsteady File=u01", "Plan File=p01",
                   "Y Axis Title=Elevation",
                   "X Axis Title(s)=Main Channel Distance",
                   "BEGIN DESCRIPTION:",
                   "Vale do Itajai -- rede em malha 2D",
                   "END DESCRIPTION:"]) + "\n")
    malha.fluxo(op, nome_area, bcs, list(mare), inicio)
    terreno = os.path.abspath(os.path.join(pasta_1d, "Terrain",
                                           "taha_ai_Terreno_v2.hdf"))
    vprojeto.rasmap(op, terreno)
    malha.geometria(op, pol, bcs, nome_area, eixo=rede)

    # ------------------------------------------- pontos de computacao
    # o escritor do Ras.exe NAO gera os pontos sozinho ("2D Flow Area was
    # not created successfully" com Points= 0); a grade vai pronta no texto
    import shapely as _shp
    g01_2d = op.caminho(f"{op.projeto}.g01")
    t = open(g01_2d, encoding="latin-1").read().replace("\r\n", "\n") \
        .split("\n")
    inset = pol.buffer(-15.0)
    x0, y0, x1, y1 = pol.bounds
    # fase deslocada ~14 m: face de celula que coincide EXATAMENTE com a
    # fronteira entre camadas do terreno reprova no avaliador de cobertura
    # do RAS (medido: 4 faces todas em y=7020153.4, a emenda das camadas)
    gx, gy = np.meshgrid(
        np.arange(x0 + op.celula / 2 + 13.7, x1, op.celula),
        np.arange(y0 + op.celula / 2 + 13.7, y1, op.celula))
    dentro = _shp.contains_xy(inset, gx.ravel(), gy.ravel())
    px, py = gx.ravel()[dentro], gy.ravel()[dentro]
    plano = np.empty(2 * len(px))
    plano[0::2] = px
    plano[1::2] = py
    corpo, lin = [], ""
    for k, v in enumerate(plano):
        lin += "%16.4f" % v
        if (k + 1) % 4 == 0:
            corpo.append(lin)
            lin = ""
    if lin:
        corpo.append(lin)
    ip = next(i for i, l in enumerate(t)
              if l.startswith("Storage Area 2D Points="))
    t = t[:ip] + ["Storage Area 2D Points= %d " % len(px)] + corpo \
        + t[ip + 1:]
    # o escritor REGENERA a grade do Point Generation Data (ignora a fase
    # dos pontos escritos); a origem deslocada tem de ir NA STRING, senao
    # as faces voltam a cair na emenda das camadas do terreno
    for k, l in enumerate(t):
        if l.startswith("Storage Area Point Generation Data="):
            t[k] = ("Storage Area Point Generation Data=%.2f,%.2f,%.1f,%.1f"
                    % (x0 + op.celula / 2 + 13.7, y0 + op.celula / 2 + 13.7,
                       op.celula, op.celula))
            break
    from ras_io import escrever as _escrever
    _escrever(g01_2d, "\n".join(t))
    print(f"      {len(px)} pontos de computacao gravados")

    from ras_commander import RasPlan, init_ras_project
    init_ras_project(op.caminho(f"{op.projeto}.prj"), op.ras_exe)
    RasPlan.set_2d_flow_options(
        op.caminho(f"{op.projeto}.p01"), mesh_name=nome_area,
        equation_set=op.equacao_2d,
        water_surface_tolerance=op.ztol_2d, volume_tolerance=op.ztol_2d,
        max_iterations=op.max_iter_2d, theta=1.0, theta_warmup=1.0,
        initial_conditions_time_hours=op.ic_horas,
        time_step_use_courant=bool(op.courant_2d),
        time_step_max_courant=op.courant_2d or None,
        time_step_min_courant=(op.courant_2d / 3.0) if op.courant_2d
        else None,
        computation_interval=op.intervalo_2d, include_default=True)
    print(f"      plano 2D: {op.equacao_2d}, dt={op.intervalo_2d}, "
          f"Courant<={op.courant_2d or 'fixo'}")

    # -------------------------------------------------- compilar via Ras.exe
    # a biblioteca nao compila texto->HDF; o PREPROCESSADOR GEOMETRICO do
    # proprio Ras.exe compila (Run HTab=-1, UNet=0) -- o mesmo truque do
    # ler_erros_geometria.preparar_hdf
    import shutil
    from ras_commander import RasCmdr
    from ras_commander.geom import GeomMesh
    p01 = op.caminho(f"{op.projeto}.p01")
    guarda = p01 + ".antes_do_htab"
    shutil.copy2(p01, guarda)
    t = open(p01, encoding="latin-1").read()
    t = t.replace("Run UNet=-1", "Run UNet=0", 1) \
         .replace("Run PostProcess=-1", "Run PostProcess=0", 1) \
         .replace("Run RASMapper=-1", "Run RASMapper=0", 1)
    open(p01, "w", encoding="ascii", errors="replace").write(t)
    try:
        ras = init_ras_project(op.caminho(f"{op.projeto}.prj"), op.ras_exe)
        RasCmdr.compute_plan("01", ras_object=ras, force_rerun=True,
                             clear_geompre=True)
    finally:
        shutil.move(guarda, p01)
    hdf = op.caminho(f"{op.projeto}.g01.hdf")
    if not os.path.exists(hdf):
        raise SystemExit("o preprocessador nao produziu o g01.hdf")
    print(f"      geometria compilada -> {os.path.basename(hdf)}")

    ras = init_ras_project(op.caminho(f"{op.projeto}.prj"), op.ras_exe)
    r = GeomMesh.generate("01", mesh_name=nome_area, cell_size=op.celula,
                          min_face_length_ratio=op.face_minima,
                          ras_object=ras)
    n = getattr(r, "cell_count", None) or getattr(r, "n_cells", None)
    print(f"      malha gerada: {n if n is not None else '?'} celulas")
    for c in GeomMesh.detect_bc_conflicts(hdf, op.celula):
        print(f"      AVISO conflito de contorno: {c}")
    GeomMesh.set_geometry_association("01", terrain_hdf_path=terreno,
                                      ras_object=ras)
    print(f"      terreno associado: {os.path.basename(terreno)}")
    if not GeomMesh.compute_property_tables("01", mesh_name=nome_area,
                                            force=True, ras_object=ras):
        raise SystemExit("compute_property_tables falhou: sem sub-grade "
                         "a celula fica plana e a calha some")
    print("      tabelas de sub-grade: ok")
    print("\npronto para rodar: "
          f"{op.caminho(op.projeto + '.prj')}")


if __name__ == "__main__":
    main(sys.argv[1:])
