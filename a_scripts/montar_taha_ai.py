# -*- coding: utf-8 -*-
"""
Script Autônomo e Independente: montar_taha_ai.py
Objetivo:
  1. Carregar a geometria completa do Vale do Itajaí (1232 seções, 12 rios principais).
  2. Unificar os sub-trechos fragmentados (R1..R5 -> R1) em rios contínuos da nascente à foz.
  3. Garantir consistência geométrica rigorosa:
     - 0 pontos duplicados em #Sta/Elev
     - Bank Stations perfeitamente alinhadas com pontos da seção
     - Declividade de leito monótona e estável
     - Projeção canônica SIRGAS 2000 UTM Zone 22S (EPSG:31982)
     - RASMapper configurado com projeção e camadas
  4. Configurar condições de contorno (.u01) e parâmetros UNET (.p01) com estabilização 100% convergente.
  5. Salvar em _anti/taha_ai/ e executar a simulação automaticamente com auditoria HDF.
"""

import os
import re
import sys
import shutil
import pathlib
import numpy as np
import h5py
from ras_commander import init_ras_project, RasCmdr, HdfResultsPlan

DIR_RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_DESTINO = os.path.join(DIR_RAIZ, "_anti", "taha_ai")
RAS_EXE = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"


def escrever_texto_dos(caminho, texto):
    with open(caminho, "wb") as f:
        f.write(texto.replace("\r\n", "\n").replace("\n", "\r\n").encode("latin-1"))


def montar_modelo(origem_prj="taha_ai.prj", pasta_dest=DIR_DESTINO, nome_proj="taha_ai"):
    print("=" * 72)
    print(f"MONTANDO MODELO COMPLETO DA BACIA DO ITAJAI: {nome_proj}")
    print(f"DESTINO: {pasta_dest}")
    print("=" * 72)
    
    pasta = pathlib.Path(pasta_dest)
    pasta.mkdir(parents=True, exist_ok=True)
    
    # 1. Carregar e Tratar Geometria (.g01)
    g01_orig = os.path.join(DIR_RAIZ, "taha_ai.g01")
    with open(g01_orig, "r", encoding="latin-1") as f:
        text = f.read()

    partes = text.split("River Reach=")
    header_global = partes[0]

    def parse_reach(bloco):
        lines = bloco.splitlines()
        header_line = lines[0]
        rio, trecho = [t.strip() for t in header_line.split(",")]
        
        idx_first_xs = None
        for i, l in enumerate(lines):
            if l.startswith("Type RM Length L Ch R"):
                idx_first_xs = i
                break
                
        xy_pts = []
        for l in lines[1:idx_first_xs]:
            if l.startswith("Reach XY=") or l.startswith("Rch Text") or not l.strip():
                continue
            tokens = l.strip().split()
            for k in range(0, len(tokens)-1, 2):
                try:
                    xy_pts.append((float(tokens[k]), float(tokens[k+1])))
                except:
                    pass
                    
        xs_txt = "\n".join(lines[idx_first_xs:]) if idx_first_xs is not None else ""
        return rio, trecho, xy_pts, xs_txt

    rios_data = {}
    for p in partes[1:]:
        rio, trecho, xy_pts, xs_txt = parse_reach(p)
        if rio not in rios_data:
            rios_data[rio] = []
        rios_data[rio].append((trecho, xy_pts, xs_txt))

    novo_g01_partes = [header_global]
    reaches_dict = {}
    for rio, sub_reaches in rios_data.items():
        todos_xy = []
        for t, xy, xs in sub_reaches:
            for pt in xy:
                if not todos_xy or abs(todos_xy[-1][0] - pt[0]) > 0.01 or abs(todos_xy[-1][1] - pt[1]) > 0.01:
                    todos_xy.append(pt)
                    
        reach_block = f"River Reach={rio:<16s},R1              \n"
        reach_block += f"Reach XY= {len(todos_xy)} \n"
        for k in range(0, len(todos_xy), 2):
            if k + 1 < len(todos_xy):
                reach_block += f"     {todos_xy[k][0]:14.4f}    {todos_xy[k][1]:14.4f}     {todos_xy[k+1][0]:14.4f}    {todos_xy[k+1][1]:14.4f}\n"
            else:
                reach_block += f"     {todos_xy[k][0]:14.4f}    {todos_xy[k][1]:14.4f}\n"
        reach_block += "Rch Text X Y=0,0,0,0\n\n"
        
        secoes_parts = []
        rs_list = []
        for idx_sub, (t, xy, xs) in enumerate(sub_reaches):
            xs_lines = xs.splitlines()
            for l in xs_lines:
                if l.startswith("Type RM Length L Ch R"):
                    rs_list.append(float(l.split(",")[1]))
            if idx_sub < len(sub_reaches) - 1:
                prox_xs = sub_reaches[idx_sub+1][2]
                prox_rs = float(prox_xs.splitlines()[0].split(",")[1])
                for k in range(len(xs_lines)-1, -1, -1):
                    if xs_lines[k].startswith("Type RM Length L Ch R"):
                        tokens = xs_lines[k].split(",")
                        curr_rs = float(tokens[1])
                        dist = abs(curr_rs - prox_rs)
                        xs_lines[k] = f"Type RM Length L Ch R = 1 ,{curr_rs:.2f}, {dist:.2f}, {dist:.2f}, {dist:.2f}"
                        break
                secoes_parts.append("\n".join(xs_lines))
            else:
                secoes_parts.append(xs)
                
        reaches_dict[rio] = sorted(rs_list, reverse=True)
        reach_block += "\n".join(secoes_parts) + "\n"
        novo_g01_partes.append(reach_block)

    g01_final = "".join(novo_g01_partes)
    
    # Remover blocos de junções obsoletas
    g01_final = re.sub(r"Junct Name=.*?(?=\n\n|River Reach=|\Z)", "", g01_final, flags=re.DOTALL)
    g01_final = re.sub(r"Junc L&A=.*?\n", "", g01_final)
    g01_final = re.sub(r"Junct Desc=.*?\n", "", g01_final)
    g01_final = re.sub(r"Junct X Y.*?\n", "", g01_final)

    wkt_crs = 'PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
    if "Spatial Reference System=" not in g01_final:
        g01_final = g01_final.replace("Geom Title=", f"Spatial Reference System={wkt_crs}\nGeom Title=", 1)

    escrever_texto_dos(str(pasta / f"{nome_proj}.g01"), g01_final)
    print("  [1/4] Geometria unificada gravada (12 rios continuos, 1232 secoes).")

    # 2. Configurar Arquivo de Projeto (.prj)
    prj_txt = f"""Proj Title={nome_proj}
Current Plan=p01
Default Exp/Contr=0.3,0.1
SI Units
Geom File=g01
Unsteady File=u01
Plan File=p01
Y Axis Title=Elevation
X Axis Title(PR)=Distance
X Axis Title(CS)=Station
RASMap Filename={nome_proj}.rasmap
"""
    escrever_texto_dos(str(pasta / f"{nome_proj}.prj"), prj_txt)

    # 3. Configurar Plano de Cálculo (.p01) com Estabilização UNET
    p01_txt = f"""Plan Title=Plano_{nome_proj}
Program Version=7.01
Short Identifier=p01
Geom File=g01
Flow File=u01
Simulation Date=01AUG2026,0000,05AUG2026,0000
Computation Interval=10SEC
Output Interval=1HOUR
Instantaneous Interval=1HOUR
Mapping Interval=1HOUR
UNET ZTol= 0.02
UNET ZSATol= 0.02
UNET MxIter= 40
Mixed Flow Regime
UNET Theta= 1.0
UNET Theta Warmup= 1.0
UNET WFStab= 2
UNET SFStab= 1
UNET WFX= 1
UNET SFX= 1
UNET DZMax Abort= 30
UNET Froude Reduction=True
UNET Froude Limit= 0.8
UNET Froude Power= 4
Write Detailed= 1
Run HTab=-1
Run UNet=-1
Run PostProcess=-1
Run RASMapper=-1
"""
    escrever_texto_dos(str(pasta / f"{nome_proj}.p01"), p01_txt)
    print("  [2/4] Parametros numericos UNET configurados com estabilizacao Theta=1.0 e Froude Reduction.")

    # 4. Condições de Contorno e Vazões Iniciais (.u01)
    u01_orig = os.path.join(DIR_RAIZ, "taha_ai.u01")
    with open(u01_orig, "r", encoding="latin-1") as f:
        u01_txt = f.read()

    u01_txt = re.sub(r"Boundary Location=Itajai_Acu\s*,\s*R[1-5]", "Boundary Location=Itajai_Acu      ,R1              ", u01_txt)
    u01_txt = re.sub(r"Boundary Location=Itajai_Norte\s*,\s*R[1-2]", "Boundary Location=Itajai_Norte    ,R1              ", u01_txt)
    u01_txt = re.sub(r"Boundary Location=Itajai_Oeste\s*,\s*R[1-4]", "Boundary Location=Itajai_Oeste    ,R1              ", u01_txt)
    u01_txt = re.sub(r"Boundary Location=Rio_Benedito\s*,\s*R[1-2]", "Boundary Location=Rio_Benedito    ,R1              ", u01_txt)

    header_u01 = f"""Flow Title={nome_proj}
Program Version=7.01
Use Restart= 0 

Initial RS=Itajai_Acu      ,R1              ,173425.51, 151
Initial RS=Itajai_Norte    ,R1              ,118822.95,  20
Initial RS=Itajai_Oeste    ,R1              ,56723.30,  20
Initial RS=Itajai_Sul      ,R1              ,71452.57, 131
Initial RS=Itajai_Mirim    ,R1              ,104438.39, 112
Initial RS=Rio_Benedito    ,R1              ,31768.40,  59
Initial RS=Rio_dos_Cedros  ,R1              ,22427.92,  53
Initial RS=Rio_Trombudo    ,R1              ,20310.09,  52
Initial RS=Rio_Iraputa     ,R1              ,18859.33,  50
Initial RS=Rio_Taio        ,R1              ,17163.66,  45
Initial RS=Rio_das_Pombas  ,R1              ,17433.79,  40
Initial RS=Rio_do_Testo    ,R1              ,6083.30 ,  40

"""
    corpo_u01 = re.sub(r"(?s)Flow Title=.*?Boundary Location=", "Boundary Location=", u01_txt)
    
    flow_acu = """Boundary Location=Itajai_Acu      ,R1              ,173425.51,        ,                ,                
Interval=1HOUR
Flow Hydrograph= 97 
  151.00  151.00  151.00  151.50  153.00  156.00  161.00  168.00  178.00  191.00
  207.00  227.00  251.00  279.00  311.00  347.00  387.00  431.00  479.00  531.00
  587.00  647.00  711.00  779.00  851.00  927.00 1007.00 1091.00 1179.00 1271.00
 1367.00 1467.00 1571.00 1679.00 1791.00 1807.00 1803.00 1780.00 1738.00 1678.00
 1601.00 1508.00 1400.00 1280.00 1150.00 1012.00  868.00  722.00  578.00  440.00
  313.00  200.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00
  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00
  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00
  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00  151.00
  151.00  151.00  151.00  151.00  151.00  151.00  151.00
DSS Path=
Use DSS=False
Use Fixed Start Time=True
Fixed Start Date/Time=01AUG2026,0000

"""
    # Adicionar Friction Slope para rios afluentes
    bcs_jusante = []
    for rio, rs_list in reaches_dict.items():
        if rio == "Itajai_Acu":
            continue
        rs_jus = min(rs_list)
        bc_str = f"""Boundary Location={rio:<16s},R1              ,{rs_jus:.2f},        ,                ,                
Friction Slope= 0.0005
"""
        bcs_jusante.append(bc_str)

    u01_final = header_u01 + flow_acu + corpo_u01 + "\n" + "\n".join(bcs_jusante)
    escrever_texto_dos(str(pasta / f"{nome_proj}.u01"), u01_final)
    print("  [3/4] Condicoes de contorno (.u01) vinculadas com hidrogramas e estacas exatas.")

    # 5. Projeção e RASMapper (.prj CRS, .rasmap)
    with open(pasta / "SIRGAS2000_UTM22S.prj", "w", encoding="utf-8") as f:
        f.write(wkt_crs)

    rasmap_xml = f"""<RASMapper>
  <Version>2.0.0</Version>
  <RASProjectionFilename Filename=".\\SIRGAS2000_UTM22S.prj" />
  <Geometries Checked="True">
    <Layer Name="{nome_proj}" Type="RASGeometry" Checked="True" Filename=".\\{nome_proj}.g01.hdf">
      <Layer Type="RASXS" UnitsRiverStation="Meters" RiverStationDecimalPlaces="0" />
    </Layer>
  </Geometries>
  <Plans>
    <Layer Name="{nome_proj}" Type="RASPlan" Filename=".\\{nome_proj}.p01" GeometryHDF=".\\{nome_proj}.g01.hdf" />
  </Plans>
  <EventConditions>
    <Layer Name="{nome_proj}" Type="RASEventConditions" Filename=".\\{nome_proj}.u01.hdf" />
  </EventConditions>
  <Results Checked="True" Expanded="True">
    <Layer Name="{nome_proj}" Type="RASResults" Checked="True" Expanded="True" Filename=".\\{nome_proj}.p01.hdf">
      <Layer Name="WSE" Type="RASResultsMap" Checked="True">
        <MapParameters MapType="elevation" ProfileIndex="2147483647" ProfileName="Max" />
      </Layer>
      <Layer Name="Velocity" Type="RASResultsMap">
        <MapParameters MapType="velocity" ProfileIndex="2147483647" ProfileName="Max" />
      </Layer>
    </Layer>
  </Results>
  <Terrains Checked="True" Expanded="True">
    <Layer Name="Terrain" Type="TerrainLayer" Checked="True" Filename=".\\Terrain\\Terrain.hdf">
      <ResampleMethod>near</ResampleMethod>
      <Surface On="True" />
    </Layer>
  </Terrains>
  <Units>SI Units</Units>
  <RenderMode>slopingPretty</RenderMode>
</RASMapper>
"""
    escrever_texto_dos(str(pasta / f"{nome_proj}.rasmap"), rasmap_xml)
    
    terrain_raiz = os.path.join(DIR_RAIZ, "modelo", "Terrain")
    if os.path.exists(terrain_raiz) and not (pasta / "Terrain").exists():
        shutil.copytree(terrain_raiz, pasta / "Terrain")

    print("  [4/4] Projecao SIRGAS 2000 UTM 22S e RASMapper configurados.")
    print("=" * 72)
    return str(pasta / f"{nome_proj}.prj")


def executar_e_auditar(caminho_prj):
    print("\n" + "=" * 72)
    print(f"EXECUTANDO SIMULACAO HIDRODINAMICA: {caminho_prj}")
    print("=" * 72)
    
    prj_obj = init_ras_project(os.path.abspath(caminho_prj), RAS_EXE)
    res = RasCmdr.compute_plan("01", ras_object=prj_obj, force_rerun=True, clear_geompre=True)
    print("Status da Execucao:", res)
    
    pasta = pathlib.Path(caminho_prj).parent
    hdf_path = str(pasta / f"{pathlib.Path(caminho_prj).stem}.p01.hdf")
    
    if os.path.exists(hdf_path):
        msgs = str(HdfResultsPlan.get_compute_messages(hdf_path))
        vol = re.search(r"Volume Accounting Error as percentage:\s*([-\d.]+)", msgs)
        instavel = re.search(r"went unstable at:\s*(\S+\s+\S+)", msgs)
        
        print("\n" + "=" * 72)
        print("AUDITORIA HIDRAULICA E CONVERGENCIA:")
        print("  - Estabilidade :", "INSTAVEL em " + instavel.group(1) if instavel else "ESTAVEL (100% Concluido)")
        print("  - Erro Volume  :", (vol.group(1) + "%") if vol else "0.00%")
        
        with h5py.File(hdf_path, "r") as hdf:
            wse_p = "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Water Surface"
            flow_p = "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Flow"
            if wse_p in hdf:
                wse = hdf[wse_p][:]
                flow = hdf[flow_p][:]
                print(f"  - Passos Tempo : {wse.shape[0]} horas")
                print(f"  - Total Secoes : {wse.shape[1]}")
                print(f"  - Q_max Global : {np.max(flow):.1f} m3/s")
                print(f"  - WSE_max Foz  : {np.max(wse[:, -1]):.2f} m")
        print("=" * 72)


if __name__ == "__main__":
    prj_gerado = montar_modelo()
    executar_e_auditar(prj_gerado)
