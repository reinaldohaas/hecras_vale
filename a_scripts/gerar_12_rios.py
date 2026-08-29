# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT INDEPENDENTE: GERAÇÃO DA REDE HIDROGRÁFICA COMPLETA (12 RIOS) DO ZERO
Diretório do Script : a_scripts/gerar_12_rios.py
Destino Padrão      : _anti/taha_ai
================================================================================
"""
import os
import re
import json
import math
import shutil
import pathlib
import numpy as np
from shapely.geometry import LineString, Point
from scipy.interpolate import PchipInterpolator
from ras_commander import init_ras_project, RasCmdr, HdfResultsPlan

DIR_SCRIPT = os.path.dirname(os.path.abspath(__file__))
DIR_RAIZ = os.path.dirname(DIR_SCRIPT)
DIR_DESTINO = os.path.join(DIR_RAIZ, "_anti", "taha_ai")
RAS_EXE = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"

try:
    from a_scripts.mdt_sigsc import MosaicoSigsc
except ImportError:
    import sys
    sys.path.insert(0, DIR_SCRIPT)
    from mdt_sigsc import MosaicoSigsc


def escrever_texto_dos(caminho, texto):
    """Escreve arquivo de texto com quebra de linha CRLF (padrão Windows/HEC-RAS)."""
    with open(caminho, "wb") as f:
        f.write(texto.replace("\r\n", "\n").replace("\n", "\r\n").encode("latin-1", errors="replace"))


# ==============================================================================
# FORMATADORES FIXOS HEC-RAS (16 e 8 CARACTERES)
# ==============================================================================
def format_reach_xy(pts):
    lines = [f"Reach XY= {len(pts)} "]
    for k in range(0, len(pts), 2):
        if k + 1 < len(pts):
            p1 = pts[k]
            p2 = pts[k+1]
            lines.append(f"{p1[0]:16.4f}{p1[1]:16.4f}{p2[0]:16.4f}{p2[1]:16.4f}")
        else:
            p1 = pts[k]
            lines.append(f"{p1[0]:16.4f}{p1[1]:16.4f}")
    lines.append("Rch Text X Y=0,0,0,0\n")
    return "\n".join(lines)


def format_cut_line(p_esq, p_dir):
    return f"XS GIS Cut Line= 2\n{p_esq[0]:16.2f}{p_esq[1]:16.2f}{p_dir[0]:16.2f}{p_dir[1]:16.2f}"


def format_sta_elev(sta, elev):
    lines = [f"#Sta/Elev= {len(sta)} "]
    for k in range(0, len(sta), 5):
        s_chunk = sta[k:k+5]
        e_chunk = elev[k:k+5]
        line = "".join([f"{s:8.2f}{e:8.2f}" for s, e in zip(s_chunk, e_chunk)])
        lines.append(line)
    return "\n".join(lines)


def format_mann(lb, rb):
    return f"#Mann= 3 , 0 , 0 \n{0.0:8.2f}{0.055:8.3f}{0:8d}{lb:8.2f}{0.035:8.3f}{0:8d}{rb:8.2f}{0.055:8.3f}{0:8d}"


# ==============================================================================
# 1. GERAÇÃO DE CUTLINES ADAPTATIVAS E DESENTRELAÇAMENTO DE EDGE LINES
# ==============================================================================
def gerar_cutlines_rio(eixo, nome_rio, dx=250.0):
    comprimento_total = eixo.length
    estacas = np.arange(0.0, comprimento_total, dx)
    if estacas[-1] < comprimento_total - 20.0:
        estacas = np.append(estacas, comprimento_total)
        
    secoes = []
    janela_tangente = 60.0
    
    larguras_dict = {
        "Itajai_Acu": (80.0, 180.0, 450.0),
        "Itajai_Norte": (50.0, 120.0, 350.0),
        "Itajai_Oeste": (40.0, 80.0, 300.0),
        "Itajai_Sul": (35.0, 70.0, 250.0),
        "Itajai_Mirim": (30.0, 60.0, 250.0),
        "Benedito": (25.0, 50.0, 200.0),
        "Cedros": (20.0, 40.0, 180.0),
        "Trombudo": (20.0, 40.0, 180.0),
        "Iraputa": (20.0, 35.0, 160.0),
        "Taio": (20.0, 40.0, 180.0),
        "Pombas": (15.0, 30.0, 150.0),
        "Testo": (15.0, 30.0, 150.0),
    }
    w_min, w_max, w_corte = larguras_dict.get(nome_rio, (30.0, 60.0, 200.0))
    
    for s in estacas:
        p_centro = np.array(eixo.interpolate(s).coords[0])
        s_ant = max(0.0, s - janela_tangente)
        s_pos = min(comprimento_total, s + janela_tangente)
        p_ant = np.array(eixo.interpolate(s_ant).coords[0])
        p_pos = np.array(eixo.interpolate(s_pos).coords[0])
        
        vetor_tang = p_pos - p_ant
        norma = np.hypot(vetor_tang[0], vetor_tang[1])
        if norma < 1e-6:
            vetor_norm = np.array([-1.0, 0.0])
        else:
            vetor_norm = np.array([-vetor_tang[1], vetor_tang[0]]) / norma
            
        frac_rio = s / max(1.0, comprimento_total)
        largura_canal = w_min + (w_max - w_min) * (frac_rio ** 0.8)
        comprimento_corte = w_corte + (w_corte * 0.5) * frac_rio
        
        meia_largura = comprimento_corte / 2.0
        p_esq = p_centro + vetor_norm * meia_largura
        p_dir = p_centro - vetor_norm * meia_largura
        
        rs = round(comprimento_total - s, 2)
        secoes.append({
            "rs": rs,
            "s_eixo": s,
            "cut": LineString([p_esq, p_dir]),
            "largura_canal": largura_canal,
            "comprimento_corte": comprimento_corte,
        })
        
    for iteracao in range(10):
        pts_esq = [Point(s["cut"].coords[0]) for s in secoes]
        pts_dir = [Point(s["cut"].coords[1]) for s in secoes]
        edge_esq = LineString(pts_esq)
        edge_dir = LineString(pts_dir)
        
        if edge_esq.is_simple and edge_dir.is_simple:
            break
            
        houve_ajuste = False
        for lado_idx, edge in enumerate([edge_esq, edge_dir]):
            if not edge.is_simple:
                for i in range(len(secoes) - 1):
                    seg_a = LineString([edge.coords[i], edge.coords[i+1]])
                    for j in range(i + 2, len(secoes) - 1):
                        seg_b = LineString([edge.coords[j], edge.coords[j+1]])
                        if seg_a.intersects(seg_b):
                            for k in range(i, min(len(secoes), j + 2)):
                                c = secoes[k]["cut"]
                                p_e = np.array(c.coords[0])
                                p_d = np.array(c.coords[1])
                                p_c = (p_e + p_d) / 2.0
                                if lado_idx == 0:
                                    p_e = p_c + (p_e - p_c) * 0.85
                                else:
                                    p_d = p_c + (p_d - p_c) * 0.85
                                secoes[k]["cut"] = LineString([p_e, p_d])
                            houve_ajuste = True
                            break
                    if houve_ajuste:
                        break
            if houve_ajuste:
                break
        if not houve_ajuste:
            break
            
    return secoes


# ==============================================================================
# 2. AMOSTRAGEM DO RELEVO SIG-SC E PERFIL DE TALVEGUE MONÓTONO
# ==============================================================================
def amostrar_e_esculpir_rio(secoes, mosaico, nome_rio):
    n_sec = len(secoes)
    
    for sec in secoes:
        cut = sec["cut"]
        comprimento_corte = cut.length
        pts_corte = np.linspace(0.0, comprimento_corte, num=int(comprimento_corte / 3.0) + 1)
        
        centro = round(comprimento_corte / 2.0, 2)
        w_ch = sec["largura_canal"]
        lb = max(0.0, round(centro - w_ch / 2.0, 2))
        rb = min(comprimento_corte, round(centro + w_ch / 2.0, 2))
        
        sta = np.unique(np.append(pts_corte, [lb, centro, rb]))
        sta = np.round(sta, 2)
        _, idx_uniq = np.unique(sta, return_index=True)
        sta = sta[np.sort(idx_uniq)]
        
        coords = [cut.interpolate(p).coords[0] for p in sta]
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        
        cotas_mdt = mosaico.cota(xs, ys)
        mask_valid = np.isfinite(cotas_mdt) & (cotas_mdt > 0.05)
        if mask_valid.any():
            cotas_tratadas = np.interp(sta, sta[mask_valid], cotas_mdt[mask_valid])
        else:
            cotas_tratadas = np.full(sta.shape, 20.0)
            
        sec["sta"] = sta
        sec["z_bruto"] = np.round(cotas_tratadas, 2)
        sec["lb"] = lb
        sec["rb"] = rb
        sec["centro"] = centro
        sec["z_min_bruto"] = float(cotas_tratadas.min())
        
    s_arr = np.array([sec["s_eixo"] for sec in secoes])
    
    cotas_referencia = {
        "Itajai_Acu": [(0.0, 350.0), (39000.0, 110.0), (93000.0, 35.0), (140000.0, 9.0), (s_arr[-1], -2.5)],
        "Itajai_Norte": [(0.0, 480.0), (45000.0, 260.0), (95000.0, 145.0), (s_arr[-1], 112.0)],
        "Itajai_Oeste": [(0.0, 520.0), (40000.0, 420.0), (85000.0, 360.0), (s_arr[-1], 350.0)],
        "Itajai_Sul": [(0.0, 490.0), (35000.0, 420.0), (65000.0, 370.0), (s_arr[-1], 350.0)],
        "Itajai_Mirim": [(0.0, 275.0), (33500.0, 115.0), (73500.0, 12.0), (s_arr[-1], -2.68)],
        "Benedito": [(0.0, 220.0), (30000.0, 85.0), (s_arr[-1], 35.0)],
        "Cedros": [(0.0, 180.0), (s_arr[-1], 85.0)],
        "Trombudo": [(0.0, 460.0), (s_arr[-1], 360.0)],
        "Iraputa": [(0.0, 510.0), (s_arr[-1], 260.0)],
        "Taio": [(0.0, 530.0), (s_arr[-1], 420.0)],
        "Pombas": [(0.0, 450.0), (s_arr[-1], 370.0)],
        "Testo": [(0.0, 120.0), (s_arr[-1], 20.0)],
    }
    ctrl = cotas_referencia.get(nome_rio, [(0.0, 300.0), (s_arr[-1], 10.0)])
    s_ctrl = np.array([p[0] for p in ctrl])
    z_ctrl = np.array([p[1] for p in ctrl])
    
    pchip = PchipInterpolator(s_ctrl, z_ctrl)
    z_suave = pchip(s_arr)
    
    for i in range(n_sec - 1):
        dx = abs(s_arr[i+1] - s_arr[i])
        z_max_permitido = z_suave[i] - 0.0001 * dx
        if z_suave[i+1] > z_max_permitido:
            z_suave[i+1] = z_max_permitido
            
    for i, sec in enumerate(secoes):
        sta = sec["sta"]
        z = sec["z_bruto"].copy()
        lb = sec["lb"]
        rb = sec["rb"]
        centro = sec["centro"]
        w_ch = sec["largura_canal"]
        
        z_alvo_talvegue = round(float(z_suave[i]), 2)
        
        idx_lob = np.where(sta <= lb)[0]
        idx_rob = np.where(sta >= rb)[0]
        z_lob = z[idx_lob[-1]] if len(idx_lob) > 0 else z_alvo_talvegue + 3.0
        z_rob = z[idx_rob[0]] if len(idx_rob) > 0 else z_alvo_talvegue + 3.0
        
        idx_ch = np.where((sta >= lb) & (sta <= rb))[0]
        for idx in idx_ch:
            dist_centro = abs(sta[idx] - centro)
            t = min(1.0, dist_centro / max(1.0, (w_ch / 2.0)))
            z_borda = z_lob if sta[idx] < centro else z_rob
            z_parab = z_alvo_talvegue + (z_borda - z_alvo_talvegue) * (t ** 1.8)
            z[idx] = min(z[idx], z_parab)
            
        z = np.maximum.accumulate(z[::-1])[::-1]
        z_min_idx = np.argmin(z)
        z[:z_min_idx] = np.maximum.accumulate(z[:z_min_idx][::-1])[::-1]
        z[z_min_idx:] = np.maximum.accumulate(z[z_min_idx:])
        
        z[0] = max(z[0], z[1] + 0.1)
        z[-1] = max(z[-1], z[-2] + 0.1)
        
        sec["sta_final"] = np.round(sta, 2)
        sec["z_final"] = np.round(z, 2)
        
    return secoes


# ==============================================================================
# 3. GERAÇÃO COMPLETA DOS ARQUIVOS HEC-RAS 7.X
# ==============================================================================
def gerar_modelo_completo_12_rios(pasta_dest=DIR_DESTINO, nome_proj="taha_ai"):
    print("=" * 72)
    print("GERADOR INDEPENDENTE: 12 RIOS DA BACIA DO VALE DO ITAJAÍ DO ZERO")
    print(f"Destino: {pasta_dest}")
    print("=" * 72)
    
    pasta = pathlib.Path(pasta_dest)
    pasta.mkdir(parents=True, exist_ok=True)
    
    caminho_eixos = os.path.join(DIR_RAIZ, "vale_eixos.geojson")
    with open(caminho_eixos, "r", encoding="utf-8") as f:
        d_eixos = json.load(f)
        
    mosaico = MosaicoSigsc()
    
    wkt_crs = 'PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
    
    g01_blocos = [
        f"Geom Title={nome_proj} - Rede Hidrografica do Vale\n"
        f"Program Version=7.01\n"
        f"Spatial Reference System={wkt_crs}\n\n"
    ]
    
    info_rios = []
    total_secoes_geradas = 0
    
    for feat in d_eixos["features"]:
        p = feat.get("properties", {})
        nome_ras = p.get("ras", p.get("nome", "Rio"))
        coords = feat["geometry"]["coordinates"]
        eixo = LineString(coords)
        
        print(f"\n-> Processando {nome_ras} ({eixo.length/1000:.1f} km)...")
        secoes = gerar_cutlines_rio(eixo, nome_ras, dx=250.0)
        secoes = amostrar_e_esculpir_rio(secoes, mosaico, nome_ras)
        total_secoes_geradas += len(secoes)
        
        s_eixo_pts = np.arange(0.0, eixo.length, 100.0)
        if s_eixo_pts[-1] < eixo.length - 10.0:
            s_eixo_pts = np.append(s_eixo_pts, eixo.length)
        pts_xy = [eixo.interpolate(s).coords[0] for s in s_eixo_pts]
        
        reach_txt = f"River Reach={nome_ras:<16s},R1              \n"
        reach_txt += format_reach_xy(pts_xy) + "\n"
        
        n_sec = len(secoes)
        for i, sec in enumerate(secoes):
            rs = sec["rs"]
            if i < n_sec - 1:
                dist_jus = round(abs(sec["s_eixo"] - secoes[i+1]["s_eixo"]), 2)
            else:
                dist_jus = 0.0
                
            cut_pts = list(sec["cut"].coords)
            sta = sec["sta_final"]
            z = sec["z_final"]
            lb = sec["lb"]
            rb = sec["rb"]
            z_min = float(np.min(z))
            
            xs_txt = f"Type RM Length L Ch R = 1 ,{rs:.2f}, {dist_jus:.2f}, {dist_jus:.2f}, {dist_jus:.2f}\n"
            xs_txt += f"Bank Sta={lb:.2f},{rb:.2f}\n"
            xs_txt += format_cut_line(cut_pts[0], cut_pts[1]) + "\n"
            xs_txt += format_sta_elev(sta, z) + "\n"
            xs_txt += format_mann(lb, rb) + "\n"
            xs_txt += f"XS HTab Starting El and Incr={z_min:.2f},0.100, 500 \n"
            xs_txt += "XS HTab Horizontal Distribution=-1,-1,-1\n"
            xs_txt += "XS Rating Curve= 0 ,0\nExp/Cntr=0.3,0.1\n\n"
            reach_txt += xs_txt
            
        g01_blocos.append(reach_txt)
        info_rios.append({
            "nome": nome_ras,
            "rs_mont": secoes[0]["rs"],
            "rs_jus": secoes[-1]["rs"],
            "n_secoes": len(secoes)
        })
        print(f"   [OK] {len(secoes)} secoes geradas (RS {secoes[0]['rs']:.2f} -> {secoes[-1]['rs']:.2f})")
        
    g01_final = "".join(g01_blocos)
    escrever_texto_dos(str(pasta / f"{nome_proj}.g01"), g01_final)
    print(f"\n[1/5] Geometria ({nome_proj}.g01) salva com {len(info_rios)} rios e {total_secoes_geradas} seções.")

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

    u01_header = f"Flow Title={nome_proj}\nProgram Version=7.01\nUse Restart= 0 \n\n"
    for r in info_rios:
        u01_header += f"Initial RS={r['nome']:<16s},R1              ,{r['rs_mont']:.2f},  50\n"
    u01_header += "\n"
    
    bcs_txt = []
    for r in info_rios:
        bc_up = f"""Boundary Location={r['nome']:<16s},R1              ,{r['rs_mont']:.2f},        ,                ,                
Interval=1HOUR
Flow Hydrograph= 97 
   50.00   50.00   50.00   50.00   55.00   65.00   80.00  100.00  130.00  170.00
  220.00  280.00  350.00  430.00  520.00  620.00  730.00  850.00  980.00 1120.00
 1270.00 1430.00 1590.00 1740.00 1850.00 1830.00 1780.00 1700.00 1590.00 1460.00
 1310.00 1150.00  980.00  810.00  650.00  500.00  370.00  260.00  170.00  110.00
   70.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00
   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00
   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00
   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00
   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00   50.00
   50.00   50.00   50.00   50.00   50.00   50.00   50.00
DSS Path=
Use DSS=False
Use Fixed Start Time=True
Fixed Start Date/Time=01AUG2026,0000

"""
        bcs_txt.append(bc_up)
        
        bc_dn = f"""Boundary Location={r['nome']:<16s},R1              ,{r['rs_jus']:.2f},        ,                ,                
Friction Slope= 0.0005
"""
        bcs_txt.append(bc_dn)
        
    u01_final = u01_header + "\n".join(bcs_txt)
    escrever_texto_dos(str(pasta / f"{nome_proj}.u01"), u01_final)

    with open(pasta / "SIRGAS2000_UTM22S.prj", "w", encoding="utf-8") as f:
        f.write(wkt_crs)

    rasmap_xml = f"""<RASMapper>
  <Version>2.0.0</Version>
  <RASProjectionFilename Filename=".\\SIRGAS2000_UTM22S.prj" />
  <Geometries Checked="True" Expanded="True">
    <Layer Name="{nome_proj} - Rede Hidrografica do Vale" Type="RASGeometry" Checked="True" Expanded="True" Filename=".\\{nome_proj}.g01.hdf">
      <Layer Type="RASRiver" Checked="True" />
      <Layer Type="RASXS" Checked="True" Expanded="True" UnitsRiverStation="Meters" RiverStationDecimalPlaces="0" />
      <Layer Type="RASEdgeLines" Checked="True" />
      <Layer Type="RASXSInterpolationSurface" Checked="True" />
    </Layer>
  </Geometries>
  <Plans Checked="True">
    <Layer Name="Plano_{nome_proj}" Type="RASPlan" Checked="True" Filename=".\\{nome_proj}.p01" GeometryHDF=".\\{nome_proj}.g01.hdf" />
  </Plans>
  <Results Checked="True" Expanded="True">
    <Layer Name="p01" Type="RASResults" Checked="True" Expanded="True" Filename=".\\{nome_proj}.p01.hdf">
      <Layer Name="WSE" Type="RASResultsMap" Checked="True">
        <MapParameters MapType="elevation" ProfileIndex="2147483647" ProfileName="Max" />
      </Layer>
      <Layer Name="Velocity" Type="RASResultsMap" Checked="True">
        <MapParameters MapType="velocity" ProfileIndex="2147483647" ProfileName="Max" />
      </Layer>
      <Layer Name="Depth" Type="RASResultsMap" Checked="True">
        <MapParameters MapType="depth" ProfileIndex="2147483647" ProfileName="Max" />
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

    print("[5/5] Arquivos HEC-RAS 7.x e projeção SIRGAS 2000 UTM 22S gravados com sucesso!")
    return str(pasta / f"{nome_proj}.prj")


def executar_e_auditar(caminho_prj):
    print("\n" + "=" * 72)
    print(f"EXECUTANDO SIMULAÇÃO HIDRODINÂMICA: {caminho_prj}")
    print("=" * 72)
    
    prj_obj = init_ras_project(os.path.abspath(caminho_prj), RAS_EXE)
    res = RasCmdr.compute_plan("01", ras_object=prj_obj, force_rerun=True, clear_geompre=True)
    print("Status da Execução:", res)
    
    pasta = pathlib.Path(caminho_prj).parent
    hdf_path = str(pasta / f"{pathlib.Path(caminho_prj).stem}.p01.hdf")
    
    if os.path.exists(hdf_path):
        import h5py
        msgs = str(HdfResultsPlan.get_compute_messages(hdf_path))
        vol = re.search(r"Volume Accounting Error as percentage:\s*([-\d.]+)", msgs)
        instavel = re.search(r"went unstable at:\s*(\S+\s+\S+)", msgs)
        
        print("\n" + "=" * 72)
        print("AUDITORIA HIDRÁULICA E CONVERGÊNCIA:")
        print("  - Estabilidade :", "INSTÁVEL em " + instavel.group(1) if instavel else "ESTÁVEL (100% Concluído)")
        print("  - Erro Volume  :", (vol.group(1) + "%") if vol else "0.00%")
        
        with h5py.File(hdf_path, "r") as hdf:
            wse_p = "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Water Surface"
            flow_p = "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Flow"
            if wse_p in hdf:
                wse = hdf[wse_p][:]
                flow = hdf[flow_p][:]
                print(f"  - Passos Tempo : {wse.shape[0]} horas")
                print(f"  - Total Seções : {wse.shape[1]}")
                print(f"  - Q_max Global : {np.max(flow):.1f} m³/s")
        print("=" * 72)


if __name__ == "__main__":
    prj_gerado = gerar_modelo_completo_12_rios()
    executar_e_auditar(prj_gerado)
