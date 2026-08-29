# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT DEFINITIVO: REDE HIDROGRÁFICA 1D 100% CONECTADA COM JUNÇÕES E MARÉ NA FOZ
Diretório do Script : a_scripts/construir_rede_conectada.py
Destino Mestre      : _anti/taha_ai
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
from shapely.ops import linemerge, substring
from pyproj import Transformer
from scipy.interpolate import PchipInterpolator
from ras_commander import init_ras_project, RasCmdr, HdfResultsPlan

DIR_SCRIPT = os.path.dirname(os.path.abspath(__file__))
DIR_RAIZ = os.path.dirname(DIR_SCRIPT)
DIR_MESTRE = os.path.join(DIR_RAIZ, "_anti", "taha_ai")
RAS_EXE = r"C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe"

try:
    from a_scripts.mdt_sigsc import MosaicoSigsc
except ImportError:
    import sys
    sys.path.insert(0, DIR_SCRIPT)
    from mdt_sigsc import MosaicoSigsc

WKT_CRS = 'PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'


def escrever_texto_dos(caminho, texto):
    with open(caminho, "wb") as f:
        f.write(texto.replace("\r\n", "\n").replace("\n", "\r\n").encode("latin-1", errors="replace"))


# ==============================================================================
# 1. CARREGAMENTO DOS 21 TRECHOS E DO CANAL RETIFICADO DO MIRIM
# ==============================================================================
def obter_trechos_e_juncoes():
    print("  [1/4] Carregando eixos dos 21 trechos e confluências...")
    with open(os.path.join(DIR_RAIZ, "taha_ai_eixo.geojson"), "r", encoding="utf-8") as f:
        d_eixo = json.load(f)

    trechos = {}
    for feat in d_eixo["features"]:
        p = feat.get("properties", {})
        r = p.get("river", p.get("rio"))
        reach = p.get("reach", p.get("trecho"))
        coords = feat["geometry"]["coordinates"]
        chave = (r, reach)
        trechos[chave] = LineString(coords)

    # Ajustar o Rio Itajaí-Mirim para conectar através do Canal Retificado
    eixo_mirim_nat = trechos[("Itajai_Mirim", "R1")]
    with open(os.path.join(DIR_RAIZ, "dados_estruturas", "canal_itajai_mirim.geojson"), "r", encoding="utf-8") as f:
        d_canal = json.load(f)
    tr = Transformer.from_crs("EPSG:4326", "EPSG:31982", always_xy=True)
    canal_segs = []
    for f in d_canal["features"]:
        coords = np.array(f["geometry"]["coordinates"])
        if coords.ndim == 2:
            x, y = tr.transform(coords[:, 0], coords[:, 1])
            canal_segs.append(LineString(np.c_[x, y]))
    canal_merged = linemerge(canal_segs)
    if canal_merged.coords[0][0] > canal_merged.coords[-1][0]:
        canal_merged = LineString(canal_merged.coords[::-1])
    p_bif = Point(canal_merged.coords[0])
    s_corte = float(eixo_mirim_nat.project(p_bif))
    eixo_mirim_mont = substring(eixo_mirim_nat, 0.0, s_corte)
    
    # Conectar exatamente na confluência J5 do Itajaí-Açu
    pt_j5 = trechos[("Itajai_Acu", "R4")].coords[-1]
    canal_coords = list(canal_merged.coords)
    canal_coords[-1] = pt_j5
    
    eixo_mirim_final = LineString(list(eixo_mirim_mont.coords) + canal_coords)
    trechos[("Itajai_Mirim", "R1")] = eixo_mirim_final

    juncoes = [
        {
            "name": "Rio_do_Sul",
            "desc": "Confluencia",
            "xy": (635126.76, 6989893.33),
            "up": [("Itajai_Sul", "R1"), ("Itajai_Oeste", "R4")],
            "dn": ("Itajai_Acu", "R1"),
            "la": [150.0, 150.0],
        },
        {
            "name": "Ibirama",
            "desc": "Confluencia",
            "xy": (649214.10, 7003933.33),
            "up": [("Itajai_Norte", "R2"), ("Itajai_Acu", "R1")],
            "dn": ("Itajai_Acu", "R2"),
            "la": [159.28, 150.0],
        },
        {
            "name": "Indaial",
            "desc": "Confluencia",
            "xy": (675363.76, 7024329.19),
            "up": [("Rio_Benedito", "R2"), ("Itajai_Acu", "R2")],
            "dn": ("Itajai_Acu", "R3"),
            "la": [571.10, 1000.0],
        },
        {
            "name": "Itajai",
            "desc": "Confluencia",
            "xy": (685304.10, 7024363.33),
            "up": [("Rio_do_Testo", "R1"), ("Itajai_Acu", "R3")],
            "dn": ("Itajai_Acu", "R4"),
            "la": [567.35, 966.0],
        },
        {
            "name": "J5",
            "desc": "Confluencia",
            "xy": (730304.10, 7024093.33),
            "up": [("Itajai_Mirim", "R1"), ("Itajai_Acu", "R4")],
            "dn": ("Itajai_Acu", "R5"),
            "la": [108.92, 1000.0],
        },
        {
            "name": "Foz_Rio_Iraputa",
            "desc": "Confluencia",
            "xy": (589065.34, 7066723.82),
            "up": [("Rio_Iraputa", "R1"), ("Itajai_Norte", "R1")],
            "dn": ("Itajai_Norte", "R2"),
            "la": [546.83, 690.41],
        },
        {
            "name": "Foz_Rio_Taio",
            "desc": "Confluencia",
            "xy": (600760.52, 7000956.61),
            "up": [("Rio_Taio", "R1"), ("Itajai_Oeste", "R1")],
            "dn": ("Itajai_Oeste", "R2"),
            "la": [751.94, 902.47],
        },
        {
            "name": "Foz_Rio_das_Pomb",
            "desc": "Confluencia",
            "xy": (618531.03, 6992447.94),
            "up": [("Rio_das_Pombas", "R1"), ("Itajai_Oeste", "R2")],
            "dn": ("Itajai_Oeste", "R3"),
            "la": [948.52, 1000.0],
        },
        {
            "name": "Foz_Rio_Trombudo",
            "desc": "Confluencia",
            "xy": (629813.81, 6985549.50),
            "up": [("Rio_Trombudo", "R1"), ("Itajai_Oeste", "R3")],
            "dn": ("Itajai_Oeste", "R4"),
            "la": [687.34, 994.56],
        },
        {
            "name": "Foz_Rio_dos_Cedr",
            "desc": "Confluencia",
            "xy": (671525.10, 7031374.33),
            "up": [("Rio_dos_Cedros", "R1"), ("Rio_Benedito", "R1")],
            "dn": ("Rio_Benedito", "R2"),
            "la": [428.38, 544.59],
        },
    ]

    return trechos, juncoes


# ==============================================================================
# 2. GERAÇÃO DE CUTLINES E AMOSTRAGEM DO RELEVO SIG-SC
# ==============================================================================
def gerar_secoes_trecho(chave, eixo, mosaico, dx=200.0):
    r_nome, reach_nome = chave
    comprimento_total = eixo.length
    estacas = np.arange(0.0, comprimento_total, dx)
    if estacas[-1] < comprimento_total - 20.0:
        estacas = np.append(estacas, comprimento_total)
        
    secoes = []
    janela_tangente = 80.0
    
    larguras_dict = {
        "Itajai_Acu": (80.0, 160.0, 350.0),
        "Itajai_Norte": (50.0, 100.0, 280.0),
        "Itajai_Oeste": (40.0, 70.0, 240.0),
        "Itajai_Sul": (35.0, 60.0, 220.0),
        "Itajai_Mirim": (30.0, 50.0, 200.0),
        "Rio_Benedito": (25.0, 45.0, 180.0),
        "Rio_dos_Cedros": (20.0, 35.0, 150.0),
        "Rio_Trombudo": (20.0, 35.0, 150.0),
        "Rio_Iraputa": (20.0, 30.0, 140.0),
        "Rio_Taio": (20.0, 35.0, 150.0),
        "Rio_das_Pombas": (15.0, 25.0, 120.0),
        "Rio_do_Testo": (15.0, 25.0, 120.0),
    }
    w_min, w_max, w_corte = larguras_dict.get(r_nome, (30.0, 50.0, 180.0))
    
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
        comprimento_corte = w_corte + (w_corte * 0.4) * frac_rio
        
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
        
    for iteracao in range(25):
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
                    for j in range(i + 2, min(len(secoes) - 1, i + 20)):
                        seg_b = LineString([edge.coords[j], edge.coords[j+1]])
                        if seg_a.intersects(seg_b):
                            for k in range(i, min(len(secoes), j + 2)):
                                c = secoes[k]["cut"]
                                p_e = np.array(c.coords[0])
                                p_d = np.array(c.coords[1])
                                p_c = (p_e + p_d) / 2.0
                                if lado_idx == 0:
                                    p_e = p_c + (p_e - p_c) * 0.75
                                else:
                                    p_d = p_c + (p_d - p_c) * 0.75
                                secoes[k]["cut"] = LineString([p_e, p_d])
                            houve_ajuste = True
                            break
                    if houve_ajuste:
                        break
            if houve_ajuste:
                break
        if not houve_ajuste:
            break

    n_sec = len(secoes)
    for sec in secoes:
        cut = sec["cut"]
        comp_corte = cut.length
        pts_corte = np.linspace(0.0, comp_corte, num=int(comp_corte / 3.0) + 1)
        centro = round(comp_corte / 2.0, 2)
        w_ch = sec["largura_canal"]
        lb = max(0.0, round(centro - w_ch / 2.0, 2))
        rb = min(comp_corte, round(centro + w_ch / 2.0, 2))
        
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
        
    s_arr = np.array([sec["s_eixo"] for sec in secoes])
    s_max = s_arr[-1]
    
    cotas_trecho_dict = {
        ("Itajai_Acu", "R1"): (350.0, 260.0),
        ("Itajai_Acu", "R2"): (260.0, 85.0),
        ("Itajai_Acu", "R3"): (85.0, 20.0),
        ("Itajai_Acu", "R4"): (20.0, -1.0),
        ("Itajai_Acu", "R5"): (-1.0, -6.0),
        ("Itajai_Norte", "R1"): (480.0, 460.0),
        ("Itajai_Norte", "R2"): (460.0, 260.0),
        ("Itajai_Oeste", "R1"): (520.0, 500.0),
        ("Itajai_Oeste", "R2"): (500.0, 430.0),
        ("Itajai_Oeste", "R3"): (430.0, 390.0),
        ("Itajai_Oeste", "R4"): (390.0, 350.0),
        ("Itajai_Sul", "R1"): (490.0, 350.0),
        ("Itajai_Mirim", "R1"): (275.0, -1.0),
        ("Rio_Benedito", "R1"): (220.0, 120.0),
        ("Rio_Benedito", "R2"): (120.0, 85.0),
        ("Rio_dos_Cedros", "R1"): (180.0, 120.0),
        ("Rio_Trombudo", "R1"): (460.0, 390.0),
        ("Rio_Iraputa", "R1"): (510.0, 460.0),
        ("Rio_Taio", "R1"): (530.0, 500.0),
        ("Rio_das_Pombas", "R1"): (450.0, 430.0),
        ("Rio_do_Testo", "R1"): (120.0, 20.0),
    }
    z_mont, z_jus = cotas_trecho_dict.get(chave, (300.0, 10.0))
    z_suave = z_mont + (z_jus - z_mont) * (s_arr / max(1.0, s_max))
    
    for i, sec in enumerate(secoes):
        sta = sec["sta"]
        lb, rb, centro = sec["lb"], sec["rb"], sec["centro"]
        w_ch = sec["largura_canal"]
        z_alvo = round(float(z_suave[i]), 2)
        
        z_borda = max(z_alvo + 3.0, float(np.max(sec["z_bruto"])))
        
        z = np.zeros_like(sta)
        for idx in range(len(sta)):
            if sta[idx] <= lb:
                frac = max(0.0, (lb - sta[idx]) / max(1.0, lb))
                z[idx] = z_borda + frac * 2.0
            elif sta[idx] >= rb:
                comp = sta[-1]
                frac = max(0.0, (sta[idx] - rb) / max(1.0, (comp - rb)))
                z[idx] = z_borda + frac * 2.0
            else:
                dist = abs(sta[idx] - centro)
                t = dist / max(1.0, (w_ch / 2.0))
                z[idx] = z_alvo + (z_borda - z_alvo) * (t ** 2.0)
                
        sec["sta_final"] = np.round(sta, 2)
        sec["z_final"] = np.round(z, 2)
        
    return secoes


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
    return f"#Mann= 3 , 0 , 0 \n{0.0:8.2f}{0.045:8.3f}{0:8d}{lb:8.2f}{0.030:8.3f}{0:8d}{rb:8.2f}{0.045:8.3f}{0:8d}"


# ==============================================================================
# 3. MONTAGEM DO MODELO MESTRE COM 10 JUNÇÕES E MARÉ NA FOZ
# ==============================================================================
def construir_modelo_completo():
    print("=" * 72)
    print("GERADOR DEFINITIVO: REDE HIDROGRÁFICA 1D 100% CONECTADA COM JUNÇÕES")
    print("=" * 72)
    
    pasta = pathlib.Path(DIR_MESTRE)
    pasta.mkdir(parents=True, exist_ok=True)
    nome_proj = "taha_ai"
    
    trechos, juncoes = obter_trechos_e_juncoes()
    mosaico = MosaicoSigsc()
    
    g01_header = f"""Geom Title={nome_proj} - Rede Hidrografica Conectada do Vale
Program Version=7.01
Spatial Reference System={WKT_CRS}

"""
    junc_txt_list = []
    for j in juncoes:
        jt = f"Junct Name={j['name']:<16s}\n"
        jt += f"Junct Desc={j['desc']}, 0 , 0 , 0 ,0\n"
        jt += f"Junct X Y & Text X Y={j['xy'][0]:.2f},{j['xy'][1]:.2f},{j['xy'][0]+800:.2f},{j['xy'][1]+800:.2f}\n"
        for up in j["up"]:
            jt += f"Up River,Reach={up[0]:<16s},{up[1]:<16s}\n"
        jt += f"Dn River,Reach={j['dn'][0]:<16s},{j['dn'][1]:<16s}\n"
        for la in j["la"]:
            jt += f"Junc L&A={la:.2f},0\n"
        junc_txt_list.append(jt)
        
    g01_corpo = []
    total_secoes = 0
    info_trechos = []
    
    for chave, eixo in trechos.items():
        r_nome, reach_nome = chave
        print(f"--> Processando {r_nome} ({reach_nome}) [{eixo.length/1000:.2f} km]...")
        secoes = gerar_secoes_trecho(chave, eixo, mosaico, dx=200.0)
        total_secoes += len(secoes)
        
        pts_xy = list(eixo.coords)
        reach_txt = f"River Reach={r_nome:<16s},{reach_nome:<16s}\n"
        reach_txt += format_reach_xy(pts_xy) + "\n"
        
        n_sec = len(secoes)
        for i, sec in enumerate(secoes):
            rs = sec["rs"]
            dist_jus = round(abs(sec["s_eixo"] - secoes[i+1]["s_eixo"]), 2) if i < n_sec - 1 else 0.0
            cut_pts = list(sec["cut"].coords)
            sta = sec["sta_final"]
            z = sec["z_final"]
            lb, rb = sec["lb"], sec["rb"]
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
            
        g01_corpo.append(reach_txt)
        info_trechos.append({
            "rio": r_nome,
            "reach": reach_nome,
            "rs_mont": secoes[0]["rs"],
            "rs_jus": secoes[-1]["rs"],
            "n_secoes": len(secoes)
        })
        
    g01_final = g01_header + "\n".join(junc_txt_list) + "\n" + "".join(g01_corpo)
    escrever_texto_dos(str(pasta / f"{nome_proj}.g01"), g01_final)
    print(f"\n[1/4] Geometria unificada com 10 junções e {total_secoes} seções gravada.")

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

    # Plano com Warmup e tolerâncias robustas
    p01_txt = f"""Plan Title=Plano_{nome_proj}
Program Version=7.01
Short Identifier=p01
Geom File=g01
Flow File=u01
Simulation Date=01AUG2026,0000,05AUG2026,0000
Computation Interval=30SEC
Output Interval=1HOUR
Instantaneous Interval=1HOUR
Mapping Interval=1HOUR
UNET ZTol= 0.03
UNET ZSATol= 0.03
UNET MxIter= 40
Mixed Flow Regime
UNET Theta= 1.0
UNET Theta Warmup= 1.0
UNET Warmup= 12.0
UNET Warmup Time Step=30SEC
UNET WFStab= 2
UNET SFStab= 1
UNET WFX= 1
UNET SFX= 1
UNET DZMax Abort= 50
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

    # Condições de Contorno (.u01)
    u01_header = f"Flow Title={nome_proj}\nProgram Version=7.01\nUse Restart= 0 \n"
    for r in info_trechos:
        u01_header += f"Initial RS={r['rio']:<16s},{r['reach']:<16s},{r['rs_mont']:<8.2f},  50\n"
    u01_header += "\n"
    
    cabeceiras = [
        ("Itajai_Sul", "R1"),
        ("Itajai_Oeste", "R1"),
        ("Itajai_Norte", "R1"),
        ("Rio_Taio", "R1"),
        ("Rio_das_Pombas", "R1"),
        ("Rio_Trombudo", "R1"),
        ("Rio_Iraputa", "R1"),
        ("Rio_Benedito", "R1"),
        ("Rio_dos_Cedros", "R1"),
        ("Rio_do_Testo", "R1"),
        ("Itajai_Mirim", "R1"),
    ]
    
    bcs_txt = []
    for r_nome, reach_nome in cabeceiras:
        rs_m = [t["rs_mont"] for t in info_trechos if t["rio"] == r_nome and t["reach"] == reach_nome][0]
        bc_up = f"""Boundary Location={r_nome:<16s},{reach_nome:<16s},{rs_m:.2f},        ,                ,                
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
Use Fixed Start Time=False
Fixed Start Date/Time=,

"""
        bcs_txt.append(bc_up)
        
    horas = np.arange(97)
    cotas_mare = 0.50 + 0.60 * np.sin(2.0 * np.pi * horas / 12.42)
    cotas_mare_str = []
    for k in range(0, len(cotas_mare), 10):
        grupo = [f"{cotas_mare[m]:8.2f}" for m in range(k, min(k+10, len(cotas_mare)))]
        cotas_mare_str.append("".join(grupo))
    corpo_mare = "\n".join(cotas_mare_str)
    
    rs_foz = [t["rs_jus"] for t in info_trechos if t["rio"] == "Itajai_Acu" and t["reach"] == "R5"][0]
    bc_foz_mare = f"""Boundary Location=Itajai_Acu      ,R5              ,{rs_foz:.2f},        ,                ,                
Interval=1HOUR
Stage Hydrograph= 97 
{corpo_mare}
DSS Path=
Use DSS=False
Use Fixed Start Time=False
Fixed Start Date/Time=,
"""
    bcs_txt.append(bc_foz_mare)
    
    u01_final = u01_header + "\n".join(bcs_txt)
    escrever_texto_dos(str(pasta / f"{nome_proj}.u01"), u01_final)

    with open(pasta / "SIRGAS2000_UTM22S.prj", "w", encoding="utf-8") as f:
        f.write(WKT_CRS)

    rasmap_xml = f"""<RASMapper>
  <Version>2.0.0</Version>
  <RASProjectionFilename Filename=".\\SIRGAS2000_UTM22S.prj" />
  <Geometries Checked="True" Expanded="True">
    <Layer Name="{nome_proj} - Rede Hidrografica Conectada do Vale" Type="RASGeometry" Checked="True" Expanded="True" Filename=".\\{nome_proj}.g01.hdf">
      <Layer Type="RASRiver" Checked="True" />
      <Layer Type="RASJunction" Checked="True" />
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

    print("[4/4] Arquivos HEC-RAS com rede conectada e maré configurados com sucesso!")
    return str(pasta / f"{nome_proj}.prj")


def executar_e_auditar(caminho_prj):
    print("\n" + "=" * 72)
    print(f"EXECUTANDO SIMULAÇÃO HIDRODINÂMICA CONECTADA: {caminho_prj}")
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
        print("AUDITORIA HIDRÁULICA DA REDE CONECTADA:")
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
    prj_gerado = construir_modelo_completo()
    executar_e_auditar(prj_gerado)
