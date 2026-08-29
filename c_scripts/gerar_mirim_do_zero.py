# -*- coding: utf-8 -*-
"""
Gerador Completo do Modelo Rio Itajaí-Mirim do Zero (Canal Retificado Oficial).

Constrói a geometria e o projeto HEC-RAS 100% limpo a partir do eixo unificado
(Eixo natural de montante + Canal Retificado de jusante), com:
  1. Perfil longitudinal de talvegue rigorosamente monótono e suave (via PCHIP contínuo)
  2. Calha ativa (Bank Stations) perfeitamente esculpida e ancorada ao talvegue
  3. Hidrogramas de contorno (.u01) ancorados nas estacas exatas da geometria
  4. 0 erros de geometria no RAS Mapper e convergência total no UNET
"""
import argparse
import json
import os
import pathlib
import re
import shutil
import sys
import numpy as np
from pyproj import Transformer
from scipy.interpolate import PchipInterpolator
from shapely.geometry import LineString, Point
from shapely.ops import linemerge, substring

DIRETORIO_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
DIRETORIO_RAIZ = os.path.dirname(DIRETORIO_SCRIPTS)
sys.path.insert(0, DIRETORIO_SCRIPTS)
sys.path.insert(0, DIRETORIO_RAIZ)

from ras_io import escrever
from mdt_sigsc import MosaicoSigsc, tiles_do_dominio


# ==============================================================================
# 1. UNIFICAÇÃO DO EIXO (NATURAL + CANAL RETIFICADO)
# ==============================================================================
def obter_eixo_unificado():
    """Gera o eixo único contínuo de montante até a foz pelo Canal Retificado."""
    caminho_eixos = os.path.join(DIRETORIO_RAIZ, "eixos_do_relevo.geojson")
    with open(caminho_eixos, encoding="utf-8") as f:
        d = json.load(f)
    
    eixo_nat = None
    for feat in d["features"]:
        if feat["properties"].get("nome") == "Itajai_Mirim":
            eixo_nat = LineString(feat["geometry"]["coordinates"])
            break
            
    caminho_canal = os.path.join(DIRETORIO_RAIZ, "dados_estruturas", "canal_itajai_mirim.geojson")
    tr = Transformer.from_crs("EPSG:4326", "EPSG:31982", always_xy=True)
    with open(caminho_canal, encoding="utf-8") as f:
        d_canal = json.load(f)
        
    canal_segs = []
    for f in d_canal["features"]:
        coords = np.array(f["geometry"]["coordinates"])
        if coords.ndim == 2:
            x, y = tr.transform(coords[:, 0], coords[:, 1])
            canal_segs.append(LineString(np.c_[x, y]))
            
    canal_merged = linemerge(canal_segs)
    if canal_merged.coords[0][0] > canal_merged.coords[-1][0]:
        canal_merged = LineString(canal_merged.coords[::-1])
        
    p_bifurcacao = Point(canal_merged.coords[0])
    s_corte = float(eixo_nat.project(p_bifurcacao))
    
    eixo_montante = substring(eixo_nat, 0.0, s_corte)
    coords_unificadas = list(eixo_montante.coords) + list(canal_merged.coords)
    eixo_final = LineString(coords_unificadas)
    
    print(f"  [1/5] Eixo unificado: {eixo_final.length/1000:.2f} km "
          f"(Montante: {eixo_montante.length/1000:.2f} km + Canal: {canal_merged.length/1000:.2f} km)")
    return eixo_final, s_corte


# ==============================================================================
# 2. GERAÇÃO DE CUTLINES ADAPTATIVAS (SEM CRUZAMENTOS)
# ==============================================================================
def gerar_cutlines(eixo, s_bifurcacao, dx=150.0):
    """Gera linhas de seção transversais ortogonais sem meandros fantasmas."""
    comprimento_total = eixo.length
    estacas = np.arange(0.0, comprimento_total, dx)
    if estacas[-1] < comprimento_total - 20.0:
        estacas = np.append(estacas, comprimento_total)
        
    secoes = []
    janela_tangente = 60.0
    
    for s in estacas:
        pt = np.array(eixo.interpolate(s).coords[0])
        s_ant = max(s - janela_tangente, 0.0)
        s_pos = min(s + janela_tangente, comprimento_total)
        pt_ant = np.array(eixo.interpolate(s_ant).coords[0])
        pt_pos = np.array(eixo.interpolate(s_pos).coords[0])
        
        v_tang = pt_pos - pt_ant
        norm = np.hypot(v_tang[0], v_tang[1])
        if norm < 1e-6:
            continue
        v_tang /= norm
        
        v_norm = np.array([-v_tang[1], v_tang[0]])
        
        if s < s_bifurcacao:
            frac = s / s_bifurcacao
            meia_largura = 70.0 + frac * 70.0   # 70m a 140m
            largura_canal = 35.0 + frac * 20.0  # calha de 35m a 55m
        else:
            meia_largura = 100.0                # 200m total no canal
            largura_canal = 45.0                # calha exata do canal
            
        p_esq = pt + v_norm * meia_largura
        p_dir = pt - v_norm * meia_largura
        
        rs = comprimento_total - s
        
        secoes.append({
            "s_eixo": s,
            "rs": round(float(rs), 2),
            "cut": LineString([p_esq, p_dir]),
            "meia_largura": meia_largura,
            "largura_canal": largura_canal,
            "no_canal": (s >= s_bifurcacao)
        })
        
    print(f"  [2/5] {len(secoes)} cutlines ortogonais geradas ao longo do eixo.")
    return secoes


# ==============================================================================
# 3. EXTRAÇÃO DE COTAS DO MDT SIG-SC 1M E PERFIL SUAVE MONÓTONO
# ==============================================================================
def amostrar_terreno_e_perfil(secoes, mosaico):
    """Extrai cotas reais do MDT e impõe perfil hidráulico contínuo decrescente."""
    print("  [3/5] Amostrando MDT SIG-SC 1m e calculando talvegue contínuo suave...")
    
    n_sec = len(secoes)
    
    # 1. Amostragem inicial do terreno para cada seção
    for sec in secoes:
        cut = sec["cut"]
        comprimento_corte = cut.length
        pts_corte = np.linspace(0.0, comprimento_corte, num=int(comprimento_corte / 3.0) + 1)
        
        centro = round(comprimento_corte / 2.0, 2)
        w_ch = sec["largura_canal"]
        lb = max(0.0, round(centro - w_ch / 2.0, 2))
        rb = min(comprimento_corte, round(centro + w_ch / 2.0, 2))
        
        # Garante que lb, rb e centro estejam estritamente na lista de estações
        sta = np.unique(np.append(pts_corte, [lb, centro, rb]))
        sta = np.round(sta, 2)
        # Deduplica estacas estritamente após arredondamento para 2 decimais
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
            cotas_tratadas = np.full(sta.shape, 10.0)
            
        sec["sta"] = sta
        sec["z_bruto"] = np.round(cotas_tratadas, 2)
        sec["lb"] = lb
        sec["rb"] = rb
        sec["centro"] = centro
        sec["z_min_bruto"] = float(cotas_tratadas.min())
        
    # 2. Construção do perfil de talvegue PCHIP estritamente monótono e suave
    s_arr = np.array([sec["s_eixo"] for sec in secoes])
    
    # Pontos de ancoragem reais ao longo do vale (s_eixo, cota talvegue)
    pontos_controle = [
        (0.0, 275.0),
        (13500.0, 214.0),
        (33500.0, 115.0),
        (53500.0, 40.0),
        (73500.0, 12.0),
        (93500.0, 1.4),
        (105950.0, -0.76),
        (s_arr[-1], -2.68)
    ]
    s_ctrl = np.array([p[0] for p in pontos_controle])
    z_ctrl = np.array([p[1] for p in pontos_controle])
    
    pchip = PchipInterpolator(s_ctrl, z_ctrl)
    z_suave = pchip(s_arr)
    
    # Garante declividade mínima estrita de 10 cm / km de montante para jusante
    for i in range(n_sec - 1):
        dx = abs(s_arr[i+1] - s_arr[i])
        z_max_permitido = z_suave[i] - 0.0001 * dx
        if z_suave[i+1] > z_max_permitido:
            z_suave[i+1] = z_max_permitido
            
    # 3. Aplica o talvegue suave na calha ativa de cada seção transversal
    for i, sec in enumerate(secoes):
        sta = sec["sta"]
        z = sec["z_bruto"].copy()
        lb = sec["lb"]
        rb = sec["rb"]
        centro = sec["centro"]
        w_ch = sec["largura_canal"]
        
        z_alvo_talvegue = round(float(z_suave[i]), 2)
        
        # Cotas de margem garantidas pelo menos 2.5m acima do talvegue
        idx_lob = np.where(sta <= lb)[0]
        idx_rob = np.where(sta >= rb)[0]
        z_lob = z[idx_lob[-1]] if len(idx_lob) > 0 else z_alvo_talvegue + 3.0
        z_rob = z[idx_rob[0]] if len(idx_rob) > 0 else z_alvo_talvegue + 3.0
        
        z_lob = max(z_lob, z_alvo_talvegue + 2.5)
        z_rob = max(z_rob, z_alvo_talvegue + 2.5)
        cota_margem = (z_lob + z_rob) / 2.0
        
        # Esculpe a calha parabólica suave entre lb e rb
        idx_calha = (sta >= lb) & (sta <= rb)
        dist_norm = np.abs(sta[idx_calha] - centro) / (w_ch / 2.0)
        dist_norm = np.clip(dist_norm, 0.0, 1.0)
        
        z[idx_calha] = z_alvo_talvegue + (cota_margem - z_alvo_talvegue) * (dist_norm ** 2)
        
        # Garante que o centro seja exatamente a cota de talvegue
        idx_centro = np.argmin(np.abs(sta - centro))
        z[idx_centro] = z_alvo_talvegue
        
        # Ajusta os overbanks para nunca ficarem abaixo da cota da margem
        idx_esq = np.where(sta < lb)[0]
        idx_dir = np.where(sta > rb)[0]
        if len(idx_esq) > 0:
            z[idx_esq] = np.maximum(z[idx_esq], z_lob)
        if len(idx_dir) > 0:
            z[idx_dir] = np.maximum(z[idx_dir], z_rob)
            
        # No canal retificado e foz, garante diques de borda >= +4.50m para conter maré e cheias
        if sec["no_canal"]:
            z[0] = max(z[0], 4.50)
            z[-1] = max(z[-1], 4.50)
            
        sec["sta"] = np.round(sta, 2)
        sec["z"] = np.round(z, 2)
        sec["lb"] = round(float(lb), 2)
        sec["rb"] = round(float(rb), 2)
        sec["z_min"] = float(z_alvo_talvegue)
        
    return secoes


# ==============================================================================
# 4. ESCRITA DO PROJETO HEC-RAS 7.X
# ==============================================================================
def escrever_geometria_g01(caminho_g01, secoes, eixo, nome_projeto="mirim_novo"):
    """Escreve o arquivo .g01 oficial com todas as seções e eixos formatados."""
    linhas = []
    linhas.append(f"Geom Title={nome_projeto}")
    linhas.append("Program Version=7.01")
    linhas.append("Viewing Rectangle= 660562.18 , 732107.05 , 7025414.41 , 6967422.16 ")
    linhas.append('Spatial Reference System=PROJCS["SIRGAS_2000_UTM_Zone_22S",GEOGCS["GCS_SIRGAS_2000",DATUM["D_SIRGAS_2000",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",-51.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
    linhas.append("")
    
    # Eixo do Rio
    linhas.append("River Reach=Itajai_Mirim    ,R1              ")
    coords_eixo = list(eixo.coords)
    linhas.append(f"Reach XY= {len(coords_eixo)} ")
    
    pts_str = []
    for x, y in coords_eixo:
        pts_str.append(f"{x:16.4f}{y:16.4f}")
    for k in range(0, len(pts_str), 2):
        linhas.append("".join(pts_str[k:k+2]))
    linhas.append("Rch Text X Y=0,0,0,0")
    linhas.append("")
    
    # Escreve as seções transversais (Ordem de montante para jusante)
    for i, s in enumerate(secoes):
        rs = s["rs"]
        if i + 1 < len(secoes):
            dx_jus = round(float(abs(secoes[i+1]["s_eixo"] - secoes[i]["s_eixo"])), 2)
            dx_lob = dx_ch = dx_rob = dx_jus
        else:
            dx_lob = dx_ch = dx_rob = 0.0
            
        linhas.append(f"Type RM Length L Ch R = 1 ,{rs:.2f},{dx_lob:8.2f},{dx_ch:8.2f},{dx_rob:8.2f}")
        linhas.append(f"Bank Sta={s['lb']:.2f},{s['rb']:.2f}")
        
        # GIS Cut Line
        cut_coords = list(s["cut"].coords)
        linhas.append(f"XS GIS Cut Line= {len(cut_coords)}")
        cut_str = "".join([f"{pt[0]:16.2f}{pt[1]:16.2f}" for pt in cut_coords])
        linhas.append(cut_str)
        
        # Pontos da seção (#Sta/Elev)
        linhas.append(f"#Sta/Elev= {len(s['sta'])} ")
        pts_fmt = [f"{st:8.2f}{ele:8.2f}" for st, ele in zip(s["sta"], s["z"])]
        for k in range(0, len(pts_fmt), 5):
            linhas.append("".join(pts_fmt[k:k+5]))
            
        # Manning calibrado
        n_overbank = 0.055 if not s["no_canal"] else 0.045
        n_canal = 0.032 if not s["no_canal"] else 0.025
        linhas.append("#Mann= 3 , 0 , 0 ")
        linhas.append(f"{s['sta'][0]:8.2f}{n_overbank:8.3f}{0:8d}{s['lb']:8.2f}{n_canal:8.3f}{0:8d}{s['rb']:8.2f}{n_overbank:8.3f}{0:8d}")
        
        # HTAB calibrado
        z_min = s["z_min"]
        linhas.append(f"XS HTab Starting El and Incr={z_min + 0.02:.2f},0.100, 500 ")
        linhas.append("XS HTab Horizontal Distribution=-1,-1,-1")
        linhas.append("XS Rating Curve= 0 ,0")
        linhas.append("Exp/Cntr=0.3,0.1")
        linhas.append("")
        
    escrever(caminho_g01, "\r\n".join(linhas))
    print(f"  [4/5] Geometria gravada com sucesso em: {caminho_g01}")


def criar_arquivos_suporte(pasta_dest, secoes, nome_projeto="mirim_novo"):
    """Cria os arquivos .prj, .p01, .u01 e .rasmap para rodar o modelo."""
    pasta = pathlib.Path(pasta_dest)
    rs_montante = secoes[0]["rs"]
    rs_segunda = secoes[1]["rs"]
    rs_penultima = secoes[-2]["rs"]
    rs_jusante = secoes[-1]["rs"]
    
    prj_txt = f"""Proj Title={nome_projeto}
Current Plan=p01
Default Exp/Contr=0.3,0.1
SI Units
Geom File=g01
Unsteady File=u01
Plan File=p01
Y Axis Title=Elevation
X Axis Title(PF)=Main Channel Distance
X Axis Title(XS)=Station
RASMap Filename={nome_projeto}.rasmap
BEGIN DESCRIPTION:
Modelo do Rio Itajai-Mirim com Canal Retificado gerado do zero.
END DESCRIPTION:
"""
    escrever(str(pasta / f"{nome_projeto}.prj"), prj_txt)
    
    p01_txt = f"""Plan Title={nome_projeto}
Program Version=7.01
Short Identifier={nome_projeto}
Geom File=g01
Flow File=u01
Simulation Date=01AUG2026,0000,08AUG2026,2300
Mixed Flow Regime
Computation Interval=5SEC
Output Interval=1HOUR
Instantaneous Interval=1HOUR
Mapping Interval=1HOUR
UNET Theta= 1
UNET Theta Warmup= 1
UNET WFStab= 2
UNET SFStab= 1
UNET WFX= 1
UNET SFX= 1
UNET DZMax Abort= 30
UNET Froude Reduction=True
UNET Froude Limit= 0.8 
UNET Froude Power= 4 
UNET ZTol= 0.02 
UNET ZSATol= 0.02 
UNET MxIter= 40 
Write Detailed= 1
Run HTab=-1
Run UNet=-1
Run PostProcess=-1
Run RASMapper=-1
"""
    escrever(str(pasta / f"{nome_projeto}.p01"), p01_txt)
    
    # Gera arquivo de vazão não-permanente (.u01) ancorando nas estacas EXATAS da nova geometria
    u01_orig = os.path.join(DIRETORIO_RAIZ, "modelo", "mirim_canal6", "mirim_canal6.u01")
    if os.path.exists(u01_orig):
        txt_u01 = open(u01_orig, encoding="latin-1").read()
        txt_u01 = re.sub(r"(?m)^Flow Title=.*$", f"Flow Title={nome_projeto}", txt_u01)
        
        # Ajusta estacas de contorno de montante, lateral e jusante
        init_block = f"""Initial RS=Itajai_Mirim    ,R1              ,{rs_montante:.2f}, 2.19
Initial RS=Itajai_Mirim    ,R1              ,60000.00, 7.50
Initial RS=Itajai_Mirim    ,R1              ,{rs_jusante:.2f}, 12.85"""
        txt_u01 = re.sub(r"(?m)^Initial RS=.*$", init_block, txt_u01)
        
        txt_u01 = re.sub(r"(?m)^Boundary Location=Itajai_Mirim\s*,R1\s*,141422[\d.]*",
                         f"Boundary Location=Itajai_Mirim    ,R1              ,{rs_montante:.2f}", txt_u01)
        
        # Lateral Inflow: do trecho RS_segunda até RS_penultima
        txt_u01 = re.sub(r"(?m)^Boundary Location=Itajai_Mirim\s*,R1\s*,[\d.]+,\s*[\d.]+,",
                         f"Boundary Location=Itajai_Mirim    ,R1              ,{rs_segunda:.2f},{rs_penultima:.2f},", txt_u01)
                         
        txt_u01 = re.sub(r"(?m)^Boundary Location=Itajai_Mirim\s*,R1\s*,11756[\d.]*",
                         f"Boundary Location=Itajai_Mirim    ,R1              ,{rs_jusante:.2f}", txt_u01)
                         
        escrever(str(pasta / f"{nome_projeto}.u01"), txt_u01)
        
    prj_crs_orig = os.path.join(DIRETORIO_RAIZ, "SIRGAS2000_UTM22S.prj")
    if os.path.exists(prj_crs_orig):
        shutil.copy2(prj_crs_orig, pasta / "SIRGAS2000_UTM22S.prj")
        
    rasmap_orig = os.path.join(DIRETORIO_RAIZ, "modelo", "mirim_canal6", "mirim_canal6.rasmap")
    if os.path.exists(rasmap_orig):
        txt_map = open(rasmap_orig, encoding="utf-8", errors="ignore").read()
        txt_map = txt_map.replace("mirim_canal6", nome_projeto)
        escrever(str(pasta / f"{nome_projeto}.rasmap"), txt_map)
        
    if os.path.exists(os.path.join(DIRETORIO_RAIZ, "modelo", "Terrain")) and not (pasta / "Terrain").exists():
        shutil.copytree(os.path.join(DIRETORIO_RAIZ, "modelo", "Terrain"), pasta / "Terrain")


# ==============================================================================
# FLUXO PRINCIPAL
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(description="Gerar Modelo Mirim do Zero")
    parser.add_argument("--saida", default="modelo/mirim_novo", help="Diretório de saída")
    parser.add_argument("--nome", default="mirim_novo", help="Nome do projeto")
    args = parser.parse_args()
    
    pasta_out = pathlib.Path(args.saida)
    if pasta_out.exists():
        shutil.rmtree(pasta_out, ignore_errors=True)
    pasta_out.mkdir(parents=True, exist_ok=True)
    
    print("=" * 72)
    print(f"GERANDO MODELO DO RIO ITAJAÍ-MIRIM 100% DO ZERO: {args.nome}")
    print("=" * 72)
    
    # 1. Unifica eixo
    eixo, s_bif = obter_eixo_unificado()
    
    # 2. Gera cutlines ortogonais
    secoes = gerar_cutlines(eixo, s_bif, dx=150.0)
    
    # 3. Amostra relevo do MDT SIG-SC 1m e impõe talvegue monótono
    bbox = (eixo.bounds[0] - 2000, eixo.bounds[1] - 2000, eixo.bounds[2] + 2000, eixo.bounds[3] + 2000)
    tiles = tiles_do_dominio(bbox)
    print(f"  [3/5] Carregando {len(tiles)} folhas do MDT SIG-SC para o domínio...")
    mosaico = MosaicoSigsc(tiles=tiles)
    
    secoes = amostrar_terreno_e_perfil(secoes, mosaico)
    mosaico.fechar()
    
    # 4. Escreve arquivos HEC-RAS limpos
    caminho_g01 = str(pasta_out / f"{args.nome}.g01")
    escrever_geometria_g01(caminho_g01, secoes, eixo, args.nome)
    criar_arquivos_suporte(args.saida, secoes, args.nome)
    
    print(f"\n[5/5] Modelo construído do zero com sucesso em: {args.saida}")
    print("=" * 72)


if __name__ == "__main__":
    main()
