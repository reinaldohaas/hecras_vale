# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT INDEPENDENTE: AUDITORIA DE QUALIDADE GEOMÉTRICA DO RIO ITAJAÍ-MIRIM
Diretório do Script : a_scripts/auditar_mirim.py
Alvo Padrão         : modelos/_anti/mirim/mirim.g01
================================================================================
"""
import argparse
import os
import re
import sys
import numpy as np

DIR_SCRIPT = os.path.dirname(os.path.abspath(__file__))
DIR_RAIZ = os.path.dirname(DIR_SCRIPT)


def ler_secoes_g01(caminho_g01):
    """Lê as seções transversais diretamente do arquivo .g01."""
    if not os.path.exists(caminho_g01):
        raise FileNotFoundError(f"Arquivo não encontrado: {caminho_g01}")
        
    with open(caminho_g01, "r", encoding="latin-1") as f:
        linhas = f.readlines()
        
    secoes = []
    secao_atual = None
    lendo_pts = False
    
    for l in linhas:
        if l.startswith("Type RM Length L Ch R"):
            if secao_atual:
                secoes.append(secao_atual)
            partes = l.strip().split(",")
            rs = float(partes[1])
            secao_atual = {"rs": rs, "sta": [], "z": [], "lb": None, "rb": None}
            lendo_pts = False
        elif l.startswith("Bank Sta="):
            partes = l.strip().replace("Bank Sta=", "").split(",")
            if secao_atual:
                secao_atual["lb"] = float(partes[0])
                secao_atual["rb"] = float(partes[1])
        elif l.startswith("#Sta/Elev="):
            lendo_pts = True
        elif l.startswith("#Mann=") or l.startswith("XS HTab") or l.startswith("Exp/Cntr="):
            lendo_pts = False
        elif lando_pts := lando_pts if False else lando_pts if False else lando_pts if False else False:
            pass
        elif lendo_pts:
            tokens = re.findall(r"[-+]?\d*\.\d+|\d+", l)
            for k in range(0, len(tokens)-1, 2):
                secao_atual["sta"].append(float(tokens[k]))
                secao_atual["z"].append(float(tokens[k+1]))
                
    if secao_atual:
        secoes.append(secao_atual)
        
    for s in secoes:
        s["sta"] = np.array(s["sta"])
        s["z"] = np.array(s["z"])
        s["z_min"] = float(np.min(s["z"])) if len(s["z"]) > 0 else np.nan
        
    return secoes


def auditar_geometria(caminho_g01):
    print("=" * 72)
    print(f"AUDITORIA GEOMÉTRICA: {caminho_g01}")
    print("=" * 72)
    
    secoes = ler_secoes_g01(caminho_g01)
    n_sec = len(secoes)
    print(f"Total de seções analisadas: {n_sec}")
    
    if n_sec == 0:
        print("Nenhuma seção encontrada.")
        return
        
    rs_arr = np.array([s["rs"] for s in secoes])
    z_min_arr = np.array([s["z_min"] for s in secoes])
    
    # 1. Declividade e Subidas de Fundo
    dz = z_min_arr[:-1] - z_min_arr[1:]
    dx = rs_arr[:-1] - rs_arr[1:]
    s0 = dz / dx
    
    subidas = np.where(s0 <= 0)[0]
    print(f"\n[1] PERFIL DE TALVEGUE:")
    print(f"  - Cota Montante  : {z_min_arr[0]:.2f} m")
    print(f"  - Cota Foz       : {z_min_arr[-1]:.2f} m")
    print(f"  - Declividade Med: {s0.mean()*100:.3f}% ({s0.mean()*1000:.2f} m/km)")
    print(f"  - Declividade Max: {s0.max()*100:.3f}% ({s0.max()*1000:.2f} m/km) [RS {rs_arr[np.argmax(s0)]:.2f}]")
    print(f"  - Declividade Min: {s0.min()*100:.3f}% ({s0.min()*1000:.2f} m/km) [RS {rs_arr[np.argmin(s0)]:.2f}]")
    print(f"  - Degraus Adversos (subidas): {len(subidas)} (ideal = 0)")
    
    # 2. Bank Stations e Pontos Duplicados
    erros_bancos = 0
    duplicatas = 0
    prof_baixa = 0
    
    for s in secoes:
        sta = s["sta"]
        z = s["z"]
        # Verifica duplicatas
        if len(sta) != len(np.unique(sta)):
            duplicatas += 1
        # Verifica se lb e rb estão estritamente nos dados
        if not np.any(np.isclose(sta, s["lb"], atol=1e-2)) or not np.any(np.isclose(sta, s["rb"], atol=1e-2)):
            erros_bancos += 1
        # Profundidade útil da calha
        idx_lb = np.argmin(np.abs(sta - s["lb"]))
        idx_rb = np.argmin(np.abs(sta - s["rb"]))
        prof = min(z[idx_lb], z[idx_rb]) - s["z_min"]
        if prof < 1.5:
            prof_baixa += 1
            
    print(f"\n[2] CONSISTÊNCIA GEOMÉTRICA:")
    print(f"  - Seções com Pontos Duplicados : {duplicatas} (ideal = 0)")
    print(f"  - Bank Stations Desalinhados   : {erros_bancos} (ideal = 0)")
    print(f"  - Seções com Profundidade <1.5m: {prof_baixa}")
    print("=" * 72)


def main():
    parser = argparse.ArgumentParser(description="Auditoria geométrica independente do Rio Itajaí-Mirim")
    parser.add_argument("--g01", default=os.path.join(DIR_RAIZ, "modelos", "_anti", "mirim", "mirim.g01"),
                        help="Caminho do arquivo .g01")
    args = parser.parse_args()
    
    auditar_geometria(args.g01)


if __name__ == "__main__":
    main()
