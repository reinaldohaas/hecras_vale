"""
Script Executável de Demonstração Completa da Bacia do Rio Itajaí:
Rede Integrada (Itajaí-Mirim + Itajaí-Açu + Rios do Oeste, Sul, Norte e Benedito com Barragens).

Execução:
    python examples/bacia_completa_demo.py
"""

import sys
import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, BASE_DIR)

from src.river import RiverNetwork, RiverReach
from src.unit_hydrograph import UnitHydrograph, scs_curve_number_excess
from src.routing import FloodRouter
from src.muskingum import MuskingumReach
from src.calibration import HydrographValidator
from src.visualization import FloodVisualizer

def run_bacia_completa_demo():
    print("\n" + "=" * 78)
    print("   MODELO HIDROLÓGICO INTEGRADO DE TODA A BACIA DO RIO ITAJAÍ")
    print("=" * 78)
    
    # 1. Transformação Chuva-Vazão em cada Sub-Bacia pelo SCS Curvilíneo
    print("\n[1] Gerando Hidrogramas de Cheia das Sub-Bacias...")
    
    # Sub-bacias e suas áreas (km2), tc (h), Chuva P (mm), CN
    subbasins = {
        'oeste': {'area': 3120.0, 'tc': 14.2, 'p': 120.0, 'cn': 76.0, 'name': 'Rio do Oeste (Taió)'},
        'sul': {'area': 2280.0, 'tc': 12.8, 'p': 130.0, 'cn': 78.0, 'name': 'Rio do Sul (Ituporanga)'},
        'norte': {'area': 3450.0, 'tc': 16.5, 'p': 140.0, 'cn': 75.0, 'name': 'Rio Hercílio (José Boiteux)'},
        'benedito': {'area': 1540.0, 'tc': 10.4, 'p': 110.0, 'cn': 77.0, 'name': 'Rio Benedito (Timbó)'},
        'mirim': {'area': 1680.0, 'tc': 11.8, 'p': 115.0, 'cn': 79.0, 'name': 'Rio Itajaí-Mirim (Brusque)'},
        'luis_alves': {'area': 580.0, 'tc': 7.2, 'p': 105.0, 'cn': 78.0, 'name': 'Rio Luís Alves'}
    }
    
    hydrographs = {}
    for k, sb in subbasins.items():
        # Chuva efetiva pelo Curve Number
        p_inc = np.zeros(25)
        # Distribuição SCS Tipo II
        for h in range(1, 25):
            frac = (0.5 * (h/12.0)**2) if h <= 12 else (1.0 - 0.5 * ((24.0 - h)/12.0)**2)
            p_inc[h] = sb['p'] * frac
        p_inc_diff = np.diff(np.insert(p_inc, 0, 0.0))
        pe_eff = scs_curve_number_excess(p_inc_diff, cn=sb['cn'])
        
        uh = UnitHydrograph(area_km2=sb['area'], tc_hours=sb['tc'], method='scs_curvilinear')
        q = uh.convolve(pe_inc_mm=pe_eff, dt_hours=1.0, base_flow=15.0, total_hours=48)
        hydrographs[k] = q
        print(f"  • {sb['name']:<30} Área: {sb['area']:>6.0f} km² | Pico: {np.max(q):>6.1f} m³/s em t={np.argmax(q):>2}h")

    # 2. Propagação e Confluências da Rede
    print("\n[2] Executando Propagação Fluvial na Rede Hidrográfica (Muskingum)...")
    
    # 2.1 Confluência em Rio do Sul
    m_oeste_riosul = MuskingumReach(1, "Oeste -> Rio do Sul", k_hours=8.0, x_param=0.20)
    m_sul_riosul = MuskingumReach(2, "Sul -> Rio do Sul", k_hours=6.0, x_param=0.20)
    
    q_oeste_jus = m_oeste_riosul.route(hydrographs['oeste'])
    q_sul_jus = m_sul_riosul.route(hydrographs['sul'])
    q_rio_do_sul = q_oeste_jus + q_sul_jus
    
    # 2.2 Confluência em Blumenau
    m_riosul_blumenau = MuskingumReach(3, "Rio do Sul -> Blumenau", k_hours=14.0, x_param=0.25)
    m_norte_blumenau = MuskingumReach(4, "José Boiteux -> Blumenau", k_hours=10.0, x_param=0.20)
    m_benedito_blumenau = MuskingumReach(5, "Benedito -> Blumenau", k_hours=4.0, x_param=0.20)
    
    q_altovale_blu = m_riosul_blumenau.route(q_rio_do_sul)
    q_norte_blu = m_norte_blumenau.route(hydrographs['norte'])
    q_benedito_blu = m_benedito_blumenau.route(hydrographs['benedito'])
    q_blumenau = q_altovale_blu + q_norte_blu + q_benedito_blu
    
    # 2.3 Rio Itajaí-Mirim
    m_mirim_brusque = MuskingumReach(6, "Mirim -> Brusque", k_hours=5.0, x_param=0.20)
    m_brusque_foz = MuskingumReach(7, "Brusque -> Canal Itajaí", k_hours=4.0, x_param=0.25)
    
    q_brusque = m_mirim_brusque.route(hydrographs['mirim'])
    q_mirim_foz = m_brusque_foz.route(q_brusque)
    
    # 2.4 Foz Total em Itajaí (Itajaí-Açu + Itajaí-Mirim + Luís Alves)
    m_blu_itajai = MuskingumReach(8, "Blumenau -> Itajaí", k_hours=8.0, x_param=0.30)
    q_acu_foz = m_blu_itajai.route(q_blumenau)
    q_itajai_total = q_acu_foz + q_mirim_foz + hydrographs['luis_alves'] * 0.7
    
    print("\n--- Resultados nas Estações Principais da Bacia ---")
    print(f"  1. 🏙️ Rio do Sul (Confluência Alto Vale):  Pico = {np.max(q_rio_do_sul):>6.1f} m³/s em t = {np.argmax(q_rio_do_sul):>2}h")
    print(f"  2. 🏙️ Blumenau Centro (Médio Vale):         Pico = {np.max(q_blumenau):>6.1f} m³/s em t = {np.argmax(q_blumenau):>2}h")
    print(f"  3. 🏙️ Brusque Centro (Itajaí-Mirim):        Pico = {np.max(q_brusque):>6.1f} m³/s em t = {np.argmax(q_brusque):>2}h")
    print(f"  4. 🌊 Itajaí Foz (Deságue Total Oceano):    Pico = {np.max(q_itajai_total):>6.1f} m³/s em t = {np.argmax(q_itajai_total):>2}h")
    
    print("\n" + "=" * 78)
    print("   MODELO INTEGRADO DA BACIA COMPLETA EXECUTADO COM SUCESSO!")
    print("=" * 78 + "\n")

if __name__ == '__main__':
    run_bacia_completa_demo()
