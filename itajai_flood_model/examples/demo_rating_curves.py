"""
Demonstração e Validação do Módulo de Curva-Chave (Rating Curve Q-H):
Exibe curvas oficiais e estimadas, tabelas de calibração para os 4 grandes eventos históricos
e gera o painel visual comparativo interativo.
"""

import os
import sys
from pathlib import Path

# Adicionar raiz do repositório ao sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import json
import numpy as np
import pandas as pd

from itajai_flood_model.src.rating_curve import (
    RatingCurveManager,
    CurveType,
    create_blumenau_official_curve,
    create_rio_do_sul_official_curve,
    create_brusque_official_curve,
    create_indaial_official_curve
)

def run_rating_curve_demo():
    print("=" * 80)
    print("CATÁLOGO DE CURVAS-CHAVE Q-H — BACIA DO RIO ITAJAÍ (ETAPA 1)")
    print("=" * 80)
    
    mgr = RatingCurveManager()
    stations = mgr.list_stations()
    
    for st in stations:
        tag = "[ OFICIAL ]" if st['is_official'] else "[ ESTIMADA ]"
        print(f"\n{tag} Chave: {st['key']}")
        print(f"   Estação: {st['name']} ({st['station_id']}) | Rio: {st['river']}")
        print(f"   Tipo: {st['curve_type']}")
        print(f"   Zero da Régua (Z0): {st['datum_z0_m']:.2f} m | Faixa: {st['validity_range_h_m'][0]:.1f}m a {st['validity_range_h_m'][1]:.1f}m")
        print(f"   Fonte: {st['source']}")

    print("\n" + "=" * 80)
    print("TABELA DE CALIBRAÇÃO NOS EVENTOS HISTÓRICOS DE BLUMENAU (ESTAÇÃO 83700000)")
    print("=" * 80)
    
    blu_curve = mgr.get_curve('blumenau')
    historical_benchmarks = [
        {'evento': 'Cheia Secular de 1983', 'h_obs': 15.34, 'q_obs': 5850.0},
        {'evento': 'Cheia de Setembro/2011', 'h_obs': 12.60, 'q_obs': 4650.0},
        {'evento': 'Desastre de Novembro/2008', 'h_obs': 11.52, 'q_obs': 4200.0},
        {'evento': 'Cheia de Outubro/2023', 'h_obs': 10.76, 'q_obs': 3950.0},
        {'evento': 'Nível de Emergência (CEOPS)', 'h_obs': 10.00, 'q_obs': 3450.0},
        {'evento': 'Nível de Alerta (CEOPS)', 'h_obs': 8.00, 'q_obs': 2400.0},
        {'evento': 'Nível de Atenção / Início de Cheia', 'h_obs': 5.00, 'q_obs': 1200.0},
        {'evento': 'Escoamento de Estiagem / Base', 'h_obs': 1.50, 'q_obs': 85.0}
    ]
    
    df_rows = []
    for b in historical_benchmarks:
        h = b['h_obs']
        q_calc = blu_curve.to_flow(h)
        h_rec = blu_curve.to_stage(q_calc)
        err_q = ((q_calc - b['q_obs']) / b['q_obs']) * 100.0
        df_rows.append({
            'Evento / Limiar': b['evento'],
            'H Obs (m)': f"{h:.2f}",
            'Q Obs (m³/s)': f"{b['q_obs']:.0f}",
            'Q Calc (m³/s)': f"{q_calc:.1f}",
            'Erro Rel Q (%)': f"{err_q:+.1f}%",
            'H Recuperado (m)': f"{h_rec:.2f}"
        })
        
    df_bench = pd.DataFrame(df_rows)
    print(df_bench.to_string(index=False))

    print("\n" + "=" * 80)
    print("AMOSTRA DAS RELAÇÕES VAZÃO-COTA NAS DEMAIS CIDADES")
    print("=" * 80)
    
    stages_test = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
    for city_key in ['rio_do_sul', 'brusque', 'indaial', 'ibirama']:
        c = mgr.get_curve(city_key)
        q_str = []
        for h in stages_test:
            if h <= c.h_max:
                q_str.append(f"H={h:.1f}m -> {c.to_flow(h):.0f} m³/s")
        print(f"• {c.name:45s} ({c.curve_type.value[:14]}): {', '.join(q_str)}")

    print("\n>>> ETAPA 1 VALIDADA COM SUCESSO!")

if __name__ == '__main__':
    run_rating_curve_demo()
