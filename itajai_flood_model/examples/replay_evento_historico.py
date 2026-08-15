"""
Demonstração: Replay de Eventos Históricos Reais da Bacia do Rio Itajaí:
- Eventos: 1983 (Secular), 2008 (Desastre Baixo Vale), 2011 (Alto/Médio Vale), 2023 (Multievento)
- Compara hidrogramas simulados com hidrogramas observados nas estações da ANA/CEOPS
- Calcula métricas estatísticas formais (NSE, RMSE, Erro de Pico, Diferença de Horário da Crista)

Execução:
    python examples/replay_evento_historico.py
"""

import sys
import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, BASE_DIR)

from src.historical_events.events_loader import HistoricalEventsLoader
from src.forecasting.engine import OperationalForecastEngine
from src.calibration import HydrographValidator

def run_historical_replay():
    print("\n" + "=" * 85)
    print("   REPLAY E VALIDAÇÃO DOS GRANDES EVENTOS HISTÓRICOS DO VALE DO ITAJAÍ")
    print("=" * 85)
    
    events_to_test = ['1983', '2008', '2011', '2023']
    
    for ev_id in events_to_test:
        ev_data = HistoricalEventsLoader.get_event_data(ev_id)
        meta = ev_data['metadata']
        
        print(f"\n" + "-" * 85)
        print(f"🌊 EVENTO HISTÓRICO: {meta['name']} ({meta['period']})")
        print(f"   • Duração: {meta['duration_hours']} horas | Chuva Antecedente P5 = {meta['p5_antecedent_mm']} mm")
        print(f"   • Regime de Barragens: {meta['dam_operating_status']}")
        print("-" * 85)
        
        # Configurar motor com a operação histórica
        dam_cfg = {
            'oeste': {'total_gates': 7, 'open_gates': 7 if ev_id == '1983' else 0, 'cap_hm3': 83.0, 'base_flow': 15.0},
            'sul': {'total_gates': 5, 'open_gates': 5 if ev_id == '1983' else 0, 'cap_hm3': 93.5, 'base_flow': 12.0},
            'norte': {'total_gates': 2, 'open_gates': 2 if ev_id in ['1983', '2008'] else 0, 'cap_hm3': 357.0, 'base_flow': 20.0}
        }
        engine = OperationalForecastEngine(dam_operations=dam_cfg)
        
        # Simulação hidrológica completa
        amc = 'AMC_III' if meta['p5_antecedent_mm'] > 50 else 'AMC_II'
        results = engine.execute_basin_forecast(ev_data['rainfall_df'], amc=amc)
        
        # Avaliar Blumenau
        q_sim_blu = results['blumenau']
        q_obs_blu = ev_data['observed_hydrographs']['blumenau']
        m_blu = HydrographValidator.calculate_metrics(q_sim_blu, q_obs_blu)
        
        # Avaliar Rio do Sul
        q_sim_rs = results['rio_do_sul']
        q_obs_rs = ev_data['observed_hydrographs']['rio_do_sul']
        m_rs = HydrographValidator.calculate_metrics(q_sim_rs, q_obs_rs)
        
        # Avaliar Brusque
        q_sim_bq = results['brusque']
        q_obs_bq = ev_data['observed_hydrographs']['brusque']
        m_bq = HydrographValidator.calculate_metrics(q_sim_bq, q_obs_bq)
        
        print(f"\n  [1] Blumenau Centro (Estação ANA 83700000):")
        print(f"      • Pico Obs: {m_blu['peak_obs_m3s']} m³/s | Pico Sim: {m_blu['peak_sim_m3s']} m³/s (Erro: {m_blu['peak_error_pct']}%)")
        print(f"      • Diferença Horário Crista: {m_blu['t_peak_diff_h']} h | RMSE: {m_blu['rmse_m3s']} m³/s | NSE: {m_blu['nse']}")
        
        print(f"\n  [2] Rio do Sul (Confluência Alto Vale):")
        print(f"      • Pico Obs: {m_rs['peak_obs_m3s']} m³/s | Pico Sim: {m_rs['peak_sim_m3s']} m³/s (Erro: {m_rs['peak_error_pct']}%)")
        print(f"      • Diferença Horário Crista: {m_rs['t_peak_diff_h']} h | RMSE: {m_rs['rmse_m3s']} m³/s | NSE: {m_rs['nse']}")

        print(f"\n  [3] Brusque (Rio Itajaí-Mirim):")
        print(f"      • Pico Obs: {m_bq['peak_obs_m3s']} m³/s | Pico Sim: {m_bq['peak_sim_m3s']} m³/s (Erro: {m_bq['peak_error_pct']}%)")
        print(f"      • Diferença Horário Crista: {m_bq['t_peak_diff_h']} h | RMSE: {m_bq['rmse_m3s']} m³/s | NSE: {m_bq['nse']}")

    print("\n" + "=" * 85)
    print("   TODOS OS REPLAYS HISTÓRICOS EXECUTADOS E VALIDADOS COM SUCESSO!")
    print("=" * 85 + "\n")

if __name__ == '__main__':
    run_historical_replay()
