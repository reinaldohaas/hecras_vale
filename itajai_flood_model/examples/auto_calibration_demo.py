"""
Demonstração: Calibração Automática Multievento e Validação Cega:
- Conjunto de Treinamento (Calibração): 1983 + 2008 + 2011
- Conjunto de Teste Independente (Validação Cega): 2023
- Avalia generalização e convergência dos parâmetros CN, Tc e K

Execução:
    python examples/auto_calibration_demo.py
"""

import sys
import os

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, BASE_DIR)

from src.historical_events.auto_calibration import BasinAutoCalibrator

def run_auto_calibration_demo():
    print("\n" + "=" * 80)
    print("   CALIBRAÇÃO AUTOMÁTICA MULTIEVENTO & VALIDAÇÃO INDEPENDENTE")
    print("=" * 80)
    
    calibrator = BasinAutoCalibrator(target_station='blumenau')
    
    print("\n[1] Iniciando Otimização Multievento nos eventos de Calibração:")
    print("    • Eventos de Treino: 1983 (Secular) + 2008 (Desastre Baixo Vale) + 2011 (Alto Vale)")
    print("    • Evento de Teste Cego: 2023 (Recente Multievento)")
    print("    • Estação Alvo: Blumenau Centro (Estação ANA 83700000)\n")
    
    results = calibrator.calibrate(
        training_event_ids=['1983', '2008', '2011'],
        validation_event_id='2023'
    )
    
    opt = results['optimal_parameters']
    print(f"--- Parâmetros Ótimos Encontrados ---")
    print(f"  • Multiplicador CN (Curve Number):  {opt['cn_multiplier']}")
    print(f"  • Multiplicador Tc (Concentração):  {opt['tc_multiplier']}")
    print(f"  • Multiplicador K (Routing Muskingum): {opt['k_multiplier']}")
    print(f"  • Score da Função Objetivo (Treino): {results['best_training_objective_score']}")
    
    print(f"\n--- Desempenho no Conjunto de Calibração (Treino) ---")
    for tr in results['training_results']:
        ev_id = tr['event_id']
        m = tr['metrics']
        print(f"  • Evento {ev_id}: NSE = {m['nse']:>6.4f} | RMSE = {m['rmse_m3s']:>6.1f} m³/s | Erro de Pico = {m['peak_error_pct']:>5.1f}%")
        
    val = results['validation_results']
    m_val = val['metrics']
    print(f"\n--- Desempenho no Conjunto de Validação Independente (Teste Cego - 2023) ---")
    print(f"  • Evento {val['event_id']}: NSE = {m_val['nse']:>6.4f} | RMSE = {m_val['rmse_m3s']:>6.1f} m³/s | Erro de Pico = {m_val['peak_error_pct']:>5.1f}%")
    print(f"  • Horário de Crista: Diferença de apenas {m_val['t_peak_diff_h']} hora(s)")
    
    print("\n" + "=" * 80)
    print("   CALIBRAÇÃO E VALIDAÇÃO CONCLUÍDAS COM SUCESSO!")
    print("=" * 80 + "\n")

if __name__ == '__main__':
    run_auto_calibration_demo()
