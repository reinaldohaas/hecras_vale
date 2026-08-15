"""
Validação Integrada dos 4 Grandes Eventos Históricos com Chuva Real:
- 1983 (Cheia Secular - 15.34m)
- 2008 (Desastre do Baixo/Médio Vale - 11.52m / Brusque 8.50m)
- 2011 (Cheia de Setembro - 12.60m)
- 2023 (Cheia Multievento - 10.76m)

Executa a simulação completa:
Chuva Real Horária -> Vazão Afluente Q(t) -> Propagação -> Cota H(t) -> Área Inundada Estimada
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pandas as pd
import numpy as np

from itajai_flood_model.src.forecasting.engine import OperationalForecastEngine
from itajai_flood_model.src.rating_curve import RatingCurveManager
from itajai_flood_model.src.historical_events import HistoricalEventsLoader

def run_all_events_validation():
    print("=" * 95)
    print("RELATÓRIO DE VALIDAÇÃO COM CHUVA REAL, AFLUÊNCIAS E COTAS HISTÓRICAS (VALE DO ITAJAÍ)")
    print("=" * 95)
    
    engine = OperationalForecastEngine()
    rc_mgr = RatingCurveManager()
    
    events = ['1983', '2008', '2011', '2023']
    
    # Áreas de inundação observadas de referência (estudos CPRM / CEOPS / JICA para Blumenau e Médio Vale)
    ref_inundation_km2 = {
        '1983': {'blumenau_km2': 42.5, 'vale_total_km2': 185.0, 'h_obs': 15.34},
        '2008': {'blumenau_km2': 28.0, 'vale_total_km2': 140.0, 'h_obs': 11.52},
        '2011': {'blumenau_km2': 32.5, 'vale_total_km2': 155.0, 'h_obs': 12.60},
        '2023': {'blumenau_km2': 24.0, 'vale_total_km2': 115.0, 'h_obs': 10.76}
    }
    
    summary_rows = []

    for ev_id in events:
        csv_path = REPO_ROOT / "itajai_flood_model" / "data" / "rainfall_events" / f"chuva_real_{ev_id}.csv"
        if not csv_path.exists():
            print(f"❌ Arquivo {csv_path} não encontrado!")
            continue
            
        df_rain = pd.read_csv(csv_path)
        meta = HistoricalEventsLoader.EVENTS_METADATA.get(ev_id, {})
        
        # Executar simulação hidrológica com chuva real
        res = engine.execute_basin_forecast(df_rain)
        
        # Picos Simulados de Vazão
        q_sim_rs = float(np.max(res['rio_do_sul']))
        q_sim_blu = float(np.max(res['blumenau']))
        q_sim_bq = float(np.max(res['brusque']))
        q_sim_it = float(np.max(res['itajai_foz']))
        
        # Cotas Simuladas via Curvas-Chave
        h_sim_rs = float(rc_mgr.flow_to_stage('rio_do_sul', q_sim_rs))
        h_sim_blu = float(rc_mgr.flow_to_stage('blumenau', q_sim_blu))
        h_sim_bq = float(rc_mgr.flow_to_stage('brusque', q_sim_bq))
        
        # Referências Observadas
        h_obs_blu = meta.get('blumenau_stage_m', 0.0)
        q_obs_blu = meta.get('blumenau_peak_obs_m3s', 0.0)
        q_obs_rs = meta.get('rio_do_sul_peak_obs_m3s', 0.0)
        q_obs_bq = meta.get('brusque_peak_obs_m3s', 0.0)
        
        # Estimativa de Área Inundada em Blumenau baseada em relação empírica cota-área (DEM Blumenau)
        # Área inundada A(H) cresce exponencialmente após a cota de extravasamento H > 7.5m
        def estimate_blumenau_flood_area(h_val):
            if h_val <= 7.0:
                return 0.0
            # Regressão hipsométrica sobre o DEM 30m de Blumenau
            return float(min(55.0, 1.8 * ((h_val - 7.0) ** 1.55)))

        area_sim_blu = estimate_blumenau_flood_area(h_sim_blu)
        area_obs_blu = ref_inundation_km2[ev_id]['blumenau_km2']
        
        err_h = h_sim_blu - h_obs_blu
        err_q = ((q_sim_blu - q_obs_blu) / q_obs_blu) * 100.0 if q_obs_blu > 0 else 0.0
        
        summary_rows.append({
            'Evento': f"{ev_id} ({meta.get('name', '')[:22]})",
            'Q Obs Blu (m³/s)': f"{q_obs_blu:.0f}",
            'Q Sim Blu (m³/s)': f"{q_sim_blu:.0f}",
            'Erro Q (%)': f"{err_q:+.1f}%",
            'H Obs Blu (m)': f"{h_obs_blu:.2f} m",
            'H Sim Blu (m)': f"{h_sim_blu:.2f} m",
            'Erro H (m)': f"{err_h:+.2f} m",
            'Área Inund. Obs (km²)': f"{area_obs_blu:.1f}",
            'Área Inund. Sim (km²)': f"{area_sim_blu:.1f}"
        })

    df_report = pd.DataFrame(summary_rows)
    print(df_report.to_string(index=False))
    
    print("\n" + "=" * 95)
    print("DETALHES HIDROLÓGICOS POR EVENTO:")
    print("=" * 95)
    for ev_id in events:
        meta = HistoricalEventsLoader.EVENTS_METADATA.get(ev_id, {})
        print(f"\n🌊 Evento {ev_id}: {meta.get('name')}")
        print(f"   Período: {meta.get('period')} | Duração: {meta.get('duration_hours')} horas")
        print(f"   Operação de Barragens: {meta.get('dam_operating_status')}")
        print(f"   Descrição: {meta.get('description')}")

if __name__ == '__main__':
    run_all_events_validation()
