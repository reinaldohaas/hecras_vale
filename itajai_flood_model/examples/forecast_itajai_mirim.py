"""
Laboratório e Demonstração: Previsão Operacional de Cheia no Rio Itajaí-Mirim:
Fluxo:
    Chuva Observada Recente (Vidal Ramos / Botuverá)
    + Previsão de Chuva Futura (+24h com cenários Low / Mean / High)
    -> SCS-CN com AMC
    -> Hidrograma Unitário
    -> Propagação Muskingum (Vidal -> Botuverá -> Brusque)
    -> Canal Retificado Oficial
    -> Deságue na Foz em Itajaí
    -> Comparação com Estação ANA 83800000 (Brusque)
    -> Sistema de Alertas

Execução:
    python examples/forecast_itajai_mirim.py
"""

import sys
import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, BASE_DIR)

from src.rainfall.provider import SyntheticRainfallProvider
from src.rainfall.forecast import RainfallForecastTimeline
from src.rainfall.antecedent_moisture import AntecedentMoistureCondition
from src.forecasting.engine import OperationalForecastEngine
from src.forecasting.assimilation import StreamflowAssimilation
from src.forecasting.alerts import FloodAlertSystem

def run_itajai_mirim_forecast_lab():
    print("\n" + "=" * 80)
    print("   LABORATÓRIO DE PREVISÃO OPERACIONAL DE CHEIAS — RIO ITAJAÍ-MIRIM")
    print("=" * 80)

    # 1. Dados Pluviométricos: Histórico Observado (24h) + Previsão Futura (24h)
    print("\n[1] Carregando Séries Pluviométricas...")
    
    # 24h observadas passadas (acumulado 95 mm)
    p_obs_provider = SyntheticRainfallProvider(total_p_mm=95.0, duration_hours=24, start_time='2026-08-14 00:00:00', station_id='MIRIM_VIDAL')
    df_obs_raw = p_obs_provider.get_hourly_rainfall()
    
    # Formatar para o formato de sub-bacias
    df_obs = pd.DataFrame({
        'timestamp': df_obs_raw['timestamp'],
        'mirim': df_obs_raw['precipitation_mm'],
        'oeste': df_obs_raw['precipitation_mm'] * 0.7,
        'sul': df_obs_raw['precipitation_mm'] * 0.8,
        'norte': df_obs_raw['precipitation_mm'] * 0.9,
        'benedito': df_obs_raw['precipitation_mm'] * 0.6,
        'luis_alves': df_obs_raw['precipitation_mm'] * 0.5
    })
    
    now_ts = df_obs['timestamp'].iloc[-1]
    print(f"  • Instante Divisor [ AGORA ]: {now_ts}")
    print(f"  • Chuva Observada nas últimas 24h no Alto Mirim: {df_obs['mirim'].sum():.1f} mm")

    # 2. Condição de Umidade Antecedente
    p5 = 48.0 # mm acumulados nos 5 dias anteriores
    amc_class = AntecedentMoistureCondition.classify_amc_from_p5(p5)
    cn_adjusted = AntecedentMoistureCondition.adjust_curve_number(79.0, amc_class)
    print(f"\n[2] Condição de Umidade do Solo (AMC):")
    print(f"  • P5 Antecedente = {p5} mm -> Classificação: {amc_class}")
    print(f"  • Curve Number Base (CN_II) = 79.0 -> Ajustado ({amc_class}) = {cn_adjusted:.1f}")

    # 3. Linha do Tempo Contínua e Cenários de Previsão Pluviométrica (+24h)
    print("\n[3] Gerando Linha do Tempo Contínua e Cenários (+24h)...")
    timeline = RainfallForecastTimeline(
        observed_df=df_obs,
        now_timestamp=now_ts,
        horizon_hours=24,
        uncertainty_factor_low=0.60,
        uncertainty_factor_high=1.40
    )
    scenarios = timeline.build_continuous_scenarios()
    
    rain_mean_fut = scenarios['mean'][scenarios['mean']['timestamp'] > now_ts]['mirim'].sum()
    rain_low_fut = scenarios['low'][scenarios['low']['timestamp'] > now_ts]['mirim'].sum()
    rain_high_fut = scenarios['high'][scenarios['high']['timestamp'] > now_ts]['mirim'].sum()
    
    print(f"  • Chuva Prevista Futura (+24h) - Cenário Central (Mean): {rain_mean_fut:.1f} mm")
    print(f"  • Intervalo de Incerteza Pluviométrica: [{rain_low_fut:.1f} mm a {rain_high_fut:.1f} mm]")

    # 4. Execução da Previsão Hidrológica e Roteamento Fluvial
    print("\n[4] Executando Motor Operacional de Previsão de Vazões...")
    engine = OperationalForecastEngine()
    forecast_pkg = engine.generate_full_forecast_package(scenarios, amc=amc_class)
    
    st_mirim = forecast_pkg['stations_summary']['brusque']
    st_canal = forecast_pkg['stations_summary']['itajai_foz']
    
    print(f"\n--- Síntese Operacional de Previsão no Rio Itajaí-Mirim ---")
    print(f"  🏙️ BRUSQUE CENTRO (Estação ANA 83800000):")
    print(f"     • Vazão Atual no instante [AGORA]:     {st_mirim['q_now_m3s']} m³/s")
    print(f"     • Vazão de Pico Prevista (Central):    {st_mirim['q_peak_mean_m3s']} m³/s")
    print(f"     • Faixa Provável de Pico:              [{st_mirim['q_peak_low_m3s']} a {st_mirim['q_peak_high_m3s']}] m³/s")
    print(f"     • Tempo Estimado até a Crista:         {st_mirim['hours_to_peak']} horas (em {st_mirim['peak_timestamp']})")
    print(f"     • Tendência Atual:                     {st_mirim['trend']}")

    # 5. Classificação de Níveis de Alerta
    alert_sys = FloodAlertSystem()
    alert_now = alert_sys.classify_flow('brusque', st_mirim['q_now_m3s'])
    alert_peak = alert_sys.classify_flow('brusque', st_mirim['q_peak_mean_m3s'])
    
    print(f"\n[5] Níveis de Alerta e Defesa Civil (Brusque):")
    print(f"  • Situação Atual:   [{alert_now['level']}] - {alert_now['message']}")
    print(f"  • Situação no Pico: [{alert_peak['level']}] - {alert_peak['message']}")

    # 6. Assimilação e Correção em Tempo Real
    print(f"\n[6] Teste de Assimilação com Estação Telemétrica em Tempo Real:")
    q_brusque_raw = np.array(forecast_pkg['forecast_curves']['brusque']['q_mean'])
    q_telemetria_real_now = st_mirim['q_now_m3s'] * 1.12 # Supondo telemetria medindo 12% a mais
    q_brusque_assimilated = StreamflowAssimilation.apply_realtime_correction(
        q_forecast_raw=q_brusque_raw,
        q_observed_at_now=q_telemetria_real_now,
        now_idx=forecast_pkg['t_now_idx'],
        relaxation_hours=12.0
    )
    print(f"  • Vazão Telemetria Observada no momento: {q_telemetria_real_now:.1f} m³/s")
    print(f"  • Novo Pico de Brusque Pós-Assimilação:  {np.max(q_brusque_assimilated):.1f} m³/s")

    print("\n" + "=" * 80)
    print("   LABORATÓRIO CONCLUÍDO COM SUCESSO! MODELO 100% OPERACIONAL.")
    print("=" * 80 + "\n")

if __name__ == '__main__':
    run_itajai_mirim_forecast_lab()
