"""
Exemplo de Download de Chuva Real e Desagregação Horária para a Bacia do Rio Itajaí:
1. Baixa chuva horária real da API de Reanálise Histórica (Open-Meteo / ERA5) para as 10 sub-bacias
2. Demonstra a desagregação temporal de totais diários (24h -> 1h)
3. Executa a simulação hidrológica com a chuva real e calcula vazões e cotas.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pandas as pd
import numpy as np

from itajai_flood_model.src.rainfall import (
    RainfallDownloader,
    RainfallDisaggregator,
    DisaggregationMethod
)
from itajai_flood_model.src.forecasting.engine import OperationalForecastEngine
from itajai_flood_model.src.rating_curve import RatingCurveManager

def demo_real_rainfall_workflow():
    print("=" * 85)
    print("DEMO: DOWNLOAD DE CHUVA REAL HORÁRIA E SIMULAÇÃO HIDROLÓGICA (VALE DO ITAJAÍ)")
    print("=" * 85)

    # 1. Download do Desastre de Novembro de 2008 (20 a 26/Nov/2008)
    print("\n1. Baixando precipitação horária real do evento de Novembro/2008...")
    save_file_2008 = REPO_ROOT / "itajai_flood_model" / "data" / "rainfall_events" / "chuva_real_2008.csv"
    df_2008 = RainfallDownloader.download_hourly_event(
        start_date="2008-11-20",
        end_date="2008-11-26",
        save_path=str(save_file_2008)
    )

    print("\nResumo da Precipitação Total Acumulada por Sub-Bacia em 2008:")
    for col in df_2008.columns:
        if col != 'timestamp':
            tot = df_2008[col].sum()
            pk = df_2008[col].max()
            print(f"   • {col:14s}: Acumulado = {tot:6.1f} mm | Pico Horário = {pk:4.1f} mm/h")

    # 2. Executar Motor Hidrológico com a Chuva Real
    print("\n2. Executando Simulação Hidrológica com a Chuva Real de 2008...")
    engine = OperationalForecastEngine()
    sim_results = engine.execute_basin_forecast(df_2008)
    
    # 3. Converter Vazões em Níveis/Cotas com o RatingCurveManager
    rc_mgr = RatingCurveManager()
    
    print("\nResultados Hidrológicos & Hidráulicos Calculados com Chuva Real:")
    stations_map = {
        'rio_do_sul': {'col': 'rio_do_sul', 'key': 'rio_do_sul', 'name': 'Rio do Sul (Confluência)'},
        'blumenau': {'col': 'blumenau', 'key': 'blumenau', 'name': 'Blumenau Centro'},
        'brusque': {'col': 'brusque', 'key': 'brusque', 'name': 'Brusque Centro (Itajaí-Mirim)'},
        'itajai_foz': {'col': 'itajai_foz', 'key': 'itajai_foz', 'name': 'Itajaí (Foz & Canal)'}
    }
    
    for st_id, st_info in stations_map.items():
        q_series = sim_results[st_info['col']]
        max_q = float(np.max(q_series))
        t_max = int(np.argmax(q_series))
        
        # Converter para cota
        curve = rc_mgr.get_curve(st_info['key'])
        stage_max = float(curve.to_stage(max_q))
        
        print(f"   🏙️ {st_info['name']:32s}: Pico Q = {max_q:6.1f} m³/s | Cota H = {stage_max:5.2f} m | Hora do Pico: t={t_max}h")

    # 4. Demonstração de Desagregação Temporal Diária -> Horária
    print("\n" + "=" * 85)
    print("3. DEMONSTRAÇÃO: DESAGREGAÇÃO TEMPORAL DE CHUVA DIÁRIA (24h -> 1h)")
    print("=" * 85)
    
    # Exemplo: 120 mm registrados em 1 dia
    daily_precip = [120.0]
    df_disagg = RainfallDisaggregator.disaggregate_daily_series(
        daily_totals=daily_precip,
        method=DisaggregationMethod.SCS_TYPE_II,
        start_datetime="2023-10-08 00:00:00"
    )
    
    print(f"Total Diário de Entrada: {daily_precip[0]:.1f} mm")
    print(f"Total Integrado Horário: {df_disagg['precipitation_mm_h'].sum():.1f} mm (Conservação de Massa = 100%)")
    print(f"Pico Horário Resultante (SCS Tipo II): {df_disagg['precipitation_mm_h'].max():.2f} mm/h (às 12:00h)")
    
    print("\nPrimeiras 12 horas do hietograma desagregado:")
    for _, row in df_disagg.head(12).iterrows():
        ts_str = str(row['timestamp'])
        print(f"   Timestamp: {ts_str} -> {row['precipitation_mm_h']:5.2f} mm/h")

    print("\n>>> PROCESSO CONCLUÍDO COM SUCESSO!")

if __name__ == '__main__':
    demo_real_rainfall_workflow()
