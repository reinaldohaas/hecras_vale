"""
Script Executável de Demonstração Completa: Modelo Hidrológico do Rio Itajaí-Mirim.

Execução:
    python examples/itajai_mirim_demo.py
"""

import sys
import os
import numpy as np
import pandas as pd

# Configurar caminhos do projeto
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, BASE_DIR)

from src.river import RiverNetwork
from src.unit_hydrograph import UnitHydrograph, scs_curve_number_excess
from src.routing import FloodRouter
from src.calibration import HydrographValidator
from src.visualization import FloodVisualizer
from src.mapping import RiverMapper

def run_itajai_mirim_demo():
    print("\n" + "=" * 75)
    print("   MODELO HIDROLÓGICO/HIDRODINÂMICO DO RIO ITAJAÍ-MIRIM (FASE 1)")
    print("=" * 75)
    
    # 1. Carregamento dos dados desacoplados do rio
    data_dir = os.path.join(BASE_DIR, 'data', 'itajai_mirim')
    reaches_csv = os.path.join(data_dir, 'reaches.csv')
    stations_csv = os.path.join(data_dir, 'stations.csv')
    rainfall_csv = os.path.join(data_dir, 'rainfall.csv')
    discharge_csv = os.path.join(data_dir, 'discharge.csv')
    
    print(f"\n[1] Carregando dados da rede fluvial em: {reaches_csv}")
    network = RiverNetwork.from_csv(reaches_csv_path=reaches_csv, stations_csv_path=stations_csv)
    
    df_summary = network.summary()
    print("\n--- Trechos Discretizados do Rio Itajaí-Mirim ---")
    print(df_summary.to_string(index=False))
    print(f"Extensão Total da Calha: {network.get_total_length_km():.1f} km")
    
    # 2. Geração do Hidrograma Unitário e Inflow em Montante
    print("\n[2] Gerando Hidrograma de Entrada em Montante via SCS Unit Hydrograph...")
    df_rain = pd.read_csv(rainfall_csv)
    pe_inc = scs_curve_number_excess(df_rain['rainfall_incremental_mm'].values, cn=79.0)
    
    # Bacia de montante (Alto Itajaí-Mirim até Vidal Ramos/Botuverá): Área ~ 650 km², tc ~ 8.5h
    uh_montante = UnitHydrograph(area_km2=650.0, tc_hours=8.5, method='scs_curvilinear')
    inflow_montante = uh_montante.convolve(pe_inc_mm=pe_inc, dt_hours=1.0, base_flow=18.0, total_hours=48)
    
    print(f"    - Pico de Entrada em Montante: {np.max(inflow_montante):.1f} m³/s no t = {np.argmax(inflow_montante)}h")
    print(f"    - Vazão de Base: 18.0 m³/s")
    
    # 3. Propagação Sequencial Trecho a Trecho por Muskingum
    print("\n[3] Executando Propagação Fluvial Sequencial (Muskingum)...")
    router = FloodRouter(network=network, method='muskingum', dt_hours=1.0)
    
    # Contribuições laterais simuladas para as microbacias intermediárias (Brusque e Guabiruba)
    lateral_inflows = {
        3: 15.0 * np.exp(-0.5 * ((np.arange(49) - 14.0) / 3.0) ** 2), # ribeirões de Brusque
        4: 10.0 * np.exp(-0.5 * ((np.arange(49) - 16.0) / 3.0) ** 2)  # ribeirões do Baixo Mirim
    }
    
    results = router.execute_routing(upstream_inflow=inflow_montante, lateral_inflows=lateral_inflows)
    metrics = results['metrics']
    
    print("\n--- Resultados Hidrológicos da Propagação ---")
    print(f"  • Pico de Entrada (Trecho 1):  {metrics['peak_inflow_m3s']:.1f} m³/s (em t = {metrics['t_peak_inflow_h']:.1f}h)")
    print(f"  • Pico de Saída (Foz Itajaí): {metrics['peak_outflow_m3s']:.1f} m³/s (em t = {metrics['t_peak_outflow_h']:.1f}h)")
    print(f"  • Atenuação (Redução de Pico): {metrics['peak_reduction_m3s']:.1f} m³/s (-{metrics['peak_reduction_pct']:.1f}%)")
    print(f"  • Atraso de Pico (Lag Time):   +{metrics['lag_time_hours']:.1f} horas")
    print(f"  • Alargamento da Onda (FWHM):  +{metrics['wave_broadening_hours']:.1f} horas")
    
    # 4. Validação Estatística contra Dados da Estação Brusque Centro
    print("\n[4] Realizando Validação Estatística contra Estação Brusque Centro...")
    df_obs = pd.read_csv(discharge_csv)
    q_obs_brusque = df_obs['discharge_observed_m3s'].values
    
    # Hidrograma simulado na estação de Brusque (Saída do Trecho 3)
    q_sim_brusque = results['reach_outflows'][3]
    
    val_metrics = HydrographValidator.calculate_metrics(q_sim=q_sim_brusque, q_obs=q_obs_brusque, dt_hours=1.0)
    HydrographValidator.print_validation_report(val_metrics)
    
    # 5. Geração de Gráficos e Saídas Visuais
    print("\n[5] Gerando Figuras e Visualizações em 'output_plots/'...")
    plots_dir = os.path.join(BASE_DIR, 'output_plots')
    visualizer = FloodVisualizer(output_dir=plots_dir)
    
    p1 = visualizer.plot_reach_propagation(results, save_name="1_propagacao_trechos.png")
    p2 = visualizer.plot_attenuation_and_lag(results, save_name="2_atenuacao_atraso_pico.png")
    reach_lens = [r.length_km for r in network.reaches]
    p3 = visualizer.plot_space_time_heatmap(results, reach_lengths_km=reach_lens, save_name="3_diagrama_espaco_tempo.png")
    p4 = visualizer.plot_validation_comparison(
        t_hours=results['time_hours'][:len(q_obs_brusque)],
        q_sim=q_sim_brusque[:len(q_obs_brusque)],
        q_obs=q_obs_brusque,
        metrics=val_metrics,
        save_name="4_validacao_simulado_vs_observado.png"
    )
    
    mapper = RiverMapper(stations_csv_path=stations_csv, reaches_csv_path=reaches_csv)
    p5 = mapper.plot_river_map(output_path=os.path.join(plots_dir, "5_mapa_trechos_itajai_mirim.png"))
    
    print(f"  ✓ [Figura 1] Hidrogramas por Trecho:        {p1}")
    print(f"  ✓ [Figura 2] Análise de Atenuação e Lag:     {p2}")
    print(f"  ✓ [Figura 3] Diagrama Espaço-Temporal Q(x,t):{p3}")
    print(f"  ✓ [Figura 4] Validação Simulado x Observado: {p4}")
    print(f"  ✓ [Figura 5] Mapa dos Trechos e Estações:    {p5}")
    
    print("\n" + "=" * 75)
    print("   DEMONSTRAÇÃO CONCLUÍDA COM SUCESSO!")
    print("=" * 75 + "\n")

if __name__ == '__main__':
    run_itajai_mirim_demo()
