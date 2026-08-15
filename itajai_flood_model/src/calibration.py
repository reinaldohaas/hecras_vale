"""
Módulo de Calibração e Validação de Modelos Hidrológicos (HydrographValidator).
Calcula métricas estatísticas formais de aderência:
- Erro no Pico (m³/s e %)
- Diferença no Horário do Pico (horas)
- RMSE (Root Mean Square Error)
- NSE (Nash-Sutcliffe Efficiency)
- Erro de Volume Total (%)
"""

import numpy as np
import pandas as pd
from typing import Dict, Any, Union, List

class HydrographValidator:
    """
    Avalia a qualidade do ajuste entre o hidrograma simulado e o observado.
    """
    @staticmethod
    def calculate_metrics(q_sim: Union[List[float], np.ndarray], q_obs: Union[List[float], np.ndarray], dt_hours: float = 1.0) -> Dict[str, Any]:
        """
        Calcula as métricas estatísticas entre Q_simulado e Q_observado.
        """
        sim = np.asarray(q_sim, dtype=float)
        obs = np.asarray(q_obs, dtype=float)
        
        # Alinhar tamanhos
        n = min(len(sim), len(obs))
        sim = sim[:n]
        obs = obs[:n]
        
        # Picos e tempos de pico
        p_sim = float(np.max(sim))
        t_sim = float(np.argmax(sim) * dt_hours)
        
        p_obs = float(np.max(obs))
        t_obs = float(np.argmax(obs) * dt_hours)
        
        # Erro de Pico
        delta_peak = p_sim - p_obs
        peak_pct_error = (delta_peak / p_obs) * 100.0 if p_obs > 0 else 0.0
        
        # Diferença no horário do pico
        delta_t_peak = t_sim - t_obs
        
        # RMSE
        rmse = float(np.sqrt(np.mean((sim - obs) ** 2)))
        
        # NSE (Nash-Sutcliffe Efficiency)
        obs_mean = float(np.mean(obs))
        denominator_nse = np.sum((obs - obs_mean) ** 2)
        numerator_nse = np.sum((obs - sim) ** 2)
        
        if denominator_nse > 0:
            nse = float(1.0 - (numerator_nse / denominator_nse))
        else:
            nse = 1.0 if numerator_nse == 0 else -999.0
            
        # Balanço Volumétrico (m³ e hm³)
        vol_sim_m3 = float(np.sum(sim) * (dt_hours * 3600.0))
        vol_obs_m3 = float(np.sum(obs) * (dt_hours * 3600.0))
        
        vol_sim_hm3 = vol_sim_m3 / 1e6
        vol_obs_hm3 = vol_obs_m3 / 1e6
        
        vol_pct_error = ((vol_sim_m3 - vol_obs_m3) / vol_obs_m3) * 100.0 if vol_obs_m3 > 0 else 0.0
        
        return {
            'num_points': n,
            'peak_sim_m3s': round(p_sim, 2),
            'peak_obs_m3s': round(p_obs, 2),
            'peak_error_m3s': round(delta_peak, 2),
            'peak_error_pct': round(peak_pct_error, 2),
            't_peak_sim_h': round(t_sim, 2),
            't_peak_obs_h': round(t_obs, 2),
            't_peak_diff_h': round(delta_t_peak, 2),
            'rmse_m3s': round(rmse, 2),
            'nse': round(nse, 4),
            'volume_sim_hm3': round(vol_sim_hm3, 2),
            'volume_obs_hm3': round(vol_obs_hm3, 2),
            'volume_error_pct': round(vol_pct_error, 2)
        }

    @staticmethod
    def print_validation_report(metrics: Dict[str, Any]):
        """
        Imprime um relatório formatado das métricas de calibração/validação.
        """
        print("=" * 60)
        print("         RELATÓRIO DE VALIDAÇÃO DO HIDROGRAMA")
        print("=" * 60)
        print(f"Número de passos avaliados:   {metrics['num_points']}")
        print(f"Pico Observado:               {metrics['peak_obs_m3s']} m³/s (em t = {metrics['t_peak_obs_h']}h)")
        print(f"Pico Simulado:                {metrics['peak_sim_m3s']} m³/s (em t = {metrics['t_peak_sim_h']}h)")
        print(f"Erro no Pico:                 {metrics['peak_error_m3s']} m³/s ({metrics['peak_error_pct']}%)")
        print(f"Diferença Horário do Pico:    {metrics['t_peak_diff_h']} horas")
        print(f"RMSE:                         {metrics['rmse_m3s']} m³/s")
        print(f"Coeficiente NSE (Nash):       {metrics['nse']}")
        print(f"Volume Observado:             {metrics['volume_obs_hm3']} hm³")
        print(f"Volume Simulado:              {metrics['volume_sim_hm3']} hm³")
        print(f"Erro Volumétrico:             {metrics['volume_error_pct']}%")
        print("=" * 60)
