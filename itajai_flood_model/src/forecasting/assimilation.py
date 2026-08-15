"""
Módulo de Assimilação e Correção de Previsão de Vazão (StreamflowAssimilation):
- Compara continuamente a vazão modelada (Q_modelo) com a vazão observada nas estações (Q_obs)
- Calcula métricas em tempo real (Erro, Viés, RMSE, NSE)
- Aplica correção transparente no instante 'AGORA' com decaimento exponencial suave para a previsão futura
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any, Union
from ..calibration import HydrographValidator

class StreamflowAssimilation:
    """
    Ajustador operacional de hidrograma com base em telemetria fluviométrica em tempo real.
    """
    
    @staticmethod
    def evaluate_model_vs_observed(q_sim: Union[List[float], np.ndarray],
                                  q_obs: Union[List[float], np.ndarray],
                                  dt_hours: float = 1.0) -> Dict[str, Any]:
        """Calcula métricas de desempenho entre simulado e observado."""
        return HydrographValidator.calculate_metrics(q_sim=q_sim, q_obs=q_obs, dt_hours=dt_hours)

    @staticmethod
    def apply_realtime_correction(q_forecast_raw: np.ndarray,
                                  q_observed_at_now: float,
                                  now_idx: int,
                                  relaxation_hours: float = 12.0) -> np.ndarray:
        """
        Ajusta a previsão futura a partir da discrepância observada no instante t_now.
        
        Equação:
            K_corr = Q_obs(t_now) / Q_mod(t_now)
            Q_corrigido(t) = Q_mod(t) * [ 1 + (K_corr - 1) * exp(-(t - t_now) / tau) ]
        """
        q_mod = np.asarray(q_forecast_raw, dtype=float).copy()
        n = len(q_mod)
        if now_idx < 0 or now_idx >= n:
            return q_mod
            
        q_mod_now = q_mod[now_idx]
        if q_mod_now <= 0 or q_observed_at_now <= 0:
            return q_mod
            
        k_corr = q_observed_at_now / q_mod_now
        # Limitar fator de correção entre 0.3 e 3.0 para estabilidade
        k_corr = float(np.clip(k_corr, 0.30, 3.00))
        
        q_corrected = q_mod.copy()
        
        # Para o futuro (t >= now_idx)
        for t in range(now_idx, n):
            dt = float(t - now_idx)
            decay_factor = np.exp(-dt / max(1.0, relaxation_hours))
            factor = 1.0 + (k_corr - 1.0) * decay_factor
            q_corrected[t] = max(0.0, q_mod[t] * factor)
            
        return q_corrected
