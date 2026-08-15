"""
Módulo de Propagação de Ondas de Cheia pelo Método Clássico de Muskingum.
Implementa:
- Cálculo automatizado de coeficientes (C0, C1, C2)
- Verificação de estabilidade numérica e conservação de massa (C0 + C1 + C2 = 1.0)
- Sub-passos automáticos quando a condição 2KX < dt < 2K(1-X) for violada
"""

import numpy as np
from typing import Tuple, Dict, Any, Union, List, Optional

class MuskingumReach:
    """
    Representa um trecho fluvial propagado pelo método de Muskingum clássico.
    """
    def __init__(self, reach_id: Union[int, str], name: str, k_hours: float, x_param: float = 0.20, dt_hours: float = 1.0):
        """
        Parâmetros:
            reach_id: Identificador único do trecho
            name: Nome descritivo do trecho
            k_hours: Constante de tempo de propagação (K) em horas
            x_param: Fator de ponderação de atenuação (X) [0.0 a 0.5]
            dt_hours: Intervalo de tempo de cálculo (dt) em horas
        """
        self.reach_id = reach_id
        self.name = name
        self.k_hours = float(k_hours)
        self.x_param = float(x_param)
        self.dt_hours = float(dt_hours)
        
        # Validar parâmetro X
        if not (0.0 <= self.x_param <= 0.5):
            raise ValueError(f"Parâmetro X deve estar no intervalo [0.0, 0.5]. Recebido: {self.x_param}")
        if self.k_hours <= 0:
            raise ValueError(f"Constante K deve ser positiva. Recebido: {self.k_hours}")
            
        self.c0, self.c1, self.c2, self.num_substeps = self._compute_coefficients()
        
    def _compute_coefficients(self) -> Tuple[float, float, float, int]:
        """
        Calcula C0, C1 e C2 com verificação de estabilidade e sub-discretização automática se necessário.
        Critério clássico: 2KX <= dt <= 2K(1-X)
        """
        k = self.k_hours
        x = self.x_param
        dt = self.dt_hours
        
        # Determinar número de sub-passos necessários para manter C0 >= 0
        # Se 2KX > dt, divide o trecho em N sub-trechos com K_sub = K/N tal que 2*(K/N)*X <= dt
        if 2.0 * k * x > dt and x > 0:
            n_sub = int(np.ceil((2.0 * k * x) / dt))
        else:
            n_sub = 1
            
        k_sub = k / n_sub
        denom = 2.0 * k_sub * (1.0 - x) + dt
        c0 = (dt - 2.0 * k_sub * x) / denom
        c1 = (dt + 2.0 * k_sub * x) / denom
        c2 = (2.0 * k_sub * (1.0 - x) - dt) / denom
        
        return c0, c1, c2, n_sub

    def get_stability_report(self) -> Dict[str, Any]:
        """
        Retorna o diagnóstico de estabilidade numérica do trecho.
        """
        k_eff = self.k_hours / self.num_substeps
        cond_inf = 2.0 * k_eff * self.x_param
        cond_sup = 2.0 * k_eff * (1.0 - self.x_param)
        is_stable = (self.c0 >= 0.0) and (self.c1 >= 0.0) and (self.c2 >= 0.0)
        sum_c = self.c0 + self.c1 + self.c2
        
        return {
            'reach_id': self.reach_id,
            'name': self.name,
            'K_total': self.k_hours,
            'X': self.x_param,
            'dt': self.dt_hours,
            'sub_reaches': self.num_substeps,
            'C0': round(self.c0, 6),
            'C1': round(self.c1, 6),
            'C2': round(self.c2, 6),
            'sum_coefficients': round(sum_c, 6),
            'limit_lower_2KX': round(cond_inf, 4),
            'limit_upper_2K_1_minus_X': round(cond_sup, 4),
            'stable_positive_coeffs': is_stable
        }

    def route(self, inflow: Union[List[float], np.ndarray], lateral_inflow: Optional[Union[List[float], np.ndarray]] = None) -> np.ndarray:
        """
        Propaga o hidrograma de entrada ao longo do trecho.
        
        Parâmetros:
            inflow: Série temporal da vazão afluente (m³/s)
            lateral_inflow: Série temporal de vazão lateral incremental (m³/s, opcional)
            
        Retorna:
            outflow: Série temporal da vazão efluente na saída do trecho (m³/s)
        """
        i_curr = np.asarray(inflow, dtype=float).copy()
        n_steps = len(i_curr)
        
        # Propagação sequencial pelos sub-passos
        for _ in range(self.num_substeps):
            q_out = np.zeros(n_steps)
            q_out[0] = i_curr[0]
            
            for t in range(0, n_steps - 1):
                val = self.c0 * i_curr[t + 1] + self.c1 * i_curr[t] + self.c2 * q_out[t]
                q_out[t + 1] = max(0.0, val)
                
            i_curr = q_out
            
        # Adicionar contribuição lateral se fornecida
        if lateral_inflow is not None:
            lat_arr = np.asarray(lateral_inflow, dtype=float)
            i_curr += lat_arr[:n_steps]
            
        return i_curr
