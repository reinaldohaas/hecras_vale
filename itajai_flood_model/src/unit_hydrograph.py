"""
Módulo de Geração de Hidrograma Unitário e Transformação Chuva-Vazão.
Desacoplado de dados externos. Suporta:
- SCS Curvilíneo (NRCS NEH-4)
- SCS Triangular
- Hidrograma Unitário Sintético Padrão
"""

import numpy as np
from typing import List, Tuple, Optional, Union

# Curva adimensional padrão do SCS (t/tp vs q/qp) - NRCS National Engineering Handbook NEH-4
SCS_NEH4_RATIOS = np.array([
    [0.0, 0.000], [0.1, 0.030], [0.2, 0.100], [0.3, 0.190], [0.4, 0.310],
    [0.5, 0.470], [0.6, 0.660], [0.7, 0.820], [0.8, 0.930], [0.9, 0.990],
    [1.0, 1.000], [1.1, 0.990], [1.2, 0.930], [1.3, 0.860], [1.4, 0.780],
    [1.5, 0.680], [1.6, 0.560], [1.7, 0.460], [1.8, 0.390], [1.9, 0.330],
    [2.0, 0.280], [2.2, 0.207], [2.4, 0.147], [2.6, 0.107], [2.8, 0.077],
    [3.0, 0.055], [3.5, 0.029], [4.0, 0.015], [4.5, 0.007], [5.0, 0.000]
])

def scs_curve_number_excess(p_series_mm: Union[List[float], np.ndarray], cn: float = 76.0) -> np.ndarray:
    """
    Calcula a chuva efetiva incremental (Pe) pelo método do SCS Curve Number (CN).
    
    Equações:
        S = (25400 / CN) - 254   (Armazenamento potencial máximo em mm)
        Ia = 0.2 * S             (Perdas iniciais em mm)
        Pe_acum = (P_acum - Ia)^2 / (P_acum - Ia + S) para P_acum > Ia, senão 0
        Pe_inc[t] = Pe_acum[t] - Pe_acum[t-1]
    """
    p_arr = np.asarray(p_series_mm, dtype=float)
    p_acum = np.cumsum(p_arr)
    
    s_val = (25400.0 / cn) - 254.0
    ia_val = 0.2 * s_val
    
    pe_acum = np.where(p_acum > ia_val, ((p_acum - ia_val) ** 2) / (p_acum - ia_val + s_val), 0.0)
    pe_inc = np.diff(np.insert(pe_acum, 0, 0.0))
    return np.maximum(0.0, pe_inc)


class UnitHydrograph:
    """
    Classe para geração de Hidrogramas Unitários e convolução com chuva efetiva.
    """
    def __init__(self, area_km2: float, tc_hours: float, method: str = 'scs_curvilinear'):
        """
        Parâmetros:
            area_km2: Área de drenagem da bacia em km²
            tc_hours: Tempo de concentração da bacia em horas
            method: 'scs_curvilinear' | 'scs_triangular' | 'synthetic_standard'
        """
        self.area_km2 = float(area_km2)
        self.tc_hours = float(tc_hours)
        self.method = method.lower()
        
    def generate_ordinates(self, dt_hours: float = 1.0, duration_hours: float = 48.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Gera as ordenadas do Hidrograma Unitário para 1 mm de chuva efetiva.
        
        Retorna:
            (tempo_horas, vazao_unit_m3s)
        """
        # Tempo de pico do hidrograma: tp = dt/2 + 0.6 * tc
        tp = 0.5 * dt_hours + 0.6 * self.tc_hours
        tp = max(0.5, tp)
        
        # Vazão de pico unitária (m³/s por mm de chuva)
        # qp = (2.08 * A * 1mm) / tp
        qp = (2.08 * self.area_km2 * 1.0) / tp
        
        t_eval = np.arange(0, duration_hours + dt_hours, dt_hours)
        
        if self.method == 'scs_curvilinear':
            t_dim = SCS_NEH4_RATIOS[:, 0]
            q_dim = SCS_NEH4_RATIOS[:, 1]
            t_over_tp = t_eval / tp
            q_unit = np.interp(t_over_tp, t_dim, q_dim, right=0.0) * qp
            
        elif self.method == 'scs_triangular':
            # Tempo de base tb = 2.67 * tp
            tb = 2.67 * tp
            q_unit = np.where(
                t_eval <= tp,
                qp * (t_eval / tp),
                np.where(t_eval <= tb, qp * (tb - t_eval) / (tb - tp), 0.0)
            )
        else: # synthetic_standard (Gaussiano assimétrico)
            sigma = tp / 2.5
            q_unit = qp * np.exp(-0.5 * ((t_eval - tp) / sigma) ** 2)
            
        # Normalização rigorosa de volume: Integral(Q) * dt = Area * 1 mm
        vol_teorico_m3 = self.area_km2 * 1e6 * 0.001 # 1 mm = 0.001 m
        vol_calc_m3 = np.sum(q_unit) * (dt_hours * 3600.0)
        if vol_calc_m3 > 0:
            q_unit = q_unit * (vol_teorico_m3 / vol_calc_m3)
            
        return t_eval, q_unit

    def convolve(self, pe_inc_mm: Union[List[float], np.ndarray], dt_hours: float = 1.0, base_flow: float = 10.0, total_hours: int = 48) -> np.ndarray:
        """
        Aplica a convolução linear: Q(t) = sum_k ( Pe[k] * U[t - k*dt] ) + Qbase
        """
        pe_arr = np.asarray(pe_inc_mm, dtype=float)
        _, u_ordinates = self.generate_ordinates(dt_hours=dt_hours, duration_hours=float(total_hours))
        
        q_direct = np.convolve(pe_arr, u_ordinates)[:total_hours + 1]
        
        # Garantir comprimento fixo
        if len(q_direct) < total_hours + 1:
            q_padded = np.zeros(total_hours + 1)
            q_padded[:len(q_direct)] = q_direct
            q_direct = q_padded
            
        return q_direct + base_flow
