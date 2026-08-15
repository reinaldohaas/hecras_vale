"""
Módulo de Propagação por Muskingum-Cunge.
Deriva os parâmetros K e X a partir da geometria da calha, declividade e rugosidade de Manning:
- Celeridade da onda (c = dQ/dA)
- Coeficiente de difusão hidráulica (D)
- K = dx / c
- X = 0.5 * (1 - Q / (B * S0 * c * dx))
"""

import numpy as np
from typing import Tuple, Dict, Any, Union, List, Optional
from .muskingum import MuskingumReach

class MuskingumCungeReach:
    """
    Trecho fluvial propagado pelo método de Muskingum-Cunge com parâmetros físicos.
    """
    def __init__(self, reach_id: Union[int, str], name: str, length_km: float, slope_m_km: float, 
                 width_m: float = 35.0, manning_n: float = 0.038, reference_q_m3s: float = 200.0, dt_hours: float = 1.0):
        """
        Parâmetros:
            reach_id: Identificador do trecho
            name: Nome descritivo
            length_km: Extensão do trecho (dx) em km
            slope_m_km: Declividade média do fundo (S0) em m/km
            width_m: Largura média da calha retangular/trapezoidal (B) em metros
            manning_n: Coeficiente de rugosidade de Manning
            reference_q_m3s: Vazão de referência para celeridade (m³/s)
            dt_hours: Intervalo de cálculo em horas
        """
        self.reach_id = reach_id
        self.name = name
        self.length_km = float(length_km)
        self.slope_m_km = float(slope_m_km)
        self.width_m = float(width_m)
        self.manning_n = float(manning_n)
        self.reference_q_m3s = float(reference_q_m3s)
        self.dt_hours = float(dt_hours)
        
        # Converter unidades
        self.dx_m = self.length_km * 1000.0
        self.s0 = max(1e-5, self.slope_m_km / 1000.0) # m/m
        
        # Calcular K e X físicos
        self.k_hours, self.x_param, self.celerity_ms = self._derive_cunge_parameters()
        
        # Instanciar o propagador de Muskingum com os K e X calculados
        self.muskingum_solver = MuskingumReach(
            reach_id=self.reach_id,
            name=self.name,
            k_hours=self.k_hours,
            x_param=self.x_param,
            dt_hours=self.dt_hours
        )
        
    def _derive_cunge_parameters(self) -> Tuple[float, float, float]:
        """
        Calcula K (horas), X (adimensional) e celeridade c (m/s) via Manning e Cunge.
        """
        b = self.width_m
        n = self.manning_n
        s0 = self.s0
        q = self.reference_q_m3s
        dx = self.dx_m
        
        # Profundidade normal (y) estimada para seção larga: Q = (1/n) * B * y^(5/3) * S0^(1/2)
        # y = [ (Q * n) / (B * S0^0.5) ]^(3/5)
        y_norm = ((q * n) / (b * np.sqrt(s0))) ** 0.6
        y_norm = max(0.2, y_norm)
        
        # Velocidade média v = Q / (B * y)
        v_mean = q / (b * y_norm)
        
        # Celeridade da onda cinemática: c = 5/3 * v para canais largos
        c = (5.0 / 3.0) * v_mean
        c = max(0.5, min(8.0, c)) # limites físicos razoáveis (m/s)
        
        # Tempo de trânsito K = dx / c (convertido para horas)
        k_sec = dx / c
        k_hours = k_sec / 3600.0
        
        # Parâmetro de difusão X = 0.5 * (1 - Q / (B * S0 * c * dx))
        num_x = q
        den_x = b * s0 * c * dx
        
        if den_x > 0:
            x_calc = 0.5 * (1.0 - (num_x / den_x))
        else:
            x_calc = 0.20
            
        x_param = max(0.05, min(0.48, x_calc))
        return k_hours, x_param, c

    def get_cunge_summary(self) -> Dict[str, Any]:
        """
        Retorna parâmetros hidráulicos derivados pelo método Cunge.
        """
        return {
            'reach_id': self.reach_id,
            'name': self.name,
            'length_km': self.length_km,
            'slope_m_km': self.slope_m_km,
            'width_m': self.width_m,
            'manning_n': self.manning_n,
            'celerity_m_s': round(self.celerity_ms, 2),
            'derived_K_hours': round(self.k_hours, 2),
            'derived_X': round(self.x_param, 3),
            'dt_hours': self.dt_hours
        }

    def route(self, inflow: Union[List[float], np.ndarray], lateral_inflow: Optional[Union[List[float], np.ndarray]] = None) -> np.ndarray:
        """
        Propaga o hidrograma de vazão pelo solver de Muskingum-Cunge.
        """
        return self.muskingum_solver.route(inflow=inflow, lateral_inflow=lateral_inflow)
