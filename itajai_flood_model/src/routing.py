"""
Módulo de Roteamento Hidrológico Multi-Trechos (FloodRouter).
Executa a propagação sequencial da onda de cheia trecho a trecho:
Inflow -> Routing Trecho 1 -> Outflow 1 (+ Contribuição Lateral) -> Trecho 2 -> ...
"""

import numpy as np
from typing import Dict, List, Optional, Any
from .river import RiverNetwork

class FloodRouter:
    """
    Controlador de propagação da onda de cheia pela rede hidrográfica.
    """
    def __init__(self, network: RiverNetwork, method: str = 'muskingum', dt_hours: float = 1.0):
        """
        Parâmetros:
            network: Objeto RiverNetwork com os trechos configurados
            method: 'muskingum' | 'muskingum_cunge'
            dt_hours: Passo de tempo de cálculo em horas
        """
        self.network = network
        self.method = method.lower()
        self.dt_hours = float(dt_hours)
        
    def execute_routing(self, upstream_inflow: np.ndarray, lateral_inflows: Optional[Dict[int, np.ndarray]] = None) -> Dict[str, Any]:
        """
        Executa a propagação sequencial ao longo de todos os trechos da rede.
        
        Parâmetros:
            upstream_inflow: Hidrograma de entrada em montante do Trecho 1 (m³/s)
            lateral_inflows: Dicionário {reach_id: array_vazao_lateral}
            
        Retorna:
            Dicionário com os hidrogramas em cada seção e métricas de propagação.
        """
        lateral_inflows = lateral_inflows or {}
        n_steps = len(upstream_inflow)
        
        reach_inflows: Dict[int, np.ndarray] = {}
        reach_outflows: Dict[int, np.ndarray] = {}
        reach_solvers: Dict[int, Any] = {}
        
        current_inflow = np.asarray(upstream_inflow, dtype=float).copy()
        
        # Iterar sobre os trechos ordenados
        for reach in self.network.reaches:
            rid = reach.reach_id
            reach_inflows[rid] = current_inflow.copy()
            
            # Criar solver configurado
            if self.method == 'muskingum_cunge':
                solver = reach.create_cunge_solver(dt_hours=self.dt_hours)
            else:
                solver = reach.create_muskingum_solver(dt_hours=self.dt_hours)
                
            reach_solvers[rid] = solver
            
            # Obter aporte lateral se houver
            lat_q = lateral_inflows.get(rid, None)
            
            # Propagar pelo trecho
            outflow = solver.route(inflow=current_inflow, lateral_inflow=lat_q)
            reach_outflows[rid] = outflow
            
            # O efluente deste trecho é o afluente do próximo
            current_inflow = outflow.copy()
            
        # Calcular Métricas de Propagação (Pico, Atraso e Atenuação)
        q_in_total = upstream_inflow
        q_out_total = current_inflow
        
        peak_in = float(np.max(q_in_total))
        t_peak_in = int(np.argmax(q_in_total)) * self.dt_hours
        
        peak_out = float(np.max(q_out_total))
        t_peak_out = int(np.argmax(q_out_total)) * self.dt_hours
        
        peak_reduction_m3s = peak_in - peak_out
        peak_reduction_pct = (peak_reduction_m3s / peak_in) * 100.0 if peak_in > 0 else 0.0
        lag_time_hours = t_peak_out - t_peak_in
        
        # Alargamento da onda (Largura da base a 50% do pico - Full Width at Half Maximum)
        def compute_fwhm(q_series):
            p = np.max(q_series)
            half = p * 0.5
            idx_above = np.where(q_series >= half)[0]
            if len(idx_above) > 1:
                return (idx_above[-1] - idx_above[0]) * self.dt_hours
            return 0.0
            
        fwhm_in = compute_fwhm(q_in_total)
        fwhm_out = compute_fwhm(q_out_total)
        wave_broadening_hours = fwhm_out - fwhm_in
        
        results = {
            'dt_hours': self.dt_hours,
            'time_steps': list(range(n_steps)),
            'time_hours': [round(t * self.dt_hours, 2) for t in range(n_steps)],
            'reach_inflows': reach_inflows,
            'reach_outflows': reach_outflows,
            'final_outflow': q_out_total,
            'metrics': {
                'peak_inflow_m3s': round(peak_in, 2),
                't_peak_inflow_h': round(t_peak_in, 2),
                'peak_outflow_m3s': round(peak_out, 2),
                't_peak_outflow_h': round(t_peak_out, 2),
                'peak_reduction_m3s': round(peak_reduction_m3s, 2),
                'peak_reduction_pct': round(peak_reduction_pct, 2),
                'lag_time_hours': round(lag_time_hours, 2),
                'fwhm_in_hours': round(fwhm_in, 2),
                'fwhm_out_hours': round(fwhm_out, 2),
                'wave_broadening_hours': round(wave_broadening_hours, 2)
            }
        }
        return results
