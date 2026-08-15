"""
Curvas-Chave Estimadas por Hidráulica Fluvial (Manning / Seção Transversal / DEM):
Utilizadas para pontos e seções sem régua linimétrica oficial ou para análise de seções transversais
naturais ao longo das 10 calhas da bacia.

Utiliza o Método dos Subtrechos / Seção Composta Dividida (Divided Channel Method - Padrão HEC-RAS / Ven Te Chow):
Q_total(H) = Q_leito_principal(H) + Q_planicie_inundacao(H)
onde cada subseção calcula sua área molhada A_i, perímetro P_i, raio hidráulico R_i = A_i/P_i e vazão Q_i.
Garante monotonicidade física estrita e continuidade na transição de extravasamento.
"""

from typing import Dict, Any, Optional, Union, List, Tuple
import numpy as np
from .base import BaseRatingCurve, CurveType

class CrossSectionGeometry:
    """
    Representação geométrica da seção transversal do rio (Trapezoidal / Composta / Perfil DEM).
    """
    def __init__(self, bottom_width_b_m: float, side_slope_z: float,
                 bankfull_depth_m: float, floodplain_width_m: float = 0.0,
                 manning_n_main: float = 0.038, manning_n_floodplain: float = 0.065):
        self.b = float(bottom_width_b_m)
        self.z = float(side_slope_z) # talude horizontal : 1 vertical (ex: 1.5)
        self.h_bankfull = float(bankfull_depth_m)
        self.b_floodplain = float(floodplain_width_m)
        self.n_main = float(manning_n_main)
        self.n_fp = float(manning_n_floodplain)

    def compute_hydraulics(self, depth_h: float, sqrt_s0: float = 0.02) -> Dict[str, float]:
        """
        Calcula os parâmetros hidráulicos e a vazão total pelo Método das Subseções Divididas (HEC-RAS).
        """
        h = max(0.0, float(depth_h))
        if h <= 0.0:
            return {
                'area_total_m2': 0.0,
                'perimeter_total_m': 0.0,
                'radius_eff_m': 0.0,
                'top_width_m': 0.0,
                'q_total_m3s': 0.0,
                'q_main_m3s': 0.0,
                'q_fp_m3s': 0.0,
                'n_eff': self.n_main
            }
            
        if h <= self.h_bankfull:
            # Escoamento confinado na calha principal
            area_main = self.b * h + self.z * (h ** 2)
            perim_main = self.b + 2.0 * h * np.sqrt(1.0 + self.z ** 2)
            top_width = self.b + 2.0 * self.z * h
            r_main = area_main / max(0.001, perim_main)
            q_main = (1.0 / self.n_main) * area_main * (r_main ** (2.0/3.0)) * sqrt_s0
            
            return {
                'area_total_m2': float(area_main),
                'perimeter_total_m': float(perim_main),
                'radius_eff_m': float(r_main),
                'top_width_m': float(top_width),
                'q_total_m3s': float(q_main),
                'q_main_m3s': float(q_main),
                'q_fp_m3s': 0.0,
                'n_eff': self.n_main
            }
        else:
            # Escoamento composto dividido: Calha principal cheia + Planície
            # 1. Calha principal (altura total h)
            area_main = self.b * h + self.z * (h ** 2)
            # No método de interface vertical suave, o perímetro de atrito sólido é a calha
            perim_main = self.b + 2.0 * self.h_bankfull * np.sqrt(1.0 + self.z ** 2)
            r_main = area_main / max(0.001, perim_main)
            q_main = (1.0 / self.n_main) * area_main * (r_main ** (2.0/3.0)) * sqrt_s0
            
            # 2. Planície de inundação (camada superior de espessura h - h_bankfull)
            h_fp = h - self.h_bankfull
            area_fp = self.b_floodplain * h_fp
            perim_fp = self.b_floodplain + 2.0 * h_fp
            r_fp = area_fp / max(0.001, perim_fp)
            q_fp = (1.0 / self.n_fp) * area_fp * (r_fp ** (2.0/3.0)) * sqrt_s0
            
            area_total = area_main + area_fp
            perim_total = perim_main + perim_fp
            w_top_main = self.b + 2.0 * self.z * h
            top_width = w_top_main + self.b_floodplain
            q_total = q_main + q_fp
            
            r_eff = area_total / max(0.001, perim_total)
            n_eff = (1.0 / q_total) * area_total * (r_eff ** (2.0/3.0)) * sqrt_s0 if q_total > 0 else self.n_main
            
            return {
                'area_total_m2': float(area_total),
                'perimeter_total_m': float(perim_total),
                'radius_eff_m': float(r_eff),
                'top_width_m': float(top_width),
                'q_total_m3s': float(q_total),
                'q_main_m3s': float(q_main),
                'q_fp_m3s': float(q_fp),
                'n_eff': float(n_eff)
            }


class HydraulicRatingCurve(BaseRatingCurve):
    """
    Curva-chave estimada a partir da equação de Manning e geometria transversal.
    """
    def __init__(self, station_id: str, name: str, river: str,
                 geometry: CrossSectionGeometry, bed_slope_s0: float,
                 datum_z0_m: float = 0.0, h_min: float = 0.1, h_max: float = 16.0,
                 source: str = "Estimativa Hidráulica Manning / DEM",
                 metadata: Optional[Dict[str, Any]] = None):
        super().__init__(station_id, name, river, CurveType.ESTIMATED_HYDRAULIC,
                         datum_z0_m, h_min, h_max, source, metadata)
        self.geometry = geometry
        self.bed_slope_s0 = max(0.00001, float(bed_slope_s0))
        self.sqrt_s0 = np.sqrt(self.bed_slope_s0)

    def to_flow(self, h: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Aplica o método de Manning com subseções compostas."""
        h_arr = np.asarray(h, dtype=float)
        is_scalar = (h_arr.ndim == 0)
        h_flat = np.atleast_1d(h_arr)
        
        q_out = np.zeros_like(h_flat)
        for i, val in enumerate(h_flat):
            if val <= 0.0:
                q_out[i] = 0.0
            else:
                hyd = self.geometry.compute_hydraulics(val, self.sqrt_s0)
                q_out[i] = max(0.0, hyd['q_total_m3s'])
                
        if is_scalar:
            return float(q_out[0])
        return q_out

    def to_stage(self, q: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Inversão numérica de vazão para nível através de busca rápida interpolada."""
        q_arr = np.asarray(q, dtype=float)
        is_scalar = (q_arr.ndim == 0)
        q_flat = np.atleast_1d(q_arr)
        
        h_grid = np.linspace(self.h_min, self.h_max, 500)
        q_grid = self.to_flow(h_grid)
        
        h_out = np.interp(q_flat, q_grid, h_grid)
        
        if is_scalar:
            return float(h_out[0])
        return h_out

    def get_section_details(self, depth_h: float) -> Dict[str, Any]:
        """Retorna os detalhes geométricos e hidráulicos para visualização gráfica."""
        hyd = self.geometry.compute_hydraulics(depth_h, self.sqrt_s0)
        return {
            'depth_h_m': float(depth_h),
            'flow_q_m3s': hyd['q_total_m3s'],
            'wet_area_m2': hyd['area_total_m2'],
            'wetted_perimeter_m': hyd['perimeter_total_m'],
            'hydraulic_radius_m': hyd['radius_eff_m'],
            'top_width_m': hyd['top_width_m'],
            'manning_n_eff': hyd['n_eff'],
            'bed_slope_s0': self.bed_slope_s0
        }
