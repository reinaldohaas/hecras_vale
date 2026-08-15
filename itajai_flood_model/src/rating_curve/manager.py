"""
Gerenciador de Curvas-Chave (RatingCurveManager):
Catálogo unificado das curvas de calibração para todas as estações e trechos da Bacia do Rio Itajaí.
Permite alternância, comparação e consulta com rastreabilidade formal da origem dos dados.
"""

from typing import Dict, Any, Optional, Union, List
import numpy as np
from .base import BaseRatingCurve, CurveType
from .official import (
    create_blumenau_official_curve,
    create_rio_do_sul_official_curve,
    create_brusque_official_curve,
    create_indaial_official_curve
)
from .hydraulic import CrossSectionGeometry, HydraulicRatingCurve

class RatingCurveManager:
    """
    Gerencia e fornece acesso a todas as curvas-chave da bacia.
    """
    def __init__(self):
        self._curves: Dict[str, BaseRatingCurve] = {}
        self._load_default_curves()

    def _load_default_curves(self):
        """Carrega o conjunto padrão de curvas oficiais e estimadas."""
        # 1. Curvas Oficiais
        self.register_curve('blumenau', create_blumenau_official_curve())
        self.register_curve('rio_do_sul', create_rio_do_sul_official_curve())
        self.register_curve('brusque', create_brusque_official_curve())
        self.register_curve('indaial', create_indaial_official_curve())

        # 2. Curvas Estimadas (Manning / Geometria de Seção do DEM)
        # Ibirama (Confluência Hercílio / Alto Açú)
        geo_ibirama = CrossSectionGeometry(
            bottom_width_b_m=65.0,
            side_slope_z=1.2,
            bankfull_depth_m=6.0,
            floodplain_width_m=80.0,
            manning_n_main=0.040,
            manning_n_floodplain=0.070
        )
        curve_ibirama = HydraulicRatingCurve(
            station_id="83300000_EST",
            name="Ibirama (Confluência Hercílio - Estimada)",
            river="Rio Itajaí-Açu",
            geometry=geo_ibirama,
            bed_slope_s0=0.00095, # Declividade DEM ~ 0.95 m/km
            datum_z0_m=118.50,
            h_min=0.30,
            h_max=16.00,
            source="Curva Estimada (Manning / Seção DEM 30m Copernicus)",
            metadata={'is_estimated': True}
        )
        self.register_curve('ibirama', curve_ibirama)

        # Itajaí Foz (Exutório & Baixo Vale)
        geo_itajai = CrossSectionGeometry(
            bottom_width_b_m=140.0,
            side_slope_z=2.0,
            bankfull_depth_m=5.0,
            floodplain_width_m=200.0,
            manning_n_main=0.032,
            manning_n_floodplain=0.055
        )
        curve_itajai = HydraulicRatingCurve(
            station_id="83900000_EST",
            name="Itajaí (Foz & Canal Oficial - Estimada)",
            river="Rio Itajaí-Açu / Foz",
            geometry=geo_itajai,
            bed_slope_s0=0.00015, # Declividade Baixo Vale ~ 0.15 m/km
            datum_z0_m=0.50,
            h_min=0.20,
            h_max=8.00,
            source="Curva Estimada (Manning / Seção Foz DEM 30m)",
            metadata={'is_estimated': True}
        )
        self.register_curve('itajai_foz', curve_itajai)

        # Também adicionamos a versão estimada de Blumenau para comparação de métodos
        geo_blumenau = CrossSectionGeometry(
            bottom_width_b_m=110.0,
            side_slope_z=1.5,
            bankfull_depth_m=7.5,
            floodplain_width_m=160.0,
            manning_n_main=0.036,
            manning_n_floodplain=0.065
        )
        curve_blumenau_est = HydraulicRatingCurve(
            station_id="83700000_EST",
            name="Blumenau Centro (Estimada por Manning)",
            river="Rio Itajaí-Açu",
            geometry=geo_blumenau,
            bed_slope_s0=0.00035, # Declividade Médio Vale ~ 0.35 m/km
            datum_z0_m=11.20,
            h_min=0.20,
            h_max=18.00,
            source="Estimativa Hidráulica por Manning (Para fins de calibração)",
            metadata={'is_estimated': True}
        )
        self.register_curve('blumenau_estimated', curve_blumenau_est)

    def register_curve(self, key: str, curve: BaseRatingCurve):
        """Registra uma curva-chave no catálogo."""
        self._curves[key] = curve

    def get_curve(self, key: str) -> BaseRatingCurve:
        """Obtém uma curva por chave."""
        if key not in self._curves:
            raise KeyError(f"Curva-chave '{key}' não encontrada no gerenciador. Disponíveis: {list(self._curves.keys())}")
        return self._curves[key]

    def flow_to_stage(self, key: str, q: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Converte vazão em nível de régua para uma estação específica."""
        curve = self.get_curve(key)
        return curve.to_stage(q)

    def stage_to_flow(self, key: str, h: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Converte nível de régua em vazão para uma estação específica."""
        curve = self.get_curve(key)
        return curve.to_flow(h)

    def list_stations(self) -> List[Dict[str, Any]]:
        """Lista todas as estações com seus metadados e tipo de curva."""
        res = []
        for k, c in self._curves.items():
            info = c.summary()
            info['key'] = k
            res.append(info)
        return res
