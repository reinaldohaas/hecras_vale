"""
Módulo de Curva-Chave (Rating Curve):
Conversão bidirecional entre Vazão Q (m³/s) e Nível/Cota H (m).
Diferenciação estrita entre Curvas Oficiais (ANA/CEOPS) e Curvas Estimadas (Manning/Seção DEM).
"""

from enum import Enum
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, Union, Tuple, List
import numpy as np

class CurveType(Enum):
    OFFICIAL_OBSERVED = "CURVA OFICIAL / OBSERVADA (ANA/CEOPS)"
    ESTIMATED_HYDRAULIC = "CURVA ESTIMADA (Manning / Seção DEM)"

class BaseRatingCurve(ABC):
    """
    Classe base abstrata para relações Vazão x Cota (Q-H).
    """
    def __init__(self, station_id: str, name: str, river: str,
                 curve_type: CurveType, datum_z0_m: float = 0.0,
                 h_min: float = 0.0, h_max: float = 20.0,
                 source: str = "", metadata: Optional[Dict[str, Any]] = None):
        self.station_id = station_id
        self.name = name
        self.river = river
        self.curve_type = curve_type
        self.datum_z0_m = datum_z0_m  # Cota absoluta do zero da régua (m no nível do mar)
        self.h_min = h_min
        self.h_max = h_max
        self.source = source
        self.metadata = metadata or {}

    @abstractmethod
    def to_stage(self, q: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calcula o nível de régua H (m) dada a vazão Q (m³/s)."""
        pass

    @abstractmethod
    def to_flow(self, h: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calcula a vazão Q (m³/s) dado o nível de régua H (m)."""
        pass

    def to_absolute_water_level(self, h: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Converte nível de régua H (m) em cota altimétrica absoluta (m acima do nível do mar)."""
        return np.asarray(h) + self.datum_z0_m

    def get_curve_table(self, n_points: int = 50) -> Dict[str, np.ndarray]:
        """Gera tabela de valores discretos H x Q ao longo da faixa válida."""
        h_vals = np.linspace(self.h_min, self.h_max, n_points)
        q_vals = self.to_flow(h_vals)
        z_vals = self.to_absolute_water_level(h_vals)
        return {
            'stage_h_m': np.round(h_vals, 3),
            'flow_q_m3s': np.round(q_vals, 2),
            'abs_level_z_m': np.round(z_vals, 3)
        }

    def summary(self) -> Dict[str, Any]:
        """Resumo informativo com identificação estrita de tipo de curva."""
        return {
            'station_id': self.station_id,
            'name': self.name,
            'river': self.river,
            'curve_type': self.curve_type.value,
            'is_official': (self.curve_type == CurveType.OFFICIAL_OBSERVED),
            'datum_z0_m': self.datum_z0_m,
            'validity_range_h_m': [self.h_min, self.h_max],
            'source': self.source
        }
