"""
Módulo de Classificação de Lâmina d'Água e Zoneamento de Risco de Inundação:
Categoriza a profundidade h(x,y) em faixas padronizadas de risco da Defesa Civil / CPRM.
"""

from enum import Enum
from typing import Dict, Any, List, Tuple
import numpy as np

class FloodDepthClass(Enum):
    DRY = "Seco (h < 0.05m)"
    LOW = "Lâmina Baixa (0.05m <= h < 0.50m - Tráfego difícil / Alagamento leve)"
    MEDIUM = "Lâmina Média (0.50m <= h < 1.50m - Invasão de residências / Risco a pedestres)"
    HIGH = "Lâmina Alta (1.50m <= h < 3.00m - Resgate de barco / Risco estrutural)"
    VERY_HIGH = "Lâmina Extrema (h >= 3.00m - Cobertura de telhados / Colapso)"

class DepthClassifier:
    """Classifica e resume os dados de lâmina d'água em zonas de risco."""

    @staticmethod
    def classify_depths(depth_grid: np.ndarray) -> Dict[str, Any]:
        """
        Calcula a distribuição percentual e em km² de cada classe de inundação.
        """
        d = np.asarray(depth_grid, dtype=float)
        
        mask_dry = (d < 0.05)
        mask_low = (d >= 0.05) & (d < 0.50)
        mask_med = (d >= 0.50) & (d < 1.50)
        mask_high = (d >= 1.50) & (d < 3.00)
        mask_very_high = (d >= 3.00)

        total_wet = np.sum(~mask_dry)
        
        return {
            'count_low': int(np.sum(mask_low)),
            'count_medium': int(np.sum(mask_med)),
            'count_high': int(np.sum(mask_high)),
            'count_very_high': int(np.sum(mask_very_high)),
            'total_wet_cells': int(total_wet),
            'pct_low': float(np.round((np.sum(mask_low) / max(1, total_wet)) * 100, 1)),
            'pct_medium': float(np.round((np.sum(mask_med) / max(1, total_wet)) * 100, 1)),
            'pct_high': float(np.round((np.sum(mask_high) / max(1, total_wet)) * 100, 1)),
            'pct_very_high': float(np.round((np.sum(mask_very_high) / max(1, total_wet)) * 100, 1)),
        }
