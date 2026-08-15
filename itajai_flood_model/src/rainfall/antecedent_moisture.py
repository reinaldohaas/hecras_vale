"""
Módulo de Condição de Umidade Antecedente do Solo (Antecedent Moisture Condition - AMC):
- Calcula o índice P5 (Chuva acumulada nos 5 dias anteriores ao evento)
- Classifica automaticamente em AMC I (Seco), AMC II (Normal) ou AMC III (Saturado/Úmido)
- Converte rigorosamente o Curve Number padrão (CN_II) para CN_I ou CN_III
"""

import numpy as np
import pandas as pd
from typing import Tuple, Dict, Any, Optional, Union

class AntecedentMoistureCondition:
    """
    Gerencia a saturação do solo e ajuste dinâmico do Curve Number (CN).
    """
    
    @staticmethod
    def classify_amc_from_p5(p5_mm: float, is_dormant_season: bool = False) -> str:
        """
        Classifica a condição de umidade segundo a tabela padrão do SCS / NRCS NEH-4:
        
        Período de Crescimento (Primavera/Verão):
          - AMC I (Seco):      P5 < 35 mm
          - AMC II (Médio):    35 mm <= P5 <= 53 mm
          - AMC III (Úmido):   P5 > 53 mm
          
        Período de Dormência (Outono/Inverno):
          - AMC I (Seco):      P5 < 13 mm
          - AMC II (Médio):    13 mm <= P5 <= 28 mm
          - AMC III (Úmido):   P5 > 28 mm
        """
        p5 = float(p5_mm)
        if is_dormant_season:
            if p5 < 13.0:
                return 'AMC_I'
            elif p5 <= 28.0:
                return 'AMC_II'
            else:
                return 'AMC_III'
        else:
            if p5 < 35.0:
                return 'AMC_I'
            elif p5 <= 53.0:
                return 'AMC_II'
            else:
                return 'AMC_III'

    @staticmethod
    def adjust_curve_number(cn_ii: float, amc_class: str) -> float:
        """
        Ajusta o Curve Number CN_II para CN_I ou CN_III segundo as fórmulas empíricas do SCS:
        
        CN_I = CN_II / (2.281 - 0.0128 * CN_II)
        CN_III = CN_II / (0.427 + 0.00573 * CN_II)
        """
        cn2 = float(cn_ii)
        if amc_class == 'AMC_I' or amc_class == 'I':
            cn1 = cn2 / (2.281 - 0.0128 * cn2)
            return float(np.clip(cn1, 10.0, 99.0))
        elif amc_class == 'AMC_III' or amc_class == 'III':
            cn3 = cn2 / (0.427 + 0.00573 * cn2)
            return float(np.clip(cn3, 10.0, 99.0))
        else: # AMC_II
            return cn2

    @staticmethod
    def compute_p5_from_series(daily_rainfall_series: Union[pd.Series, np.ndarray, list]) -> float:
        """
        Soma as últimas 5 diárias da série anterior ao início do evento.
        """
        arr = np.asarray(daily_rainfall_series, dtype=float)
        if len(arr) == 0:
            return 20.0 # Padrão neutro AMC II
        p5 = float(np.sum(arr[-5:]))
        return max(0.0, p5)
