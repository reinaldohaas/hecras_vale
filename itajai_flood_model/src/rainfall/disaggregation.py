"""
Módulo de Desagregação Temporal de Precipitação (Chuva Diária 24h -> Chuva Horária 1h):
Transforma totais pluviométricos diários (ex: de estações convencionais da ANA/EPAGRI/INMET)
em hietogramas horários detalhados P(t) [mm/h].

Métodos suportados:
1. SCS Tipo II (Padrão para eventos convectivos/frontais intensos no Sul do Brasil)
2. Método Triangular / Huff (Quartis de pico antecipado, central ou tardio)
3. Método de Sifalda (Distribuição assimétrica com pico no terço central)
4. Distribuição Uniforme
"""

from enum import Enum
from typing import Union, List, Dict, Optional
import numpy as np
import pandas as pd

class DisaggregationMethod(Enum):
    SCS_TYPE_II = "SCS Tipo II (Pico Central Acentuado)"
    SIFALDA = "Sifalda (Pico no Terço Central)"
    HUFF_Q2 = "Huff 2º Quartil (Pico em 30-50% da duração)"
    TRIANGULAR = "Triangular Assimétrico"
    UNIFORM = "Uniforme / Constante"

class RainfallDisaggregator:
    """
    Desagrega acumulados diários de chuva (24 horas) em taxas horárias (mm/h).
    """

    @staticmethod
    def get_distribution_weights(method: DisaggregationMethod = DisaggregationMethod.SCS_TYPE_II,
                                 n_hours: int = 24) -> np.ndarray:
        """
        Retorna vetor de pesos normalizados (soma = 1.0) para cada uma das n_hours.
        """
        t = np.arange(1, n_hours + 1)
        
        if method == DisaggregationMethod.SCS_TYPE_II:
            # Curva acumulada adimensional SCS Tipo II
            t_frac = t / float(n_hours)
            p_acum = np.zeros(n_hours)
            for i, tf in enumerate(t_frac):
                if tf <= 0.5:
                    p_acum[i] = 2.0 * (tf ** 2)
                else:
                    p_acum[i] = 1.0 - 2.0 * ((1.0 - tf) ** 2)
            # Diferenças incrementais horárias
            p_inc = np.diff(np.insert(p_acum, 0, 0.0))
            weights = p_inc / np.sum(p_inc)

        elif method == DisaggregationMethod.SIFALDA:
            # Pico no centro com cauda suave
            center = n_hours * 0.45
            sigma = n_hours * 0.20
            weights = np.exp(-0.5 * ((t - center) / sigma) ** 2)
            weights = weights / np.sum(weights)

        elif method == DisaggregationMethod.HUFF_Q2:
            # Pico no 2º quartil (t ~ 35-40% da duração)
            center = n_hours * 0.38
            sigma = n_hours * 0.18
            weights = np.exp(-0.5 * ((t - center) / sigma) ** 2)
            weights = weights / np.sum(weights)

        elif method == DisaggregationMethod.TRIANGULAR:
            # Hietograma triangular com pico em t_p = 0.4 * n_hours
            tp = 0.4 * n_hours
            weights = np.where(t <= tp, t / tp, (n_hours - t) / max(0.1, n_hours - tp))
            weights = np.maximum(0.0, weights)
            weights = weights / np.sum(weights)

        else: # UNIFORM
            weights = np.ones(n_hours) / float(n_hours)

        return np.round(weights, 5)

    @classmethod
    def disaggregate_daily_series(cls, daily_totals: Union[List[float], np.ndarray, pd.Series],
                                  method: DisaggregationMethod = DisaggregationMethod.SCS_TYPE_II,
                                  start_datetime: Optional[str] = None) -> pd.DataFrame:
        """
        Recebe uma série de totais diários (ex: 7 dias) e gera a série horária contínua (7 * 24 = 168h).
        """
        daily_arr = np.asarray(daily_totals, dtype=float)
        weights_24h = cls.get_distribution_weights(method=method, n_hours=24)
        
        hourly_precip = []
        for day_val in daily_arr:
            day_precip = float(day_val) * weights_24h
            hourly_precip.extend(day_precip)
            
        hourly_arr = np.round(np.array(hourly_precip), 2)
        n_total_hours = len(hourly_arr)
        
        if start_datetime:
            time_index = pd.date_range(start=start_datetime, periods=n_total_hours, freq='h')
        else:
            time_index = [f"Hora {h:03d}h" for h in range(n_total_hours)]
            
        return pd.DataFrame({
            'timestamp': time_index,
            'precipitation_mm_h': hourly_arr
        })
