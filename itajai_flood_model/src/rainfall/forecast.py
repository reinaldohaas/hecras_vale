"""
Módulo de Previsão Pluviométrica Temporal Contínua (RainfallForecastTimeline):
- Junta o passado observado (t <= t_agora) com o futuro previsto (t > t_agora)
- Gera 3 cenários de previsão:
    * P_low:   Cenário Seco / Inferior (ex: 0.6x da previsão central)
    * P_mean:  Cenário Determinístico Central (Modelo Meteorológico)
    * P_high:  Cenário Chuvoso / Superior (ex: 1.4x da previsão central)
- Suporta horizontes configuráveis: 6h, 12h, 24h, 48h, 72h
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any, Union

class RainfallForecastTimeline:
    """
    Linha do tempo contínua combinando observação e previsão com quantificação de incerteza.
    """
    def __init__(self, observed_df: pd.DataFrame,
                 forecast_df: Optional[pd.DataFrame] = None,
                 now_timestamp: Optional[Union[str, pd.Timestamp]] = None,
                 horizon_hours: int = 24,
                 uncertainty_factor_low: float = 0.60,
                 uncertainty_factor_high: float = 1.40):
        """
        Parâmetros:
            observed_df: DataFrame com colunas ['timestamp', sub-bacias...]
            forecast_df: DataFrame com previsão para o futuro (opcional, pode ser gerado sinteticamente)
            now_timestamp: Instante divisor 'AGORA' (se None, assume o último timestamp observado)
            horizon_hours: Horizonte futuro de previsão em horas (6, 12, 24, 48, 72)
            uncertainty_factor_low: Multiplicador do cenário inferior
            uncertainty_factor_high: Multiplicador do cenário superior
        """
        self.observed_df = observed_df.sort_values(by='timestamp').reset_index(drop=True)
        self.horizon_hours = int(horizon_hours)
        self.factor_low = float(uncertainty_factor_low)
        self.factor_high = float(uncertainty_factor_high)
        
        if now_timestamp is not None:
            self.now_timestamp = pd.to_datetime(now_timestamp)
        else:
            self.now_timestamp = self.observed_df['timestamp'].iloc[-1]
            
        self.forecast_df = forecast_df

    def build_continuous_scenarios(self) -> Dict[str, pd.DataFrame]:
        """
        Constrói os 3 DataFrames contínuos (P_low, P_mean, P_high) para cada sub-bacia.
        
        Retorna:
            {
                'mean': df_continuous_mean,
                'low': df_continuous_low,
                'high': df_continuous_high,
                'metadata': {
                    'now_timestamp': self.now_timestamp,
                    'horizon_hours': self.horizon_hours,
                    'total_hours': N
                }
            }
        """
        # 1. Separar histórico observado até 'now_timestamp'
        hist_df = self.observed_df[self.observed_df['timestamp'] <= self.now_timestamp].copy()
        subbasin_cols = [c for c in self.observed_df.columns if c != 'timestamp']
        
        # 2. Obter ou gerar previsão futura
        fut_timestamps = [self.now_timestamp + pd.Timedelta(hours=i+1) for i in range(self.horizon_hours)]
        
        if self.forecast_df is not None and not self.forecast_df.empty:
            # Usar previsão fornecida
            fut_raw = self.forecast_df[self.forecast_df['timestamp'] > self.now_timestamp].copy()
            fut_raw = fut_raw[fut_raw['timestamp'] <= fut_timestamps[-1]]
            
            # Reindexar para garantir cobertura do horizonte
            fut_mean_df = pd.DataFrame({'timestamp': fut_timestamps})
            for col in subbasin_cols:
                if col in fut_raw.columns:
                    fut_mean_df = pd.merge(fut_mean_df, fut_raw[['timestamp', col]], on='timestamp', how='left')
                else:
                    fut_mean_df[col] = 0.0
            fut_mean_df = fut_mean_df.fillna(0.0)
        else:
            # Decaimento exponencial da última chuva observada se não houver arquivo externo
            fut_mean_df = pd.DataFrame({'timestamp': fut_timestamps})
            for col in subbasin_cols:
                last_val = hist_df[col].iloc[-1] if not hist_df.empty else 5.0
                # Curva de decaimento suave
                decay = np.array([last_val * np.exp(-0.15 * i) for i in range(1, self.horizon_hours + 1)])
                fut_mean_df[col] = np.round(decay, 2)

        # 3. Criar Cenários Futuros (Low, Mean, High)
        fut_low_df = fut_mean_df.copy()
        fut_high_df = fut_mean_df.copy()
        
        for col in subbasin_cols:
            fut_low_df[col] = np.round(fut_mean_df[col] * self.factor_low, 2)
            fut_high_df[col] = np.round(fut_mean_df[col] * self.factor_high, 2)

        # 4. Concatenar Passado + Futuro para os 3 cenários
        df_mean = pd.concat([hist_df, fut_mean_df], ignore_index=True).sort_values('timestamp').reset_index(drop=True)
        df_low = pd.concat([hist_df, fut_low_df], ignore_index=True).sort_values('timestamp').reset_index(drop=True)
        df_high = pd.concat([hist_df, fut_high_df], ignore_index=True).sort_values('timestamp').reset_index(drop=True)
        
        # Adicionar coluna indicativa 'phase' (OBSERVED | FORECASTED)
        for df in [df_mean, df_low, df_high]:
            df['phase'] = np.where(df['timestamp'] <= self.now_timestamp, 'OBSERVED', 'FORECASTED')

        return {
            'mean': df_mean,
            'low': df_low,
            'high': df_high,
            'now_timestamp': self.now_timestamp,
            'horizon_hours': self.horizon_hours,
            'subbasin_cols': subbasin_cols
        }
