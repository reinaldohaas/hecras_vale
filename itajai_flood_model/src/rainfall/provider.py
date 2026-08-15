"""
Provedores de Dados Pluviométricos (Rainfall Providers):
- Interface abstrata genérica para integração com múltiplas fontes (ANA, CEMADEN, EPAGRI, INMET, CSV, Sintético)
- Padronização em DataFrame com flags de consistência e controle de qualidade
- Detecção e preenchimento de pequenas lacunas
"""

import os
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Dict, Any, Union

@dataclass
class RainfallObservation:
    """Estrutura padronizada para uma medição pluviométrica."""
    timestamp: pd.Timestamp
    station_id: str
    latitude: float
    longitude: float
    precipitation_mm: float
    source: str = "GENERIC"
    quality_flag: str = "VALID" # VALID | ESTIMATED | SUSPECT | MISSING


class RainfallProvider(ABC):
    """Interface abstrata genérica para provedores de chuva."""
    
    @abstractmethod
    def get_hourly_rainfall(self, start_time: Optional[Union[str, pd.Timestamp]] = None,
                            end_time: Optional[Union[str, pd.Timestamp]] = None,
                            station_ids: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Retorna DataFrame com colunas:
        ['timestamp', 'station_id', 'latitude', 'longitude', 'precipitation_mm', 'source', 'quality_flag']
        """
        pass

    @abstractmethod
    def get_stations(self) -> pd.DataFrame:
        """
        Retorna metadados das estações:
        ['station_id', 'name', 'river_subbasin', 'latitude', 'longitude', 'altitude_m']
        """
        pass


class CSVRainfallProvider(RainfallProvider):
    """
    Provedor baseado em arquivos CSV locais com verificação de consistência e preenchimento de lacunas.
    """
    def __init__(self, rainfall_csv_path: str, stations_csv_path: Optional[str] = None):
        self.rainfall_csv_path = rainfall_csv_path
        self.stations_csv_path = stations_csv_path
        self._df_rainfall: Optional[pd.DataFrame] = None
        self._df_stations: Optional[pd.DataFrame] = None
        self._load_and_validate()

    def _load_and_validate(self):
        if not os.path.exists(self.rainfall_csv_path):
            raise FileNotFoundError(f"Arquivo pluviométrico não encontrado: {self.rainfall_csv_path}")

        df = pd.read_csv(self.rainfall_csv_path)
        
        # Normalização de colunas
        col_map = {
            'DataHora': 'timestamp', 'datahora': 'timestamp', 'time': 'timestamp', 'hora': 'timestamp',
            'Chuva': 'precipitation_mm', 'chuva': 'precipitation_mm', 'precipitacao_mm': 'precipitation_mm', 'p_mm': 'precipitation_mm',
            'Estacao': 'station_id', 'estacao': 'station_id', 'cod_estacao': 'station_id', 'id': 'station_id',
            'Lat': 'latitude', 'lat': 'latitude', 'Lon': 'longitude', 'lon': 'longitude'
        }
        df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})
        
        if 'timestamp' in df.columns:
            df['timestamp'] = pd.to_datetime(df['timestamp'])
        else:
            # Caso seja numérico sequencial (ex: hora 0, 1, 2...)
            if 'hour' in df.columns or 't' in df.columns:
                h_col = 'hour' if 'hour' in df.columns else 't'
                base_dt = pd.Timestamp('2023-10-01 00:00:00')
                df['timestamp'] = [base_dt + pd.Timedelta(hours=int(h)) for h in df[h_col]]
            else:
                df['timestamp'] = pd.date_range(start='2023-10-01 00:00:00', periods=len(df), freq='h')

        if 'quality_flag' not in df.columns:
            df['quality_flag'] = 'VALID'
        if 'source' not in df.columns:
            df['source'] = 'CSV_LOCAL'
            
        # Filtro de controle de qualidade
        # 1. Valores negativos -> 0.0 com flag ESTIMATED
        neg_mask = df['precipitation_mm'] < 0.0
        if neg_mask.any():
            df.loc[neg_mask, 'quality_flag'] = 'ESTIMATED'
            df.loc[neg_mask, 'precipitation_mm'] = 0.0
            
        # 2. Valores absurdos (> 250 mm em 1 hora para a região do Itajaí) -> flag SUSPECT
        suspect_mask = df['precipitation_mm'] > 250.0
        if suspect_mask.any():
            df.loc[suspect_mask, 'quality_flag'] = 'SUSPECT'

        # 3. Tratamento de NaN / NaT
        nan_mask = df['precipitation_mm'].isna()
        if nan_mask.any():
            df.loc[nan_mask, 'quality_flag'] = 'ESTIMATED'
            # Interpolação linear para até 3 horas consecutivas de falha
            df['precipitation_mm'] = df.groupby('station_id')['precipitation_mm'].transform(
                lambda s: s.interpolate(method='linear', limit=3).fillna(0.0)
            )

        self._df_rainfall = df

        # Carregar estações
        if self.stations_csv_path and os.path.exists(self.stations_csv_path):
            df_st = pd.read_csv(self.stations_csv_path)
            self._df_stations = df_st
        else:
            # Gerar resumo a partir do pluviômetro
            st_list = df['station_id'].unique()
            st_data = []
            for st in st_list:
                sub = df[df['station_id'] == st]
                lat = sub['latitude'].iloc[0] if 'latitude' in sub.columns and not sub['latitude'].isna().all() else -27.0
                lon = sub['longitude'].iloc[0] if 'longitude' in sub.columns and not sub['longitude'].isna().all() else -49.5
                st_data.append({
                    'station_id': st,
                    'name': f"Estação Pluviométrica {st}",
                    'river_subbasin': 'Itajai',
                    'latitude': lat,
                    'longitude': lon,
                    'altitude_m': 100.0
                })
            self._df_stations = pd.DataFrame(st_data)

    def get_hourly_rainfall(self, start_time: Optional[Union[str, pd.Timestamp]] = None,
                            end_time: Optional[Union[str, pd.Timestamp]] = None,
                            station_ids: Optional[List[str]] = None) -> pd.DataFrame:
        df = self._df_rainfall.copy()
        if start_time is not None:
            df = df[df['timestamp'] >= pd.to_datetime(start_time)]
        if end_time is not None:
            df = df[df['timestamp'] <= pd.to_datetime(end_time)]
        if station_ids is not None:
            df = df[df['station_id'].isin(station_ids)]
        return df.sort_values(by=['timestamp', 'station_id']).reset_index(drop=True)

    def get_stations(self) -> pd.DataFrame:
        return self._df_stations.copy()


class SyntheticRainfallProvider(RainfallProvider):
    """
    Provedor de Chuva Sintética (NRCS Tipo II ou Huff) para compatibilidade retroativa e cenários de projeto.
    """
    def __init__(self, total_p_mm: float = 120.0, duration_hours: int = 24, distribution: str = 'scs_type_2',
                 start_time: str = '2026-08-15 00:00:00', station_id: str = 'SINTETICA_VALE'):
        self.total_p_mm = float(total_p_mm)
        self.duration_hours = int(duration_hours)
        self.distribution = distribution.lower()
        self.start_time = pd.to_datetime(start_time)
        self.station_id = station_id

    def get_hourly_rainfall(self, start_time: Optional[Union[str, pd.Timestamp]] = None,
                            end_time: Optional[Union[str, pd.Timestamp]] = None,
                            station_ids: Optional[List[str]] = None) -> pd.DataFrame:
        n_hours = self.duration_hours
        p_inc = np.zeros(n_hours + 1)
        
        for h in range(1, n_hours + 1):
            if self.distribution == 'scs_type_2':
                mid = n_hours / 2.0
                frac = (0.5 * (h / mid)**2) if h <= mid else (1.0 - 0.5 * ((n_hours - h) / mid)**2)
            else: # Uniforme
                frac = h / float(n_hours)
            p_inc[h] = self.total_p_mm * frac
            
        p_hourly = np.diff(p_inc) # Comprimento n_hours
        
        timestamps = [self.start_time + pd.Timedelta(hours=i) for i in range(n_hours)]
        
        df = pd.DataFrame({
            'timestamp': timestamps,
            'station_id': self.station_id,
            'latitude': -27.10,
            'longitude': -49.35,
            'precipitation_mm': np.round(p_hourly, 2),
            'source': f'SYNTHETIC_{self.distribution.upper()}',
            'quality_flag': 'VALID'
        })
        
        if start_time is not None:
            df = df[df['timestamp'] >= pd.to_datetime(start_time)]
        if end_time is not None:
            df = df[df['timestamp'] <= pd.to_datetime(end_time)]
            
        return df.reset_index(drop=True)

    def get_stations(self) -> pd.DataFrame:
        return pd.DataFrame([{
            'station_id': self.station_id,
            'name': 'Posto de Chuva de Projeto Sintética',
            'river_subbasin': 'Bacia do Rio Itajaí',
            'latitude': -27.10,
            'longitude': -49.35,
            'altitude_m': 150.0
        }])
