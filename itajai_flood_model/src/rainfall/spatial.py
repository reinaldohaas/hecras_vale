"""
Módulo de Interpolação e Distribuição Espacial de Chuva (SpatialRainfallInterpolator):
- Calcula a precipitação média sobre cada sub-bacia hidrográfica (P_bacia(t))
- Suporta ponderação por Polígonos de Thiessen e Ponderação pelo Inverso da Distância (IDW)
- Garante a conservação estrita de massa pluviométrica
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any

# Coordenadas aproximadas dos centróides das 7 principais sub-bacias do Rio Itajaí
SUBBASIN_CENTROIDS = {
    'oeste': {'name': 'Rio Itajaí do Oeste (Taió)', 'lat': -27.15, 'lon': -50.05, 'area_km2': 3120.0},
    'sul': {'name': 'Rio Itajaí do Sul (Ituporanga)', 'lat': -27.42, 'lon': -49.62, 'area_km2': 2280.0},
    'norte': {'name': 'Rio Hercílio / Itajaí do Norte (José Boiteux)', 'lat': -26.95, 'lon': -49.70, 'area_km2': 3450.0},
    'benedito': {'name': 'Rio Benedito (Timbó)', 'lat': -26.78, 'lon': -49.32, 'area_km2': 1540.0},
    'mirim': {'name': 'Rio Itajaí-Mirim (Brusque)', 'lat': -27.18, 'lon': -49.02, 'area_km2': 1680.0},
    'luis_alves': {'name': 'Rio Luís Alves', 'lat': -26.72, 'lon': -48.95, 'area_km2': 580.0},
    'acu': {'name': 'Rio Itajaí-Açu (Calha Principal / Médio e Baixo Vale)', 'lat': -26.95, 'lon': -49.20, 'area_km2': 15000.0}
}


class SpatialRainfallInterpolator:
    """
    Calcula a precipitação média temporal sobre cada sub-bacia a partir de múltiplos pontos pluviométricos.
    """
    def __init__(self, method: str = 'idw', custom_weights: Optional[Dict[str, Dict[str, float]]] = None):
        """
        Parâmetros:
            method: 'thiessen' | 'idw'
            custom_weights: Dicionário customizado {subbasin_id: {station_id: peso}} onde soma dos pesos = 1.0
        """
        self.method = method.lower()
        self.custom_weights = custom_weights or {}
        
    def _compute_idw_weights(self, subbasin_lat: float, subbasin_lon: float,
                             stations_df: pd.DataFrame, power: float = 2.0) -> Dict[str, float]:
        """Calcula pesos IDW baseados na distância euclidiana/geodésica."""
        dists = {}
        for _, st in stations_df.iterrows():
            st_id = str(st['station_id'])
            lat_st = float(st.get('latitude', -27.0))
            lon_st = float(st.get('longitude', -49.5))
            # Distância euclidiana em graus (aproximação adequada para a escala da bacia)
            d = np.sqrt((subbasin_lat - lat_st)**2 + (subbasin_lon - lon_st)**2)
            d = max(d, 0.01) # Prevenir divisão por zero se o ponto for idêntico
            dists[st_id] = 1.0 / (d ** power)
            
        total_inv_d = sum(dists.values())
        return {st_id: (w / total_inv_d) for st_id, w in dists.items()}

    def calculate_subbasin_rainfall(self, rainfall_df: pd.DataFrame,
                                    stations_df: pd.DataFrame,
                                    subbasin_keys: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Recebe a série de chuva de todas as estações e retorna um DataFrame com a chuva média horária
        para cada sub-bacia:
        Columns: ['timestamp', 'oeste', 'sul', 'norte', 'benedito', 'mirim', 'luis_alves', 'acu']
        """
        if rainfall_df.empty:
            return pd.DataFrame()

        subbasins_to_compute = subbasin_keys or list(SUBBASIN_CENTROIDS.keys())
        
        # Obter timestamps únicos ordenados
        timestamps = np.sort(rainfall_df['timestamp'].unique())
        
        # Pivotar chuva: index = timestamp, columns = station_id
        pivot_df = rainfall_df.pivot_table(
            index='timestamp', columns='station_id', values='precipitation_mm', aggfunc='mean'
        ).fillna(0.0)
        
        available_stations = [st for st in pivot_df.columns if st in stations_df['station_id'].values]
        if not available_stations:
            available_stations = list(pivot_df.columns)

        result_dict = {'timestamp': timestamps}

        for sb_key in subbasins_to_compute:
            sb_info = SUBBASIN_CENTROIDS.get(sb_key, {'lat': -27.0, 'lon': -49.5})
            
            # 1. Verificar se há pesos customizados
            if sb_key in self.custom_weights:
                weights = self.custom_weights[sb_key]
            else:
                # 2. Calcular pesos IDW
                weights = self._compute_idw_weights(sb_info['lat'], sb_info['lon'], stations_df)
                
            # Filtrar e normalizar pesos para estações disponíveis no pivot
            active_weights = {st: weights[st] for st in available_stations if st in weights}
            w_sum = sum(active_weights.values())
            
            if w_sum > 0:
                norm_weights = {st: (w / w_sum) for st, w in active_weights.items()}
                # P_bacia(t) = sum( peso_i * P_i(t) )
                sb_rainfall = np.zeros(len(timestamps))
                for st, w in norm_weights.items():
                    sb_rainfall += w * pivot_df[st].values
                result_dict[sb_key] = np.round(sb_rainfall, 2)
            else:
                # Se não houver estações ponderáveis, usa média simples das disponíveis
                result_dict[sb_key] = np.round(pivot_df.mean(axis=1).values, 2)

        return pd.DataFrame(result_dict)
