"""
Baixador Automático de Chuva Histórica Real e Desagregador (RainfallDownloader):
Acessa bases de dados históricas e reanálise (Open-Meteo / ERA5 / ANA) para obter
precipitação horária real em todas as 10 sub-bacias do Vale do Itajaí.
"""

import os
import json
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd

from .disaggregation import RainfallDisaggregator, DisaggregationMethod

class RainfallDownloader:
    """
    Realiza o download de séries históricas de precipitação horária e diária
    para as sub-bacias da Bacia Hidrográfica do Rio Itajaí.
    """

    # Coordenadas geográficas dos centróides pluviométricos de cada sub-bacia
    SUB_BASIN_CENTROIDS = {
        'oeste': {'name': 'Alto Vale - Rio do Oeste (Taió)', 'lat': -27.115, 'lon': -49.998},
        'mirim_doce': {'name': 'Alto Vale - Rio Mirim Doce', 'lat': -27.195, 'lon': -50.075},
        'sul': {'name': 'Alto Vale - Rio do Sul (Ituporanga)', 'lat': -27.414, 'lon': -49.605},
        'perimbo': {'name': 'Alto Vale - Rio Perimbó (Petrolândia)', 'lat': -27.535, 'lon': -49.705},
        'trombudo': {'name': 'Alto Vale - Rio Trombudo (Agrolândia)', 'lat': -27.300, 'lon': -49.792},
        'norte': {'name': 'Médio Vale - Rio Hercílio (Boiteux / Ibirama)', 'lat': -26.960, 'lon': -49.628},
        'benedito': {'name': 'Médio Vale - Rio Benedito (Timbó / Pomerode)', 'lat': -26.820, 'lon': -49.270},
        'mirim': {'name': 'Médio/Baixo - Rio Itajaí-Mirim (Brusque / Botuverá)', 'lat': -27.098, 'lon': -48.912},
        'luis_alves': {'name': 'Baixo Vale - Rio Luís Alves', 'lat': -26.720, 'lon': -48.930},
        'acu': {'name': 'Médio Vale - Tronco Principal (Blumenau / Gaspar)', 'lat': -26.918, 'lon': -49.066}
    }

    @classmethod
    def download_hourly_event(cls, start_date: str, end_date: str,
                              timezone: str = "America/Sao_Paulo",
                              save_path: Optional[str] = None) -> pd.DataFrame:
        """
        Baixa chuva horária real para todas as 10 sub-bacias no intervalo de datas (YYYY-MM-DD).
        Retorna DataFrame com coluna 'timestamp' e colunas para cada sub-bacia (mm/h).
        """
        print(f"📡 Baixando precipitação horária real para 10 sub-bacias ({start_date} a {end_date})...")
        
        results = {}
        time_index = None

        for sb_key, info in cls.SUB_BASIN_CENTROIDS.items():
            lat = info['lat']
            lon = info['lon']
            
            url = (
                f"https://archive-api.open-meteo.com/v1/archive?"
                f"latitude={lat}&longitude={lon}&start_date={start_date}&end_date={end_date}"
                f"&hourly=precipitation&timezone={urllib.parse.quote(timezone)}"
            )
            
            req = urllib.request.Request(url, headers={'User-Agent': 'ItajaiFloodModel/2.0'})
            try:
                with urllib.request.urlopen(req, timeout=20) as resp:
                    data = json.loads(resp.read().decode('utf-8'))
                    
                hourly = data.get('hourly', {})
                if time_index is None:
                    time_index = hourly.get('time', [])
                    
                p_arr = hourly.get('precipitation', [])
                results[sb_key] = [max(0.0, float(v)) if v is not None else 0.0 for v in p_arr]
                total_mm = sum(results[sb_key])
                max_h = max(results[sb_key]) if results[sb_key] else 0.0
                print(f"   ✓ {info['name']:55s}: Total = {total_mm:6.1f} mm | Pico = {max_h:4.1f} mm/h")
                
            except Exception as e:
                print(f"   ❌ Erro ao baixar {sb_key}: {e}")
                if time_index:
                    results[sb_key] = [0.0] * len(time_index)

        df_out = pd.DataFrame({'timestamp': time_index})
        for k in cls.SUB_BASIN_CENTROIDS.keys():
            if k in results:
                df_out[k] = results[k]
            else:
                df_out[k] = 0.0
                
        if save_path:
            p = Path(save_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            df_out.to_csv(save_path, index=False)
            print(f"💾 Dados salvos com sucesso em: {save_path}")

        return df_out

    @classmethod
    def download_and_disaggregate_daily_event(cls, daily_data_dict: Dict[str, List[float]],
                                              start_date: str,
                                              method: DisaggregationMethod = DisaggregationMethod.SCS_TYPE_II) -> pd.DataFrame:
        """
        Recebe totais diários por sub-bacia e aplica a desagregação horária.
        """
        out_dict = {}
        for sb_key, daily_vals in daily_data_dict.items():
            df_sb = RainfallDisaggregator.disaggregate_daily_series(
                daily_totals=daily_vals,
                method=method,
                start_datetime=f"{start_date} 00:00:00"
            )
            if 'timestamp' not in out_dict:
                out_dict['timestamp'] = df_sb['timestamp'].tolist()
            out_dict[sb_key] = df_sb['precipitation_mm_h'].tolist()
            
        return pd.DataFrame(out_dict)
