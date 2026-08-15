"""
Carregador e Processador de Perfis Altimétricos dos Rios (DEM Copernicus 30m):
Extrai cotas de fundo Z_bed(x), distâncias acumuladas x (km) e declividades locais
para as 10 calhas fluviais da Bacia do Rio Itajaí.
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

class DEMProfileLoader:
    """
    Carrega e processa perfis topográficos e batimétricos dos rios a partir do modelo DEM.
    """
    def __init__(self, dem_json_path: Optional[str] = None):
        if dem_json_path is None:
            # Caminho padrão
            root = Path(__file__).resolve().parent.parent.parent.parent
            dem_json_path = str(root / "app" / "itajai_real_dem_model.json")
            
        self.dem_json_path = dem_json_path
        self.profiles: Dict[str, Dict[str, Any]] = {}
        self._load_profiles()

    def _load_profiles(self):
        p = Path(self.dem_json_path)
        if not p.exists():
            raise FileNotFoundError(f"Arquivo DEM JSON não encontrado: {self.dem_json_path}")
            
        with open(self.dem_json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            
        raw_profiles = data.get('river_profiles', {})
        for r_key, p_data in raw_profiles.items():
            # Extrair distâncias e cotas
            dists = p_data.get('dists_km') or p_data.get('distances_km') or []
            z_raw = p_data.get('z_dem') or p_data.get('elevations') or []
            coords = p_data.get('coords') or []
            
            n_pts = len(coords)
            if not dists and n_pts > 0:
                length_km = float(p_data.get('length_km', 50.0))
                dists = np.linspace(0.0, length_km, n_pts).tolist()
                
            if not z_raw and n_pts > 0:
                # Estimar perfil com declividade média
                h_drop = float(p_data.get('h_drop_m', 100.0))
                z_start = 350.0
                z_raw = np.linspace(z_start, z_start - h_drop, n_pts).tolist()

            dists_arr = np.array(dists, dtype=float)
            z_arr = np.array(z_raw, dtype=float)

            # Suavizar ruído de elevação do DEM (filtro de média móvel + monotonicidade decrescente)
            z_smooth = self._smooth_and_enforce_downstream_slope(z_arr)

            self.profiles[r_key] = {
                'name': p_data.get('name', r_key),
                'n_points': n_pts,
                'length_km': float(dists_arr[-1]) if len(dists_arr) > 0 else 0.0,
                'distances_km': dists_arr,
                'z_bed_raw': z_arr,
                'z_bed_smooth': z_smooth,
                'coords': coords,
                'slope_m_km': float(p_data.get('slope_m_km', 1.0))
            }

    @staticmethod
    def _smooth_and_enforce_downstream_slope(z: np.ndarray, window: int = 5) -> np.ndarray:
        """
        Aplica suavização e garante que o fundo do rio flua estritamente para jusante (sem depressões artificiais).
        """
        if len(z) == 0:
            return z
            
        # 1. Média móvel
        z_pad = np.pad(z, (window//2, window//2), mode='edge')
        z_conv = np.convolve(z_pad, np.ones(window)/window, mode='valid')
        z_conv = z_conv[:len(z)]
        
        # 2. Monotonicidade decrescente montante -> jusante
        z_mono = np.zeros_like(z_conv)
        z_mono[0] = z_conv[0]
        for i in range(1, len(z_conv)):
            z_mono[i] = min(z_mono[i-1] - 0.01, z_conv[i])
            
        return np.round(z_mono, 2)

    def get_river_profile(self, river_key: str) -> Dict[str, Any]:
        """Obtém o perfil processado de uma calha fluvial."""
        if river_key not in self.profiles:
            raise KeyError(f"Rio '{river_key}' não encontrado. Disponíveis: {list(self.profiles.keys())}")
        return self.profiles[river_key]
