"""
Modelo HAND Espaço-Temporal Dinâmico (DynamicSpatiotemporalHAND):
Calcula a inundação 2D levando em conta a cota exata do rio H(x, t) EM CADA LUGAR e EM CADA HORA:
1. Mapeia cada pixel (r, c) da bacia para seu rio específico e quilometragem/estação de drenagem x_reach.
2. Extrai H_drain(r, c, t) a partir da matriz hidrodinâmica espaço-temporal H(x, t) no instante t.
3. Calcula a profundidade instantânea: depth(r, c, t) = max(0, H_drain(r, c, t) - HAND(r, c)).
4. Aplica conectividade hidráulica com o canal.
5. Exporta as camadas horárias da mancha de inundação para reprodução dinâmica.
"""

import sys
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import scipy.ndimage as ndimage
import rasterio
from rasterio.features import rasterize, shapes
from shapely.geometry import shape, mapping

class DynamicSpatiotemporalHAND:
    """
    Motor HAND Espaço-Temporal: calcula a mancha de inundação 2D variando hora a hora e ponto a ponto.
    """
    def __init__(self, dem_tif_path: str, dem_model_json_path: str):
        self.dem_tif_path = Path(dem_tif_path)
        self.dem_model_json_path = Path(dem_model_json_path)
        
        self.dem: Optional[np.ndarray] = None
        self.hand: Optional[np.ndarray] = None
        self.river_id_grid: Optional[np.ndarray] = None
        self.station_frac_grid: Optional[np.ndarray] = None
        self.stream_mask: Optional[np.ndarray] = None
        self.transform: Optional[rasterio.Affine] = None
        self.crs = None
        self.bounds = None
        
        self.river_keys = ['acu', 'oeste', 'sul', 'norte', 'benedito', 'mirim', 'luis_alves', 'trombudo', 'mirim_doce', 'perimbo']
        self.river_key_to_id = {k: i + 1 for i, k in enumerate(self.river_keys)}
        self.id_to_river_key = {i + 1: k for i, k in enumerate(self.river_keys)}
        
        self._initialize_hand_and_river_mapping()

    def _initialize_hand_and_river_mapping(self):
        if not self.dem_tif_path.exists():
            raise FileNotFoundError(f"GeoTIFF não encontrado: {self.dem_tif_path}")
            
        with rasterio.open(self.dem_tif_path) as src:
            self.dem = src.read(1).astype(float)
            self.transform = src.transform
            self.crs = src.crs
            self.bounds = src.bounds
            
        with open(self.dem_model_json_path, 'r', encoding='utf-8') as f:
            dem_model = json.load(f)
            
        profiles = dem_model.get('river_profiles', {})
        nrows, ncols = self.dem.shape
        
        # Grade com ID do rio e fração de distância ao longo do rio
        river_channel_id = np.zeros((nrows, ncols), dtype=np.int32)
        river_channel_frac = np.zeros((nrows, ncols), dtype=np.float32)
        channel_seeds = np.zeros((nrows, ncols), dtype=bool)

        inv_transform = ~self.transform

        # Densificar coordenadas de cada rio e mapear para a grade
        for r_key, prof in profiles.items():
            if r_key not in self.river_key_to_id:
                continue
            r_id = self.river_key_to_id[r_key]
            coords = prof.get('coords', [])
            if len(coords) < 2:
                continue
                
            # Interpolar pontos a cada ~50m para garantir continuidade na grade raster
            dense_coords = []
            dense_fracs = []
            n_orig = len(coords)
            for i in range(n_orig - 1):
                p1 = coords[i]
                p2 = coords[i+1]
                f1 = i / float(n_orig - 1)
                f2 = (i + 1) / float(n_orig - 1)
                sub_steps = 6
                for s in range(sub_steps):
                    w = s / float(sub_steps)
                    lon = (1.0 - w) * p1[0] + w * p2[0]
                    lat = (1.0 - w) * p1[1] + w * p2[1]
                    f_val = (1.0 - w) * f1 + w * f2
                    dense_coords.append((lon, lat))
                    dense_fracs.append(f_val)
                    
            for (lon, lat), f_val in zip(dense_coords, dense_fracs):
                c, r = inv_transform * (lon, lat)
                c, r = int(round(c)), int(round(r))
                if 0 <= r < nrows and 0 <= c < ncols:
                    river_channel_id[r, c] = r_id
                    river_channel_frac[r, c] = f_val
                    channel_seeds[r, c] = True

        self.stream_mask = channel_seeds
        
        # Calcular distância euclidiana até a célula de drenagem mais próxima
        dist, indices = ndimage.distance_transform_edt(~self.stream_mask, return_indices=True)
        nearest_r = indices[0]
        nearest_c = indices[1]
        
        # Atribuir HAND, ID do rio e fração de estação para TODAS as células da bacia
        z_drainage = self.dem[nearest_r, nearest_c]
        self.hand = np.maximum(0.0, self.dem - z_drainage)
        
        self.river_id_grid = river_channel_id[nearest_r, nearest_c]
        self.station_frac_grid = river_channel_frac[nearest_r, nearest_c]

    def compute_hourly_inundation(self, river_spatial_h_matrices: Dict[str, np.ndarray],
                                  t_hour: int,
                                  max_dist_cells: int = 150,
                                  min_depth_m: float = 0.15) -> Dict[str, Any]:
        """
        Calcula a lâmina de inundação 2D para a hora t específica:
        - river_spatial_h_matrices: {r_key: matriz shape (n_sections, 49) com cotas locais H(x, t) em metros}
        """
        nrows, ncols = self.dem.shape
        h_drain_grid = np.zeros((nrows, ncols), dtype=float)

        # Atribuir o nível d'água H(x, t) daquela hora exata e daquele local exato
        for r_key, r_id in self.river_key_to_id.items():
            mask_r = (self.river_id_grid == r_id)
            if not np.any(mask_r):
                continue
                
            h_matrix = river_spatial_h_matrices.get(r_key)
            if h_matrix is not None:
                n_sec = h_matrix.shape[0]
                t_idx = min(t_hour, h_matrix.shape[1] - 1)
                # Perfil de cota do rio naquele instante t
                h_profile_t = h_matrix[:, t_idx]
                
                # Mapear a fração do trecho [0, 1] para a cota correspondente
                fracs_in_reach = self.station_frac_grid[mask_r]
                sec_indices = np.clip((fracs_in_reach * (n_sec - 1)).astype(int), 0, n_sec - 1)
                h_drain_grid[mask_r] = h_profile_t[sec_indices]
            else:
                # Nível padrão de calha se matriz não informada
                h_drain_grid[mask_r] = 1.50

        # Lâmina d'água instantânea: depth(r, c, t) = max(0, H_drain(r, c, t) - HAND(r, c))
        depth_raw = np.maximum(0.0, h_drain_grid - self.hand)
        
        # Filtro de alcance lateral do vale
        dist_to_stream = ndimage.distance_transform_edt(~self.stream_mask)
        depth_raw[dist_to_stream > max_dist_cells] = 0.0

        # Filtro de Conectividade Hidráulica 2D com a calha (LUT ultrarrápido)
        is_wet = (depth_raw >= min_depth_m)
        labeled, num_features = ndimage.label(is_wet, structure=np.ones((3, 3)))
        
        stream_labels = np.unique(labeled[self.stream_mask])
        stream_labels = stream_labels[stream_labels > 0]
        
        lut = np.zeros(num_features + 1, dtype=bool)
        lut[stream_labels] = True
        connected_mask = lut[labeled]
        depth_connected = np.where(connected_mask, depth_raw, 0.0)

        # Cálculo de métricas
        dx_m = self.transform[0] * 111320.0 * np.cos(np.radians(np.mean([self.bounds.bottom, self.bounds.top])))
        dy_m = abs(self.transform[4]) * 110570.0
        cell_area_km2 = (dx_m * dy_m) / 1e6
        
        flooded_cells = int(np.sum(depth_connected >= min_depth_m))
        area_km2 = float(flooded_cells * cell_area_km2)
        vol_hm3 = float(np.sum(depth_connected) * (dx_m * dy_m) / 1e6)

        return {
            'time_hour': t_hour,
            'depth_raster': depth_connected,
            'area_km2': round(area_km2, 2),
            'volume_hm3': round(vol_hm3, 2),
            'max_depth_m': round(float(np.max(depth_connected)), 2) if flooded_cells > 0 else 0.0
        }
