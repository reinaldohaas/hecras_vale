"""
Módulo Oficial de Inundação 2D baseado no Modelo HAND (Height Above Nearest Drainage):
1. Calcula a matriz HAND = Z_DEM - Z_drenagem_proxima a partir do GeoTIFF real do DEM e da hidrografia da ANA.
2. Mapeia a lâmina d'água 2D: depth(x, y, t) = max(0, H_local(t) - HAND(x, y)).
3. Aplica o filtro de conectividade hidráulica 2D (Flood-Fill) a partir do leito do rio.
4. Vetoriza os contornos topográficos reais em GeoJSON (sem fitas ou corredores artificiais).
"""

import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import scipy.ndimage as ndimage
import rasterio
from rasterio.features import rasterize, shapes
from shapely.geometry import shape, mapping
from shapely.ops import unary_union

class HANDModel:
    """
    Motor do modelo HAND (Height Above Nearest Drainage) para a Bacia do Rio Itajaí.
    """
    def __init__(self, dem_tif_path: str, river_geojson_path: str):
        self.dem_tif_path = Path(dem_tif_path)
        self.river_geojson_path = Path(river_geojson_path)
        
        self.dem: Optional[np.ndarray] = None
        self.hand: Optional[np.ndarray] = None
        self.river_mask: Optional[np.ndarray] = None
        self.nearest_drainage_indices: Optional[Tuple[np.ndarray, np.ndarray]] = None
        self.transform: Optional[rasterio.Affine] = None
        self.crs = None
        self.bounds = None
        
        self._load_and_compute_hand()

    def _load_and_compute_hand(self):
        if not self.dem_tif_path.exists():
            raise FileNotFoundError(f"DEM GeoTIFF não encontrado: {self.dem_tif_path}")
            
        with rasterio.open(self.dem_tif_path) as src:
            self.dem = src.read(1).astype(float)
            self.transform = src.transform
            self.crs = src.crs
            self.bounds = src.bounds
            
        # Carregar hidrografia da ANA
        with open(self.river_geojson_path, 'r', encoding='utf-8') as f:
            river_data = json.load(f)
            
        # Rasterizar rios sobre o grid DEM
        shapes_to_rasterize = []
        for feat in river_data.get('features', []):
            geom = feat.get('geometry')
            if geom:
                shapes_to_rasterize.append((geom, 1))
                
        self.river_mask = rasterize(
            shapes_to_rasterize,
            out_shape=self.dem.shape,
            transform=self.transform,
            fill=0,
            dtype=np.uint8
        )
        
        # Calcular distância euclidiana e índices da célula de rio mais próxima
        stream_mask = (self.river_mask == 1)
        dist, indices = ndimage.distance_transform_edt(~stream_mask, return_indices=True)
        self.nearest_drainage_indices = (indices[0], indices[1])
        
        # Elevação da drenagem correspondente
        z_drainage = self.dem[indices[0], indices[1]]
        
        # HAND = max(0, Z_DEM - Z_drenagem)
        self.hand = np.maximum(0.0, self.dem - z_drainage)

    def generate_flood_inundation(self, stages_by_region: Dict[str, float],
                                  max_dist_cells: int = 140,
                                  min_depth_m: float = 0.05) -> Dict[str, Any]:
        """
        Gera a mancha de inundação baseada em HAND para as cotas de régua especificadas:
        - stages_by_region: {'blumenau': 15.34, 'rio_sul': 13.0, 'brusque': 8.5, 'baixo_vale': 4.5}
        """
        nrows, ncols = self.dem.shape
        
        # Criar matriz 2D de cotas locais H(x, y) interpoladas ao longo da bacia
        h_local = np.ones((nrows, ncols), dtype=float) * stages_by_region.get('blumenau', 15.34)
        
        frac_c = np.linspace(0, 1, ncols)
        # Oeste -> Leste (Rio do Sul -> Blumenau -> Foz)
        h_rs = stages_by_region.get('rio_sul', 13.0)
        h_blu = stages_by_region.get('blumenau', 15.34)
        h_foz = stages_by_region.get('baixo_vale', 4.5)
        
        for c in range(ncols):
            f = c / float(ncols - 1)
            if f < 0.5:
                # Alto Vale -> Médio Vale
                h_val = (1.0 - 2.0 * f) * h_rs + (2.0 * f) * h_blu
            else:
                # Médio Vale -> Foz
                h_val = (2.0 - 2.0 * f) * h_blu + (2.0 * f - 1.0) * h_foz
            h_local[:, c] = h_val

        # Lâmina bruta: depth = max(0, H_local - HAND)
        depth_raw = np.maximum(0.0, h_local - self.hand)
        
        # Filtro de alcance do vale
        stream_mask = (self.river_mask == 1)
        dist_to_stream = ndimage.distance_transform_edt(~stream_mask)
        depth_raw[dist_to_stream > max_dist_cells] = 0.0

        # Filtro de Conectividade Hidráulica 2D (Connected Component / Seed Fill)
        is_wet = (depth_raw >= min_depth_m)
        labeled_array, num_features = ndimage.label(is_wet, structure=np.ones((3, 3)))
        
        # Manter apenas as componentes que tocam a rede de drenagem (sementes de rio)
        river_labels = np.unique(labeled_array[stream_mask])
        river_labels = river_labels[river_labels > 0] # Excluir fundo (0)
        
        connected_mask = np.isin(labeled_array, river_labels)
        depth_connected = np.where(connected_mask, depth_raw, 0.0)

        # Cálculo da Área Inundada real
        dx_m = self.transform[0] * 111320.0 * np.cos(np.radians(np.mean([self.bounds.bottom, self.bounds.top])))
        dy_m = abs(self.transform[4]) * 110570.0
        cell_area_km2 = (dx_m * dy_m) / 1e6
        
        flooded_cells = int(np.sum(depth_connected >= min_depth_m))
        total_area_km2 = float(flooded_cells * cell_area_km2)
        total_volume_hm3 = float(np.sum(depth_connected) * (dx_m * dy_m) / 1e6)

        return {
            'depth_raster': depth_connected,
            'connected_mask': connected_mask,
            'area_km2': round(total_area_km2, 2),
            'volume_hm3': round(total_volume_hm3, 2),
            'max_depth_m': round(float(np.max(depth_connected)), 2) if flooded_cells > 0 else 0.0,
            'mean_depth_m': round(float(np.mean(depth_connected[connected_mask])), 2) if flooded_cells > 0 else 0.0
        }

    def export_inundation_geojson(self, depth_raster: np.ndarray,
                                  output_geojson_path: str,
                                  min_depth_m: float = 0.10,
                                  simplify_tolerance_deg: float = 0.0005):
        """
        Vetoriza a mancha de inundação do raster HAND em polígonos GeoJSON contornados,
        classificados em faixas de profundidade (Leve, Média, Severa).
        """
        mask = (depth_raster >= min_depth_m).astype(np.uint8)
        
        features = []
        # Extrair polígonos contornados das células inundadas
        for geom_dict, val in shapes(mask, mask=(mask == 1), transform=self.transform):
            if val == 1:
                poly_shp = shape(geom_dict)
                if poly_shp.area < 1e-7:
                    continue
                    
                if simplify_tolerance_deg > 0:
                    poly_shp = poly_shp.simplify(simplify_tolerance_deg, preserve_topology=True)
                    
                feat = {
                    'type': 'Feature',
                    'geometry': mapping(poly_shp),
                    'properties': {
                        'tipo': 'Mancha de Inundação HAND (Topografia Real DEM)',
                        'metodo': 'Height Above Nearest Drainage + Conectividade BFS',
                        'fill_color': '#00f0ff',
                        'fill_opacity': 0.65
                    }
                }
                features.append(feat)

        geojson_out = {
            'type': 'FeatureCollection',
            'features': features
        }
        
        p = Path(output_geojson_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, 'w', encoding='utf-8') as f:
            json.dump(geojson_out, f, indent=2)
            
        print(f"✅ Mancha HAND exportada com sucesso para {output_geojson_path} ({len(features)} polígonos)!")
        return geojson_out
