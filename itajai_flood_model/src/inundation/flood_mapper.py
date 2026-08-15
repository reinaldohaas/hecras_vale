"""
Gerador e Exportador de Manchas de Inundação GeoJSON (FloodplainMapper):
Gera geometrias vetoriais georreferenciadas (GeoJSON FeatureCollection)
com classificação de lâmina d'água e métricas espaciais para visualização em Leaflet e GIS.
"""

import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

class FloodplainMapper:
    """
    Constrói polígonos de mancha de inundação 2D a partir dos perfis de linha d'água
    e do modelo digital de elevação.
    """

    @staticmethod
    def generate_flood_corridor_polygons(river_coords: List[Tuple[float, float]],
                                         reach_depths: np.ndarray,
                                         h_bank: float = 7.5,
                                         max_corridor_width_deg: float = 0.025) -> Dict[str, Any]:
        """
        Gera polígonos de faixa de inundação ao longo de uma calha fluvial com base na profundidade
        e no extravasamento das margens.
        """
        features = []
        n_pts = len(river_coords)
        if n_pts < 2 or len(reach_depths) < n_pts:
            return {'type': 'FeatureCollection', 'features': []}

        for i in range(n_pts - 1):
            lon1, lat1 = river_coords[i]
            lon2, lat2 = river_coords[i+1]
            
            d_val = float(0.5 * (reach_depths[i] + reach_depths[i+1]))
            
            if d_val > h_bank:
                # Transbordo ativo: largura da mancha proporcional ao excesso de cota (H - H_bank)
                excess_h = d_val - h_bank
                # Largura angular (aproximadamente 1km a cada 3m de sobre-elevação)
                w = min(max_corridor_width_deg, 0.003 + 0.004 * (excess_h ** 0.8))
                
                # Vetor normal ao segmento
                dx = lon2 - lon1
                dy = lat2 - lat1
                norm = np.sqrt(dx**2 + dy**2) + 1e-9
                nx = -dy / norm
                ny = dx / norm
                
                # Vértices do trapézio de inundação
                p1 = [lon1 + nx * w, lat1 + ny * w]
                p2 = [lon2 + nx * w, lat2 + ny * w]
                p3 = [lon2 - nx * w, lat2 - ny * w]
                p4 = [lon1 - nx * w, lat1 - ny * w]
                poly_coords = [[p1, p2, p3, p4, p1]]

                # Classificação de lâmina
                depth_class = "Severa" if excess_h > 2.0 else ("Moderada" if excess_h > 0.8 else "Leve")
                color_hex = "#ff0055" if excess_h > 2.0 else ("#f59e0b" if excess_h > 0.8 else "#00f0ff")

                feat = {
                    'type': 'Feature',
                    'geometry': {
                        'type': 'Polygon',
                        'coordinates': poly_coords
                    },
                    'properties': {
                        'segment_index': i,
                        'depth_m': round(d_val, 2),
                        'excess_stage_m': round(excess_h, 2),
                        'depth_class': depth_class,
                        'fill_color': color_hex,
                        'fill_opacity': 0.65
                    }
                }
                features.append(feat)

        return {
            'type': 'FeatureCollection',
            'features': features
        }

    @classmethod
    def export_flood_geojson(cls, geojson_data: Dict[str, Any], output_path: str):
        """Salva a mancha em arquivo GeoJSON."""
        p = Path(output_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(geojson_data, f, indent=2, ensure_ascii=False)
