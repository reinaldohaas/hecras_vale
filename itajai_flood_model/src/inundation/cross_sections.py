"""
Módulo de Seções Transversais Hidráulicas e Delimitação da Mancha 2D (HEC-RAS Style):
Define seções transversais com batimetria de fundo, cotas de margem plena e declividades de várzea.
Calcula a largura de inundação W_left e W_right para cada cota de água Z_water(x, t) e delimita o polígono da mancha.
"""

import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

class RiverCrossSection:
    """Representa uma seção transversal do vale com geometria de calha e várzeas."""
    def __init__(self, station_id: int, river_name: str,
                 lon: float, lat: float, dist_km: float,
                 z_bed: float, h_bank: float,
                 channel_width_m: float,
                 slope_left_m_per_m: float,
                 slope_right_m_per_m: float,
                 max_floodplain_width_m: float = 3500.0):
        self.station_id = station_id
        self.river_name = river_name
        self.lon = lon
        self.lat = lat
        self.dist_km = dist_km
        self.z_bed = z_bed
        self.h_bank = h_bank
        self.z_bank = z_bed + h_bank
        self.channel_width_m = channel_width_m
        self.slope_left = slope_left_m_per_m
        self.slope_right = slope_right_m_per_m
        self.max_floodplain_width_m = max_floodplain_width_m

    def compute_top_width(self, z_water: float) -> Dict[str, float]:
        """Calcula a largura de inundação à esquerda, direita e total."""
        if z_water <= self.z_bed:
            return {'w_left_m': 0.0, 'w_right_m': 0.0, 'w_total_m': 0.0, 'h_channel_m': 0.0, 'h_overbank_m': 0.0}

        h_channel = z_water - self.z_bed
        
        if z_water <= self.z_bank:
            # Confinado no leito menor
            w_left = self.channel_width_m / 2.0
            w_right = self.channel_width_m / 2.0
            h_overbank = 0.0
        else:
            # Extravasamento na planície de inundação
            h_overbank = z_water - self.z_bank
            w_left = (self.channel_width_m / 2.0) + min(self.max_floodplain_width_m, h_overbank / max(1e-4, self.slope_left))
            w_right = (self.channel_width_m / 2.0) + min(self.max_floodplain_width_m, h_overbank / max(1e-4, self.slope_right))

        return {
            'w_left_m': float(w_left),
            'w_right_m': float(w_right),
            'w_total_m': float(w_left + w_right),
            'h_channel_m': float(h_channel),
            'h_overbank_m': float(h_overbank)
        }


class CrossSectionDelineator:
    """Gera a mancha 2D conectando as larguras de inundação das seções transversais."""

    @staticmethod
    def build_cross_sections_for_river(coords: List[Tuple[float, float]],
                                       z_bed_arr: List[float],
                                       river_key: str) -> List[RiverCrossSection]:
        n_pts = len(coords)
        sections = []
        
        # Parâmetros morfológicos específicos por vale
        specs = {
            'acu': {'h_bank': 7.5, 'b0': 150.0, 's_l': 0.0035, 's_r': 0.0040, 'max_w': 3200.0},
            'mirim': {'h_bank': 5.0, 'b0': 45.0, 's_l': 0.0030, 's_r': 0.0035, 'max_w': 1800.0},
            'oeste': {'h_bank': 5.5, 'b0': 50.0, 's_l': 0.0045, 's_r': 0.0050, 'max_w': 1500.0},
            'sul': {'h_bank': 5.0, 'b0': 45.0, 's_l': 0.0040, 's_r': 0.0045, 'max_w': 1600.0},
            'norte': {'h_bank': 6.0, 'b0': 80.0, 's_l': 0.0150, 's_r': 0.0180, 'max_w': 800.0}, # Vale mais encaixado/garganta
            'benedito': {'h_bank': 5.0, 'b0': 40.0, 's_l': 0.0060, 's_r': 0.0065, 'max_w': 1000.0},
            'luis_alves': {'h_bank': 4.5, 'b0': 30.0, 's_l': 0.0040, 's_r': 0.0045, 'max_w': 1400.0},
            'trombudo': {'h_bank': 4.5, 'b0': 30.0, 's_l': 0.0050, 's_r': 0.0055, 'max_w': 1100.0},
            'mirim_doce': {'h_bank': 4.0, 'b0': 25.0, 's_l': 0.0080, 's_r': 0.0085, 'max_w': 900.0},
            'perimbo': {'h_bank': 4.0, 'b0': 25.0, 's_l': 0.0070, 's_r': 0.0075, 'max_w': 900.0}
        }
        sp = specs.get(river_key, {'h_bank': 5.0, 'b0': 50.0, 's_l': 0.005, 's_r': 0.005, 'max_w': 1500.0})

        for i in range(n_pts):
            lon, lat = coords[i]
            z_b = float(z_bed_arr[i]) if i < len(z_bed_arr) else 0.0
            
            # Ajuste de garganta em trechos montanhosos (ex: Subida/Ibirama no Rio Itajaí-Açu)
            s_l = sp['s_l']
            s_r = sp['s_r']
            if river_key == 'acu' and 20 <= i <= 45: # Garganta entre Rio do Sul e Indaial
                s_l = 0.018 # Encosta íngreme
                s_r = 0.022
            elif river_key == 'acu' and i > 50: # Blumenau, Gaspar, Ilhota, Itajaí (Planície Ampla)
                s_l = 0.0025
                s_r = 0.0028

            sec = RiverCrossSection(
                station_id=i,
                river_name=river_key,
                lon=lon, lat=lat,
                dist_km=float(i),
                z_bed=z_b,
                h_bank=sp['h_bank'],
                channel_width_m=sp['b0'],
                slope_left_m_per_m=s_l,
                slope_right_m_per_m=s_r,
                max_floodplain_width_m=sp['max_w']
            )
            sections.append(sec)
            
        return sections

    @staticmethod
    def delineate_flood_polygon(sections: List[RiverCrossSection],
                                z_water_profile: List[float]) -> Optional[Dict[str, Any]]:
        """
        Delineia o polígono vetorial da mancha conectando os limites de inundação
        esquerdo e direito calculados para cada seção transversal.
        """
        n_sec = len(sections)
        if n_sec < 2 or len(z_water_profile) < n_sec:
            return None

        meters_per_deg_lat = 111320.0
        
        left_pts = []
        right_pts = []
        is_flooding = False

        for i in range(n_sec):
            sec = sections[i]
            z_w = float(z_water_profile[i])
            w_info = sec.compute_top_width(z_w)
            
            if w_info['h_overbank_m'] > 0.05:
                is_flooding = True

            # Vetor normal à seção transversal
            if i < n_sec - 1:
                dx = sections[i+1].lon - sec.lon
                dy = sections[i+1].lat - sec.lat
            else:
                dx = sec.lon - sections[i-1].lon
                dy = sec.lat - sections[i-1].lat

            norm = np.sqrt(dx**2 + dy**2) + 1e-9
            nx = -dy / norm
            ny = dx / norm

            meters_per_deg_lon = meters_per_deg_lat * np.cos(np.radians(sec.lat))
            
            w_left_deg_x = (w_info['w_left_m'] * nx) / meters_per_deg_lon
            w_left_deg_y = (w_info['w_left_m'] * ny) / meters_per_deg_lat
            
            w_right_deg_x = (w_info['w_right_m'] * nx) / meters_per_deg_lon
            w_right_deg_y = (w_info['w_right_m'] * ny) / meters_per_deg_lat

            left_pts.append([sec.lon + w_left_deg_x, sec.lat + w_left_deg_y])
            right_pts.append([sec.lon - w_right_deg_x, sec.lat - w_right_deg_y])

        if not is_flooding:
            return None

        # Fechar polígono de mancha
        poly_coords = left_pts + list(reversed(right_pts)) + [left_pts[0]]
        
        return {
            'type': 'Feature',
            'geometry': {
                'type': 'Polygon',
                'coordinates': [poly_coords]
            },
            'properties': {
                'river': sections[0].river_name,
                'max_stage_m': float(np.max(z_water_profile)),
                'min_stage_m': float(np.min(z_water_profile))
            }
        }
