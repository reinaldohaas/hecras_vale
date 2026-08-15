"""
Motor de Cálculo de Cota Espacial e Perfil Longitudinal de Linha d'Água (SpatialStageEngine):
Converte a matriz de vazões espacial Q(x, t) [m³/s] em matrizes contínuas de:
1. Profundidade da água H(x, t) [m]
2. Cota absoluta da linha d'água Z_water(x, t) [m acima do nível do mar]
3. Identificação dos trechos com extravasamento da calha para a várzea (H > H_margem).
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np
from pathlib import Path

from .dem_profile_loader import DEMProfileLoader
from itajai_flood_model.src.rating_curve.manager import RatingCurveManager

class SpatialStageEngine:
    """
    Calcula perfis longitudinais contínuos de linha d'água ao longo de toda a rede hidrográfica.
    """
    def __init__(self, dem_loader: Optional[DEMProfileLoader] = None,
                 rating_manager: Optional[RatingCurveManager] = None):
        self.dem_loader = dem_loader or DEMProfileLoader()
        self.rc_mgr = rating_manager or RatingCurveManager()
        
        # Parâmetros hidráulicos médios de calha principal por rio (Largura b, Talude z, Rugosidade n, Altura da Margem)
        self.river_hydraulic_params = {
            'acu': {'b_start': 55.0, 'b_end': 150.0, 'z': 1.5, 'n': 0.036, 'h_bank': 7.5},
            'oeste': {'b_start': 25.0, 'b_end': 50.0, 'z': 1.2, 'n': 0.038, 'h_bank': 5.5},
            'mirim_doce': {'b_start': 15.0, 'b_end': 28.0, 'z': 1.1, 'n': 0.040, 'h_bank': 4.5},
            'sul': {'b_start': 22.0, 'b_end': 48.0, 'z': 1.2, 'n': 0.038, 'h_bank': 5.0},
            'perimbo': {'b_start': 12.0, 'b_end': 25.0, 'z': 1.1, 'n': 0.040, 'h_bank': 4.0},
            'trombudo': {'b_start': 15.0, 'b_end': 30.0, 'z': 1.2, 'n': 0.039, 'h_bank': 4.5},
            'norte': {'b_start': 35.0, 'b_end': 70.0, 'z': 1.3, 'n': 0.038, 'h_bank': 6.0},
            'benedito': {'b_start': 20.0, 'b_end': 45.0, 'z': 1.2, 'n': 0.038, 'h_bank': 5.0},
            'mirim': {'b_start': 20.0, 'b_end': 60.0, 'z': 1.4, 'n': 0.035, 'h_bank': 5.5},
            'luis_alves': {'b_start': 12.0, 'b_end': 28.0, 'z': 1.2, 'n': 0.038, 'h_bank': 4.5}
        }

    def compute_reach_depth_and_stage(self, river_key: str,
                                      q_spatial_matrix: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calcula as matrizes H(x, t) e Z_water(x, t) para uma calha fluvial.
        
        Parâmetros:
        - river_key: Identificador do rio (ex: 'acu', 'oeste', 'mirim')
        - q_spatial_matrix: Matriz shape (n_sections, n_timesteps) com vazões em m³/s
        
        Retorna dicionário com:
        - 'distances_km': Vetor 1D das distâncias ao longo do rio
        - 'z_bed_m': Vetor 1D da cota de fundo
        - 'depth_h_m': Matriz 2D (n_sections, n_timesteps) de profundidade da lâmina d'água
        - 'z_water_m': Matriz 2D (n_sections, n_timesteps) da cota absoluta da linha d'água
        - 'is_overtopping': Matriz booleana 2D (n_sections, n_timesteps) indicando transbordo
        """
        prof = self.dem_loader.get_river_profile(river_key)
        dists = prof['distances_km']
        z_bed = prof['z_bed_smooth']
        n_sections = len(dists)
        
        q_mat = np.asarray(q_spatial_matrix, dtype=float)
        n_steps, n_times = q_mat.shape
        
        if n_steps != n_sections:
            # Reamostrar q_mat para casar com as seções
            q_resampled = np.zeros((n_sections, n_times))
            orig_idx = np.linspace(0, n_sections - 1, n_steps)
            for t in range(n_times):
                q_resampled[:, t] = np.interp(np.arange(n_sections), orig_idx, q_mat[:, t])
            q_mat = q_resampled

        params = self.river_hydraulic_params.get(river_key, {'b_start': 30.0, 'b_end': 60.0, 'z': 1.3, 'n': 0.038, 'h_bank': 5.0})
        b_vec = np.linspace(params['b_start'], params['b_end'], n_sections)
        z_slope = params['z']
        manning_n = params['n']
        h_bank = params['h_bank']

        depth_mat = np.zeros((n_sections, n_times))
        z_water_mat = np.zeros((n_sections, n_times))
        overtop_mat = np.zeros((n_sections, n_times), dtype=bool)

        # Declividade local do fundo S0 (m/m)
        dz = np.diff(z_bed)
        dx = np.diff(dists) * 1000.0 # em metros
        s0_local = np.abs(dz / np.maximum(10.0, dx))
        s0_local = np.pad(s0_local, (0, 1), mode='edge')
        s0_local = np.maximum(0.00005, s0_local) # Evitar declividade zero

        # Resolver Manning local para cada ponto e tempo
        for i in range(n_sections):
            b_i = b_vec[i]
            sqrt_s0 = np.sqrt(s0_local[i])
            
            # Grade de busca rápida para inversão de Manning
            h_search = np.linspace(0.1, 20.0, 300)
            area_search = b_i * h_search + z_slope * (h_search ** 2)
            perim_search = b_i + 2.0 * h_search * np.sqrt(1.0 + z_slope ** 2)
            r_search = area_search / np.maximum(0.01, perim_search)
            q_search = (1.0 / manning_n) * area_search * (r_search ** (2.0/3.0)) * sqrt_s0

            for t in range(n_times):
                q_val = max(0.1, q_mat[i, t])
                h_val = float(np.interp(q_val, q_search, h_search))
                depth_mat[i, t] = np.round(h_val, 2)
                z_water_mat[i, t] = z_bed[i] + h_val
                overtop_mat[i, t] = (h_val > h_bank)

        # Garantir monotonicidade decrescente estrita da linha d'água no sentido do escoamento
        for t in range(n_times):
            for i in range(1, n_sections):
                if z_water_mat[i, t] >= z_water_mat[i-1, t]:
                    # Forçar declividade suave de jusante
                    z_water_mat[i, t] = z_water_mat[i-1, t] - 0.02
                    depth_mat[i, t] = max(0.2, z_water_mat[i, t] - z_bed[i])

        return {
            'river_key': river_key,
            'name': prof['name'],
            'distances_km': dists,
            'z_bed_m': z_bed,
            'depth_h_m': depth_mat,
            'z_water_m': np.round(z_water_mat, 2),
            'is_overtopping': overtop_mat,
            'h_bank_m': h_bank
        }
