"""
Testes Unitários para o Módulo de Cota Espacial e Perfil Longitudinal (ETAPA 2):
Verifica:
1. Carregamento correto dos perfis do DEM e garantia de monotonicidade de fundo
2. Cálculo das matrizes H(x, t) e Z_water(x, t)
3. Monotonicidade decrescente estrita da linha d'água de montante para jusante
4. Detecção precisa de extravasamento da calha (overtopping)
"""

import unittest
import numpy as np
from pathlib import Path

from itajai_flood_model.src.spatial_stage import DEMProfileLoader, SpatialStageEngine

class TestSpatialStageProfile(unittest.TestCase):

    def setUp(self):
        self.loader = DEMProfileLoader()
        self.engine = SpatialStageEngine(dem_loader=self.loader)

    def test_dem_profile_loader(self):
        """Verifica se os 10 rios foram carregados com cotas e distâncias consistentes."""
        rivers = ['acu', 'oeste', 'mirim_doce', 'sul', 'perimbo', 'trombudo', 'norte', 'benedito', 'mirim', 'luis_alves']
        for r_key in rivers:
            prof = self.loader.get_river_profile(r_key)
            self.assertGreater(prof['n_points'], 0, f"Perfil {r_key} vazio")
            self.assertGreater(prof['length_km'], 0.0, f"Comprimento do rio {r_key} inválido")
            
            # Verificar se o fundo suavizado não sobe de montante para jusante
            z_bed = prof['z_bed_smooth']
            diffs = np.diff(z_bed)
            self.assertTrue(np.all(diffs <= 0.0), f"Fundo do rio {r_key} possui subida de jusante")

    def test_spatial_stage_engine_calculation(self):
        """Verifica o cálculo de profundidade e cota absoluta da água."""
        prof = self.loader.get_river_profile('acu')
        n_pts = prof['n_points']
        n_times = 49
        
        # Criar matriz sintética de cheia (passagem de onda de cheia de 100 a 4000 m³/s)
        q_mat = np.zeros((n_pts, n_times))
        for t in range(n_times):
            peak_factor = np.sin(np.pi * (t / 48.0)) ** 2
            q_mat[:, t] = 100.0 + 3900.0 * peak_factor

        result = self.engine.compute_reach_depth_and_stage('acu', q_mat)
        
        self.assertEqual(result['river_key'], 'acu')
        self.assertEqual(result['depth_h_m'].shape, (n_pts, n_times))
        self.assertEqual(result['z_water_m'].shape, (n_pts, n_times))
        
        # Verificar se a profundidade é sempre positiva
        self.assertTrue(np.all(result['depth_h_m'] >= 0.1))
        
        # Verificar se Z_water = Z_bed + H
        for t in range(n_times):
            diff = result['z_water_m'][:, t] - (result['z_bed_m'] + result['depth_h_m'][:, t])
            self.assertTrue(np.all(np.abs(diff) < 0.5))

    def test_water_surface_downstream_monotonicity(self):
        """A linha d'água deve descer estritamente em direção à foz."""
        prof = self.loader.get_river_profile('acu')
        n_pts = prof['n_points']
        q_mat = np.ones((n_pts, 24)) * 2500.0 # Cheia de 2500 m³/s
        
        result = self.engine.compute_reach_depth_and_stage('acu', q_mat)
        z_water = result['z_water_m']
        
        for t in range(24):
            diffs = np.diff(z_water[:, t])
            self.assertTrue(np.all(diffs <= 0.0), f"Linha d'água no tempo {t}h sobe de montante para jusante!")

    def test_flood_overtopping_detection(self):
        """Verifica a detecção de transbordo (overtopping) nas margens."""
        prof = self.loader.get_river_profile('acu')
        n_pts = prof['n_points']
        
        # Vazão muito baixa (50 m³/s) -> não deve transbordar
        q_low = np.ones((n_pts, 10)) * 50.0
        res_low = self.engine.compute_reach_depth_and_stage('acu', q_low)
        self.assertFalse(np.any(res_low['is_overtopping'][:, 0]), "Vazão de estiagem falsamente gerou transbordo")
        
        # Vazão extrema (5500 m³/s - Cheia de 1983) -> deve transbordar a calha
        q_flood = np.ones((n_pts, 10)) * 5500.0
        res_flood = self.engine.compute_reach_depth_and_stage('acu', q_flood)
        self.assertTrue(np.any(res_flood['is_overtopping'][:, 0]), "Cheia extrema não detectou transbordo")

if __name__ == '__main__':
    unittest.main()
