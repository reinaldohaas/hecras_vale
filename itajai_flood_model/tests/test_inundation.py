"""
Testes Unitários para o Módulo de Manchas de Inundação e Conectividade Hidráulica (ETAPA 3):
Verifica:
1. Filtro de conectividade hidráulica (eliminação de depressões fechadas isoladas)
2. Cálculo de área e volume na grade topográfica
3. Classificação de profundidade e faixas de risco
4. Geração de polígonos GeoJSON válidos para visualização
"""

import unittest
import numpy as np
from pathlib import Path

from itajai_flood_model.src.inundation import (
    InundationGrid,
    HydraulicConnectivityFilter,
    DepthClassifier,
    FloodplainMapper
)

class TestInundationMapping(unittest.TestCase):

    def test_connectivity_filter_eliminates_isolated_pits(self):
        """
        Garante que uma depressão isolada sem conexão com o rio seja eliminada (h=0),
        enquanto a várzea conectada ao rio permaneça inundada.
        """
        # Grade 10x10
        # Coluna 1: Rio (conectado)
        # Coluna 2-3: Várzea conectada
        # Coluna 5: Dique / Elevação natural (seco)
        # Coluna 7-8: Depressão fechada isolada
        
        depth_raw = np.zeros((10, 10))
        depth_raw[:, 1:4] = 1.5   # Várzea conectada (1.5m de água)
        depth_raw[:, 7:9] = 2.0   # Depressão isolada (2.0m de água aparente)
        
        river_mask = np.zeros((10, 10), dtype=bool)
        river_mask[:, 1] = True   # Rio na coluna 1
        
        filtered = HydraulicConnectivityFilter.apply_connectivity_filter(depth_raw, river_mask)
        
        # A várzea conectada (colunas 1 a 3) DEVE estar inundada
        self.assertTrue(np.all(filtered[:, 1:4] == 1.5))
        
        # A depressão isolada (colunas 7 a 8) DEVE ser ZERADA (eliminada pelo filtro)
        self.assertTrue(np.all(filtered[:, 7:9] == 0.0))

    def test_inundation_grid_area_volume(self):
        """Verifica o cálculo de área inundada (km²) e volume de cheia (hm³)."""
        grid = InundationGrid(bounds=(-49.5, -27.5, -48.5, -26.5), cell_size_m=30.0, nrows=50, ncols=50)
        
        z_dem = np.ones((50, 50)) * 10.0
        z_dem[20:30, 20:30] = 5.0 # Várzea rebaixada a 5m
        grid.set_topography(z_dem)
        
        # Semente do rio
        river_coords = [(-49.0, -27.0)]
        grid.rasterize_river_line(river_coords)
        
        # Superfície da água em 8.0m
        z_water = np.ones((50, 50)) * 8.0
        
        res = grid.compute_flood_depths(z_water, apply_connectivity=False)
        
        # Apenas as células rebaixadas (10x10 = 100 células) devem inundar com 3m
        self.assertEqual(res['wet_cells_count'], 100)
        expected_area_km2 = 100 * (30.0 * 30.0) / 1e6 # 0.09 km²
        self.assertAlmostEqual(res['flooded_area_km2'], expected_area_km2, places=2)
        self.assertAlmostEqual(res['max_depth_m'], 3.0, places=1)

    def test_depth_classifier(self):
        """Verifica a contagem e percentuais das classes de risco."""
        depth_mat = np.array([
            [0.0, 0.2, 0.8],
            [1.8, 3.5, 0.0]
        ])
        res = DepthClassifier.classify_depths(depth_mat)
        
        self.assertEqual(res['total_wet_cells'], 4)
        self.assertEqual(res['count_low'], 1)       # 0.2m
        self.assertEqual(res['count_medium'], 1)    # 0.8m
        self.assertEqual(res['count_high'], 1)      # 1.8m
        self.assertEqual(res['count_very_high'], 1) # 3.5m

    def test_flood_mapper_corridor_polygons(self):
        """Verifica a geração de polígonos GeoJSON."""
        coords = [(-49.066, -26.918), (-49.060, -26.915), (-49.055, -26.910)]
        depths = np.array([9.5, 10.2, 8.8]) # Acima de h_bank = 7.5
        
        geojson = FloodplainMapper.generate_flood_corridor_polygons(coords, depths, h_bank=7.5)
        
        self.assertEqual(geojson['type'], 'FeatureCollection')
        self.assertGreater(len(geojson['features']), 0)
        self.assertEqual(geojson['features'][0]['geometry']['type'], 'Polygon')
        self.assertEqual(len(geojson['features'][0]['geometry']['coordinates'][0]), 5)

    def test_river_cross_section_top_width(self):
        """Verifica a largura de inundação calculada a partir do perfil transversal."""
        from itajai_flood_model.src.inundation import RiverCrossSection, CrossSectionDelineator
        
        sec = RiverCrossSection(
            station_id=1, river_name="acu",
            lon=-49.066, lat=-26.918, dist_km=105.0,
            z_bed=1.30, h_bank=7.50, # z_bank = 8.80m
            channel_width_m=150.0,
            slope_left_m_per_m=0.0035, # 1m a cada ~285m
            slope_right_m_per_m=0.0040 # 1m a cada 250m
        )
        
        # 1. Dentro do leito menor (Z = 5.0m < 8.80m)
        w_in = sec.compute_top_width(z_water=5.0)
        self.assertEqual(w_in['w_left_m'], 75.0)
        self.assertEqual(w_in['w_right_m'], 75.0)
        self.assertEqual(w_in['h_overbank_m'], 0.0)
        
        # 2. Extravasamento de cheia (Z = 15.34m > 8.80m -> h_overbank = 6.54m)
        w_over = sec.compute_top_width(z_water=15.34)
        self.assertGreater(w_over['h_overbank_m'], 6.0)
        self.assertGreater(w_over['w_left_m'], 1500.0) # Várzea esquerda inundada >1.5km
        self.assertGreater(w_over['w_right_m'], 1400.0) # Várzea direita inundada >1.4km

if __name__ == '__main__':
    unittest.main()
