"""
Suíte de Testes Hidráulicos e de Consistência 2D (ETAPA CRÍTICA):
Valida os 8 requisitos fundamentais da hidráulica de inundação:
1. TESTE 1: Q = 0 -> Área inundada na várzea = 0
2. TESTE 2: Nível baixo -> Área confinada / pequena
3. TESTE 3: Nível elevado -> Área maior de extravasamento
4. TESTE 4: Nível crescente -> Área não pode diminuir arbitrariamente (Monotonicidade)
5. TESTE 5: Conectividade -> Depressão isolada não é inundada no mapa conectado
6. TESTE 6: Confluência -> Conservação rigorosa de vazão (sum Q_in = Q_out)
7. TESTE 7: Bifurcação do Itajaí-Mirim -> Q_total = Q_canal + Q_braço
8. TESTE 8: Foz e Remanso -> Variação de H_ocean(t) altera cota para montante (sem restrição dZ/dx <= 0)
"""

import unittest
import numpy as np
from pathlib import Path

from itajai_flood_model.src.inundation.hydraulic_network import HydraulicNode, RiverBranch, ItajaiHydraulicNetwork
from itajai_flood_model.src.inundation.water_surface_raster import WaterSurfaceRasterEngine
from itajai_flood_model.src.inundation.validation_metrics import FloodValidationMetrics

class TestHydraulicInundation2D(unittest.TestCase):

    def setUp(self):
        self.raster_engine = WaterSurfaceRasterEngine(grid_shape=(50, 50))
        self.raster_engine.build_synthetic_valley_dem()
        
        # Rio sintético central
        river_coords = [(-50.0 + 0.025 * i, -27.05) for i in range(40)]
        z_bed = [300.0 - 6.0 * i for i in range(40)]
        self.sample_branch = RiverBranch("teste", "Rio Teste", river_coords, z_bed, "no_up", "no_down", b_start=40, b_end=80)
        self.raster_engine.rasterize_river_corridor({'teste': self.sample_branch})

    def test_1_zero_flow_zero_floodplain_area(self):
        """TESTE 1: Q = 0 -> Área inundada na várzea = 0."""
        # Se cota da água estiver abaixo da cota do terreno
        z_water_dry = np.ones((50, 50)) * -10.0
        res = self.raster_engine.compute_2d_inundation(z_water_dry)
        self.assertEqual(res['area_connected_km2'], 0.0)
        self.assertEqual(res['volume_connected_hm3'], 0.0)

    def test_2_and_3_stage_area_progression(self):
        """TESTE 2 e 3: Nível baixo -> área pequena; Nível alto -> área maior."""
        z_water_low = np.ones((50, 50)) * 120.0
        res_low = self.raster_engine.compute_2d_inundation(z_water_low)
        
        z_water_high = np.ones((50, 50)) * 250.0
        res_high = self.raster_engine.compute_2d_inundation(z_water_high)
        
        self.assertGreater(res_high['area_connected_km2'], res_low['area_connected_km2'])
        self.assertGreater(res_high['volume_connected_hm3'], res_low['volume_connected_hm3'])

    def test_4_monotonicity_stage_area(self):
        """TESTE 4: Nível crescente -> Área não pode diminuir arbitrariamente."""
        stages = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0]
        areas = []
        for s in stages:
            z_w = np.ones((50, 50)) * s
            res = self.raster_engine.compute_2d_inundation(z_w)
            areas.append(res['area_connected_km2'])
            
        for i in range(1, len(areas)):
            self.assertGreaterEqual(areas[i], areas[i-1], f"Área diminuiu de stage {stages[i-1]} para {stages[i]}")

    def test_5_connectivity_eliminates_isolated_pits(self):
        """TESTE 5: Depressão isolada não é inundada no mapa conectado."""
        # Injetar uma depressão isolada longe do corredor do rio
        dem_with_pit = np.copy(self.raster_engine.z_dem)
        dem_with_pit[5:10, 5:10] = 5.0 # Depressão rebaixada no canto
        self.raster_engine.z_dem = dem_with_pit
        
        # Permitir corredor cobrir a área para testar o filtro
        self.raster_engine.river_corridor_mask[5:10, 5:10] = True
        # As sementes do canal estão longe (linha central)
        
        z_w = np.ones((50, 50)) * 50.0
        res = self.raster_engine.compute_2d_inundation(z_w)
        
        # O mapa geométrico bruto DEVE detectar a poça
        self.assertGreater(res['area_geometric_km2'], res['area_connected_km2'])
        self.assertGreater(res['isolated_pits_area_km2'], 0.0)
        # As células da depressão desconectada (5:10, 5:10) DEVEM ter profundidade zero no mapa conectado
        self.assertTrue(np.all(res['connected_depth_m'][5:10, 5:10] == 0.0))

    def test_6_confluence_mass_balance(self):
        """TESTE 6: Confluência -> Conservação rigorosa de vazão (sum Q_in = Q_out)."""
        node_conf = HydraulicNode("conf_rio_sul", "Rio do Sul", -49.64, -27.21, 330.0)
        node_conf.inflow_rivers = ["oeste", "sul", "trombudo"]
        
        inflows = {"oeste": 1450.0, "sul": 820.0, "trombudo": 350.0}
        q_out = node_conf.compute_confluence_balance(inflows)
        
        self.assertEqual(q_out, 1450.0 + 820.0 + 350.0) # 2620.0 m³/s
        self.assertEqual(node_conf.q_in_total, node_conf.q_out_total)

    def test_7_itajai_mirim_bifurcation_split(self):
        """TESTE 7: Bifurcação do Itajaí-Mirim -> Q_total = Q_canal + Q_braco."""
        node_bif = HydraulicNode("bif_mirim", "Bifurcação Mirim", -48.74, -26.96, 8.0, node_type="bifurcation")
        
        q_total_mirim = 1650.0 # Cheia de 2008
        splits = node_bif.compute_bifurcation_split(q_total_mirim, {"canal_retificado": 0.70, "braco_velho": 0.30})
        
        self.assertAlmostEqual(splits['canal_retificado'], 1155.0) # 70%
        self.assertAlmostEqual(splits['braco_velho'], 495.0)       # 30%
        self.assertAlmostEqual(splits['canal_retificado'] + splits['braco_velho'], q_total_mirim)

    def test_8_ocean_boundary_backwater_effect(self):
        """TESTE 8: Foz e Remanso -> Variação de H_ocean altera cota para montante sem bloqueio dZ/dx <= 0."""
        # Criar trecho de estuário (Gaspar -> Itajaí Foz)
        coords = [(-49.0 + 0.05 * i, -26.9) for i in range(10)]
        z_bed = [1.3 - 0.5 * i for i in range(10)] # Fundo cai de 1.3m para -3.2m
        estuary = RiverBranch("estuario", "Baixo Vale", coords, z_bed, "no_gaspar", "no_foz", b_start=150, b_end=250)
        
        q_flow = np.ones(10) * 1200.0
        
        # 1. Maré normal (Z_ocean = 0.0m)
        res_normal_tide = estuary.compute_backwater_water_surface(q_flow, downstream_z_water=0.0)
        
        # 2. Maré meteorológica alta (Storm surge: Z_ocean = 3.5m)
        res_high_tide = estuary.compute_backwater_water_surface(q_flow, downstream_z_water=3.5)
        
        # A cota da água na foz DEVE subir com a maré
        self.assertGreater(res_high_tide['z_water_m'][-1], res_normal_tide['z_water_m'][-1])
        # O remanso DEVE se propagar para as seções imediatamente a montante da foz
        self.assertGreater(res_high_tide['z_water_m'][-2], res_normal_tide['z_water_m'][-2])

if __name__ == '__main__':
    unittest.main()
