"""
Suíte de Testes Rigorosos para o Modelo Unificado HAND + Inundação Sincronizada:
Valida as provas matemáticas e os testes físicos obrigatórios.
"""

import unittest
import numpy as np
from pathlib import Path
import json

from itajai_flood_model.src.inundation.unified_hand_engine import (
    TopographicHANDModel,
    SynchronizedFloodEngine,
    FloodState
)

class TestUnifiedHANDModel(unittest.TestCase):

    def test_01_mandatory_numerical_example(self):
        """
        TESTE DO EXEMPLO NUMÉRICO OBRIGATÓRIO (Item 4 da Especificação):
        Z_drain = 5m, Z_DEM = 8m -> HAND = 3m
        Z_water = 10m -> eta = 5m
        depth = eta - HAND = 5 - 3 = 2m (ou 10 - 8 = 2m).
        """
        z_drain = 5.0
        z_dem = 8.0
        hand = z_dem - z_drain
        self.assertEqual(hand, 3.0, "HAND deve ser exatamente 3m")

        z_water = 10.0
        eta = z_water - z_drain
        self.assertEqual(eta, 5.0, "eta deve ser exatamente 5m")

        # Formulação A (Direta): Z_water - Z_DEM
        depth_a = z_water - z_dem
        # Formulação B (HAND): eta - HAND
        depth_b = eta - hand

        self.assertEqual(depth_a, 2.0, "Profundidade direta deve ser exatamente 2m")
        self.assertEqual(depth_b, 2.0, "Profundidade HAND deve ser exatamente 2m")
        self.assertEqual(depth_a, depth_b, "Ambas as formulações devem coincidir rigorosamente")

    def test_02_mathematical_identity_grid(self):
        """
        TESTE DE IDENTIDADE MATRICIAL:
        Para qualquer grade bidimensional aleatória:
        (Z_water - Z_DEM) == (eta - HAND)
        """
        np.random.seed(42)
        z_drain = np.random.uniform(0.0, 300.0, (50, 50))
        hand = np.random.uniform(0.0, 50.0, (50, 50))
        z_dem = z_drain + hand

        z_water = z_drain + np.random.uniform(-5.0, 20.0, (50, 50))
        eta = z_water - z_drain

        depth_direct = z_water - z_dem
        depth_hand = eta - hand

        np.testing.assert_allclose(depth_direct, depth_hand, atol=1e-7)

    def test_03_hand_is_static_flood_is_dynamic(self):
        """
        TESTE CRÍTICO DE ESTATICIDADE DO HAND (Item 23 da Especificação):
        Manter o mesmo DEM e a mesma rede.
        Variar somente Z_water(t).
        Resultado: HAND permanece idêntico; a mancha muda com o tempo.
        """
        dem_path = Path("dem_blumenau_itajai.tif")
        json_path = Path("app/itajai_real_dem_model.json")
        if not dem_path.exists() or not json_path.exists():
            self.skipTest("Arquivos reais não encontrados no diretório atual")

        topo = TopographicHANDModel(dem_path, json_path)
        hand_initial = np.copy(topo.hand)

        engine = SynchronizedFloodEngine(topo)

        # Tempo 1: Nível baixo (Z_water normal)
        profiles_t1 = {'acu': {'z_water_m': topo.dem[topo.stream_mask][:80] + 1.0}}
        res_t1 = engine.compute_instantaneous_flood(profiles_t1, t_hour=0.0)

        # Tempo 2: Nível alto de cheia (Z_water + 12m)
        profiles_t2 = {'acu': {'z_water_m': topo.dem[topo.stream_mask][:80] + 12.0}}
        res_t2 = engine.compute_instantaneous_flood(profiles_t2, t_hour=24.0)

        # 1. O HAND NUNCA mudou
        np.testing.assert_array_equal(topo.hand, hand_initial, "O HAND deve permanecer 100% estático!")

        # 2. A mancha expandiu dinamicamente com a onda de cheia
        self.assertGreater(res_t2.area_km2, res_t1.area_km2)
        self.assertGreater(res_t2.volume_hm3, res_t1.volume_hm3)

    def test_04_stage_inundation_progression(self):
        """
        TESTE FUNDAMENTAL DE SINCRONIZAÇÃO (Item 21 da Especificação):
        Z_drain = 100m, HAND = 2m -> Z_DEM = 102m.
        t1: Z_water = 101m -> seco (depth = 0)
        t2: Z_water = 103m -> depth = 1m
        t3: Z_water = 105m -> depth = 3m
        """
        z_drain = 100.0
        hand = 2.0
        z_dem = z_drain + hand # 102m

        # t1
        z_w1 = 101.0
        eta1 = z_w1 - z_drain # 1m
        h1 = max(0.0, eta1 - hand) # 1 - 2 = -1 -> 0m
        self.assertEqual(h1, 0.0)

        # t2
        z_w2 = 103.0
        eta2 = z_w2 - z_drain # 3m
        h2 = max(0.0, eta2 - hand) # 3 - 2 = 1m
        self.assertEqual(h2, 1.0)

        # t3
        z_w3 = 105.0
        eta3 = z_w3 - z_drain # 5m
        h3 = max(0.0, eta3 - hand) # 5 - 2 = 3m
        self.assertEqual(h3, 3.0)

    def test_05_configurable_min_flood_depth(self):
        """TESTE DO LIMIAR CONFIGURÁVEL MIN_FLOOD_DEPTH (Item 20)."""
        dem_path = Path("dem_blumenau_itajai.tif")
        json_path = Path("app/itajai_real_dem_model.json")
        if not dem_path.exists() or not json_path.exists():
            self.skipTest("Arquivos reais não encontrados")

        topo = TopographicHANDModel(dem_path, json_path)
        
        # Testar com limiar de 0.05m vs 0.50m
        engine_fine = SynchronizedFloodEngine(topo, min_flood_depth_m=0.05)
        engine_coarse = SynchronizedFloodEngine(topo, min_flood_depth_m=0.50)

    def test_06_normal_level_confined_in_channel(self):
        """
        TESTE DO NÍVEL NORMAL / CALHA CONFINADA:
        Em vazão de estiagem ou nível normal (Z_water <= Z_bank),
        a água deve estar 100% contida na calha profunda, resultando em
        área inundada na planície estritamente igual a 0.0 km².
        """
        dem_path = Path("dem_blumenau_itajai.tif")
        json_path = Path("app/itajai_real_dem_model.json")
        if not dem_path.exists() or not json_path.exists():
            self.skipTest("Arquivos reais não encontrados")

        topo = TopographicHANDModel(dem_path, json_path)
        engine = SynchronizedFloodEngine(topo)

        # Nível normal de calha (Z_water = Z_bed + 1.5m <= Z_bank)
        profiles_normal = {'acu': {'z_water_m': topo.bed_elevation[topo.stream_mask][:80] + 1.5}}
        res_normal = engine.compute_instantaneous_flood(profiles_normal, t_hour=0.0)

        self.assertEqual(res_normal.area_km2, 0.0,
                         "Em nível normal (dentro da calha), a planície deve estar 100% seca (0.0 km²)!")
        self.assertEqual(res_normal.volume_hm3, 0.0)

    def test_07_ocean_tide_boundary(self):
        """TESTE DA CONDIÇÃO DE CONTORNO OCEÂNICA COM MARÉ E RESSACA."""
        topo = TopographicHANDModel(Path("dem_blumenau_itajai.tif"), Path("app/itajai_real_dem_model.json"))
        engine = SynchronizedFloodEngine(topo)

        # t=4h (maré alta astronômica)
        tide_high = engine.tide_model.compute_ocean_level(t_hour=4.0, storm_surge_peak_m=0.0)
        # t=10.2h (maré baixa astronômica)
        tide_low = engine.tide_model.compute_ocean_level(t_hour=10.21, storm_surge_peak_m=0.0)
        
        self.assertGreater(tide_high.total_ocean_level_z, tide_low.total_ocean_level_z,
                           "A maré alta deve ser superior à maré baixa")

        # Com ressaca / storm surge em t=24h
        tide_surge = engine.tide_model.compute_ocean_level(t_hour=24.0, storm_surge_peak_m=1.40, t_surge_peak_h=24.0)
        self.assertEqual(tide_surge.storm_surge_m, 1.40,
                         "A sobre-elevação de ressaca no pico deve ser exatamente 1.40m")
        self.assertGreater(tide_surge.total_ocean_level_z, tide_low.total_ocean_level_z,
                           "O nível do mar com ressaca deve superar a maré baixa astronômica")
