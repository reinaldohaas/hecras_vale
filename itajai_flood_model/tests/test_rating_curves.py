"""
Suíte de Testes Unitários do Módulo de Curva-Chave (Rating Curve Q-H).
Executável tanto com pytest quanto diretamente via Python nativo.
"""

import sys
from pathlib import Path

# Adicionar raiz do repositório ao sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import unittest
import numpy as np

from itajai_flood_model.src.rating_curve import (
    RatingCurveManager,
    CurveType,
    create_blumenau_official_curve,
    create_rio_do_sul_official_curve,
    create_brusque_official_curve,
    create_indaial_official_curve,
    CrossSectionGeometry,
    HydraulicRatingCurve
)

class TestRatingCurves(unittest.TestCase):
    
    def test_blumenau_official_benchmarks(self):
        """Valida a curva oficial de Blumenau contra os grandes registros históricos."""
        curve = create_blumenau_official_curve()
        self.assertEqual(curve.curve_type, CurveType.OFFICIAL_OBSERVED)
        self.assertEqual(curve.station_id, "83700000")

        # 1983: 15.34m -> ~5850 m³/s
        q_1983 = curve.to_flow(15.34)
        self.assertTrue(5600 <= q_1983 <= 6100, f"Vazão 1983 fora da faixa: {q_1983}")
        h_back_1983 = curve.to_stage(q_1983)
        self.assertAlmostEqual(h_back_1983, 15.34, delta=0.05)

        # 2011: 12.60m -> ~4650 m³/s
        q_2011 = curve.to_flow(12.60)
        self.assertTrue(4400 <= q_2011 <= 4900, f"Vazão 2011 fora da faixa: {q_2011}")
        h_back_2011 = curve.to_stage(q_2011)
        self.assertAlmostEqual(h_back_2011, 12.60, delta=0.05)

        # 2008: 11.52m -> ~4200 m³/s
        q_2008 = curve.to_flow(11.52)
        self.assertTrue(3900 <= q_2008 <= 4400, f"Vazão 2008 fora da faixa: {q_2008}")

        # 2023: 10.76m -> ~3950 m³/s
        q_2023 = curve.to_flow(10.76)
        self.assertTrue(3600 <= q_2023 <= 4100, f"Vazão 2023 fora da faixa: {q_2023}")

    def test_rio_do_sul_official_curve(self):
        """Valida curva oficial de Rio do Sul (Confluência)."""
        curve = create_rio_do_sul_official_curve()
        self.assertEqual(curve.curve_type, CurveType.OFFICIAL_OBSERVED)
        
        # 13.0m (1983) -> ~3900 m³/s
        q_rs = curve.to_flow(13.0)
        self.assertTrue(3600 <= q_rs <= 4200)
        h_back = curve.to_stage(q_rs)
        self.assertAlmostEqual(h_back, 13.0, delta=0.05)

    def test_brusque_official_curve(self):
        """Valida curva oficial de Brusque (Rio Itajaí-Mirim)."""
        curve = create_brusque_official_curve()
        self.assertEqual(curve.curve_type, CurveType.OFFICIAL_OBSERVED)
        
        # 8.5m (2008) -> ~1650 m³/s
        q_bq = curve.to_flow(8.5)
        self.assertTrue(1450 <= q_bq <= 1800)
        h_back = curve.to_stage(q_bq)
        self.assertAlmostEqual(h_back, 8.5, delta=0.05)

    def test_hydraulic_manning_curve(self):
        """Valida a curva estimada por Manning e geometria da seção transversal."""
        geo = CrossSectionGeometry(
            bottom_width_b_m=80.0,
            side_slope_z=1.5,
            bankfull_depth_m=6.0,
            floodplain_width_m=100.0,
            manning_n_main=0.035,
            manning_n_floodplain=0.065
        )
        curve = HydraulicRatingCurve(
            station_id="TEST_EST",
            name="Seção Teste Manning",
            river="Rio Teste",
            geometry=geo,
            bed_slope_s0=0.0005,
            datum_z0_m=20.0
        )
        self.assertEqual(curve.curve_type, CurveType.ESTIMATED_HYDRAULIC)
        
        # Teste de monotonicidade
        h_series = np.linspace(1.0, 12.0, 10)
        q_series = curve.to_flow(h_series)
        self.assertTrue(np.all(np.diff(q_series) > 0), "A curva de vazão deve ser estritamente crescente com a cota")

        # Teste de reversibilidade Q -> H -> Q
        h_orig = 7.5
        q_val = curve.to_flow(h_orig)
        h_calc = curve.to_stage(q_val)
        self.assertAlmostEqual(h_calc, h_orig, delta=0.05)

        # Teste de detalhes da seção
        det = curve.get_section_details(h_orig)
        self.assertGreater(det['wet_area_m2'], 0)
        self.assertGreater(det['hydraulic_radius_m'], 0)
        self.assertEqual(det['flow_q_m3s'], q_val)

    def test_manager_catalog(self):
        """Valida o catálogo integrado do RatingCurveManager."""
        mgr = RatingCurveManager()
        stations = mgr.list_stations()
        self.assertGreaterEqual(len(stations), 5)
        
        # Verificar distinção estrita oficial vs estimada
        blu = mgr.get_curve('blumenau')
        self.assertEqual(blu.curve_type, CurveType.OFFICIAL_OBSERVED)
        
        ibi = mgr.get_curve('ibirama')
        self.assertEqual(ibi.curve_type, CurveType.ESTIMATED_HYDRAULIC)

        # Conversões via manager
        h_blu = mgr.flow_to_stage('blumenau', 2400.0)
        self.assertTrue(7.0 <= h_blu <= 9.0) # Nível de alerta Blumenau (~8.0m)
        
        q_blu = mgr.stage_to_flow('blumenau', 8.0)
        self.assertTrue(2200 <= q_blu <= 2600)

if __name__ == '__main__':
    unittest.main()
