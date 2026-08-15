"""
Testes unitários para o modelo de propagação de Muskingum e Hidrograma Unitário.
"""

import unittest
import numpy as np
import sys
import os

# Adicionar src ao path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.muskingum import MuskingumReach
from src.muskingum_cunge import MuskingumCungeReach
from src.unit_hydrograph import UnitHydrograph, scs_curve_number_excess
from src.calibration import HydrographValidator

class TestMuskingumHydrology(unittest.TestCase):

    def test_muskingum_mass_conservation(self):
        """Verifica se a conservação de massa é rigorosamente mantida (Integral(Qin) = Integral(Qout))."""
        reach = MuskingumReach(reach_id=1, name="Trecho Teste", k_hours=6.0, x_param=0.20, dt_hours=1.0)
        
        # Hidrograma sintético de entrada
        t = np.arange(0, 50, 1.0)
        inflow = 15.0 + 300.0 * np.exp(-0.5 * ((t - 15.0) / 4.0) ** 2)
        
        outflow = reach.route(inflow)
        
        vol_in = np.sum(inflow) * 3600.0
        vol_out = np.sum(outflow) * 3600.0
        
        # Erro relativo de volume menor que 0.5%
        rel_error = abs(vol_in - vol_out) / vol_in
        self.assertLess(rel_error, 0.005)

    def test_muskingum_attenuation_and_lag(self):
        """Verifica se há atenuação física do pico e retardo temporal."""
        reach = MuskingumReach(reach_id=1, name="Trecho Teste", k_hours=8.0, x_param=0.20, dt_hours=1.0)
        
        t = np.arange(0, 50, 1.0)
        inflow = 10.0 + 500.0 * np.exp(-0.5 * ((t - 12.0) / 3.0) ** 2)
        outflow = reach.route(inflow)
        
        peak_in = np.max(inflow)
        t_peak_in = np.argmax(inflow)
        
        peak_out = np.max(outflow)
        t_peak_out = np.argmax(outflow)
        
        # O pico de saída deve ser menor ou igual ao de entrada
        self.assertLessEqual(peak_out, peak_in)
        # O horário do pico de saída deve ser posterior ao de entrada
        self.assertGreaterEqual(t_peak_out, t_peak_in)

    def test_coefficient_sum(self):
        """Garante que a soma C0 + C1 + C2 = 1.0 com alta precisão."""
        reach = MuskingumReach(reach_id=1, name="Trecho Teste", k_hours=12.0, x_param=0.25, dt_hours=1.0)
        sum_c = reach.c0 + reach.c1 + reach.c2
        self.assertAlmostEqual(sum_c, 1.0, places=6)

    def test_muskingum_cunge_derivation(self):
        """Verifica derivação física de K e X pelo método Cunge."""
        cunge = MuskingumCungeReach(
            reach_id=1, name="Cunge Teste", length_km=25.0, slope_m_km=1.5,
            width_m=40.0, manning_n=0.035, reference_q_m3s=250.0, dt_hours=1.0
        )
        self.assertGreater(cunge.k_hours, 0)
        self.assertTrue(0.0 <= cunge.x_param <= 0.5)
        self.assertGreater(cunge.celerity_ms, 0)

    def test_validation_metrics(self):
        """Testa o cálculo do NSE e RMSE."""
        obs = np.array([10, 20, 50, 100, 70, 30, 15], dtype=float)
        sim = np.array([10, 22, 48, 95, 68, 32, 16], dtype=float)
        
        metrics = HydrographValidator.calculate_metrics(sim, obs, dt_hours=1.0)
        self.assertGreater(metrics['nse'], 0.95)
        self.assertLess(metrics['rmse_m3s'], 5.0)

if __name__ == '__main__':
    unittest.main()
