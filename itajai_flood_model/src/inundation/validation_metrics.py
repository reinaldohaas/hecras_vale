"""
Módulo de Métricas de Validação Espacial e Comparação com Eventos Históricos (ValidationMetrics):
Calcula:
1. Intersection over Union (IoU / Jaccard Index)
2. Precision (Precisão), Recall (Sensibilidade) e F1-Score
3. Erro de Cota Máxima: Delta H_max = |H_sim - H_obs|
4. Erro de Tempo da Crista: Delta t_peak = |t_sim_peak - t_obs_peak|
5. Erro Relativo de Área Inundada: |A_sim - A_obs| / A_obs
"""

from typing import Dict, Any, Optional
import numpy as np

class FloodValidationMetrics:
    """Calcula métricas estatísticas e de validação espacial para manchas de inundação."""

    @staticmethod
    def compute_spatial_metrics(simulated_mask: np.ndarray,
                                observed_mask: np.ndarray,
                                cell_area_km2: float = 0.0009) -> Dict[str, float]:
        """
        Calcula as métricas de sobreposição espacial binária (True Positive, False Positive, False Negative).
        """
        sim = np.asarray(simulated_mask, dtype=bool)
        obs = np.asarray(observed_mask, dtype=bool)

        tp = int(np.sum(sim & obs))       # Ambos inundados (Verdadeiro Positivo)
        fp = int(np.sum(sim & ~obs))      # Apenas simulado (Falso Positivo / Sobre-estimativa)
        fn = int(np.sum(~sim & obs))      # Apenas observado (Falso Negativo / Sub-estimativa)
        tn = int(np.sum(~sim & ~obs))     # Ambos secos

        area_sim_km2 = float(np.sum(sim) * cell_area_km2)
        area_obs_km2 = float(np.sum(obs) * cell_area_km2)
        area_overlap_km2 = float(tp * cell_area_km2)

        # 1. Intersection over Union (IoU / Critical Success Index CSI)
        union = tp + fp + fn
        iou = float(tp / max(1, union))

        # 2. Precision
        precision = float(tp / max(1, tp + fp))

        # 3. Recall (Hit Rate)
        recall = float(tp / max(1, tp + fn))

        # 4. F1-Score
        f1 = float((2 * precision * recall) / max(1e-6, precision + recall))

        # 5. Erro Relativo de Área
        area_err_pct = float(abs(area_sim_km2 - area_obs_km2) / max(1e-3, area_obs_km2) * 100.0)

        return {
            'iou_csi': round(iou, 4),
            'precision': round(precision, 4),
            'recall': round(recall, 4),
            'f1_score': round(f1, 4),
            'area_simulated_km2': round(area_sim_km2, 2),
            'area_observed_km2': round(area_obs_km2, 2),
            'area_overlap_km2': round(area_overlap_km2, 2),
            'area_error_pct': round(area_err_pct, 2),
            'true_positives_count': tp,
            'false_positives_count': fp,
            'false_negatives_count': fn
        }

    @staticmethod
    def compute_hydrograph_metrics(t_series: np.ndarray,
                                   q_sim: np.ndarray, q_obs: np.ndarray,
                                   h_sim: np.ndarray, h_obs: np.ndarray) -> Dict[str, float]:
        """Calcula erros de pico, volume e Nash-Sutcliffe Efficiency (NSE)."""
        q_s = np.asarray(q_sim, dtype=float)
        q_o = np.asarray(q_obs, dtype=float)
        h_s = np.asarray(h_sim, dtype=float)
        h_o = np.asarray(h_obs, dtype=float)

        peak_q_sim = float(np.max(q_s))
        peak_q_obs = float(np.max(q_o))
        peak_h_sim = float(np.max(h_s))
        peak_h_obs = float(np.max(h_o))

        t_peak_sim = float(t_series[np.argmax(q_s)])
        t_peak_obs = float(t_series[np.argmax(q_o)])

        err_peak_q = float(abs(peak_q_sim - peak_q_obs))
        err_peak_h = float(abs(peak_h_sim - peak_h_obs))
        err_t_peak = float(abs(t_peak_sim - t_peak_obs))

        # Nash-Sutcliffe Efficiency (NSE)
        denom = np.sum((q_o - np.mean(q_o)) ** 2)
        nse = 1.0 - (np.sum((q_s - q_o) ** 2) / max(1e-6, denom)) if denom > 0 else 1.0

        return {
            'nse_efficiency': round(float(nse), 4),
            'peak_flow_sim_m3s': round(peak_q_sim, 1),
            'peak_flow_obs_m3s': round(peak_q_obs, 1),
            'peak_flow_error_m3s': round(err_peak_q, 1),
            'peak_stage_sim_m': round(peak_h_sim, 2),
            'peak_stage_obs_m': round(peak_h_obs, 2),
            'peak_stage_error_m': round(err_peak_h, 2),
            'peak_time_sim_h': round(t_peak_sim, 1),
            'peak_time_obs_h': round(t_peak_obs, 1),
            'peak_timing_error_h': round(err_t_peak, 1)
        }
