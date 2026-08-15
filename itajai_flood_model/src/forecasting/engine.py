"""
Motor Operacional de Previsão de Vazões (OperationalForecastEngine):
- Transforma precipitação distribuída em vazões para cada sub-bacia
- Inclui todos os afluentes do Alto Vale: Rio Perimbó, Rio Trombudo e Rio Mirim Doce
- Realiza propagação por Muskingum estável e acoplamento com represas
- Produz previsões com incerteza [Q_low, Q_mean, Q_high] para os pontos críticos
- Calcula crista, tempo de chegada, tendência e horizonte de previsão
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any, Union

from ..unit_hydrograph import UnitHydrograph, scs_curve_number_excess
from ..muskingum import MuskingumReach
from ..rainfall.antecedent_moisture import AntecedentMoistureCondition

class OperationalForecastEngine:
    """
    Motor hidrológico para previsão operacional contínua na Bacia do Rio Itajaí.
    """
    def __init__(self, subbasin_params: Optional[Dict[str, Dict[str, Any]]] = None,
                 dam_operations: Optional[Dict[str, Dict[str, Any]]] = None):
        """
        Parâmetros das sub-bacias e operação das barragens.
        """
        self.subbasin_params = subbasin_params or {
            'oeste': {'area': 2480.0, 'tc': 14.2, 'cn_ii': 76.0, 'base_flow': 12.0, 'name': 'Rio do Oeste (Taió)'},
            'mirim_doce': {'area': 640.0, 'tc': 7.8, 'cn_ii': 77.0, 'base_flow': 7.0, 'name': 'Rio Mirim Doce (Afluente Oeste)'},
            'sul': {'area': 1250.0, 'tc': 12.8, 'cn_ii': 78.0, 'base_flow': 7.0, 'name': 'Rio do Sul (Ituporanga)'},
            'perimbo': {'area': 510.0, 'tc': 6.8, 'cn_ii': 78.0, 'base_flow': 5.0, 'name': 'Rio Perimbó (Petrolândia/Ituporanga)'},
            'trombudo': {'area': 720.0, 'tc': 8.5, 'cn_ii': 78.0, 'base_flow': 8.0, 'name': 'Rio Trombudo (Afluente Sul)'},
            'norte': {'area': 3450.0, 'tc': 16.5, 'cn_ii': 75.0, 'base_flow': 20.0, 'name': 'Rio Hercílio (José Boiteux)'},
            'benedito': {'area': 1540.0, 'tc': 10.4, 'cn_ii': 77.0, 'base_flow': 10.0, 'name': 'Rio Benedito (Timbó)'},
            'mirim': {'area': 1680.0, 'tc': 11.8, 'cn_ii': 79.0, 'base_flow': 18.0, 'name': 'Rio Itajaí-Mirim (Brusque)'},
            'luis_alves': {'area': 580.0, 'tc': 7.2, 'cn_ii': 78.0, 'base_flow': 8.0, 'name': 'Rio Luís Alves'}
        }
        
        self.dam_operations = dam_operations or {
            'oeste': {'total_gates': 7, 'open_gates': 7, 'cap_hm3': 83.0, 'base_flow': 12.0},
            'sul': {'total_gates': 5, 'open_gates': 5, 'cap_hm3': 93.5, 'base_flow': 7.0},
            'norte': {'total_gates': 2, 'open_gates': 2, 'cap_hm3': 357.0, 'base_flow': 20.0}
        }

    def _apply_dam(self, q_in: np.ndarray, dam_cfg: Dict[str, Any]) -> np.ndarray:
        """Simula amortecimento e operação de comportas de fundo."""
        tot = dam_cfg.get('total_gates', 1)
        open_g = dam_cfg.get('open_gates', tot)
        cap_hm3 = dam_cfg.get('cap_hm3', 100.0)
        base_q = dam_cfg.get('base_flow', 10.0)
        
        if open_g == tot:
            return q_in.copy()
            
        frac = open_g / float(tot)
        max_cap_m3 = cap_hm3 * 1e6
        dt_sec = 3600.0
        
        q_out = np.zeros_like(q_in)
        storage = 0.0
        
        for t in range(len(q_in)):
            if storage < max_cap_m3:
                q_allowed = base_q + frac * max(0.0, q_in[t] - base_q)
                q_out[t] = q_allowed
                storage += max(0.0, q_in[t] - q_allowed) * dt_sec
            else:
                # Vertimento livre
                q_out[t] = q_in[t]
        return q_out

    def run_subbasin_runoff(self, p_series: np.ndarray, sb_info: Dict[str, Any], amc: str = 'AMC_II') -> np.ndarray:
        """Calcula o hidrograma afluente de uma sub-bacia."""
        cn_eff = AntecedentMoistureCondition.adjust_curve_number(sb_info['cn_ii'], amc)
        pe_inc = scs_curve_number_excess(p_series, cn=cn_eff)
        
        uh = UnitHydrograph(area_km2=sb_info['area'], tc_hours=sb_info['tc'], method='scs_curvilinear')
        q = uh.convolve(pe_inc_mm=pe_inc, dt_hours=1.0, base_flow=sb_info['base_flow'], total_hours=len(p_series) - 1)
        return q

    def execute_basin_forecast(self, continuous_rain_df: pd.DataFrame, amc: str = 'AMC_II') -> Dict[str, Any]:
        """
        Executa a propagação determinística da bacia completa com os afluentes Perimbó, Trombudo e Mirim Doce.
        """
        n_steps = len(continuous_rain_df)
        sub_q = {}
        
        for k, info in self.subbasin_params.items():
            if k in continuous_rain_df.columns:
                p_arr = continuous_rain_df[k].values
            elif 'oeste' in continuous_rain_df.columns and k == 'mirim_doce':
                p_arr = continuous_rain_df['oeste'].values
            elif 'sul' in continuous_rain_df.columns and k in ('trombudo', 'perimbo'):
                p_arr = continuous_rain_df['sul'].values
            else:
                p_arr = np.zeros(n_steps)
            sub_q[k] = self.run_subbasin_runoff(p_arr, info, amc=amc)
            
        # 1. Aplicar Barragens
        q_oeste_dam = self._apply_dam(sub_q['oeste'], self.dam_operations['oeste'])
        q_sul_dam = self._apply_dam(sub_q['sul'], self.dam_operations['sul'])
        q_norte_dam = self._apply_dam(sub_q['norte'], self.dam_operations['norte'])
        
        # 2. Propagação dos Tributários Independentes (Muskingum)
        m_oeste = MuskingumReach('OESTE_RS', 'Taió a Rio do Sul', k_hours=8.0, x_param=0.20)
        m_mirim_doce = MuskingumReach('MIRIM_DOCE', 'Mirim Doce a Rio do Oeste', k_hours=4.0, x_param=0.20)
        m_sul = MuskingumReach('SUL_RS', 'Ituporanga a Rio do Sul', k_hours=6.0, x_param=0.20)
        m_perimbo = MuskingumReach('PERIMBO', 'Petrolândia a Ituporanga', k_hours=3.5, x_param=0.20)
        m_trombudo = MuskingumReach('TROMBUDO', 'Agrolândia a Rio do Sul', k_hours=5.0, x_param=0.20)
        m_norte = MuskingumReach('NORTE_IB', 'Boiteux a Ibirama', k_hours=10.0, x_param=0.20)
        m_benedito = MuskingumReach('BEN_IND', 'Timbó a Indaial', k_hours=4.0, x_param=0.20)
        m_mirim_bq = MuskingumReach('MIRIM_BQ', 'Vidal a Brusque', k_hours=10.0, x_param=0.25)
        m_mirim_foz = MuskingumReach('MIRIM_CANAL', 'Brusque a Canal Itajaí', k_hours=4.0, x_param=0.25)
        m_luis = MuskingumReach('LUIS_FOZ', 'Luís Alves a Foz', k_hours=3.0, x_param=0.20)
        
        q_oeste_jus = m_oeste.route(q_oeste_dam) + m_mirim_doce.route(sub_q['mirim_doce'])
        q_sul_jus = m_sul.route(q_sul_dam) + m_perimbo.route(sub_q['perimbo']) + m_trombudo.route(sub_q['trombudo'])
        q_norte_jus = m_norte.route(q_norte_dam)
        q_benedito_jus = m_benedito.route(sub_q['benedito'])
        q_brusque = m_mirim_bq.route(sub_q['mirim'])
        q_mirim_foz = m_mirim_foz.route(q_brusque)
        q_luis_jus = m_luis.route(sub_q['luis_alves'])
        
        # 3. Confluências Principais
        # Ponto 1: Rio do Sul (Encontro de Oeste + Mirim Doce + Sul + Perimbó + Trombudo)
        q_rio_sul = q_oeste_jus + q_sul_jus
        
        # Ponto 2: Ibirama (Rio do Sul propagado até Ibirama + Hercílio)
        m_rs_ibirama = MuskingumReach('RS_IBIRAMA', 'Rio do Sul a Ibirama', k_hours=4.0, x_param=0.20)
        q_rs_em_ibirama = m_rs_ibirama.route(q_rio_sul)
        q_ibirama = q_rs_em_ibirama + q_norte_jus
        
        # Ponto 3: Indaial (Ibirama propagado + Benedito)
        m_ib_indaial = MuskingumReach('IB_INDAIAL', 'Ibirama a Indaial', k_hours=7.0, x_param=0.22)
        q_ib_em_indaial = m_ib_indaial.route(q_ibirama)
        q_indaial = q_ib_em_indaial + q_benedito_jus
        
        # Ponto 4: Blumenau Centro (Indaial a Blumenau)
        m_ind_blu = MuskingumReach('IND_BLU', 'Indaial a Blumenau', k_hours=3.0, x_param=0.25)
        q_blumenau = m_ind_blu.route(q_indaial)
        
        # Ponto 5: Itajaí Foz (Blumenau a Itajaí + Canal Mirim + Luís Alves)
        m_blu_it = MuskingumReach('BLU_ITAJAI', 'Blumenau a Itajaí Foz', k_hours=8.0, x_param=0.30)
        q_blu_em_it = m_blu_it.route(q_blumenau)
        q_itajai_foz = q_blu_em_it + q_mirim_foz + q_luis_jus * 0.7
        
        return {
            'rio_do_sul': q_rio_sul,
            'ibirama': q_ibirama,
            'indaial': q_indaial,
            'blumenau': q_blumenau,
            'brusque': q_brusque,
            'canal_mirim': q_mirim_foz,
            'itajai_foz': q_itajai_foz
        }

    def generate_full_forecast_package(self, scenarios_dict: Dict[str, pd.DataFrame], amc: str = 'AMC_II') -> Dict[str, Any]:
        """
        Executa os cenários Low, Mean e High e gera sumário operacional com métricas de crista e alerta.
        """
        now_ts = scenarios_dict['now_timestamp']
        horizon_h = scenarios_dict['horizon_hours']
        timestamps = scenarios_dict['mean']['timestamp'].values
        
        now_idx = np.where(scenarios_dict['mean']['timestamp'] == now_ts)[0]
        t_now_idx = int(now_idx[0]) if len(now_idx) > 0 else 0
        
        res_mean = self.execute_basin_forecast(scenarios_dict['mean'], amc=amc)
        res_low = self.execute_basin_forecast(scenarios_dict['low'], amc=amc)
        res_high = self.execute_basin_forecast(scenarios_dict['high'], amc=amc)
        
        station_names = {
            'rio_do_sul': '1. Rio do Sul (Confluência Alto Vale)',
            'ibirama': '2. Ibirama (Confluência Rio Hercílio)',
            'indaial': '3. Indaial / Timbó (Confluência Benedito)',
            'blumenau': '4. Blumenau Centro (Médio Vale)',
            'brusque': '5. Brusque Centro (Rio Itajaí-Mirim)',
            'itajai_foz': '6. Itajaí (Foz Geral & Canal Retificado)'
        }
        
        stations_summary = {}
        forecast_curves = {}
        
        for st_key, st_label in station_names.items():
            q_m = res_mean[st_key]
            q_l = res_low[st_key]
            q_h = res_high[st_key]
            
            q_now = float(q_m[t_now_idx])
            
            future_slice = q_m[t_now_idx:]
            if len(future_slice) > 0:
                rel_peak_idx = int(np.argmax(future_slice))
                peak_idx = t_now_idx + rel_peak_idx
                q_peak_mean = float(q_m[peak_idx])
                q_peak_low = float(q_l[peak_idx])
                q_peak_high = float(q_h[peak_idx])
                hours_to_peak = float(rel_peak_idx)
            else:
                peak_idx = t_now_idx
                q_peak_mean = q_now
                q_peak_low = float(q_l[t_now_idx])
                q_peak_high = float(q_h[t_now_idx])
                hours_to_peak = 0.0
                
            if hours_to_peak > 1.0 and (q_peak_mean - q_now) > 50.0:
                trend = "SUBINDO ⬆️"
            elif hours_to_peak <= 1.0 and q_now >= (q_peak_mean * 0.95):
                trend = "NO TOPO 🔴"
            else:
                trend = "DESCENDO / RECESSÃO ⬇️"
                
            stations_summary[st_key] = {
                'name': st_label,
                'q_now_m3s': round(q_now, 1),
                'q_peak_mean_m3s': round(q_peak_mean, 1),
                'q_peak_low_m3s': round(q_peak_low, 1),
                'q_peak_high_m3s': round(q_peak_high, 1),
                'hours_to_peak': round(hours_to_peak, 1),
                'peak_timestamp': str(timestamps[peak_idx])[:16],
                'trend': trend
            }
            
            forecast_curves[st_key] = {
                'timestamps': [str(t)[:16] for t in timestamps],
                'q_mean': np.round(q_m, 1).tolist(),
                'q_low': np.round(q_l, 1).tolist(),
                'q_high': np.round(q_h, 1).tolist()
            }
            
        return {
            'now_timestamp': str(now_ts)[:16],
            't_now_idx': t_now_idx,
            'horizon_hours': horizon_h,
            'amc_condition': amc,
            'stations_summary': stations_summary,
            'forecast_curves': forecast_curves
        }
