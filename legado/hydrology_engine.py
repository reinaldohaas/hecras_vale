"""
Engine Hidrológico Completo para o Vale do Itajaí:
- Hidrograma Unitário Sintético do SCS / NRCS (NEH-4 Curvilíneo)
- Distribuição Temporal de Chuva de Projeto (Blocos Alternados / SCS Tipo II)
- Propagação Hidrológica de Cheias no Rio via Muskingum (Trecho a Trecho)
- Amortecimento em Barragens / Reservatórios via Modified Puls (Level-Pool Routing)
"""

import numpy as np
import json

# 1. Tabela Adimensional Padrão do Hidrograma Unitário Curvilíneo do SCS (NRCS NEH-4)
SCS_DIMENSIONLESS_RATIOS = [
    (0.0, 0.0), (0.1, 0.03), (0.2, 0.10), (0.3, 0.19), (0.4, 0.31), (0.5, 0.47),
    (0.6, 0.66), (0.7, 0.82), (0.8, 0.93), (0.9, 0.99), (1.0, 1.00), (1.1, 0.99),
    (1.2, 0.93), (1.3, 0.86), (1.4, 0.78), (1.5, 0.68), (1.6, 0.56), (1.7, 0.46),
    (1.8, 0.39), (1.9, 0.33), (2.0, 0.28), (2.2, 0.207), (2.4, 0.147), (2.6, 0.107),
    (2.8, 0.077), (3.0, 0.055), (3.5, 0.029), (4.0, 0.015), (4.5, 0.007), (5.0, 0.0)
]

def get_scs_unit_hydrograph(area_km2, tc_hours, dt_hours=1.0, total_hours=48):
    """
    Gera o Hidrograma Unitário do SCS para 1 mm de chuva efetiva.
    """
    # Tempo de pico (tp = dt/2 + 0.6*tc)
    tp = 0.5 * dt_hours + 0.6 * tc_hours
    tp = max(1.0, tp)
    
    # Vazão de pico unitária (qp = 2.08 * A * 1mm / tp) em m³/s/mm
    qp = (2.08 * area_km2 * 1.0) / tp
    
    # Interpolação da curva adimensional
    t_dim = [p[0] for p in SCS_DIMENSIONLESS_RATIOS]
    q_dim = [p[1] for p in SCS_DIMENSIONLESS_RATIOS]
    
    t_eval = np.arange(0, total_hours + 1, dt_hours)
    t_over_tp = t_eval / tp
    
    q_unit = np.interp(t_over_tp, t_dim, q_dim, right=0.0) * qp
    
    # Normalização para garantir conservação rigorosa de volume (Integral = Area * 1mm)
    vol_esperado_m3 = area_km2 * 1e6 * 0.001 # 1mm = 0.001m
    vol_calculado_m3 = np.sum(q_unit) * (dt_hours * 3600.0)
    if vol_calculado_m3 > 0:
        q_unit = q_unit * (vol_esperado_m3 / vol_calculado_m3)
        
    return t_eval, q_unit

def scs_effective_rain(p_total_mm, cn=75):
    """
    Calcula a chuva efetiva acumulada e incremental pelo método do Curve Number (SCS).
    """
    s_mm = (25400.0 / cn) - 254.0
    ia_mm = 0.2 * s_mm
    
    # Distribuição temporal da chuva em 24h pelo método dos Blocos Alternados (hietograma)
    t_24 = np.arange(1, 25)
    p_acum_rel = np.where(t_24 <= 12, 0.5 * (t_24 / 12.0)**2, 1.0 - 0.5 * ((24.0 - t_24) / 12.0)**2)
    p_acum = p_total_mm * p_acum_rel
    
    # Chuva efetiva acumulada
    pe_acum = np.where(p_acum > ia_mm, ((p_acum - ia_mm)**2) / (p_acum - ia_mm + s_mm), 0.0)
    
    # Chuva efetiva incremental horária
    pe_inc = np.diff(np.insert(pe_acum, 0, 0.0))
    
    # Expandir para 48h (zero após 24h)
    pe_48 = np.zeros(49)
    pe_48[0:24] = pe_inc
    return pe_48

def convolve_hydrograph(pe_inc, q_unit, base_flow=15.0):
    """
    Convolução discreta: Q(t) = sum( Pe[k] * U[t - k] ) + Qbase
    """
    q_dir = np.convolve(pe_inc, q_unit)[:len(pe_inc)]
    return q_dir + base_flow

def muskingum_routing(q_in, k_hours, x_weight=0.20, dt_hours=1.0):
    """
    Propagação de onda de cheia em canal fluvial pelo método de Muskingum:
    O(t+dt) = C0*I(t+dt) + C1*I(t) + C2*O(t)
    """
    denom = 2.0 * k_hours * (1.0 - x_weight) + dt_hours
    c0 = (dt_hours - 2.0 * k_hours * x_weight) / denom
    c1 = (dt_hours + 2.0 * k_hours * x_weight) / denom
    c2 = (2.0 * k_hours * (1.0 - x_weight) - dt_hours) / denom
    
    n = len(q_in)
    q_out = np.zeros(n)
    q_out[0] = q_in[0]
    
    for t in range(0, n - 1):
        q_val = c0 * q_in[t + 1] + c1 * q_in[t] + c2 * q_out[t]
        q_out[t + 1] = max(0.0, q_val)
        
    return q_out

def reservoir_puls_routing(q_in, dam_capacity_hm3, dam_active=True, dt_hours=1.0):
    """
    Amortecimento em reservatório de contenção de cheias (Modified Puls Method).
    """
    if not dam_active:
        return q_in.copy()
        
    max_cap_m3 = dam_capacity_hm3 * 1e6
    dt_sec = dt_hours * 3600.0
    
    n = len(q_in)
    q_out = np.zeros(n)
    storage = np.zeros(n)
    
    # Vazão controlada máxima das comportas de fundo
    q_gate_max = 50.0 # m³/s
    
    for t in range(0, n - 1):
        i_avg = 0.5 * (q_in[t] + q_in[t + 1])
        if storage[t] < max_cap_m3:
            q_out[t + 1] = min(q_gate_max, q_in[t + 1])
        else:
            q_out[t + 1] = q_in[t + 1]
            
        o_avg = 0.5 * (q_out[t] + q_out[t + 1])
        storage[t + 1] = max(0.0, min(max_cap_m3 * 1.2, storage[t] + (i_avg - o_avg) * dt_sec))
        
    return q_out

def simulate_itajai_basin(p_oeste=120, p_sul=130, p_norte=140, p_benedito=110, p_mirim=115,
                          dam_oeste=True, dam_sul=True, dam_norte=True, split_canal_pct=70):
    """
    Simulação hidrológica completa da Bacia do Rio Itajaí com propagação Muskingum e barragens.
    """
    subbasins = {
        'oeste': {'area': 3120.0, 'tc': 14.2, 'p': p_oeste, 'cn': 76, 'dam_cap': 83.0, 'dam_on': dam_oeste},
        'sul': {'area': 2280.0, 'tc': 12.8, 'p': p_sul, 'cn': 78, 'dam_cap': 93.5, 'dam_on': dam_sul},
        'norte': {'area': 3450.0, 'tc': 16.5, 'p': p_norte, 'cn': 75, 'dam_cap': 357.0, 'dam_on': dam_norte},
        'benedito': {'area': 1540.0, 'tc': 10.4, 'p': p_benedito, 'cn': 77, 'dam_cap': 0.0, 'dam_on': False},
        'mirim': {'area': 1680.0, 'tc': 11.8, 'p': p_mirim, 'cn': 79, 'dam_cap': 0.0, 'dam_on': False}
    }
    
    q_sub = {}
    for k, sb in subbasins.items():
        _, q_unit = get_scs_unit_hydrograph(sb['area'], sb['tc'])
        pe = scs_effective_rain(sb['p'], sb['cn'])
        q_raw = convolve_hydrograph(pe, q_unit, base_flow=15.0)
        
        if sb['dam_cap'] > 0:
            q_damped = reservoir_puls_routing(q_raw, sb['dam_cap'], sb['dam_on'])
        else:
            q_damped = q_raw
        q_sub[k] = q_damped
        
    # Trecho 1: Taió -> Rio do Sul (K = 8h, X = 0.20)
    q_oeste_riosul = muskingum_routing(q_sub['oeste'], k_hours=8.0, x_weight=0.20)
    
    # Trecho 2: Ituporanga -> Rio do Sul (K = 6h, X = 0.20)
    q_sul_riosul = muskingum_routing(q_sub['sul'], k_hours=6.0, x_weight=0.20)
    
    # 📍 CIDADE 4: RIO DO SUL (Confluência)
    q_rio_do_sul = q_oeste_riosul + q_sul_riosul
    
    # Trecho 3: Rio do Sul -> Blumenau (K = 14h, X = 0.25)
    q_alto_vale_blumenau = muskingum_routing(q_rio_do_sul, k_hours=14.0, x_weight=0.25)
    
    # Trecho 4: José Boiteux (Rio Hercílio) -> Blumenau (K = 10h, X = 0.20)
    q_norte_blumenau = muskingum_routing(q_sub['norte'], k_hours=10.0, x_weight=0.20)
    
    # Trecho 5: Rio Benedito -> Blumenau (K = 4h, X = 0.20)
    q_benedito_blumenau = muskingum_routing(q_sub['benedito'], k_hours=4.0, x_weight=0.20)
    
    # 📍 CIDADE 3: BLUMENAU (Médio Vale)
    q_blumenau = q_alto_vale_blumenau + q_norte_blumenau + q_benedito_blumenau
    
    # 📍 CIDADE 2: BRUSQUE (Rio Itajaí-Mirim)
    # Trecho 6: Alto Mirim -> Brusque (K = 5h, X = 0.20)
    q_brusque = muskingum_routing(q_sub['mirim'], k_hours=5.0, x_weight=0.20)
    
    # Trecho 7: Blumenau -> Itajaí (K = 8h, X = 0.30)
    q_acu_foz = muskingum_routing(q_blumenau, k_hours=8.0, x_weight=0.30)
    
    # Trecho 8: Brusque -> Itajaí (K = 4h, X = 0.25)
    q_mirim_foz = muskingum_routing(q_brusque, k_hours=4.0, x_weight=0.25)
    
    f_canal = split_canal_pct / 100.0
    f_braco = 1.0 - f_canal
    q_mirim_canal = q_mirim_foz * f_canal
    q_mirim_braco = q_mirim_foz * f_braco
    
    # 📍 CIDADE 1: ITAJAÍ FOZ (Deságue Total)
    q_itajai = q_acu_foz + q_mirim_canal + q_mirim_braco
    
    results = {
        'horas': list(range(49)),
        'q_rio_do_sul': [round(float(v), 1) for v in q_rio_do_sul],
        'q_blumenau': [round(float(v), 1) for v in q_blumenau],
        'q_brusque': [round(float(v), 1) for v in q_brusque],
        'q_itajai': [round(float(v), 1) for v in q_itajai],
        'q_mirim_canal': [round(float(v), 1) for v in q_mirim_canal],
        'q_mirim_braco': [round(float(v), 1) for v in q_mirim_braco],
        'q_oeste_afluente': [round(float(v), 1) for v in q_sub['oeste']],
        'q_sul_afluente': [round(float(v), 1) for v in q_sub['sul']],
        'q_norte_afluente': [round(float(v), 1) for v in q_sub['norte']]
    }
    return results

if __name__ == '__main__':
    res = simulate_itajai_basin()
    print('Simulacao Hidrologica SCS + Muskingum concluida com sucesso!')
    p_rs = max(res['q_rio_do_sul'])
    t_rs = res['q_rio_do_sul'].index(p_rs)
    p_bl = max(res['q_blumenau'])
    t_bl = res['q_blumenau'].index(p_bl)
    p_br = max(res['q_brusque'])
    t_br = res['q_brusque'].index(p_br)
    p_ij = max(res['q_itajai'])
    t_ij = res['q_itajai'].index(p_ij)
    print(f'  - Rio do Sul: {p_rs} m3/s no t={t_rs}h')
    print(f'  - Blumenau:   {p_bl} m3/s no t={t_bl}h')
    print(f'  - Brusque:    {p_br} m3/s no t={t_br}h')
    print(f'  - Itajai Foz: {p_ij} m3/s no t={t_ij}h')
