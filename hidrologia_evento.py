# -*- coding: utf-8 -*-
"""
Chuva REAL observada -> hidrogramas de entrada para o HEC-RAS.

Faz a metade hidrologica da cadeia; a hidraulica fica toda no HEC-RAS:

    chuva horaria observada
      -> perda por SCS-CN (sobre o hietograma REAL)
      -> hidrograma unitario do SCS (NEH-4)
      -> convolucao
      -> amortecimento nas barragens (Modified Puls)   [opcional]
      -> Q(t) por sub-bacia  ->  .u01

Diferenca para hydrology_engine.scs_effective_rain(): aquela funcao recebe a
chuva TOTAL e a distribui por blocos alternados sinteticos de 24 h. Aqui a
perda e aplicada sobre a serie horaria observada, preservando a forma real do
evento (que e o que diferencia 1983 de 2008 de 2011 de 2023).

O roteamento Muskingum do hydrology_engine NAO e usado de proposito: quem
propaga a onda passa a ser o solver do HEC-RAS.

Uso:
    from hidrologia_evento import hidrogramas
    q, horas = hidrogramas("2008", barragens=True)
"""
import csv
import os

import numpy as np

from hydrology_engine import get_scs_unit_hydrograph, reservoir_puls_routing

DIR_CHUVA = os.path.join("itajai_flood_model", "data", "rainfall_events")

# Parametros das sub-bacias (hydrology_engine.simulate_itajai_basin) mais a
# parcela incremental do proprio Acu, que la nao existe.
#   col        = coluna do CSV de chuva
#   dam_hm3    = capacidade da barragem de contencao (0 = sem barragem)
# Rendimento especifico usado para a vazao de base quando qbase = None.
# O engine original usava 15-40 m3/s fixos por sub-bacia, o que deixa o Acu
# com 33 m3/s para 5.033 km2 -- lamina fina demais numa secao de 1.400 m, e o
# solver instabiliza nas dezenas de horas antes de a cheia chegar.
RENDIMENTO_BASE = 0.022        # m3/s por km2 (22 L/s.km2)

SUBBACIAS = {
    "sul":      {"area": 2280.0, "tc": 12.8, "cn": 78, "col": "sul",
                 "dam_hm3": 93.5,  "qbase": None},
    "oeste":    {"area": 3120.0, "tc": 14.2, "cn": 76, "col": "oeste",
                 "dam_hm3": 83.0,  "qbase": None},
    "norte":    {"area": 3450.0, "tc": 16.5, "cn": 75, "col": "norte",
                 "dam_hm3": 357.0, "qbase": None},
    "benedito": {"area": 1540.0, "tc": 10.4, "cn": 77, "col": "benedito",
                 "dam_hm3": 0.0,   "qbase": None},
    "mirim":    {"area": 1680.0, "tc": 11.8, "cn": 79, "col": "mirim",
                 "dam_hm3": 0.0,   "qbase": None},
    # incremental do Acu: 14.871 km2 da bacia menos os 11.557 km2 nomeados
    "acu_incr": {"area": 3314.0, "tc": 12.0, "cn": 76, "col": "acu",
                 "dam_hm3": 0.0,   "qbase": None},
}


def ler_chuva(evento):
    """Serie horaria de chuva por sub-bacia (mm)."""
    caminho = os.path.join(DIR_CHUVA, f"chuva_real_{evento}.csv")
    if not os.path.exists(caminho):
        raise FileNotFoundError(caminho)
    with open(caminho, newline="", encoding="utf-8") as f:
        linhas = list(csv.DictReader(f))
    cols = [c for c in linhas[0] if c != "timestamp"]
    serie = {c: np.array([float(l[c]) for l in linhas]) for c in cols}
    return serie, len(linhas)


def chuva_efetiva_scs(p_horaria, cn):
    """Perda SCS-CN aplicada ao hietograma REAL (nao a um bloco sintetico).

    S  = 25400/CN - 254           (mm)
    Ia = 0,2 S
    Pe = (P - Ia)^2 / (P - Ia + S) para P acumulado > Ia
    """
    s_mm = (25400.0 / cn) - 254.0
    ia_mm = 0.2 * s_mm
    p_acum = np.cumsum(np.asarray(p_horaria, dtype=float))
    excedente = p_acum - ia_mm
    pe_acum = np.where(excedente > 0.0,
                       (excedente ** 2) / (excedente + s_mm), 0.0)
    pe_inc = np.diff(np.insert(pe_acum, 0, 0.0))
    return np.clip(pe_inc, 0.0, None)


def hidrogramas(evento, barragens=True, horas=None):
    """Q(t) por sub-bacia, em m3/s, para o evento pedido.

    barragens=False reproduz o cenario 'sem obras' (comportas totalmente
    abertas), que e a comparacao usada para medir o efeito das 3 barragens.
    """
    serie, n = ler_chuva(evento)
    n = horas or n
    saida = {}
    for chave, sb in SUBBACIAS.items():
        col = sb["col"]
        if col not in serie:
            continue
        p = serie[col][:n]
        if len(p) < n:                       # completa com zeros se preciso
            p = np.pad(p, (0, n - len(p)))
        pe = chuva_efetiva_scs(p, sb["cn"])
        _, u = get_scs_unit_hydrograph(sb["area"], sb["tc"],
                                       dt_hours=1.0, total_hours=n)
        qbase = sb["qbase"]
        if qbase is None:
            qbase = RENDIMENTO_BASE * sb["area"]
        q = np.convolve(pe, u)[:n] + qbase
        if barragens and sb["dam_hm3"] > 0:
            q = reservoir_puls_routing(q, sb["dam_hm3"], True, dt_hours=1.0)
            # BUG do hydrology_engine.reservoir_puls_routing: ele cria
            # q_out = np.zeros(n) e so preenche a partir do indice 1, entao
            # q_out[0] fica ZERO. Como a cabeceira do Acu e Sul+Oeste (ambos
            # com barragem), o HEC-RAS partia com vazao nula e instabilizava
            # ja no primeiro minuto. Restaura a continuidade no instante 0.
            if q[0] <= 0.0 and len(q) > 1:
                q[0] = q[1]
        saida[chave] = q
    return saida, n


if __name__ == "__main__":
    for ev in ("1983", "2008", "2011", "2023"):
        try:
            com, n = hidrogramas(ev, barragens=True)
            sem, _ = hidrogramas(ev, barragens=False)
        except FileNotFoundError:
            print(f"{ev}: sem dados de chuva"); continue
        print(f"\n=== {ev} ({n} h) ===")
        print(f"{'sub-bacia':<12}{'Q pico c/ barragem':>20}{'sem barragem':>15}"
              f"{'amortecimento':>15}")
        for k in com:
            a, b = com[k].max(), sem[k].max()
            red = f"{100*(1-a/b):5.1f}%" if b > 0 and abs(a - b) > 1e-6 else "   --"
            print(f"{k:<12}{a:>20.1f}{b:>15.1f}{red:>15}")
        tot_c = sum(com[k].max() for k in com)
        tot_s = sum(sem[k].max() for k in sem)
        print(f"{'SOMA':<12}{tot_c:>20.1f}{tot_s:>15.1f}"
              f"{100*(1-tot_c/tot_s):>14.1f}%")
