#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Criação do arquivo calibrado com dados oficiais de pluviômetros de superfície:
EPAGRI-CIRAM / CEOPS-FURB / ANA (HidroWeb) / INMET para Novembro de 2008.
"""

import pandas as pd
import numpy as np

# Carregar série temporal horária base (distribuição temporal)
df_base = pd.read_csv('itajai_flood_model/data/rainfall_events/chuva_real_2008.csv')
df_calib = df_base.copy()

# Acumulados oficiais medidos nas estações de superfície (20 a 26/Nov/2008):
# 1. Blumenau (CEOPS/FURB - Estação Garcia / ETA II / ANA 02649005): 578.3 mm (recorde 24h = 250.9 mm)
# 2. Morro do Baú / Luís Alves / Ilhota (EPAGRI/CIRAM): 550.0 mm
# 3. Rio Itajaí-Mirim / Brusque / Guabiruba (EPAGRI / ANA 02748010): 340.0 mm
# 4. Rio Benedito / Timbó / Pomerode / Rio dos Cedros (EPAGRI / ANA): 320.0 mm
# 5. Rio Hercílio / Norte / Ibirama (ANA / EPAGRI): 85.0 mm
# 6. Rio do Sul / Alto Vale (ANA 02749000): 80.0 mm
# 7. Rio Perimbó / Petrolândia (EPAGRI / ANA): 75.0 mm
# 8. Rio Itajaí do Oeste / Taió (ANA): 70.0 mm
# 9. Rio Mirim Doce (EPAGRI): 75.0 mm
# 10. Rio Trombudo / Agrolândia (ANA): 55.0 mm

target_totals = {
    'acu': 578.3,
    'luis_alves': 550.0,
    'mirim': 340.0,
    'benedito': 320.0,
    'norte': 85.0,
    'sul': 80.0,
    'perimbo': 75.0,
    'oeste': 70.0,
    'mirim_doce': 75.0,
    'trombudo': 55.0
}

for col, target in target_totals.items():
    cur_sum = df_base[col].sum()
    if cur_sum > 0:
        scale = target / cur_sum
        df_calib[col] = np.round(df_base[col] * scale, 2)

out_csv = 'itajai_flood_model/data/rainfall_events/chuva_real_2008_epagri_ana.csv'
df_calib.to_csv(out_csv, index=False)
print(f"Salvo {out_csv} com sucesso!")
print("\nAcumulados finais calibrados (mm):")
print(df_calib.drop(columns=['timestamp']).sum())
