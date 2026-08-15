"""
Ponte de Exportação para Modelagem Hidráulica HEC-RAS 1D/2D (HECRASBridge):
- Gera arquivos de hidrograma de contorno (Flow Hydrograph Boundary Conditions)
- Formata séries temporais Q(t) no padrão de datas e horas do HEC-RAS (DDMMMYYYY HHMM)
- Prepara condições de entrada em montante e afluentes laterais para o plano não permanente (.u01 / DSS)
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional

class HECRASBridge:
    """
    Gera arquivos e tabelas compatíveis com as condições de contorno do HEC-RAS.
    """
    @staticmethod
    def export_flow_hydrograph_table(timestamps: List[str], flow_series_m3s: np.ndarray,
                                    reach_name: str = "Rio Itajai-Acu",
                                    node_name: str = "Rio do Sul") -> pd.DataFrame:
        """
        Formata a tabela de vazão para o padrão do HEC-RAS Unsteady Flow Editor:
        Colunas: ['Date', 'Time', 'Flow_m3s', 'River', 'Reach', 'Node']
        """
        df = pd.DataFrame()
        ts_parsed = pd.to_datetime(timestamps)
        
        # Formato de data HEC-RAS: '01OCT2023'
        df['Date'] = [t.strftime('%d%b%Y').upper() for t in ts_parsed]
        # Formato de hora: '0100', '1200'
        df['Time'] = [t.strftime('%H%M') for t in ts_parsed]
        df['Flow_m3s'] = np.round(flow_series_m3s, 2)
        df['River'] = "Itajai"
        df['Reach'] = reach_name
        df['Node'] = node_name
        return df

    @staticmethod
    def save_unsteady_boundary_csv(output_csv_path: str, hydrographs_dict: Dict[str, Dict[str, Any]]):
        """
        Salva arquivo CSV multi-estações pronto para importação no HEC-RAS.
        """
        all_dfs = []
        for node_key, data in hydrographs_dict.items():
            ts = data['timestamps']
            q = np.asarray(data['q_mean'], dtype=float)
            df_node = HECRASBridge.export_flow_hydrograph_table(
                timestamps=ts,
                flow_series_m3s=q,
                reach_name="Calha Principal",
                node_name=node_key
            )
            all_dfs.append(df_node)
            
        if all_dfs:
            final_df = pd.concat(all_dfs, ignore_index=True)
            os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
            final_df.to_csv(output_csv_path, index=False)
            return output_csv_path
        return None
