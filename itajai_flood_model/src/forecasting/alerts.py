"""
Módulo de Classificação de Níveis de Alerta de Cheia (FloodAlertSystem):
- Carrega limiares configuráveis de arquivo externo (thresholds.json)
- Classifica o estado hidrológico: NORMAL | ATENÇÃO | ALERTA | EMERGÊNCIA
- Fornece cores e mensagens padronizadas de defesa civil
"""

import os
import json
from typing import Dict, Any, Optional

DEFAULT_THRESHOLDS_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'thresholds.json')
)

class FloodAlertSystem:
    """
    Classificador de risco e níveis de emergência.
    """
    def __init__(self, thresholds_path: Optional[str] = None):
        self.thresholds_path = thresholds_path or DEFAULT_THRESHOLDS_PATH
        self.thresholds = self._load_thresholds()

    def _load_thresholds(self) -> Dict[str, Any]:
        if os.path.exists(self.thresholds_path):
            try:
                with open(self.thresholds_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    return data.get('stations', {})
            except Exception as e:
                print(f"Aviso: Erro ao carregar {self.thresholds_path}: {e}")
        return {}

    def classify_flow(self, station_key: str, flow_m3s: float) -> Dict[str, Any]:
        """
        Retorna a classificação do nível de alerta para a vazão fornecida.
        """
        st_cfg = self.thresholds.get(station_key, {})
        q_atencao = st_cfg.get('flow_atencao_m3s', 800.0)
        q_alerta = st_cfg.get('flow_alerta_m3s', 1600.0)
        q_emergencia = st_cfg.get('flow_emergencia_m3s', 2800.0)
        
        q = float(flow_m3s)
        
        if q >= q_emergencia:
            return {
                'level': 'EMERGÊNCIA',
                'code': 'EMERGENCIA',
                'color': '#ff0055',
                'bg_color': 'rgba(255, 0, 85, 0.2)',
                'badge_class': 'badge-emergencia',
                'message': 'Cota de inundação crítica atingida na área urbana'
            }
        elif q >= q_alerta:
            return {
                'level': 'ALERTA',
                'code': 'ALERTA',
                'color': '#ef4444',
                'bg_color': 'rgba(239, 68, 68, 0.2)',
                'badge_class': 'badge-alerta',
                'message': 'Extravasamento iminente das calhas principais e arroios'
            }
        elif q >= q_atencao:
            return {
                'level': 'ATENÇÃO',
                'code': 'ATENCAO',
                'color': '#f59e0b',
                'bg_color': 'rgba(245, 158, 11, 0.2)',
                'badge_class': 'badge-atencao',
                'message': 'Elevação acelerada dos níveis fluviais'
            }
        else:
            return {
                'level': 'NORMAL',
                'code': 'NORMAL',
                'color': '#10b981',
                'bg_color': 'rgba(16, 185, 129, 0.15)',
                'badge_class': 'badge-normal',
                'message': 'Escoamento dentro da calha natural'
            }
