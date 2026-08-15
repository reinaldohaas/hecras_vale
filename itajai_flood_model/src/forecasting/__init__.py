"""
Pacote do Motor Operacional de Previsão de Vazões e Alertas de Cheia:
- OperationalForecastEngine: Execução integrada de previsão com cenários de incerteza
- StreamflowAssimilation: Comparação com estações fluviométricas e correção transparente
- FloodAlertSystem: Classificação de alerta e emergência configurável
"""

from .engine import OperationalForecastEngine
from .assimilation import StreamflowAssimilation
from .alerts import FloodAlertSystem

__all__ = [
    'OperationalForecastEngine',
    'StreamflowAssimilation',
    'FloodAlertSystem'
]
