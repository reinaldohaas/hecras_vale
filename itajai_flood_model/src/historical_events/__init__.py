"""
Pacote de Eventos Históricos Reais, Replay e Calibração Automática:
- HistoricalEventsLoader: Carregador padronizado dos eventos de 1983, 2008, 2011 e 2023
- BasinAutoCalibrator: Otimizador paramétrico mono e multi-evento com separação treino/validação
"""

from .events_loader import HistoricalEventsLoader
from .auto_calibration import BasinAutoCalibrator

__all__ = [
    'HistoricalEventsLoader',
    'BasinAutoCalibrator'
]
