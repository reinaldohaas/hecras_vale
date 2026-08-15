"""
Pacote de Gestão e Processamento de Dados Pluviométricos:
- Provedores de Chuva Observada e Prevista
- Interpolação Espacial de Chuva por Sub-bacia (Thiessen / IDW)
- Condição de Umidade Antecedente (P5 e AMC I, II, III)
- Séries Temporais Contínuas de Previsão com Cenários de Incerteza
"""

from .provider import RainfallProvider, CSVRainfallProvider, SyntheticRainfallProvider, RainfallObservation
from .spatial import SpatialRainfallInterpolator
from .antecedent_moisture import AntecedentMoistureCondition
from .forecast import RainfallForecastTimeline

__all__ = [
    'RainfallProvider',
    'CSVRainfallProvider',
    'SyntheticRainfallProvider',
    'RainfallObservation',
    'SpatialRainfallInterpolator',
    'AntecedentMoistureCondition',
    'RainfallForecastTimeline'
]
