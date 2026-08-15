"""
Pacote rainfall:
Módulo de gestão, download de chuva real, interpolação espacial e desagregação temporal.
"""

from .provider import RainfallProvider
from .spatial import SpatialRainfallInterpolator
from .antecedent_moisture import AntecedentMoistureCondition
from .forecast import RainfallForecastTimeline
from .disaggregation import RainfallDisaggregator, DisaggregationMethod
from .downloader import RainfallDownloader

__all__ = [
    'RainfallProvider',
    'SpatialRainfallInterpolator',
    'AntecedentMoistureCondition',
    'RainfallForecastTimeline',
    'RainfallDisaggregator',
    'DisaggregationMethod',
    'RainfallDownloader'
]
