"""
Pacote inundation:
Módulo de cálculo de lâmina d'água 2D, filtro de conectividade hidráulica
e mapeamento vetorial de manchas de inundação no Vale do Itajaí.
"""

from .flood_grid import InundationGrid, HydraulicConnectivityFilter
from .depth_raster import DepthClassifier, FloodDepthClass
from .flood_mapper import FloodplainMapper
from .cross_sections import RiverCrossSection, CrossSectionDelineator

__all__ = [
    'InundationGrid',
    'HydraulicConnectivityFilter',
    'DepthClassifier',
    'FloodDepthClass',
    'FloodplainMapper',
    'RiverCrossSection',
    'CrossSectionDelineator'
]
