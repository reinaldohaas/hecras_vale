"""
Pacote rating_curve:
Módulo de Curva-Chave (Rating Curve) Q-H para a Bacia do Rio Itajaí.
"""

from .base import BaseRatingCurve, CurveType
from .official import (
    SegmentedPowerCurve,
    create_blumenau_official_curve,
    create_rio_do_sul_official_curve,
    create_brusque_official_curve,
    create_indaial_official_curve
)
from .hydraulic import CrossSectionGeometry, HydraulicRatingCurve
from .manager import RatingCurveManager

__all__ = [
    'BaseRatingCurve',
    'CurveType',
    'SegmentedPowerCurve',
    'create_blumenau_official_curve',
    'create_rio_do_sul_official_curve',
    'create_brusque_official_curve',
    'create_indaial_official_curve',
    'CrossSectionGeometry',
    'HydraulicRatingCurve',
    'RatingCurveManager'
]
