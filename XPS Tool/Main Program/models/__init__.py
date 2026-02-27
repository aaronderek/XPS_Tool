"""
2DEG potential well models.
"""

from .triangular import TriangularModel
from .fang_howard import FangHowardModel
from .parabolic import ParabolicModel

__all__ = ['TriangularModel', 'FangHowardModel', 'ParabolicModel']
