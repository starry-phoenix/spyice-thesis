from enum import Enum

from .debug_constants import DebugConstants
from .real_constants import RealConstants


class Constants(Enum):
    """Enumeration class for constants."""

    REAL = RealConstants()
    DEBUG = DebugConstants()


# Example usage
# my_constant = Constants.REAL.value
# rho_br_val = my_constant.rho_br
# print(rho_br_val)
