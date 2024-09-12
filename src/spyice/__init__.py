# SPDX-License-Identifier: Apache-2.0
from .main_process import MainProcess
from .utils import SpyiceLogger
from .update_physical_values import (
    update_state_variables,
    update_enthalpy,
    update_enthalpy_solid_state,
)
from .rhs import apply_boundary_condition, correct_for_brine_movement
