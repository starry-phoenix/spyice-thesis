from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from omegaconf import DictConfig
from src.spyice.parameters.constants import Constants
from src.spyice.parameters.real_constants import RealConstants
from src.spyice.parameters.debug_constants import DebugConstants
from src.spyice.parameters.algae_constants import nutrient_cn_dsi_ice, nutrient_cn_dsi_water, carbon_cc_ice_initial, carbon_cc_water_initial
from src.spyice.utils.config_sort import read_omegaconfig
from src.spyice.utils.create_directory import create_output_directory
from src.spyice.utils.initial_userinput import (
    calculate_initial_melt_temperature,)

# TODO: Add docstrings to the functions and classes
# TODO: DEBUG this script for new enum classes.
# TODO: const parameters and varied parameters difference show
# TODO: liquidus relationship status is always "normal" in script update physical values in function calculate_melting_temperature_from_salinity

def _dt_stability_validator(dz: float, dt: float) -> None:
    """Validates the time-step (dt) based on the Fourier stability criteria.

    Args:
        dz (float): The spatial step size.
        dt (float): The time-step to be validated.
    Raises:
        ValueError: If the time-step does not follow the Fourier stability criteria.
    Returns:
        None
    """

    fourier_criteria = int(121 * (dz / 0.01) ** 2)
    if dt > fourier_criteria:
        raise ValueError(
            "Time-step not following Fourier stability criteria, choose dt < "
            + str(fourier_criteria)
        )


def fourier_number_timestep(constants):
    """Calculates the Fourier number for the given timestep.

    Args:
        constants (object): An object containing the required constants.

    Returns:
        float: The calculated Fourier number.
    """
    return 0.5 * constants.rho_br * constants.c_br * constants.dz**2 / constants.k_br

class LiquidusRelation(str, Enum):
    """Represents the liquidus relation type.

    Attributes:
        NORMAL (str): Normal liquidus relation.
        FREZCHEM (str): Frezchem liquidus relation.
    """
    NORMAL = "Normal"
    FREZCHEM = "Frezchem"

class BoundaryConditionType(str, Enum):
    """Represents the type of boundary condition.

    Attributes:
        NEUMANN (str): Neumann boundary condition.
        DIRICHLET (str): Dirichlet boundary condition.
    """
    NEUMANN = "Neumann"
    DIRICHLET = "Dirichlet"

class InitialSalinity(str, Enum):
    """Represents the type of initial salinity. User can add more salinity option SX where X is a number.
    The salinity is in parts per thousand (ppt).

    Attributes:
        S34 (str): Initial salinity of 34 ppt.
        S0 (str): Initial salinity of 0 ppt.
        S1 (str): Initial salinity of 1 ppt.
        S2 (str): Initial salinity of 2 ppt.
        S3 (str): Initial salinity of 3 ppt.
        S_LINEAR (str): Linear initial salinity.

    """
    S34 = "S34"
    S0 = "S0"
    S1 = "S1"
    S2 = "S2"
    S3 = "S3"
    S_LINEAR = "S_linear"


class InitialTemperature(str, Enum):
    """Represents the type of initial temperature.

    Attributes:
        T_Stefan (str): Initial temperature based on Stefan condition.
        T271p25 (str): Initial temperature of 271.25 K.
        T250 (str): Initial temperature of 250 K.
        T_MELT (str): Initial temperature at which the material melts.
        T_S (str): Initial temperature based on salinity.
    """
    T_STEFAN = "T_Stefan"
    T271p25 = "T271.25"
    T250 = "T250"
    T_MELT = "Tm_w"
    T_S = "T(S)"

class InitialMeltTemperature(str, Enum):
    """Represents the type of initial melt temperature.

    Attributes:
        T_MELT (str): Initial melt temperature.
        T_S (str): Initial temperature based on salinity.
    """
    ONEPHASE = "onephase"
    FREZCHEM = "Frezchem"
    TWOPHASE = "twophase"

class InitialLiquidFraction(str, Enum):
    """Represents the type of initial liquid fraction.

    Attributes:
        P0 (str): Initial liquid fraction of 0.
        P1 (str): Initial liquid fraction of 1.
        P_Stefan (str): Initial liquid fraction based on Stefan condition.
        PX (str): Initial liquid fraction based on a custom profile.
    """
    P0 = "P0"
    P1 = "P1"
    P_Stefan = "P_Stefan"
    PX = "PX"  # where X is a number

class FileNameSuffix(str, Enum):
    """Represents the suffix for the output file name.

    Attributes:
        NON_CONST_DENS_MUSHFIX (str): Non-constant density mush-fix.
        NON_CONST_DENS (str): Non-constant density.
        CONST_DENS (str): Constant density.
    """
    NON_CONST_DENS_MUSHFIX = "nonconst_dens-mushfix"
    NON_CONST_DENS = "nonconst_dens"
    CONST_DENS = "const_dens"

class TopTemperatureType(str, Enum):
    """Represents the type of top temperature condition.

    Attributes:
        STEFAN (str): Stefan condition.
        DIRICHLET (str): Dirichlet condition.
    """
    STEFAN = "Stefan"
    DIRICHLET = "Dirichlet"

@dataclass
class UserInput:
    """Represents the user input parameters for the model.

    Attributes:
        constants (RealConstants | DebugConstants): The type of constants to use.
        max_iterations (int): The maximum number of iterations.
        is_stefan (bool): Flag indicating whether Stefan condition is applied.
        is_buffo (bool): Flag indicating whether Buffo condition is applied.
        liquidus_relation_type (str): The type of liquidus relation to use.
        grid_resolution_dz (float): The grid resolution in the z-direction.
        boundary_condition_type (str): The type of boundary condition to use.
        temperature_tolerance (float): The temperature tolerance.
        salinity_tolerance (float): The salinity tolerance.
        liquid_fraction_tolerance (float): The liquid fraction tolerance.
        initial_temperature (str): The initial temperature profile.
        initial_salinity (str): The initial salinity profile.
        initial_liquid_fraction (str): The initial liquid fraction profile.
        output_suffix (str): The suffix to be added to the output files.
        temperature_top_type (str): The type of temperature condition at the top boundary.
        phase_type (int): The type of phase to consider.
        grid_timestep_dt (float): The grid timestep.
        dir_output_name (str): The name of the output directory.
        critical_liquid_fraction (float): The critical liquid fraction.
        boundary_salinity (float): The boundary salinity (automatically calculated).
        temperature_melt (float): The temperature at which the material melts (automatically calculated).
        boundary_top_temperature (float): The temperature at the top boundary (automatically calculated).
        geometry_type (int): The type of geometry.
        counter_limit (int): The counter limit.
    Methods:
        __post_init__(): Performs post-initialization tasks.
    """

    # grid_timestep_dt=config_data.time_step,
    # initial_salinity=config_data.initial_salinity,
    # dir_output_name=output_dir,
    # max_iterations=config_data.max_iterations

    # self.constants_type = self.read_omegaconfig("constants")
    # self.time_step = self.read_omegaconfig("dt")
    # self.initial_salinity = self.read_omegaconfig("S_IC")
    # self.max_iterations = self.read_omegaconfig("iter_max")
    # self.is_salinity_equation = self.read_omegaconfig("salinity")
    ...

    # --- Constants and Config ---
    constants: RealConstants | DebugConstants = Constants.REAL.value
    config_data: DictConfig = field(default_factory=dict)

    # --- Model Switches ---
    is_stefan: bool = True
    is_buffo: bool = True
    is_voller: bool = False
    is_salinity_equation: bool = False
    is_diffusiononly_equation: bool = True
    is_algae_equation: bool = False
    is_radiation_equation: bool = False

    # --- Iteration and Limits ---
    max_iterations: int = 25000
    counter_limit: int = 100000

    # --- Grid and Time Step ---
    grid_resolution_dz: float = 0.01
    grid_timestep_dt: float = 47  # in seconds

    # --- Boundary and Geometry ---
    boundary_condition_type: BoundaryConditionType = BoundaryConditionType.DIRICHLET.value
    geometry_type: int = field(init=False)
    boundary_salinity: float = field(init=False)
    boundary_top_temperature: float = field(init=False)
    temperature_top_type: TopTemperatureType = TopTemperatureType.STEFAN.value

    # --- Tolerances ---
    temperature_tolerance: float = 0.01
    salinity_tolerance: float = 0.01
    liquid_fraction_tolerance: float = 0.01

    # --- Initial Conditions ---
    initial_temperature: InitialTemperature = InitialTemperature.T_S.value
    initial_salinity: InitialSalinity = InitialSalinity.S34.value
    initial_liquid_fraction: InitialLiquidFraction = InitialLiquidFraction.P1.value
    critical_liquid_fraction: float = 0.1
    temperature_melt: float = field(init=False)

    # --- Phase and Liquidus ---
    phase_type: int = 1
    liquidus_relation_type: LiquidusRelation = LiquidusRelation.NORMAL.value

    # --- Output and Directory ---
    output_suffix: FileNameSuffix = FileNameSuffix.NON_CONST_DENS_MUSHFIX.value
    dir_output_name_hydra: str = (
        "Temperature_{S_IC}_{bc_condition}_{dz}_{dt}_{iter_max}_{cap_dens}"
    )
    dir_output_name: str = (
        "Temperature_{S_IC}_{bc_condition}_{dz}_{dt}_{iter_max}_{cap_dens}"
    )

    # --- Algae Model Parameters ---
    nutrient_cn_dsi_water: float = nutrient_cn_dsi_water
    nutrient_cn_dsi_ice: float = nutrient_cn_dsi_ice
    carbon_cc_ice_initial: float = carbon_cc_ice_initial
    carbon_cc_water_initial: float = carbon_cc_water_initial

    def __post_init__(self):
        _dt_stability_validator(self.grid_resolution_dz, self.grid_timestep_dt)

        if isinstance(self.constants, RealConstants):
            self.boundary_salinity = 34.0
            self.boundary_top_temperature = 265.0

            # melt temperature affects the liquid relation: Frezchem or Normal in src/update_physical_values.py script
            method = InitialMeltTemperature.ONEPHASE.value
            self.temperature_melt = calculate_initial_melt_temperature(self.boundary_salinity, method)

            self.geometry_type = 2
            if self.config_data:
                self.grid_timestep_dt = read_omegaconfig(self.config_data, "dt")
                self.initial_salinity = read_omegaconfig(self.config_data, "S_IC")
                self.max_iterations = read_omegaconfig(self.config_data, "iter_max")
                self.grid_resolution_dz = read_omegaconfig(self.config_data, "dz")
                self.dir_output_name = create_output_directory(
                    self.dir_output_name_hydra,
                    self.initial_salinity,
                    self.boundary_condition_type,
                    self.grid_resolution_dz,
                    self.grid_timestep_dt,
                    self.max_iterations,
                    self.output_suffix,
                )

        elif isinstance(self.constants, DebugConstants):
            self.boundary_salinity = 1.0
            self.boundary_top_temperature = -1.0
            self.temperature_melt = 0.0
            self.geometry_type = 1