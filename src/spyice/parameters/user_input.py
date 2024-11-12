from __future__ import annotations

from dataclasses import dataclass, field

from omegaconf import DictConfig
from src.spyice.parameters.constants import Constants
from src.spyice.parameters.real_constants import RealConstants
from src.spyice.parameters.debug_constants import DebugConstants
from src.spyice.utils.config_sort import read_omegaconfig
from src.spyice.utils.create_directory import create_output_directory


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

    fourier_criteria = int(50 * (dz / 0.01) ** 2)
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

    constants: RealConstants | DebugConstants = Constants.REAL.value
    config_data: DictConfig = field(default_factory=dict)
    max_iterations: int = 500
    is_stefan: bool = True
    is_buffo: bool = True
    is_voller: bool = False
    is_salinity_equation: bool = True
    liquidus_relation_type: str = "Normal"  # Normal or Frezchem
    grid_resolution_dz: float = 0.01
    boundary_condition_type: str = "Dirichlet"  # Neumann or Dirichlet
    temperature_tolerance: float = 0.01
    salinity_tolerance: float = 0.01
    liquid_fraction_tolerance: float = 0.001
    initial_temperature: str = "T(S)"  # "T_Stefan" or  "T271.25" or "T250" or "Tm_w"
    initial_salinity: str = (
        "S34"  # "S_linear" or "S34" or "S0" or "SX" where X is a number
    )
    initial_liquid_fraction: str = (
        "P1"  # "P_Stefan" or "P0" or "P1" or "PX" where X is a number
    )
    output_suffix: str = "nonconst_dens-mushfix"
    temperature_top_type: str = "Stefan"  # "Stefan" or "Dirichlet"
    phase_type: int = 1
    grid_timestep_dt: float = 10
    dir_output_name_hydra: str = (
        "Temperature_{S_IC}_{bc_condition}_{dz}_{dt}_{iter_max}_{cap_dens}"
    )
    dir_output_name: str = (
        "Temperature_{S_IC}_{bc_condition}_{dz}_{dt}_{iter_max}_{cap_dens}"
    )
    critical_liquid_fraction: float = 0.1

    boundary_salinity: float = field(init=False)
    temperature_melt: float = field(init=False)
    boundary_top_temperature: float = field(init=False)
    geometry_type: int = field(init=False)
    counter_limit: int = 100000

    def __post_init__(self):
        _dt_stability_validator(self.grid_resolution_dz, self.grid_timestep_dt)

        if isinstance(self.constants, RealConstants):
            self.boundary_salinity = 34.0
            self.boundary_top_temperature = 265.0

            # melt temperature affects the liquid relation: Frezchem or Normal in src/update_physical_values.py script
            # self.temperature_melt = 273.15 - 1.853 * self.boundary_salinity / 28.0 
            T_m = 273.15 # melt temperature as solidus temperature
            T_s = 252.05 # eutectic temperature for Sbr = 233ppt
            S_br = 233.0  # brine salinity in ppt
            self.temperature_melt = 273.15 + (T_s - T_m)*self.boundary_salinity/S_br
            # self.temperature_melt = (
            #     -(9.1969758 * (1e-05) * self.boundary_salinity**2)
            #     - 0.03942059 * self.boundary_salinity
            #     + 272.63617665
            # )

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

                # self.is_salinity_equation = read_omegaconfig(
                #     self.config_data, "salinity"
                # )

        elif isinstance(self.constants, DebugConstants):
            self.boundary_salinity = 0.0
            self.boundary_top_temperature = -1.0
            self.temperature_melt = 0.0
            self.geometry_type = 1
