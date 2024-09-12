from __future__ import annotations

from dataclasses import dataclass

from ..utils.config_sort import read_omegaconfig

from ..update_physical_values import (
    update_enthalpy,
    update_enthalpy_solid_state,
)
from ..parameters.constants import Constants
from ..parameters.results_params import ResultsParams
from ..parameters.user_input import UserInput
from .geometry_settings import GeometrySettings
from .initial_boundary_conditions import set_initial_conditions


@dataclass
class PreprocessData:
    """Class representing the preprocessing of data.

    Attributes:
        is_preprocessing (bool): Flag indicating if preprocessing is enabled.
    """

    is_preprocessing: bool = True


class PreProcess(UserInput, GeometrySettings, ResultsParams):
    """Class for preprocessing data before modeling."""

    def __init__(self, constants_type, config_data, output_dir):
        """Initialize the PreProcess class.
        Args:
            config_data (ConfigData): The configuration data object.
            output_dir (str): The output directory.
        Raises:
            None
        Returns:
            None
        """

        if constants_type == "real":
            UserInput.__init__(
                self,
                constants=Constants.REAL.value,
                config_data=config_data,
                dir_output_name=output_dir,
            )
            print("User Configuration Data Setup Complete...")
        else:
            super(UserInput, self).__init__(
                Constants.DEBUG.value, config_data=config_data
            )

        GeometrySettings.__init__(
            self,
            geom=self.geometry_type,
            dz=self.grid_resolution_dz,
        )
        print("Geometry Data Setup Complete...")
        ResultsParams.__init__(self, iter_max=self.max_iterations, nz=self.nz)
        print("Results Data Setup Complete...")

    def set_preprocess(self, config_data, output_dir):
        pass

    def preprocess(self):
        """Preprocesses the data before running the simulation.
        This method sets up the initial conditions and boundary conditions for the simulation.
        It calculates the solid enthalpy and updates the enthalpy based on temperature, salinity, and liquid fraction.
        Finally, it initializes the ice thickness and prints a message indicating that the initial and boundary conditions have been applied.
            Args:
                None
        """

        self.time_passed = set_up_iter(self.max_iterations, self.grid_timestep_dt)
        [
            self.temperature,
            self.salinity,
            self.liquid_fraction,
            self.upwind_velocity,
        ] = set_initial_conditions(
            self.nz,
            self.boundary_salinity,
            self.initial_temperature,
            self.initial_salinity,
            self.initial_liquid_fraction,
            self.boundary_top_temperature,
        )
        self.solid_enthalpy = update_enthalpy_solid_state(
            self.salinity,
            self.nz,
            self.liquidus_relation_type,
            self.temperature_melt,
        )
        self.enthalpy = update_enthalpy(
            self.temperature, self.salinity, self.liquid_fraction, self.nz
        )
        self.ice_thickness, _ = 0, 0

        print("Applied Initial & Boundary Conditions...")

    @classmethod
    def get_variables(
        cls, config, out_dir_final: str
    ) -> tuple[PreprocessData, UserInput]:
        """Retrieves variables and user input data after preprocessing.

        Args:
            cls: The class object.
            config: The configuration object.
            out_dir_final: The output directory.
        Returns:
            A tuple containing the filtered variables and user input data.
        """

        print("Preprocessing...")
        constants_type = read_omegaconfig(config, "constants")
        preprocess_obj = cls(constants_type, config, out_dir_final)
        preprocess_obj.preprocess()
        filtered_vars = dict(vars(preprocess_obj))
        userinput_data = preprocess_obj.get_userinput()
        print("Preprocessing done.")
        return cls.set_dataclass(filtered_vars, PreprocessData), userinput_data

    def get_userinput(self):
        """Returns a UserInput object with the following attributes:

        Args:
            self (PreProcess): The PreProcess instance.
        Returns:
            UserInput: A UserInput object with the following attributes:
                - constants: The constants attribute.
                - grid_timestep_dt: The grid_timestep_dt attribute.
                - initial_salinity: The initial_salinity attribute.
                - dir_output_name: The dir_output_name attribute.
                - max_iterations: The max_iterations attribute.
        """

        return UserInput(
            constants=self.constants,
            grid_timestep_dt=self.grid_timestep_dt,
            initial_salinity=self.initial_salinity,
            dir_output_name=self.dir_output_name,
            max_iterations=self.max_iterations,
        )

    @staticmethod
    def set_dataclass(data_to_be_converted: dict, dataclass: dataclass) -> dataclass:
        """Sets the attributes of a dataclass object using a dictionary.

        Args:
            data_to_be_converted (dict): A dictionary containing the attribute names and values to be set.
            dataclass (dataclass): The dataclass object to be modified.
        Returns:
            dataclass: The modified dataclass object with the attributes set.
        """

        for key, value in data_to_be_converted.items():
            setattr(dataclass, key, value)

        return dataclass


### Computes maximum number of iterations based on start and end time and time step


def set_up_iter(iter_max, grid_timestep_dt):
    """
    Sets up the iteration parameters for the simulation.
    Args:
        iter_max (int): The maximum number of iterations.
        grid_timestep_dt (float): The time step for the grid.
    Returns:
        int: Always returns 0.
    """

    dt = grid_timestep_dt
    iter_max = iter_max
    print(f"Time step set to: {dt}s")
    return 0
