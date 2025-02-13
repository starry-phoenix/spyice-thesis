from __future__ import annotations

from dataclasses import asdict, dataclass

import matplotlib.pyplot as plt
import numpy as np

from src.spyice.models.algae_model import (
    biogeochemical_model,
    biogeochemical_model_at_alldepths,
)
from src.spyice.models.radiative_model import (
    calculate_radiative_terms,
    get_salinity_source_term,
)
from src.spyice.models.stefan_problem import StefanProblem
from src.spyice.parameters.results_params import ResultsParams
from src.spyice.parameters.user_input import UserInput
from src.spyice.preprocess.initial_boundary_conditions import temperature_gradient
from src.spyice.statevariables import (
    compute_error_for_convergence,
    overwrite_statevariables,
    reset_error_for_while_loop,
    set_statevariables,
)
from src.spyice.update_physical_values import update_state_variables
from src.spyice.utils.helpers import export_residuals, t_total
from src.spyice.utils.spyice_exceptions import ConvergenceError, InvalidPhaseError

plt.style.use("src.spyice.utils.custom")
# plt.rcParams.update(
#     {
#         "text.usetex": True,
#     }
# )
plt.rcParams["text.latex.preamble"].join(
    [
        r"\usepackage{dashbox}",
        r"\setmainfont{xcolor}",
    ]
)


def locate_ice_ocean_interface(phi, dz, nz, **kwargs):
    """Locate ice ocean interface, based on liquid fraction equivalent ice thickness

    Args:
        phi (array-like): Liquid fraction [-]
        dz (float): Spatial discretization [m]
        nz (int): Number of computational nodes
        **kwargs: Additional keyword arguments
            Stefan (bool): Validation with Stefan problem (default: True)
    Returns:
        tuple: A tuple containing:
            - if_depth (float): Location of the ice-water interface/sea ice total thickness [m]
            - if_depth_index (int): Index of the 'transition cell' from ice to ocean (freezing) or water to ice (melting)
    """

    # initialize variable for for-loop

    is_stefan = kwargs.get("Stefan", True)
    Z = (nz - 1) * dz
    depth = np.linspace(dz, Z + dz, nz)
    ice_depth_index = 0
    ice_depth = 0

    for i in range(nz):
        if (
            is_stefan and phi[i] < 0.95 or not is_stefan and phi[i] < 1
        ):  # Simplified condition for Stefan problem
            ice_depth = depth[i]
            ice_depth_index = i + 1
        else:
            break
    return ice_depth, ice_depth_index


class SeaIceModel:
    """SeaIceModelClass represents a class that models the behavior of sea ice."""

    def __init__(self, dataclass: dataclass, user_dataclass: UserInput):
        """
        Args:
            dataclass (PreprocessData): The preprocessed data for the model.
            user_dataclass (UserInput): The user input data for the model.
        """

        self.mush_lowerbound = 0.1
        self.mush_upperbound = 0.9
        self.preprocess_data = dataclass
        self.ui_data = dataclass
        self.results = ResultsParams(dataclass.max_iterations, dataclass.nz)
        self.results.store_results(
            dataclass,
            dataclass.temperature[0],
            dataclass.salinity[0],
            dataclass.liquid_fraction[0],
            dataclass.enthalpy[0],
            dataclass.solid_enthalpy[0],
            dataclass.upwind_velocity[0],
            dataclass.ice_thickness,
            dataclass.time_passed,
            0,
        )

    def bc_neumann(self, phi_k, nz, bc_condition=None):
        """Apply Neumann boundary condition to the sea ice model.

        Args:
            phi_k (float): The value of phi at the k-th layer.
            nz (int): The number of layers in the sea ice model.
            bc_condition (str, optional): The type of boundary condition to apply. Defaults to None.
        Returns:
            None
        """

        if bc_condition == "Neumann":
            self.preprocess_data.temp_grad = temperature_gradient(
                phi_k, nz
            )  # T_bottom - T_top / depth according to Buffo
        else:
            self.preprocess_data.temp_grad = None

    def set_boundary_condition_type(self, critical_depth, bc_type="Neumann"):
        """Sets the boundary condition type for the model. This method sets the boundary condition type for the model. It calculates the temperature gradient based on the given critical depth and boundary condition type. If the boundary condition type is "Neumann", the temperature gradient is calculated using the formula:
        temp_grad = dz * (temperature_melt - boundary_top_temperature) / critical_depth
        If the boundary condition type is not "Neumann", the temperature gradient is set to None.

        Args:
            critical_depth (float): The critical depth value.
            bc_type (str, optional): The type of boundary condition. Defaults to "Neumann".

        Example:
            model.set_boundary_condition_type(10.0, "Neumann")
        """

        """Sets the boundary condition type for the model.
        Parameters:
        - critical_depth (float): The critical depth value.
        - bc_type (str, optional): The type of boundary condition. Defaults to "Neumann".
        Returns:
        None
        Description:
        This method sets the boundary condition type for the model. It calculates the temperature gradient based on the given critical depth and boundary condition type. If the boundary condition type is "Neumann", the temperature gradient is calculated using the formula:
        
        """

        if bc_type == "Neumann":
            self.preprocess_data.temp_grad = (
                self.preprocess_data.dz
                * (
                    self.preprocess_data.temperature_melt
                    - self.preprocess_data.boundary_top_temperature
                )
                / critical_depth
            )
        else:
            self.preprocess_data.temp_grad = None

    def t_running(self, fig, ax1, t_stefan, t_k, t_k_buffo=None, count=0):
        """Plot the temperature profile against depth.

        Args:
            fig (matplotlib.figure.Figure): The figure object to plot on.
            ax1 (matplotlib.axes.Axes): The axes object to plot on.
            t_stefan (numpy.ndarray): The temperature profile obtained analytically.
            t_k (numpy.ndarray): The temperature profile obtained numerically.
            t_k_buffo (numpy.ndarray, optional): The temperature profile obtained using Buffo method. Defaults to None.
            count (int, optional): The count value. Defaults to 0.
        Returns:
            None
        """
        depth_arr = np.linspace(0, -self.preprocess_data.Z, self.preprocess_data.nz)
        ax1.plot(t_k, depth_arr, "r--")
        ax1.plot(t_stefan, depth_arr, "k")
        if t_k_buffo is not None:
            ax1.plot(t_k_buffo, depth_arr, "b-.", alpha=0.5)
        ax1.set_ylabel("Depth [m]")
        ax1.set_xlabel("Temperature [K]")
        # ax1.set_yscale("log")
        # fig.tight_layout()
        if count == 1:
            plt.legend(["Numerical", "Analytical", "Buffo"])
        # display(fig)
        # clear_output(wait=True)

    def track_mush_for_parameter(self, phi_k_, param, param_iterlist):
        """Track the mush for a given parameter.

        Args:
            phi_k_: numpy array representing the values of phi_k
            param: numpy array representing the parameter values
            param_iterlist: list to store the tracked mush values for the parameter
            Updated list with the tracked mush values for the parameter
        """
        # check for values of phi_k_ that fall between the mush_lowerbound and mush_upperbound range
        mush_cond = (phi_k_ >= self.mush_lowerbound) & (phi_k_ <= self.mush_upperbound)
        # if there are values that fall within the range, get the index of the middle value
        if mush_cond.any():
            mush_indx = np.where(mush_cond)[0][len(np.where(mush_cond)[0]) // 2 - 1]
        else:
            mush_cond = phi_k_ > self.mush_upperbound
            mush_indx = 0

        if (
            param[phi_k_ < self.mush_lowerbound].size
            & param[phi_k_ > self.mush_upperbound].size
        ):
            param_iterlist.append(
                [
                    param[phi_k_ < self.mush_lowerbound][0],
                    param[phi_k_ == phi_k_][mush_indx],
                    param[phi_k_ > self.mush_upperbound][-1],
                ]
            )
        elif param[phi_k_ > self.mush_upperbound].size:
            param_iterlist.append(
                [
                    param[phi_k_ == phi_k_][mush_indx],
                    param[phi_k_ == phi_k_][mush_indx],
                    param[phi_k_ > self.mush_upperbound][-1],
                ]
            )
        else:
            param_iterlist.append(
                [
                    param[phi_k_ == phi_k_][mush_indx],
                    param[phi_k_ == phi_k_][mush_indx],
                    param[phi_k_ == phi_k_][mush_indx],
                ]
            )

        return param_iterlist, mush_indx

    def phi_all_mush_list(self, phi_k_, phi_all_mush_list):
        """Calculates the number of elements in phi_k_ that fall within the mush_lowerbound and mush_upperbound range.

        Args:
            phi_k_ (numpy.ndarray): The input array containing the values to be checked.
            phi_all_mush_list (list): The list to which the count of elements within the range will be appended.
        Returns:
            list: The updated phi_all_mush_list with the count of elements within the range appended.
        """

        mush_cond = (phi_k_ >= self.mush_lowerbound) & (phi_k_ <= self.mush_upperbound)
        phi_all_mush_list.append(np.count_nonzero(mush_cond))
        return phi_all_mush_list

    def convergence_loop_iteration(
        self,
        t,
        t_km1,
        s_km1,
        brine_velocity_km1,
        nutrient_cn_km1,
        phi_km1,
        thickness_index,
        source_term_temperature,
        source_term_salinity,
        x_wind_temperature,
        x_wind_salinity,
        buffo=False,
        stefan=False,
        voller=False,
        temp_grad=None,
        salinity_equation=False,
    ):
        """Performs a single iteration of the convergence loop by first resetting the parameters to their previous time step values and then running the convergence loop until convergence is reached.
        Args:
            t (float): Current temperature.
            t_km1 (float): Temperature at the previous time step.
            s_km1 (float): Salinity at the previous time step.
            phi_km1 (float): Porosity at the previous time step.
            buffo (bool, optional): Flag indicating whether to use the buffo method. Defaults to False.
            stefan (bool, optional): Flag indicating whether to use the Stefan method. Defaults to False.
            temp_grad (float, optional): Temperature gradient. Defaults to None.

        Methods called:
            - reset_iteration_parameters: Resets the iteration parameters for the sea ice model with values at previous time steps.
            - run_while_convergence_iteration: Runs the convergence loop until convergence is reached.

        Returns:
            tuple: A tuple containing the following values:
                - t_k (float): Current temperature.
                - t_prev (float): Temperature at the previous time step.
                - s_k (float): Current salinity.
                - s_prev (float): Salinity at the previous time step.
                - phi_k (float): Current porosity.
                - phi_prev (float): Porosity at the previous time step.
                - h_k (float): Current heat flux.
                - h_solid (float): Heat flux at the solid-liquid interface.
                - thickness (float): Current thickness.
                - thickness_index (int): Index of the thickness.
                - t_km1 (float): Temperature at the previous time step.
                - s_km1 (float): Salinity at the previous time step.
                - phi_km1 (float): Porosity at the previous time step.
                - temperature_liquidus (float): Liquidus temperature.
                - temperature_solidus (float): Solidus temperature.
        """

        (
            t_km1,
            s_km1,
            brine_velocity_km1,
            nutrient_cn_km1,
            phi_km1,
            temp_grad,
            t_err,
            s_err,
            phi_err,
            t_initial,
            s_initial,
            phi_initial,
            source_term_temperature,
            source_term_salinity,
            counter,
        ) = self.reset_iteration_parameters(
            t,
            t_km1,
            s_km1,
            brine_velocity_km1,
            nutrient_cn_km1,
            phi_km1,
            source_term_temperature,
            source_term_salinity,
        )
        (
            t_km1,
            s_km1,
            brine_velocity_k,
            nutrient_cn_km1,
            phi_km1,
            t_prev,
            s_prev,
            phi_prev,
            h_k,
            h_solid,
            phi_k,
            t_k,
            s_k,
            nutrient_cn_k,
            thickness,
            thickness_index,
            temperature_liquidus,
            temperature_solidus,
            x_wind_temperature,
            x_wind_salinity,
        ) = self.run_while_convergence_iteration(
            t,
            t_km1,
            s_km1,
            brine_velocity_km1,
            nutrient_cn_km1,
            phi_km1,
            thickness_index,
            buffo,
            stefan,
            voller,
            t_err,
            s_err,
            phi_err,
            source_term_temperature,
            source_term_salinity,
            x_wind_temperature,
            x_wind_salinity,
            counter,
            salinity_equation,
        )
        return (
            t_k,
            t_prev,
            s_k,
            s_prev,
            brine_velocity_k,
            nutrient_cn_k,
            nutrient_cn_km1,
            phi_k,
            phi_prev,
            h_k,
            h_solid,
            thickness,
            thickness_index,
            t_km1,
            s_km1,
            phi_km1,
            temperature_liquidus,
            temperature_solidus,
            x_wind_temperature,
            x_wind_salinity,
        )

    def run_while_convergence_iteration(
        self,
        t,
        t_initial,
        s_initial,
        brine_velocity_initial,
        nutrient_cn_initial,
        phi_initial,
        thickness_index,
        buffo,
        stefan,
        voller,
        t_err,
        s_err,
        phi_err,
        source_term_temperature,
        source_term_salinity,
        x_wind_temperature,
        x_wind_salinity,
        counter,
        _is_salinity_equation=False,
    ):
        """Runs the convergence loop until convergence is reached.

        Args:
            t: Time step
            t_km1: Previous temperature array
            s_km1: Previous salinity array
            phi_km1: Previous phi array
            buffo: Flag for Buffo
            stefan: Flag for Stefan
            t_err: Temperature error
            s_err: Salinity errorhatch
            phi_err: Phi error
            t_initial: Initial temperature
            phi_initial: Initial phi
            source_term: RHS Source term of PDE equation
            counter: Iteration counter

        Methods called:
            - initialise_previous_statevariables: Defines the previous state variables for convergence iteration (temperature,salinity, liquid_fraction).
            - update_state_variables: Update the state variables (Enthalpy, Enthalpy Solid, Liquid Fraction, Temperature, Salinity).
            - locate_ice_ocean_interface: Locate the ice-ocean interface based on the liquid fraction.
            - overwrite_statevariables:
            - track_mush_for_parameter
            - phi_all_mush_list
            - check_convergence: computes error for convergence

        Returns:
            tuple: A tuple containing the following values:
                - t_km1 (float): Temperature at the previous time step.
                - s_km1 (float): Salinity at the previous time step.
                - phi_km1 (float): Liquid fraction at the previous time step.
                - t_prev (float): Temperature at the previous time step.
                - s_prev (float): Salinity at the previous time step.
                - phi_prev (float): Liquid fraction at the previous time step.
                - h_k (float): Current heat flux.
                - h_solid (float): Heat flux at the solid-liquid interface.
                - phi_k (float): Liquid fraction at the current time step.
                - t_k (float): Current temperature.
                - s_k (float): Current salinity.
                - thickness (float): Current thickness.
                - thickness_index (int): Index of the thickness
        """

        # Set previous state variables temperature, salinity, liquid fraction respectively
        (
            t_prev,
            s_prev,
            brine_velocity_prev,
            nutrient_cn_prev,
            phi_prev,
            t_k_melt,
            t_k_A_LHS_matrix,
            temp_factor3,
            x_wind_temperature_prev,
            x_wind_salinity_prev,
        ) = overwrite_statevariables(
            t_initial,
            s_initial,
            brine_velocity_initial,
            nutrient_cn_initial,
            phi_initial,
            x_wind_temperature,
            x_wind_salinity,
        )
        residual_voller = 1
        # Run the while loop until convergence is reached
        # FIXME: tolerance check and change for new enthalpy update for inhomognoeus physical values
        # TODO: Document Voller method on overleaf
        while (
            ((residual_voller > self.preprocess_data.temperature_tolerance) & (stefan))
            or (t_err > self.preprocess_data.temperature_tolerance)
            or (phi_err > self.preprocess_data.liquid_fraction_tolerance)
            or (s_err > self.preprocess_data.salinity_tolerance)
        ):
            # Update state variables Enthalpy, Enthalpy Solid, Liquid Fraction, Temperature, Salinity respectively
            (
                h_k,
                h_solid,
                phi_k,
                t_k,
                s_k,
                brine_velocity,
                nutrient_cn,
                t_k_A_LHS_matrix,
                temp_factor3,
                t_k_melt,
                temperature_liquidus,
                temperature_solidus,
                x_wind_temperature,
                x_wind_salinity,
            ) = update_state_variables(
                self.preprocess_data,
                t_prev,
                s_prev,
                brine_velocity_prev,
                nutrient_cn_prev,
                phi_prev,
                thickness_index,
                buffo,
                stefan,
                voller,
                t_initial,
                s_initial,
                nutrient_cn_initial,
                phi_initial,
                t_k_melt,
                t_k_A_LHS_matrix,
                temp_factor3,
                source_term_salinity,
                source_term_temperature,
                x_wind_temperature_prev,
                x_wind_salinity_prev,
                _is_salinity_equation=_is_salinity_equation,
            )
            # Locate ice-ocean interface based on liquid fraction
            thickness, thickness_index = locate_ice_ocean_interface(
                phi_k,
                self.preprocess_data.grid_resolution_dz,
                self.preprocess_data.nz,
                Stefan=self.preprocess_data.is_stefan,
            )

            # Check for convergence
            t_err, s_err, phi_err, counter, residual_voller = self.check_convergence(
                t,
                counter,
                t_prev,
                s_prev,
                phi_prev,
                phi_k,
                t_k,
                s_k,
                t_initial,
                phi_initial,
                s_initial,
                voller=voller,
                A_matrix=t_k_A_LHS_matrix,
                temp_factor3=temp_factor3,
            )
            # Update state variables temperature, salinity, liquid fraction respectively
            (
                t_prev,
                s_prev,
                brine_velocity_prev,
                nutrient_cn_prev,
                phi_prev,
                t_k_melt,
                t_k_A_LHS_matrix,
                temp_factor3,
                x_wind_temperature_prev,
                x_wind_salinity_prev,
            ) = overwrite_statevariables(
                t_k,
                s_k,
                brine_velocity,
                nutrient_cn,
                phi_k,
                t_k_melt,
                x_wind_temperature,
                x_wind_salinity,
                t_k_A_LHS_matrix,
                temp_factor3,
            )
            # Track mushy layer using liquid fraction for temperature and phi values
            self.record_mushy_layer_data(
                t, t_prev, stefan, phi_prev, residual_voller, s_prev
            )
            # record thickness index at the given timestep t
            self.preprocess_data.thickness_index_total[t - 1] = thickness_index

            if counter >= self.preprocess_data.counter_limit:
                export_residuals(
                    self,
                    residuals=self.preprocess_data.residual_voller_all,
                    temperature_mushy=self.preprocess_data.t_k_iter_all,
                    phi_mushy=self.preprocess_data.all_phi_iter_all,
                    salinity_mushy=self.preprocess_data.s_k_iter_all,
                    output_dir=self.preprocess_data.dir_output_name,
                )
                msg = f"Convergence not reached at time t = {t}"
                raise ConvergenceError(msg)

        return (
            t_prev,
            s_prev,
            brine_velocity_prev,
            nutrient_cn_prev,
            phi_prev,
            t_initial,
            s_initial,
            phi_initial,
            h_k,
            h_solid,
            phi_k,
            t_k,
            s_k,
            nutrient_cn,
            thickness,
            thickness_index,
            temperature_liquidus,
            temperature_solidus,
            x_wind_temperature_prev,
            x_wind_salinity_prev,
        )

    def check_convergence(
        self,
        t,
        counter,
        t_prev,
        s_prev,
        phi_prev,
        phi_k,
        t_k,
        s_k,
        t_initial,
        phi_initial,
        s_initial,
        voller=False,
        **kwargs,
    ):
        if counter > 0:
            (
                t_err,
                t_err_full,
                s_err,
                s_err_full,
                phi_err,
                phi_err_full,
                residual_voller,
            ) = compute_error_for_convergence(
                t_k,
                t_prev,
                s_k,
                s_prev,
                phi_k,
                phi_prev,
                voller,
                temp_factor3=kwargs.get("temp_factor3"),
                t_initial=t_initial,
                phi_initial=phi_initial,
                s_initial=s_initial,
                A_matrix=kwargs.get("A_matrix"),
            )

        counter += 1
        return t_err, s_err, phi_err, counter, residual_voller

    def record_mushy_layer_data(self, t, t_km1, stefan, phi_k, residual_voller, s_prev):
        """Records the mushy layer data for temperature and phi values at time iterations t at specific time steps corresponding to initial stages, middle and final stages of the process."""
        if stefan:
            # if t in [7, 36, 720, 7200, 14400, 21600] and stefan:
            self.preprocess_data.t_k_iter, t_mush_indx = self.track_mush_for_parameter(
                phi_k, t_km1, self.preprocess_data.t_k_iter
            )
            self.preprocess_data.phi_k_iter, phi_mush_indx = (
                self.track_mush_for_parameter(
                    phi_k, phi_k, self.preprocess_data.phi_k_iter
                )
            )
            self.preprocess_data.all_phi_iter = self.phi_all_mush_list(
                phi_k, self.preprocess_data.all_phi_iter
            )
            self.preprocess_data.mush_indx_list.append([t_mush_indx, phi_mush_indx])
            self.preprocess_data.t_k_before_convergence.append(t_km1)
            self.preprocess_data.residual_voller.append(residual_voller)
            self.preprocess_data.s_k_iter.append(
                [s_prev[0], s_prev[t_mush_indx], s_prev[-1]]
            )

    def reset_iteration_parameters(
        self,
        t,
        tkm1,
        s_km1,
        brine_velocity_km1,
        nutrient_cn_km1,
        phi_km1,
        source_term_temperature,
        source_term_salinity,
    ):
        """Reset the iteration parameters for the sea ice model.

        Args:
            t (float): Current temperature.
            tkm1 (float): Temperature at the previous time step.
            s_km1 (float): Salinity at the previous time step.
            phi_km1 (float): Liquid fraction at the previous time step.
        Returns:
            tuple: A tuple containing the following iteration parameters:
                - t_km1 (float): Temperature at the previous time step.
                - s_km1 (float): Salinity at the previous time step.
                - phi_km1 (float): Liquid fraction at the previous time step.
                - temp_grad (float): Temperature gradient.
                - t_err (float): Temperature error.
                - s_err (float): Salinity error.
                - phi_err (float): Liquid fraction error.
                - t_initial (float): Initial temperature.
                - s_initial (float): Initial salinity.
                - phi_initial (float): Initial liquid fraction.
                - source_term_array (ndarray): Array of source-term values.
                - counter (int): Iteration counter.
        """
        t_err, s_err, phi_err = reset_error_for_while_loop(
            self.preprocess_data.temperature_tolerance,
            self.preprocess_data.salinity_tolerance,
            self.preprocess_data.liquid_fraction_tolerance,
        )
        t_initial, t_km1, s_initial, s_km1, phi_initial, phi_km1 = (
            self.initialize_state_variables(t, tkm1, s_km1, phi_km1)
        )
        nutrient_cn_initial = nutrient_cn_km1
        temp_grad = self.preprocess_data.temp_grad
        source_term_temperature = source_term_temperature
        source_term_salinity = source_term_salinity
        counter = 1
        self.record_iteration_data()
        return (
            t_km1,
            s_km1,
            brine_velocity_km1,
            nutrient_cn_initial,
            phi_km1,
            temp_grad,
            t_err,
            s_err,
            phi_err,
            t_initial,
            s_initial,
            phi_initial,
            source_term_temperature,
            source_term_salinity,
            counter,
        )

    def check_and_reset_any_iteration_data(
        self, parameter_list, parameter_all_list
    ) -> None:
        """Check the iteration data for temperature and phi values.
        This method appends the temperature and phi values from the current iteration to the respective arrays.
        The arrays are used to store the iteration data for further analysis.

        Args:
            parameter_list (list): The list of parameter values.
            parameter_all_list (list): The list of all parameter values.
        Returns:
            None
        """

        if np.array(parameter_list).any():
            parameter_all_list.append(parameter_list)
        parameter_list = []

    def record_iteration_data(self):
        """Records the iteration data for temperature and phi values.
        This method appends the temperature and phi values from the current iteration to the respective arrays.
        The arrays are used to store the iteration data for further analysis.

        Args:
            None
        Returns:
            None
        """

        if np.array(self.preprocess_data.t_k_iter).any():
            self.preprocess_data.t_k_iter_all.append(self.preprocess_data.t_k_iter)
        if np.array(self.preprocess_data.phi_k_iter).any():
            self.preprocess_data.phi_k_iter_all.append(self.preprocess_data.phi_k_iter)
        if np.array(self.preprocess_data.all_phi_iter).any():
            self.preprocess_data.all_phi_iter_all.append(
                self.preprocess_data.all_phi_iter
            )
        if np.array(self.preprocess_data.t_k_before_convergence).any():
            self.preprocess_data.t_k_before_convergence_all.append(
                self.preprocess_data.t_k_before_convergence
            )
        if np.array(self.preprocess_data.mush_indx_list).any():
            self.preprocess_data.mush_indx_list_all.append(
                self.preprocess_data.mush_indx_list
            )
        if np.array(self.preprocess_data.residual_voller).any():
            self.preprocess_data.residual_voller_all.append(
                self.preprocess_data.residual_voller
            )
        if np.array(self.preprocess_data.s_k_iter).any():
            self.preprocess_data.s_k_iter_all.append(self.preprocess_data.s_k_iter)
        (
            self.preprocess_data.t_k_iter,
            self.preprocess_data.phi_k_iter,
            self.preprocess_data.all_phi_iter,
            self.preprocess_data.t_k_before_convergence,
            self.preprocess_data.mush_indx_list,
            self.preprocess_data.residual_voller,
            self.preprocess_data.s_k_iter,
        ) = [], [], [], [], [], [], []

        # parameter_to_record_and_reset_list = [
        #     [self.preprocess_data.t_k_iter, self.preprocess_data.t_k_iter_all],
        #     [self.preprocess_data.phi_k_iter, self.preprocess_data.phi_k_iter_all],
        #     [self.preprocess_data.all_phi_iter, self.preprocess_data.all_phi_iter_all],
        #     [
        #         self.preprocess_data.t_k_before_convergence,
        #         self.preprocess_data.t_k_before_convergence_all,
        #     ],
        #     [
        #         self.preprocess_data.mush_indx_list,
        #         self.preprocess_data.mush_indx_list_all,
        #     ],
        # ]

        # for parameter_list, parameter_all_list in parameter_to_record_and_reset_list:
        #     self.check_and_reset_any_iteration_data(parameter_list, parameter_all_list)

    def initialize_state_variables(self, t, t_km1, s_km1, phi_km1):
        """Initializes the state variables for the sea ice model.

        Args:
            t (int): The current time step.
            t_km1 (float): The temperature at the previous time step.
            s_km1 (float): The salinity at the previous time step.
            phi_km1 (float): The liquid fraction at the previous time step.
        Returns:
            tuple: A tuple containing the initialized state variables:
                - t_initial (float): The initial temperature.
                - t_km1 (float): The temperature at the previous time step.
                - s_km1 (float): The salinity at the previous time step.
                - phi_initial (float): The initial liquid fraction.
                - phi_km1 (float): The liquid fraction at the previous time step.
                - temp_grad (float): The temperature gradient.
        """

        if t == 1:
            (
                t_initial,
                t_km1,
                s_initial,
                s_km1,
                phi_initial,
                phi_km1,
            ) = set_statevariables(
                self.preprocess_data.temperature,
                self.preprocess_data.salinity,
                self.preprocess_data.liquid_fraction,
            )
        else:
            t_initial, t_km1, s_initial, s_km1, phi_initial, phi_km1 = (
                set_statevariables(t_km1, s_km1, phi_km1)
            )
        return (
            t_initial,
            t_km1,
            s_initial,
            s_km1,
            phi_initial,
            phi_km1,
        )

    def choose_phase_type_iteration(self, t):
        """Choose the phase type iteration based on the one-phase and two-phase generalised Stefan Probem.

        Args:
            t (int): The time index.
        Returns:
            tuple: A tuple containing the following values:
                - t_stefan (float): The Stefan temperature.
                - error_depth_t (float): The error in depth.
                - depth_stefan_t (float): The depth at time t.
        Raises:
            InvalidPhaseError: If the phase type is invalid (not 1 or 2).
        """

        if self.preprocess_data.phase_type == 1:
            self.preprocess_data.depth_stefan_all = StefanProblem.stefan_problem(
                self.results.all_t_passed, self.ui_data
            )
            depth_stefan_t = self.preprocess_data.depth_stefan_all[t - 1]
            t_stefan = StefanProblem.calculate_temperature_profile(
                depth_stefan_t,
                self.preprocess_data.time_passed,
                self.preprocess_data.grid_resolution_dz,
                self.preprocess_data.nz,
                self.ui_data,
            )
        elif self.preprocess_data.phase_type == 2:
            self.preprocess_data.depth_stefan_all = (
                StefanProblem.stefan_problem_twophase(
                    self.preprocess_data.all_t_passed, self.ui_data
                )
            )
            depth_stefan_t = self.preprocess_data.depth_stefan_all[t - 1]
            t_stefan, c_stefan = StefanProblem.calculate_temperature_twophase_profiles(
                depth_stefan_t,
                self.preprocess_data.preprocess_data.time_passed,
                self.preprocess_data.grid_resolution_dz,
                self.preprocess_data.nz,
                self.ui_data,
            )
        else:
            msg = "Invalid phase type. Choose 1 or 2."
            raise InvalidPhaseError(msg)

        error_depth_t = np.abs(depth_stefan_t) - np.abs(self.results.all_thick[t - 1])
        return t_stefan, error_depth_t, depth_stefan_t

    def run_sea_ice_model(self):
        """Runs the sea ice model.

        This function iterates over a specified number of time steps and performs calculations
        to simulate the behavior of sea ice. It updates the results and saves a temperature profile
        plot at the end.

        Args:
            None

        Returns:
            None
        """

        fig1, ax1 = plt.subplots(1, 1)
        # set nutrient concentration initial values
        # TODO: reorganize the init to another script
        nutrient_cn_km1 = (
            self.results.nutrient_concentration_multiplelayers[0]
            * self.preprocess_data.nutrient_cn_dsi_water
        )
        nutrient_cn_km1[0] = self.preprocess_data.nutrient_cn_dsi_ice
        nutrient_cn_km1_buffo = nutrient_cn_km1
        carbon_cc = (
            self.results.carbon_concentration_multiplelayers[0]
            * self.preprocess_data.carbon_cc_water_initial
        )
        carbon_cc[0] = self.preprocess_data.carbon_cc_ice_initial
        # thickness
        thickness_index, thickness_index_buffo = 0, 0
        # carbon_concentration, nutrient_concentration, photosynthetic_rate_mu, radiation_algae = biogeochemical_model(self.preprocess_data.boundary_top_temperature, self.preprocess_data.boundary_salinity, 1.0, self.results.nutrient_concentration[0], self.results.carbon_concentration[0], self.preprocess_data.grid_timestep_dt, 0, 0)
        radiative_source_term = 0.0
        salinity_source_term = 0.0
        count = 0
        # set temperature salinity and liquid fraction values
        (
            t_km1,
            s_km1,
            phi_km1,
            x_wind_temperature,
            x_wind_salinity,
            brine_velocity_km1,
        ) = (
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            self.preprocess_data.upwind_velocity,
        )
        (
            t_km1_buffo,
            s_km1_buffo,
            phi_km1_buffo,
            x_wind_temperature_buffo,
            x_wind_salinity_buffo,
            brine_velocity_km1_buffo,
        ) = (
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            self.preprocess_data.upwind_velocity,
        )
        # with alive_bar(self.preprocess_data.max_iterations, force_tty=True) as bar:
        for t in range(1, self.preprocess_data.max_iterations):
            # time.sleep(0.005)
            if self.preprocess_data.is_buffo:
                (
                    t_k_buffo,
                    t_prev_buffo,
                    s_k_buffo,
                    s_prev_buffo,
                    brine_velocity_km1_buffo,
                    nutrient_cn_km1_buffo,
                    nutrient_cn_prev_buffo,
                    phi_k_buffo,
                    phi_prev_buffo,
                    h_k_buffo,
                    h_solid_buffo,
                    thickness_buffo,
                    thickness_index_buffo,
                    t_km1_buffo,
                    s_km1_buffo,
                    phi_km1_buffo,
                    temperature_liquidus_buffo,
                    temperature_solidus_buffo,
                    x_wind_temperature_buffo,
                    x_wind_salinity_buffo,
                ) = self.convergence_loop_iteration(
                    t,
                    t_km1_buffo,
                    s_km1_buffo,
                    brine_velocity_km1_buffo,
                    nutrient_cn_km1_buffo,
                    phi_km1_buffo,
                    thickness_index_buffo,
                    radiative_source_term,
                    salinity_source_term,
                    x_wind_temperature_buffo,
                    x_wind_salinity_buffo,
                    stefan=True,
                    buffo=False,
                    temp_grad=self.preprocess_data.temp_grad,
                    salinity_equation=self.preprocess_data.is_salinity_equation,
                )
            elif not self.preprocess_data.is_buffo:
                (t_k_buffo, s_k_buffo, phi_k_buffo, thickness_buffo) = (
                    self.preprocess_data.t_k_buffo_list[0],
                    self.preprocess_data.s_buffo_list[0],
                    self.preprocess_data.phi_buffo_list[0],
                    self.preprocess_data.thickness_list_buffo[0],
                )
            (
                t_k,
                t_prev,
                s_k,
                s_prev,
                brine_velocity_km1,
                nutrient_cn_k,
                nutrient_cn_km1,
                phi_k,
                phi_prev,
                h_k,
                h_solid,
                thickness,
                thickness_index,
                t_km1,
                s_km1,
                phi_km1,
                temperature_liquidus,
                temperature_solidus,
                x_wind_temperature,
                x_wind_salinity,
            ) = self.convergence_loop_iteration(
                t,
                t_km1,
                s_km1,
                brine_velocity_km1,
                nutrient_cn_km1,
                phi_km1,
                thickness_index,
                radiative_source_term,
                salinity_source_term,
                x_wind_temperature,
                x_wind_salinity,
                stefan=True,
                buffo=False,
                voller=self.preprocess_data.is_voller,
                temp_grad=self.preprocess_data.temp_grad,
                salinity_equation=self.preprocess_data.is_salinity_equation,
            )
            self.preprocess_data.time_passed = t_total(
                self.preprocess_data.time_passed,
                self.preprocess_data.grid_timestep_dt,
            )
            self.bc_neumann(
                phi_k,
                self.preprocess_data.boundary_condition_type,
                self.preprocess_data.nz,
            )
            t_stefan, error_depth_t, thickness_stefan = (
                self.choose_phase_type_iteration(t)
            )

            # TODO: Add radiation and algae
            # TODO: create a parameter script for algae model
            (
                radiative_source_term,
                salinity_source_term,
                carbon_cc,
                nutrient_cn_k,
                photosynthetic_rate_mu,
                radiation_algae,
                chla_bulk,
            ) = self.calculate_source_terms(
                carbon_cc,
                thickness_index,
                t,
                t_k,
                s_k,
                nutrient_cn_k,
                phi_k,
                thickness,
                algae_model_depth_type="all",
            )

            nutrient_cn_km1 = nutrient_cn_k
            
            self.results = self.results.store_results(
                self.results,
                t_k,
                s_k,
                phi_k,
                h_k,
                h_solid,
                self.preprocess_data.upwind_velocity[0],
                thickness,
                self.preprocess_data.time_passed,
                t,
            )
            self.results = self.results.store_results_for_iter_t(
                self.results,
                t,
                thickness_index,
                t_k,
                t_stefan,
                s_k,
                s_k_buffo,
                brine_velocity_km1,
                nutrient_cn_k,
                carbon_cc,
                phi_k,
                phi_k_buffo,
                h_k,
                h_solid,
                thickness,
                thickness_buffo,
                thickness_stefan,
                t_k_buffo,
                temperature_liquidus,
                temperature_solidus,
                carbon_cc[thickness_index],
                nutrient_cn_k[thickness_index],
                photosynthetic_rate_mu,
                radiation_algae,
                chla_bulk,
                radiative_source_term,
                salinity_source_term,
                buffo=self.preprocess_data.is_buffo,
            )
            # bar()

            if t % 1000 == 0:
                count += 1
                self.t_running(
                    fig1, ax1, t_stefan, t_k=t_k, t_k_buffo=t_k_buffo, count=count
                )

        fig1.savefig(
            f"{self.preprocess_data.dir_output_name}/TemperatureProfile.pdf",
            backend="pgf",
        )

    def calculate_source_terms(
        self,
        carbon_cc,
        thickness_index,
        t,
        t_k,
        s_k,
        nutrient_cn_k,
        phi_k,
        thickness,
        **kwargs,
    ):
        Z = (self.preprocess_data.nz - 1) * self.preprocess_data.grid_resolution_dz
        depth = np.linspace(
            self.preprocess_data.grid_resolution_dz,
            Z + self.preprocess_data.grid_resolution_dz,
            self.preprocess_data.nz,
        )

        if kwargs.get("algae_model_depth_type") == "single":
            (
                carbon_concentration,
                nutrient_concentration_atinterface,
                photosynthetic_rate_mu,
                radiation_algae,
                chla_bulk,
            ) = biogeochemical_model(
                t_k,
                s_k,
                phi_k,
                nutrient_cn_k[thickness_index],
                carbon_cc[thickness_index],
                self.preprocess_data.grid_timestep_dt,
                thickness_index,
                thickness,
            )
            nutrient_concentration = nutrient_concentration_atinterface
        elif kwargs.get("algae_model_depth_type") == "all":
            (
                carbon_concentration,
                nutrient_concentration,
                photosynthetic_rate_mu,
                radiation_algae,
                chla_bulk,
            ) = biogeochemical_model_at_alldepths(
                t_k,
                s_k,
                phi_k,
                nutrient_cn_k,
                carbon_cc,
                self.preprocess_data.grid_timestep_dt,
                depth,
            )
        else:
            raise ValueError("Invalid algae model depth type.")

        self.results.thickness_list[t] = thickness
        # salinity source term using gravity drainage
        salinity_source_term = get_salinity_source_term(
            int(thickness_index + 1),
            self.results.thickness_list,
            s_k,
            phi_k,
            self.preprocess_data.grid_timestep_dt,
            self.preprocess_data.grid_resolution_dz,
        )

        # get radiative terms all
        radiative_source_term = calculate_radiative_terms(
            depth, thickness_index, radiation_algae, kwargs.get("algae_model_depth_type")
        )
        return (
            radiative_source_term,
            salinity_source_term,
            carbon_concentration,
            nutrient_concentration,
            photosynthetic_rate_mu,
            radiation_algae,
            chla_bulk,
        )

    def add_new_parameter_results(self) -> None:
        self.results.t_k_iter_all = self.preprocess_data.t_k_iter_all
        self.results.phi_k_iter_all = self.preprocess_data.phi_k_iter_all
        self.results.all_phi_iter_all = self.preprocess_data.all_phi_iter_all
        self.results.t_k_before_convergence_all = (
            self.preprocess_data.t_k_before_convergence_all
        )
        self.results.mush_indx_list_all = self.preprocess_data.mush_indx_list_all
        self.results.residual_voller_all = self.preprocess_data.residual_voller_all
        self.results.s_k_iter_all = self.preprocess_data.s_k_iter_all

        return self.results

    @classmethod
    def get_results(
        cls, dataclass: dataclass, user_dataclass: UserInput
    ) -> ResultsParams:
        """Runs the sea ice model and returns the results.

        Args:
            cls (class): The class object.
            dataclass (PreprocessData): The dataclass containing preprocessed data.
            user_dataclass (UserInput): The dataclass containing user input.
        Returns:
            Results: The results dataclass object generated by running the sea ice model.
        """
        print("Running model...")
        results_obj = cls(dataclass, user_dataclass)
        results_obj.set_dataclass(dataclass)
        results_obj.run_sea_ice_model()
        final_results = results_obj.add_new_parameter_results()
        print("Model run complete and Ready for Analysis.")

        return final_results

    def set_dataclass(self, _dataclass):
        """Sets the dataclass attributes of the object.

        Args:
            _dataclass: An instance of the dataclass.
        Returns:
            None
        """

        data_class_obj = _dataclass()
        for key, value in asdict(data_class_obj).items():
            setattr(self, key, value)
