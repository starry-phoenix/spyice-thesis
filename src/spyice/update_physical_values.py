from __future__ import annotations

import numpy as np
from scipy.optimize import bisect, newton
from typing import TYPE_CHECKING

from src.spyice.parameters.user_input import UserInput
from src.spyice.models.advection_diffusion import AdvectionDiffusion

if TYPE_CHECKING:
    from preprocess.pre_process import PreprocessData

_ui = UserInput()
_temperature_melt, _boundary_top_temperature = (
    _ui.temperature_melt,
    _ui.boundary_top_temperature,
)
_specific_heat_ice, _specific_heat_brine, _latent_heat_water = _ui.constants.c_i, _ui.constants.c_br ,_ui.constants.L
_rho_ice, _rho_brine = _ui.constants.rho_i, _ui.constants.rho_br

def calculate_melting_temperature_from_salinity(
    _salinity, _temperature_melt=_temperature_melt, _liquid_relation="Normal"
):
    """Calculates the melting temperature of seawater based on salinity.

    Args:
        _salinity (numpy.ndarray): Array of salinity values.
        _temperature_melt (float, optional): Melting temperature. Defaults to _temperature_melt.
        _liquid_relation (str, optional): Liquid relation type. Must be either "Normal" or "Frezchem". Defaults to "Normal".
    Returns:
        numpy.ndarray: Array of melting temperature values.
    Raises:
        TypeError: If _liquid_relation is not "Normal" or "Frezchem".
    """

    if _liquid_relation not in ["Normal", "Frezchem"]:
        raise TypeError("Liquid relation not available")
    if _liquid_relation == "Normal":
        T_m = _temperature_melt # melt temperature as solidus temperature
        T_s = 252.05 # eutectic temperature for Sbr = 233ppt
        S_br = 233.0  # brine salinity in ppt
        _melting_temperature_seawater = _temperature_melt * np.ones(_salinity.shape) + (T_s - T_m)*_salinity/S_br
    elif _liquid_relation == "Frezchem":
        # _melting_temperature_seawater = _temperature_melt + (
        #     -(9.1969758 * (1e-05) * _salinity**2) - 0.03942059 * _salinity
        # )
        _melting_temperature_seawater = 272.63617665 + (
            -(9.1969758 * (1e-05) * _salinity**2) - 0.03942059 * _salinity
        )
    return _melting_temperature_seawater

def update_liquid_fraction_buffo(
    _temperature,
    _salinity,
    _liquid_fraction,
    _enthalpy,
    _enthalpy_solid,
    _nz,
    _is_stefan=False,
    _method="likebuffo",
    _pt2_system=False,
):
    """Updates the liquid fraction based on temperature, salinity, enthalpy, and other parameters.

    Args:
        _temperature (float): The temperature value.
        _salinity (float): The salinity value.
        _liquid_fraction (float): The liquid fraction value.
        _enthalpy (float): The enthalpy value.
        _enthalpy_solid (float): The solid enthalpy value.
        _nz (int): The number of vertical grid points.
        _is_stefan (bool, optional): Whether to use Stefan condition. Defaults to False.
        _method (str, optional): The method to use. Defaults to "likebuffo".
    Returns:
        tuple: A tuple containing the updated liquid fraction values.
    Raises:
        AssertionError: If the liquid fraction has a non-physical value.
    """
    _phi = np.ones(_nz, dtype=np.float64)
    _alpha = 1.853 / 28.0
    _phi_km1 = _liquid_fraction

    if _pt2_system is False:
        _temperature_difference = (
            _temperature - calculate_melting_temperature_from_salinity(_salinity)
        )
        _enthalpy = _specific_heat_ice * _temperature + _latent_heat_water * _phi_km1
        _enthalpy_solid = (
            _specific_heat_ice * calculate_melting_temperature_from_salinity(_salinity)
        )
        _enthalpy_difference = _enthalpy - _enthalpy_solid
    else:
        _enthalpy_difference = _enthalpy - _enthalpy_solid

    if _is_stefan:
        _phi = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _enthalpy_difference / _latent_heat_water,
            ),
        )
    else:
        _phi = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _enthalpy_difference / _latent_heat_water,
            ),
        )
        assert np.all(
            (_phi >= 0) & (_phi <= 1)
        ), "liquid fraction has non-physical value"

    _phi = np.where(_phi> 1.0, 1.0, np.where(_phi < 0.0, 0.0, _phi))
    _phi = np.round(_phi, 2)
    return _phi * phi_control_for_infinite_values(_phi)


def update_liquid_fraction_voller_under_relaxation(
    _temperature,
    _salinity,
    _temperature_initial,
    _liquid_fraction_previous,
    _temperature_previous,
    _salinity_previous,
    _nz,
    temp_factor_3,
    t_k_A_LHS_matrix,
    t_k_A_LHS_matrix_previous,
    under_relaxation_factor,
    _is_stefan=False,
    _pt2_system=False,
    _nonconstant_physical_properties=False,
):
    """
    Update the liquid fraction using the Voller under-relaxation method.
    Args:
        - _temperature (float): The temperature value.
        - _salinity (float): The salinity value.
        - _liquid_fraction (float): The liquid fraction value.
        - _nz (int): The number of vertical grid points.
        - temp_factor_3 (float): The temperature factor.
        - under_relaxation_factor (float): The under-relaxation factor.
        - _is_stefan (bool, optional): Whether to use Stefan condition. Defaults to False.
    Returns:
        - numpy.ndarray: The updated liquid fraction values.
    Raises:
        - AssertionError: If the liquid fraction has a non-physical value

    """
    _phi = np.ones(_nz, dtype=np.float64)
    if _pt2_system is False:
        _phi_km1 = _liquid_fraction_previous
        _temperature_difference = (
            _temperature - calculate_melting_temperature_from_salinity(_salinity)
        )
        _temperature_difference_previous = (
            _temperature_initial
            - calculate_melting_temperature_from_salinity(_salinity_previous)
        )
        _enthalpy = _specific_heat_ice * _temperature + _latent_heat_water * _phi_km1
        _enthalpy_solid = (
            _specific_heat_ice * calculate_melting_temperature_from_salinity(_salinity)
        )
        _enthalpy_difference = _enthalpy - _enthalpy_solid
    else:
        _phi_km1 = _liquid_fraction_previous
        _temperature_difference = (
            _temperature_previous
            - calculate_melting_temperature_from_salinity(_salinity_previous)
        )
        _enthalpy = (
            _specific_heat_ice * _temperature_previous + _latent_heat_water * _phi_km1
        )
        _enthalpy_solid = (
            _specific_heat_ice
            * calculate_melting_temperature_from_salinity(_salinity_previous)
        )
        _enthalpy_difference = _enthalpy - _enthalpy_solid

    a_w_temperature, a_p_temperature, a_e_temperature = (
        t_k_A_LHS_matrix.diagonal(-1),
        t_k_A_LHS_matrix.diagonal(),
        t_k_A_LHS_matrix.diagonal(1),
    )

    if _is_stefan and not _nonconstant_physical_properties:
        _phi = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _phi_km1
                + under_relaxation_factor * _temperature_difference / temp_factor_3,
            ),
        )
    elif _is_stefan and _nonconstant_physical_properties:
        phi_difference = (
            _temperature_difference_previous * (1 + a_p_temperature) / temp_factor_3
        ) - (t_k_A_LHS_matrix @ (_temperature_difference / temp_factor_3))
        _phi = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _phi_km1 + phi_difference,
            ),
        )
    else:
        _phi = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _phi_km1
                + under_relaxation_factor * _temperature_difference / temp_factor_3,
            ),
        )
        assert np.all(
            (_phi >= 0) & (_phi <= 1)
        ), "liquid fraction has non-physical value"

    return _phi * phi_control_for_infinite_values(_phi)

def update_liquid_fraction_mixture_with_under_relaxation(
    _temperature,
    _salinity,
    _liquid_fraction_previous,
    under_relaxation_factor,
    _is_stefan=False,) :

    """
    Update the liquid fraction mixture using under-relaxation.
    Parameters:
        _temperature (float): The current temperature.
        _salinity (float): The salinity of the mixture.
        _liquid_fraction_previous (float): The previous liquid fraction.
        under_relaxation_factor (float): The under-relaxation factor to be applied.
        _is_stefan (bool, optional): Flag to indicate if the Stefan condition should be applied. Defaults to False.
    Returns:
        float: The updated liquid fraction mixture.
    """

    _temperature_liquidus , _temperature_solidus = calculate_liquidus_temperature(_salinity)
    _phi_km1 = _liquid_fraction_previous
    _temperature_effective = (_temperature_liquidus - _temperature_solidus)* _phi_km1*_rho_brine/_rho_ice + _temperature_solidus
    _temperature_difference = _temperature - _temperature_effective
    _specific_heat_effective = _specific_heat_ice * (1 - _phi_km1) + _specific_heat_brine * _phi_km1    # TODO: if not this use averaged specific heat capacities 
    _specific_heat_effective = (_specific_heat_ice + _specific_heat_brine)/2 * np.ones(len(_salinity))
    _specific_heat_effective_by_L = _specific_heat_effective / _latent_heat_water

    _enthalpy = _specific_heat_ice * _temperature + _latent_heat_water * _phi_km1
    _enthalpy_solid = (
        _specific_heat_ice * calculate_melting_temperature_from_salinity(_salinity)
    )

    # TODO: original phi condition
    if _is_stefan:
        _phi = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _phi_km1
                + under_relaxation_factor * _temperature_difference * _specific_heat_effective_by_L,
            ),
        )
    # TODO: modified with temperatures of previous step
    # if _is_stefan:
    #     _phi = np.where(
    #         _temperature <= _temperature_liquidus,
    #         0,
    #         np.where(
    #             _enthalpy > (_enthalpy_solid + _latent_heat_water),
    #             1,
    #             _phi_km1
    #             + under_relaxation_factor * _temperature_difference * _specific_heat_effective_by_L,
    #         ),
    #     )
    
    return _phi * phi_control_for_infinite_values(_phi), _temperature_liquidus, _temperature_solidus

def update_liquid_fraction_mixture_with_enthalpy_equation(
    _temperature,
    _salinity,
    _liquid_fraction_previous,
    _enthalpy,
    _enthalpy_solid,
    under_relaxation_factor,
    _is_stefan=False,) :

    """
    Update the liquid fraction mixture using under-relaxation.
    Parameters:
        _temperature (float): The current temperature.
        _salinity (float): The salinity of the mixture.
        _liquid_fraction_previous (float): The previous liquid fraction.
        under_relaxation_factor (float): The under-relaxation factor to be applied.
        _is_stefan (bool, optional): Flag to indicate if the Stefan condition should be applied. Defaults to False.
    Returns:
        float: The updated liquid fraction mixture.
    """

    _temperature_liquidus , _temperature_solidus = calculate_liquidus_temperature(np.ones(len(_salinity))*34.0)
    _phi_km1 = _liquid_fraction_previous
    _alpha = _rho_brine*_specific_heat_brine/(_rho_ice*_specific_heat_ice)
    if (_alpha - 1) == 0:
        _coef = 0
    else:
        _coef = 1/(_alpha - 1)
    
    _temperature_difference = _temperature_liquidus - _temperature_solidus
    _temperature_numerator = _temperature_difference
    _temperature_denominator = _temperature - _coef * _temperature_difference
    # _enthalpy = _specific_heat_ice * _temperature + _latent_heat_water * _phi_km1
    # _enthalpy_solid = (
    #     _specific_heat_ice * calculate_melting_temperature_from_salinity(_salinity)
    # )

    if _is_stefan:
        _phi = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _phi_km1 - _phi_km1 * _temperature_numerator / _temperature_denominator,
            ),
        )
    
    return _phi * phi_control_for_infinite_values(_phi)

def calculate_liquidus_temperature(_salinity, _liquid_relation="Normal"):
    S_br = 233.0  # brine salinity in ppt
    
    if _liquid_relation == "Normal":
        T_m = _temperature_melt # melt temperature as solidus temperature
        T_s = 252.05 # eutectic temperature for Sbr = 233ppt
    elif _liquid_relation == "Frezchem":
        T_m = _temperature_melt
        T_s = 252.05 # eutectic temperature for Sbr = 233ppt

    T_l = lambda S: (T_m + (T_s - T_m)/S_br*S)  # liquidus temperature at salinity S  # noqa: E731
    T_l_S = T_l(_salinity)  # liquidus temperature at salinity S

    return T_l_S, T_s*np.ones(len(_salinity))


def update_liquid_fraction_voller_continuous_thermal_properties(
    _temperature,
    _salinity,
    _phi,
    _enthalpy,
    _enthalpy_solid,
    under_relaxation_factor,
    temp_factor_3,
    _is_stefan=False,
    _method="likebuffo",
):
    """Updates the liquid fraction based on temperature, salinity, enthalpy, and other parameters.

    Args:
        _temperature (float): The temperature value.
        _salinity (float): The salinity value.
        _liquid_fraction (float): The liquid fraction value.
        _enthalpy (float): The enthalpy value.
        _enthalpy_solid (float): The solid enthalpy value.
        _nz (int): The number of vertical grid points.
        _is_stefan (bool, optional): Whether to use Stefan condition. Defaults to False.
        _method (str, optional): The method to use. Defaults to "likebuffo".
    Returns:
        tuple: A tuple containing the updated liquid fraction values.
    Raises:
        AssertionError: If the liquid fraction has a non-physical value.
    """

    _alpha = 1.853 / 28.0

    _temperature_difference = (
        _temperature - calculate_melting_temperature_from_salinity(_salinity)
    )
    if _is_stefan:
        _phi_kp1 = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _phi + under_relaxation_factor * _temperature_difference / (temp_factor_3),
            ),
        )
    else:
        _phi_kp1 = np.where(
            _enthalpy <= _enthalpy_solid,
            0,
            np.where(
                _enthalpy > (_enthalpy_solid + _latent_heat_water),
                1,
                _phi + under_relaxation_factor * _temperature_difference / (temp_factor_3),
            ),
        )
        assert np.all(
            (_phi_kp1 >= 0) & (_phi_kp1 <= 1)
        ), "liquid fraction has non-physical value"

    return _phi_kp1 * phi_control_for_infinite_values(_phi_kp1)


def update_enthalpy_solid_state(
    _salinity, _nz, _liq_rel="Normal", _temperature_melt=_temperature_melt
):
    """
    Updates the enthalpy in the solid state based on the given parameters.
    Args:
        _salinity (float): The salinity value.
        _nz (int): The nz value.
        _liq_rel (str, optional): The liquid relation. Defaults to "Normal".
        _temperature_melt (float, optional): The melting temperature. Defaults to _temperature_melt.
    Returns:
        float: The updated enthalpy in the solid state.
    """

    _temperature_fromsalinity = calculate_melting_temperature_from_salinity(
        _salinity, _liquid_relation=_liq_rel, _temperature_melt=_temperature_melt
    )
    return _specific_heat_ice * _temperature_fromsalinity


def update_enthalpy(
    _temperature, _salinity, _liquid_fraction, _nz, _method="likebuffo"
):
    """Updates the enthalpy based on the given parameters.

    Args:
        _temperature (float): The temperature value.
        _salinity (float): The salinity value.
        _liquid_fraction (float): The liquid fraction value.
        _nz (int): The nz value.
        _method (str, optional): The method used for calculating enthalpy. Defaults to "likebuffo".
    Returns:
        float: The updated enthalpy value.
    Raises:
        TypeError: If the given method is not available.
    """

    if _method not in ["likebuffo"]:
        raise TypeError("method for enthalpy not available")
    if _method == "likebuffo":
        _enthalpy = (
            _latent_heat_water * _liquid_fraction + _specific_heat_ice * _temperature
        )
    return _enthalpy

def calculate_brine_velocity_upwindscheme(phi, phi_initial, dz, dt, nz, rho_i, rho_br, thickness_index_prev):
        '''
        compute brine velocity (no gravity) according to upwind scheme

        Arguments------------------------------------------------------------------
            phi_initial initial liquid fraction, computed in previous time step [-]
            phi
            lastly updated liquid fraction in convergence scheme of
            the present time step [-]
            dz
            dt
            nz
            spatial discretization [m]
            time discretization [s]
            number of computational nodes
        Results---------------
            wkm1
            brine velocity updated [ms-1]
        '''
        w_km1 = np.zeros(nz)
        # initial velocity zero constant in space
        e = (rho_i/rho_br)- 1
        nu = e*dz/dt
        
        for i in range(1,nz,1):
            w_km1[i] = w_km1[i-1] + nu*(phi[i] - phi_initial[i])

        return w_km1

def calculate_brine_velocity_darcyscheme(phi, phi_initial, dz, dt, nz, rho_i, rho_br, thickness_index_prev):
        '''
        compute brine velocity (no gravity)

        Arguments------------------------------------------------------------------
            phi_initial initial liquid fraction, computed in previous time step [-]
            phi
            lastly updated liquid fraction in convergence scheme of
            the present time step [-]
            dz
            dt
            nz
            spatial discretization [m]
            time discretization [s]
            number of computational nodes
        Results---------------
            wkm1
            brine velocity updated [ms-1]
        '''
        dt_phi = np.zeros(nz)
        w_km1 = np.zeros(nz)
        # initial velocity zero constant in space
        dt_phi= (phi - phi_initial)/dt
        # time derivative of phi
        dt_phi_dz = dt_phi[:] * dz
        dt_phi_integrated = np.cumsum(dt_phi_dz)
        e = (rho_i/rho_br)- 1
        w_km1 = e * dt_phi_integrated  # * 0.25 #*phi_zero

        return w_km1

def update_temperature_and_salinity(
    preprocess_data_object: PreprocessData,
    t_prev: np.ndarray,
    s_prev: np.ndarray,
    phi_k: np.ndarray,
    brine_velocity_prev: np.ndarray,
    thickness_index_prev: float,
    t_initial: np.ndarray,
    s_initial: np.ndarray,
    phi_initial: np.ndarray,
    t_melt: np.ndarray,
    source_term_temperature: np.ndarray,
    source_term_salinity: np.ndarray,
    x_wind_temperature: np.ndarray,
    x_wind_salinity: np.ndarray,
    buffo: bool,
    stefan: bool,
    voller: bool,
    _is_salinity_equation: bool,
    _nonconstant_physical_properties: bool = False,
    t_k_A_LHS_matrix_prev: np.ndarray = None,
):
    """
    Update the temperature and salinity based on the given parameters.

    Args:
        preprocess_data_object (PreprocessData): The preprocess data object.
        t_prev (numpy.ndarray): The previous temperature values.
        s_prev (numpy.ndarray): The previous salinity values.
        phi_k (numpy.ndarray): The liquid fraction values.
        t_initial (numpy.ndarray): The initial temperature values.
        s_initial (numpy.ndarray): The initial salinity values.
        phi_initial (numpy.ndarray): The initial liquid fraction values.
        source_term (numpy.ndarray): The source term values.
        buffo (bool): The buffo flag.
        stefan (bool): The stefan flag.
        _is_salinity_equation (bool, optional): Whether to consider the salinity equation. Defaults to False.
    Returns:
        tuple: A tuple containing the updated temperature and salinity values.
        a_p: The main diagonal of the matrix A for the Voller scheme liquid fraction update

    """

    advection_diffusion_temp = AdvectionDiffusion(
        "temperature",
        t_prev,
        source_term_temperature,
        t_initial,
        phi_k,
        phi_initial,
        thickness_index_prev,
        brine_velocity_prev,
        preprocess_data_object.grid_timestep_dt,
        preprocess_data_object.grid_resolution_dz,
        preprocess_data_object.nz,
        preprocess_data_object.time_passed,
        preprocess_data_object.initial_salinity,
        Stefan=stefan,
        Buffo=buffo,
        Voller=voller,
        bc_neumann=preprocess_data_object.temp_grad,
    )
    t_k, x_wind_t, A_before_correction = advection_diffusion_temp.unknowns_matrix(
        t_melt, _nonconstant_physical_properties
    )
    
    if _is_salinity_equation is True:
        advection_diffusion_salinity = AdvectionDiffusion(
            "salinity",
            s_prev,
            source_term_salinity,
            s_initial,
            phi_k,
            phi_initial,
            thickness_index_prev,
            brine_velocity_prev,
            preprocess_data_object.grid_timestep_dt,
            preprocess_data_object.grid_resolution_dz,
            preprocess_data_object.nz,
            preprocess_data_object.time_passed,
            preprocess_data_object.initial_salinity,
            Stefan=stefan,
            Buffo=buffo,
            Voller=voller,
            bc_neumann=preprocess_data_object.temp_grad,
        )
        s_k, x_wind_s, S_before_correction = (
            advection_diffusion_salinity.unknowns_matrix(
                t_melt, _nonconstant_physical_properties
            )
        )
    else:
        s_k = s_prev
    return t_k, s_k, A_before_correction, advection_diffusion_temp.factor3, x_wind_t, x_wind_s

def update_algae_transport(    
    preprocess_data_object: PreprocessData,
    nutrient_cn_prev: np.ndarray,
    phi_k: np.ndarray,
    brine_velocity_prev: np.ndarray,
    thickness_index_prev: float,
    nutrient_cn_initial: np.ndarray,
    phi_initial: np.ndarray,
    t_melt: np.ndarray,
    source_term: np.ndarray,
    buffo: bool,
    stefan: bool,
    voller: bool,
    _nonconstant_physical_properties: bool = False,):

    advection_diffusion_algae = AdvectionDiffusion(
        "salinity",
        nutrient_cn_prev,
        source_term,
        nutrient_cn_initial,
        phi_k,
        phi_initial,
        thickness_index_prev,
        brine_velocity_prev,
        preprocess_data_object.grid_timestep_dt,
        preprocess_data_object.grid_resolution_dz,
        preprocess_data_object.nz,
        preprocess_data_object.time_passed,
        preprocess_data_object.initial_salinity,
        Stefan=stefan,
        Buffo=buffo,
        Voller=voller,
        bc_neumann=preprocess_data_object.temp_grad,
    )
    nutrient_cn, x_wind_cn, _ = (
        advection_diffusion_algae.unknowns_matrix(
            t_melt, _nonconstant_physical_properties
        )
    )

    return nutrient_cn


def update_state_variables(
    preprocess_data_object: PreprocessData,
    t_prev,
    s_prev,
    brine_velocity_prev,
    nutrient_cn_prev,
    phi_prev,
    thickness_index_prev,
    buffo,
    stefan,
    voller,
    t_initial,
    s_initial,
    nutrient_cn_initial,
    phi_initial,
    t_k_melt,
    t_k_A_LHS_matrix_prev,
    temp_factor_3,
    source_term_temperature,
    source_term_salinity,
    x_wind_temperature,
    x_wind_salinity,
    _is_salinity_equation,
):
    """
    Update the state variables Temperature, Salinity, Liquid Fraction, Enthalpy and Enthalpy of solid based on the given parameters.

    Args:
        preprocess_data_object (PreprocessData): The preprocess data object.
        t_km1 (numpy.ndarray): The previous temperature values.
        s_km1 (numpy.ndarray): The previous salinity values.
        nutrient_cn (numpy.ndarray): The previous nutrient concentration values
        phi_km1 (numpy.ndarray): The previous liquid fraction values.
        buffo (bool): The buffo flag.
        stefan (bool): The stefan flag.
        t_initial (numpy.ndarray): The initial temperature values.
        s_initial (numpy.ndarray): The initial salinity values.
        nutrient_cn_initial (numpy.ndarray): The initial nutrient concentration values
        phi_initial (numpy.ndarray): The initial liquid fraction values.
        source_term (numpy.ndarray): The source term values like radiation ocean flux, algae, etc.

    Methods called:
        Updates the state variables based on the given parameters.
            - update_enthalpy: Update the enthalpy of the system with previous temperature, salinity, liquid fraction.
            - update_enthalpy_solid_state: Update the enthalpy of the solid state with previous salinity using the liquidus relation.
            - update_liquid_fraction: Update the liquid fraction of the system with previous temperature, salinity, enthalpy, solid enthalpy.
            - update_temperature_and_salinity: Update the temperature and salinity of the system with updated liquid fraction and previous temperature, salinity, initial temperature, initial salinity, initial liquid fraction, source term, buffo, stefan.
            - update_algae_transport: Update the nutrient concentration of the system 
    Returns:
        tuple: A tuple containing the following updated values:
            - h_k (numpy.ndarray): The updated enthalpy values.
            - h_solid (numpy.ndarray): The updated solid enthalpy values.
            - phi_k (numpy.ndarray): The updated liquid fraction values.
            - t_k (numpy.ndarray): The updated temperature values.
            - s_k (numpy.ndarray): The updated salinity values.
            - nutrient_cn (numpy.ndarray): The updated nutrient values 
    """

    h_k = update_enthalpy(t_prev, s_prev, phi_prev, preprocess_data_object.nz)
    h_solid = update_enthalpy_solid_state(
        s_prev,
        preprocess_data_object.nz,
        preprocess_data_object.liquidus_relation_type,
    )
    if voller:
        t_k, s_k, t_k_A_LHS_matrix, temp_factor_3, x_wind_temperature, x_wind_salinity = update_temperature_and_salinity(
            preprocess_data_object,
            t_prev,
            s_prev,
            phi_prev,
            brine_velocity_prev,
            thickness_index_prev,
            t_initial,
            s_initial,
            phi_initial,
            t_k_melt,
            source_term_temperature,
            source_term_salinity,
            x_wind_temperature,
            x_wind_salinity,
            buffo,
            stefan,
            voller,
            _is_salinity_equation,
        )
        a_p_temperature = t_k_A_LHS_matrix.diagonal()
        phi_k = update_liquid_fraction_voller_continuous_thermal_properties(
            t_k,
            s_k,
            phi_prev,
            a_p_temperature,
            temp_factor_3,
            _is_stefan=preprocess_data_object.is_stefan,
        )
        t_k_melt = calculate_melting_temperature_from_salinity(s_k)
        temperature_solidus = 0.0
        temperature_liquidus = 0.0
    elif stefan:

        # phi update for a mixture using enthalpy equation
        # FIXME: This method is not working as expected. Requires debugging
        # phi_k = update_liquid_fraction_mixture_with_enthalpy_equation(
        #     t_prev,
        #     s_prev,
        #     phi_initial,
        #     h_k,
        #     h_solid,
        #     1.0,
        #     _is_stefan=preprocess_data_object.is_stefan,
        # )

        # phi update for a mixture using Faden method from paper "An optimum Enthalpy Approach for Melting and Solidification with Volume Change"
        phi_k, temperature_liquidus, temperature_solidus = update_liquid_fraction_mixture_with_under_relaxation(
            t_prev,
            s_prev,
            phi_prev,
            1.0,
            _is_stefan=preprocess_data_object.is_stefan,)  # returns an array of [phi, T_l, T_s]

        # # phi update for a mixture using Voller method with under-relaxation
        # phi_k = update_liquid_fraction_voller_under_relaxation(
        #     t_prev,
        #     s_prev,
        #     t_initial,
        #     phi_prev,
        #     t_prev,
        #     s_prev,
        #     preprocess_data_object.nz,
        #     temp_factor_3,
        #     t_k_A_LHS_matrix_prev,
        #     t_k_A_LHS_matrix_prev,
        #     1.4,
        #     _is_stefan=preprocess_data_object.is_stefan,
        #     _pt2_system=False,
        #     _nonconstant_physical_properties=False,
        # )
        
        t_k, s_k, t_k_A_LHS_matrix, temp_factor_3, x_wind_temperature, x_wind_salinity = update_temperature_and_salinity(
            preprocess_data_object,
            t_prev,
            s_prev,
            phi_k,
            brine_velocity_prev,
            thickness_index_prev,
            t_initial,
            s_initial,
            phi_initial,
            t_k_melt,
            source_term_temperature,
            source_term_salinity,
            x_wind_temperature,
            x_wind_salinity,
            buffo,
            stefan,
            voller,
            _is_salinity_equation,
            t_k_A_LHS_matrix_prev=t_k_A_LHS_matrix_prev,
        )

        nutrient_cn = update_algae_transport(
            preprocess_data_object,
            nutrient_cn_prev,
            phi_k,
            brine_velocity_prev,
            thickness_index_prev,
            nutrient_cn_initial,
            phi_initial,
            t_k_melt,
            source_term_salinity,
            buffo,
            stefan,
            voller,
        )

        # switch algae transport off
        # nutrient_cn = 0.0*nutrient_cn

        t_k_melt = calculate_melting_temperature_from_salinity(s_k)

    elif buffo:
        phi_k = update_liquid_fraction_buffo(
            t_prev,
            s_prev,
            phi_prev,
            h_k,
            h_solid,
            preprocess_data_object.nz,
            _is_stefan=preprocess_data_object.is_stefan,
            _method="likebuffo",
        )
        t_k, s_k, t_k_A_LHS_matrix, temp_factor_3, x_wind_temperature, x_wind_salinity = update_temperature_and_salinity(
            preprocess_data_object,
            t_prev,
            s_prev,
            phi_k,
            brine_velocity_prev,
            thickness_index_prev,
            t_initial,
            s_initial,
            phi_initial,
            t_k_melt,
            source_term_temperature,
            source_term_salinity,
            x_wind_temperature,
            x_wind_salinity,
            buffo,
            stefan,
            voller,
            _is_salinity_equation=_is_salinity_equation,
        )
        t_k_melt = calculate_melting_temperature_from_salinity(s_k)
        temperature_liquidus = 0.0
        temperature_solidus = 0.0
        nutrient_cn = 0.0
    else:
        AssertionError("No method selected for liquid fraction update")

    # Voller scheme: lower, main, upper diagonal of the matrix A of Ax = b

    # calculate brine velocity
    brine_velocity = calculate_brine_velocity_upwindscheme(phi_k, phi_initial, preprocess_data_object.grid_resolution_dz, preprocess_data_object.grid_timestep_dt, preprocess_data_object.nz, _rho_ice, _rho_brine, thickness_index_prev)

    # TODO: investigate values of T_k, s_k, brine_velocity, nutrient_cn, t_k_A_LHS_matrix, temp_factor_3, t_k_melt, temperature_liquidus, temperature_solidus, x_wind_temperature, x_wind_salinity
    # TODO: if T_k, s_k change then the initial values need to change 
    return h_k, h_solid, phi_k, t_k, s_k, brine_velocity, nutrient_cn, t_k_A_LHS_matrix, temp_factor_3, t_k_melt, temperature_liquidus, temperature_solidus, x_wind_temperature, x_wind_salinity


def phi_control_for_infinite_values(_phi):
    """Calculates the control values for infinite phi values.

    Args:
        _phi (numpy.ndarray): The input array of phi values.
    Returns:
        numpy.ndarray: The control values for the given phi values.
    """

    _phi_not_low = np.where(_phi > 0.000001, 1 / _phi, np.inf)
    _phi_control = np.ones_like(_phi, dtype=np.float64)
    _phi_control[np.isinf(_phi_not_low)] = 0
    return _phi_control


def H_function(_self, _temperature, _enthalpy_s1):
    """Calculates the value of H function.

    Args:
        _self (float): The value of self.
        _temperature (float): The temperature.
        _enthalpy_s1 (float): The enthalpy.
    Returns:
        float: The calculated value of H function.
    """

    _a = 0.0000059
    _b = 1 / _a
    return (
        _temperature * _specific_heat_ice
        + _latent_heat_water
        * 0.5
        * (np.tanh(_a * _self - (_enthalpy_s1 + _latent_heat_water / 2) / _b) + 1)
        - _self
    )


def H_function_derivate(_x, _enthalpy_s1):
    """Calculates the derivative of the H function.

    Args:
        _x (float): The value of x.
        _enthalpy_s1 (float): The value of enthalpy_s1.
    Returns:
        float: The derivative of the H function.
    """

    _a = 0.0000059
    _b = 1 / _a
    return (0.5 * _a * _latent_heat_water) / (
        np.cosh(_a * _x - (_latent_heat_water / 2 + _enthalpy_s1) / _b) ** 2
    ) - 1


def H_newton_iteration(_temperature, _enthalpy_s1):
    """Performs Newton iteration to find the root of the H_function.

    Args:
        _temperature (float): The temperature value.
        _enthalpy_s1 (float): The enthalpy value.
    Returns:
        float: The root of the H_function.
    Raises:
        None
    """

    _x = 1
    return newton(
        H_function,
        x0=_x,
        fprime=H_function_derivate,
        args=(_temperature, _enthalpy_s1),
        maxiter=100,
    )


def phi_func(_enthalpy_k1, _enthalpy_s1):
    """Calculates the phi value based on the given enthalpy values.

    Args:
        _enthalpy_k1 (float): The enthalpy value for k1.
        _enthalpy_s1 (float): The enthalpy value for s1.
    Returns:
        float: The calculated phi value.
    """

    _a = 0.0000059
    _b = 1 / _a
    return 0.5 * (
        np.tanh(_a * _enthalpy_k1 - (_enthalpy_s1 + _latent_heat_water / 2) / _b) + 1
    )


def _H_update(_temperature, _enthalpy_solid, _nz):
    """Updates the enthalpy values using Newton's iteration method.

    Args:
        _temperature (list): A list of temperatures.
        _enthalpy_solid (list): A list of solid enthalpy values.
        _nz (int): The number of elements.
    Returns:
        numpy.ndarray: An array of updated enthalpy values.
    """

    return np.array(
        [
            H_newton_iteration(float(_temperature[k]), _enthalpy_solid[k])
            for k in range(_nz)
        ]
    )


def _H_newton_manual(_temperature, _enthalpy_s1):
    """Calculates the new value of _x using Newton's method.

    Args:
        _temperature (float): The temperature value.
        _enthalpy_s1 (float): The enthalpy value.
    Returns:
        float: The new value of _x.
    Raises:
        None
    """
    # code implementation
    ...

    _x = _temperature * _specific_heat_ice + _latent_heat_water
    _x_new = _x + 1000
    while np.absolute(_x - _x_new) > 0.0003:
        _x = _x_new
        _fx = (
            _temperature * _specific_heat_ice
            + _latent_heat_water
            * np.piecewise(
                _x,
                [
                    _x < _enthalpy_s1,
                    _enthalpy_s1 + _latent_heat_water >= _x >= _enthalpy_s1,
                    _x > _enthalpy_s1 + _latent_heat_water,
                ],
                [0, lambda _x: (_x - _enthalpy_s1) / _latent_heat_water, 1],
            )
            - _x
        )
        _dfxdx = np.piecewise(
            _x,
            [
                _x < _enthalpy_s1,
                _enthalpy_s1 + _latent_heat_water >= _x >= _enthalpy_s1,
                _x > _enthalpy_s1 + _latent_heat_water,
            ],
            [-1, 0, -1],
        )
        _x_new = _x - _fx / _dfxdx
    return _x_new


def _H_bisect_search(_temperature, _enthalpy_s1):
    """Performs a bisection search to find the root of the H_function within the given temperature range.

    Args:
        _temperature (float): The temperature value.
        _enthalpy_s1 (float): The enthalpy value.
    Returns:
        float: The root of the H_function within the given temperature range.
    """
    # code implementation

    _x1 = _temperature * _specific_heat_ice + _latent_heat_water * 0
    _x2 = _temperature * _specific_heat_ice + _latent_heat_water * 1
    try:
        _root = bisect(
            H_function, a=_x1, b=_x2, args=(_temperature, _enthalpy_s1), xtol=0.1
        )
    except RuntimeError:
        _low = H_function(_x1, _temperature, _enthalpy_s1)
        _high = H_function(_x2, _temperature, _enthalpy_s1)
        if _low < _high:
            _root = _x1
        elif _high < _low:
            _root = _x2
    print(_root)
