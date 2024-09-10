from __future__ import annotations

import numpy as np
from scipy.optimize import bisect, newton
from typing import TYPE_CHECKING
from parameters.user_input import UserInput
from models.advection_diffusion import AdvectionDiffusion

if TYPE_CHECKING:
    from preprocess.pre_process import PreprocessData

_ui = UserInput()
_temperature_melt, _boundary_top_temperature = (
    _ui.temperature_melt,
    _ui.boundary_top_temperature,
)
_specific_heat_ice, _latent_heat_water = _ui.constants.c_i, _ui.constants.L


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
        _melting_temperature_seawater = _temperature_melt * np.ones(_salinity.shape)
    elif _liquid_relation == "Frezchem":
        _melting_temperature_seawater = _temperature_melt - (
            -(9.1969758 * (1e-05) * _salinity**2)
            - 0.03942059 * _salinity
            + 272.63617665
        )
    return _melting_temperature_seawater


def update_liquid_fraction(
    _temperature,
    _salinity,
    _enthalpy,
    _enthalpy_solid,
    _nz,
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
    _phi = np.ones(_nz, dtype=np.float64)
    _alpha = 1.853 / 28.0

    _temperature_difference = (
        _temperature - calculate_melting_temperature_from_salinity(_salinity)
    )
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

    return _phi * phi_control_for_infinite_values(_phi)


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


def update_temperature_and_salinity(
    preprocess_data_object: PreprocessData,
    t_prev: np.ndarray,
    s_prev: np.ndarray,
    phi_k: np.ndarray,
    t_initial: np.ndarray,
    s_initial: np.ndarray,
    phi_initial: np.ndarray,
    source_term: np.ndarray,
    buffo: bool,
    stefan: bool,
    _is_salinity_equation: bool = False,
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

    """

    advection_diffusion_temp = AdvectionDiffusion(
        "temperature",
        t_prev,
        source_term,
        t_initial,
        phi_k,
        phi_initial,
        preprocess_data_object.upwind_velocity,
        preprocess_data_object.grid_timestep_dt,
        preprocess_data_object.grid_resolution_dz,
        preprocess_data_object.nz,
        preprocess_data_object.time_passed,
        preprocess_data_object.initial_salinity,
        Stefan=stefan,
        Buffo=buffo,
        bc_neumann=preprocess_data_object.temp_grad,
    )
    t_k, x_wind_t = advection_diffusion_temp.unknowns_matrix()
    # TODO: Add algae growth model here to salinity
    if _is_salinity_equation:
        advection_diffusion_salinity = AdvectionDiffusion(
            "salinity",
            s_prev,
            source_term,
            s_initial,
            phi_k,
            phi_initial,
            preprocess_data_object.upwind_velocity,
            preprocess_data_object.grid_timestep_dt,
            preprocess_data_object.grid_resolution_dz,
            preprocess_data_object.nz,
            preprocess_data_object.time_passed,
            preprocess_data_object.initial_salinity,
            Stefan=stefan,
            Buffo=buffo,
            bc_neumann=preprocess_data_object.temp_grad,
        )
        s_k, x_wind_s = advection_diffusion_salinity.unknowns_matrix()
    else:
        s_k = s_prev
    return t_k, s_k


def update_state_variables(
    preprocess_data_object: PreprocessData,
    t_prev,
    s_prev,
    phi_prev,
    buffo,
    stefan,
    t_initial,
    s_initial,
    phi_initial,
    source_term,
    _is_salinity_equation=False,
):
    """
    Update the state variables Temperature, Salinity, Liquid Fraction, Enthalpy and Enthalpy of solid based on the given parameters.

    Args:
        preprocess_data_object (PreprocessData): The preprocess data object.
        t_km1 (numpy.ndarray): The previous temperature values.
        s_km1 (numpy.ndarray): The previous salinity values.
        phi_km1 (numpy.ndarray): The previous liquid fraction values.
        buffo (bool): The buffo flag.
        stefan (bool): The stefan flag.
        t_initial (numpy.ndarray): The initial temperature values.
        s_initial (numpy.ndarray): The initial salinity values.
        phi_initial (numpy.ndarray): The initial liquid fraction values.
        source_term (numpy.ndarray): The source term values.

    Methods called:
        Updates the state variables based on the given parameters.
            - update_enthalpy: Update the enthalpy of the system with previous temperature, salinity, liquid fraction.
            - update_enthalpy_solid_state: Update the enthalpy of the solid state with previous salinity using the liquidus relation.
            - update_liquid_fraction: Update the liquid fraction of the system with previous temperature, salinity, enthalpy, solid enthalpy.
            - update_temperature_and_salinity: Update the temperature and salinity of the system with updated liquid fraction and previous temperature, salinity, initial temperature, initial salinity, initial liquid fraction, source term, buffo, stefan.
    Returns:
        tuple: A tuple containing the following updated values:
            - h_k (numpy.ndarray): The updated enthalpy values.
            - h_solid (numpy.ndarray): The updated solid enthalpy values.
            - phi_k (numpy.ndarray): The updated liquid fraction values.
            - t_k (numpy.ndarray): The updated temperature values.
            - s_k (numpy.ndarray): The updated salinity values.
    """

    h_k = update_enthalpy(t_prev, s_prev, phi_prev, preprocess_data_object.nz)
    h_solid = update_enthalpy_solid_state(
        s_prev,
        preprocess_data_object.nz,
        preprocess_data_object.liquidus_relation_type,
    )
    phi_k = update_liquid_fraction(
        t_prev,
        s_prev,
        h_k,
        h_solid,
        preprocess_data_object.nz,
        _is_stefan=preprocess_data_object.is_stefan,
    )
    t_k, s_k = update_temperature_and_salinity(
        preprocess_data_object,
        t_prev,
        s_prev,
        phi_k,
        t_initial,
        s_initial,
        phi_initial,
        source_term,
        buffo,
        stefan,
        _is_salinity_equation=_is_salinity_equation,
    )

    return h_k, h_solid, phi_k, t_k, s_k


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
