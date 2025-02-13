from __future__ import annotations

import re

import numpy as np

from src.spyice.parameters.user_input import UserInput

ui = UserInput()
(
    temperature_melt,
    boundary_top_temperature,
    boundary_salinity,
    boundary_top_temperature,
    dz,
    phi_c,
) = (
    ui.temperature_melt,
    ui.boundary_top_temperature,
    ui.boundary_salinity,
    ui.boundary_top_temperature,
    ui.grid_resolution_dz,
    ui.critical_liquid_fraction,
)


class SalinityUnavailableError(Exception):
    """Exception raised when the S_IC option is not available in initial conditions."""


def temperature_gradient(phi, nz):
    """Calculates the temperature gradient based on the given potential temperature profile and number of vertical levels.

    Args:
        phi (list): The potential temperature profile.
        nz (int): The number of vertical levels.

    Returns:
        float: The calculated temperature gradient.

    Raises:
        None

    """
    critical_depth = next((i for i in range(nz - 1) if phi[i] <= phi_c), 0)
    critical_depth = critical_depth + (critical_depth + 1) * dz
    return dz * (temperature_melt - boundary_top_temperature) / (critical_depth)


def set_inital_temperature(
    initial_temperature, nz, boundary_salinity, boundary_top_temperature
):
    """Sets the initial temperature profile based on the given parameters.

    Args:
        initial_temperature (str): The type of initial temperature profile to set.
        nz (int): The number of vertical grid points.
        boundary_salinity (float): The salinity at the boundary.
        boundary_top_temperature (float): The temperature at the top boundary.

    Returns:
        numpy.ndarray: The initial temperature profile as a 1D numpy array.

    Raises:
        None
    """

    temperature = np.zeros(nz, dtype=np.float64)
    temperature_melt = compute_melting_temperature_from_salinity(boundary_salinity)
    if (
        initial_temperature == "T(S)"
        or initial_temperature != "T_Stefan"
        and initial_temperature != "T271.25"
        and initial_temperature != "T250"
        and initial_temperature == "Tm_w"
    ):
        temperature = np.ones(nz, dtype=np.float64) * temperature_melt
    elif initial_temperature == "T_Stefan":
        temperature = np.ones(nz, dtype=np.float64) * boundary_top_temperature
    elif initial_temperature == "T271.25":
        temperature = np.ones(nz, dtype=np.float64) * 271.25
    elif initial_temperature == "T250":
        temperature = np.ones(nz, dtype=np.float64) * 250.0
    return temperature


def set_inital_salinity(initial_salinity, nz, boundary_salinity):
    """Sets the initial salinity values for each layer in the model.

    Args:
        initial_salinity (str): The type of initial salinity distribution.
        nz (int): The number of layers in the model.
        boundary_salinity (float): The salinity value at the boundary.
    Returns:
        numpy.ndarray: An array of initial salinity values for each layer.
    Raises:
        SalinityException: If the initial_salinity value is not recognized.
    """

    # sourcery skip: lift-return-into-if, remove-unnecessary-else
    salinity = np.zeros(nz, dtype=np.float64)
    if initial_salinity == "S_linear":
        salinity = np.linspace(200, boundary_salinity, nz)
    else:
        salinity_value = re.findall("[0-9]+$", initial_salinity)
        raise_salinity_exception(salinity_value)
        salinity = np.ones(nz, dtype=np.float64) * float(salinity_value[0])

    return salinity


def set_initial_liquidfraction(initial_liquid_fraction, nz):
    """Sets the initial liquid fraction based on the given input.

    Args:
        initial_liquid_fraction (str): The initial liquid fraction type.
        nz (int): The number of grid points.
    Returns:
        numpy.ndarray: The array representing the initial liquid fraction.
    Raises:
        None
    Examples:
        >>> set_initial_liquidfraction("P1", 10)
        array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
        >>> set_initial_liquidfraction("P_Stefan", 5)
        array([0., 0., 0., 0., 0.])
        >>> set_initial_liquidfraction("P0", 8)
        array([0., 0., 0., 0., 0., 0., 0., 0.])
    """

    phi = np.zeros(nz, dtype=np.float64)
    if initial_liquid_fraction == "P1":
        phi = np.ones(nz, dtype=np.float64) * 1
        phi[0] = 0.0
    elif initial_liquid_fraction in ["P_Stefan", "P0"]:
        phi = np.zeros(nz, dtype=np.float64)
        phi[0] = 1.0

    return phi


def set_initial_conditions(
    nz,
    boundary_salinity,
    initial_temperature="T0",
    initial_salinity="S1",
    initial_liquid_fraction="P1",
    boundary_top_temperature=265.0,
):
    """Sets the initial conditions for the simulation.

    Args:
        nz (int): Number of vertical grid points.
        boundary_salinity (float): Salinity value at the boundary.
        initial_temperature (str, optional): Initial temperature profile. Defaults to "T0".
        initial_salinity (str, optional): Initial salinity profile. Defaults to "S1".
        initial_liquid_fraction (str, optional): Initial liquid fraction profile. Defaults to "P1".
        boundary_top_temperature (float, optional): Temperature value at the top boundary. Defaults to 265.0.
    Returns:
        tuple: A tuple containing the temperature, salinity, liquid fraction, and upwind velocity arrays.
    """

    salinity = set_inital_salinity(initial_salinity, nz, boundary_salinity)
    temperature = set_inital_temperature(
        initial_temperature, nz, boundary_salinity, boundary_top_temperature
    )
    liquid_fraction = set_initial_liquidfraction(initial_liquid_fraction, nz)
    upwind_velocity = np.zeros(nz, dtype=np.float64)

    return temperature, salinity, liquid_fraction, upwind_velocity


def t_w3(dt):
    """Calculates the top boundary temperature and freeze date.

    Args:
        dt (float): The time step.
    Returns:
        tuple: A tuple containing the top boundary temperature array and the freeze date.
    """
    # code implementation...

    time_temp = np.arange(-2000000, 315576000, dt)
    boundary_top_temperature = (
        -9 * np.sin((2 * np.pi / 31557600) * time_temp) - 18.15 + 273.15
    )
    for i in range(len(boundary_top_temperature)):
        boundary_top_temperature[i] = max(boundary_top_temperature[i], 250)
    freeze_date = (86400 * 149) / dt

    return boundary_top_temperature, freeze_date


def set_boundary_temperature(t_passed, temperature_bottom, **kwargs):
    """Sets the boundary temperature based on the given parameters.

    Args:
        t_passed (float): The time passed.
        temperature_bottom (float): The bottom temperature.
        **kwargs: Additional keyword arguments.
    Returns:
        tuple: A tuple containing the top temperature and the bottom temperature.
    """

    top_temp = kwargs.get("top_temp")
    if top_temp == "T_const_250":
        temperature_top = 250.0
    if top_temp == "Stefan":
        temperature_top = boundary_top_temperature
        temperature_bottom = temperature_bottom
    elif top_temp == "T_const_260":
        temperature_top = 260.0
    elif top_temp == "T_const_265":
        temperature_top = 265.0
    elif top_temp == "T_W3":
        freeze_date = 86400 * 149
        temperature_top = (
            -9 * np.sin((2 * np.pi / 31557600) * (t_passed + freeze_date))
            - 18.15
            + 273.15
        )
        temperature_top = max(temperature_top, 250)

    return temperature_top, temperature_bottom


def compute_melting_temperature_from_salinity(initial_salinity):
    """Computes the melting temperature from the given initial salinity.

    Args:
        initial_salinity (float): The initial salinity value.
    Returns:
        float: The computed melting temperature.
    """

    return temperature_melt


def boundary_condition(argument, t_passed, initial_salinity, **kwargs):
    """Calculates the boundary conditions for temperature or salinity.

    Args:
        argument (str): The argument specifying whether to calculate temperature or salinity.
        t_passed (float): The time passed.
        initial_salinity (float): The initial salinity value.
        **kwargs: Additional keyword arguments.
    Returns:
        tuple: A tuple containing the boundary conditions for temperature or salinity.
    Raises:
        None
    """

    if argument == "temperature":
        temperature_bottom, temperature_top = calculate_boundary_temperature(
            t_passed, initial_salinity, kwargs
        )
        return temperature_top, temperature_bottom

    elif argument == "salinity":
        salinity_bottom, salinity_top = calculate_boundary_salinity(initial_salinity)
        return salinity_top, salinity_bottom
    return None


def calculate_boundary_salinity(initial_salinity):
    """Calculates the boundary salinity values based on the initial salinity.

    Args:
        initial_salinity (str): The initial salinity value.
    Returns:
        tuple: A tuple containing the bottom and top salinity values.
    Raises:
        SalinityException: If the initial salinity value is not valid.
    """
    # Rest of the code...

    if initial_salinity == "S_linear":
        salinity_bottom = 34.0
        salinity_top = 100.0
    else:
        salinity_value = re.findall("[0-9]+$", initial_salinity)
        raise_salinity_exception(salinity_value)
        salinity_bottom = float(salinity_value[0])
        salinity_top = float(salinity_value[0])
    return salinity_bottom, salinity_top


def calculate_boundary_temperature(t_passed, initial_salinity, kwargs):
    """Calculates the boundary temperature based on the given parameters.

    Args:
        t_passed (float): The time passed.
        initial_salinity (str): The initial salinity value.
        kwargs (dict): Additional keyword arguments.
    Returns:
        tuple: A tuple containing the bottom temperature and the top temperature.
    """
    # Rest of the code...

    salinity_value = re.findall("[0-9]+$", initial_salinity)
    if len(salinity_value) == 0:
        temperature_bottom = 271.25
    else:
        temperature_bottom = compute_melting_temperature_from_salinity(
            float(salinity_value[0])
        )

    temperature_top, temperature_bottom = set_boundary_temperature(
        t_passed, temperature_bottom, **kwargs
    )
    return temperature_bottom, temperature_top


def raise_salinity_exception(salinity_value):
    """Raises a custom exception if the salinity value is empty.

    Args:
        salinity_value (str): The salinity value.
    Raises:
        SalinityUnavailableError: If the salinity value is empty.
    """

    if len(salinity_value) == 0:
        # Replace the selected code with the custom exception
        msg = "S_IC option not available in initial condition"
        raise SalinityUnavailableError(msg)
