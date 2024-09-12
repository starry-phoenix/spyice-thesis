from __future__ import annotations

import numpy as np

from src.spyice.parameters.user_input import UserInput
from src.spyice.preprocess.initial_boundary_conditions import boundary_condition

ui = UserInput()
temperature_melt, temperature_top = ui.temperature_melt, ui.boundary_top_temperature


def apply_boundary_condition(
    argument,
    x_initial,
    source,
    factor1,
    factor3,
    a,
    delta_upwind,
    w,
    nz,
    t_passed,
    salinity_initial,
    _temperature_top,
    is_stefan,
    is_buffo=False,
    bc_neumann=None,
):
    """Creates the right hand side of the matrix equation considering source terms.

    Args:
        argument (str): Either 'salinity' for salt equation or 'temperature' for temperature equation.
        x_initial (float): Value of X at the last time step.
        source (float): Source term.
        factor1 (float): Factor 1.
        factor3 (float): Factor 3.
        a (float): A parameter.
        delta_upwind (float): Difference of ice volume fraction between this and the last time step.
        w (float): W parameter.
        nz (int): Number of computational nodes.
        t_passed (float): Time passed in seconds.
        salinity_initial (float): Initial salinity value.
        _temperature_top (float): Top temperature value.
        is_stefan (bool): Indicates if Stefan condition is used.
        is_buffo (bool, optional): Indicates if Buffo condition is used. Defaults to False.
        bc_neumann (float, optional): Neumann boundary condition. Defaults to None.
        float: The right hand side of the equation.
    """
    rhs_matrix = x_initial - factor3 * delta_upwind

    if argument == "salinity":
        salinity_bc_top, salinity_bc_bottom = boundary_condition(
            argument, t_passed, salinity_initial
        )
        rhs_matrix[-1] += factor1[-1] * salinity_bc_bottom  # Dirichlet

    elif argument == "temperature":
        temperature_bc_top, temperature_bc_bottom = boundary_condition(
            argument, t_passed, salinity_initial, top_temp=_temperature_top
        )

        if is_stefan:
            if bc_neumann is not None:
                rhs_matrix[0] -= 2 * factor1[0] * bc_neumann
            else:
                rhs_matrix[0] = temperature_top
            rhs_matrix[-1] = temperature_melt  # Dirichlet
        elif is_buffo:
            if bc_neumann is not None:
                rhs_matrix[0] -= 2 * factor1[0] * bc_neumann
            else:
                rhs_matrix[0] = temperature_top
            rhs_matrix[-1] += factor1[-1] * temperature_melt
        else:
            rhs_matrix[0] += factor1[0] * temperature_top
            rhs_matrix[-1] += factor1[-1] * temperature_melt
    else:
        print("argument for BC not available")

    return rhs_matrix


# not relevant for our topic
def correct_for_brine_movement(
    argument, x_initial, w, t_passed, nz, salinity_initial, top_temp
):
    """Corrects for brine movement based on the given arguments.

    Args:
        argument (str): The argument for correction, either "salinity" or "temperature".
        x_initial (numpy.ndarray): The initial values of x.
        w (numpy.ndarray): The values of w.
        t_passed (float): The time passed.
        nz (int): The number of elements.
        salinity_initial (float): The initial salinity value.
        top_temp (float): The top temperature value.
    Returns:
        numpy.ndarray: The corrected values of x.
    Raises:
        None
    """

    x_upwind = np.zeros(nz)
    for i in range(1, nz - 1):
        if w[i] != 0:
            sign_w_i = np.sign(w[i])
            w_neg = ((1 - abs(sign_w_i)) / (2 * sign_w_i)) * (
                x_initial[i + 1] - x_initial[i]
            )  # w < 0
            w_pos = ((abs(sign_w_i) + 1) / (2 * sign_w_i)) * (
                x_initial[i] - x_initial[i - 1]
            )  # w > 0
            x_upwind[i] = w_neg + w_pos
        else:
            x_upwind[i] = 0

    if argument == "salinity":
        salinity_bc_bottom, _ = boundary_condition(argument, t_passed, salinity_initial)
        x_upwind[0] = 0 if w[0] >= 0 else x_initial[1] - x_initial[0]
        if w[-1] > 0:
            x_upwind[-1] = x_initial[-1] - x_initial[-2]
        elif w[-1] < 0:
            x_upwind[-1] = salinity_bc_bottom - x_initial[-1]
        else:
            x_upwind[-1] = 0

    elif argument == "temperature":
        temperature_bc_top, _ = boundary_condition(
            argument, t_passed, salinity_initial, top_temp=top_temp
        )
        x_upwind[0] = (
            x_initial[0] - temperature_bc_top
            if w[0] > 0
            else (x_initial[1] - x_initial[0] if w[0] < 0 else 0)
        )
        if w[-1] > 0:
            x_upwind[-1] = x_initial[-1] - x_initial[-2]
        elif w[-1] < 0:
            x_upwind[-1] = temperature_bc_top - x_initial[-1]
        else:
            x_upwind[-1] = 0
    else:
        print("argument for BC not available")
    return x_upwind
