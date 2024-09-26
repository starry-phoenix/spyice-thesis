import numpy as np


# Function to set state variables
def set_statevariables(
    temperature_calculated, salinity_calculated, liquid_fraction_calculated
):
    """Set the state variables for temperature, salinity, and liquid fraction.

    Args:
        temperature_calculated (float): The calculated temperature.
        salinity_calculated (float): The calculated salinity.
        liquid_fraction_calculated (float): The calculated liquid fraction.
    Returns:
        tuple: A tuple containing the following state variables:
            - temperature_initial (float): The initial temperature.
            - temperature_previous (float): The previous temperature.
            - salinity_initial (float): The initial salinity.
            - salinity_previous (float): The previous salinity.
            - phi_initial (float): The initial liquid fraction.
            - phi_prev (float): The previous liquid fraction.

    """

    temperature_previous, temperature_initial = (
        temperature_calculated,
        temperature_calculated,
    )  # Set previous T as current T for initialization, # Store initial T for later use
    salinity_previous, salinity_initial = (
        salinity_calculated,
        salinity_calculated,
    )  # Set previous S as current S for initialization, # Store initial S for later use
    phi_prev, phi_initial = (
        liquid_fraction_calculated,
        liquid_fraction_calculated,
    )  # Set previous phi as current phi for initialization, # Store initial phi for later use

    return (
        temperature_initial,
        temperature_previous,
        salinity_initial,
        salinity_previous,
        phi_initial,
        phi_prev,
    )


# Function to overwrite state variables
def overwrite_statevariables(
    temperature_calculated,
    salinity_calculated,
    liquid_fraction_calculated,
    a_p_temperature=None,
    temp_factor3=None,
):
    """Overwrites the state variables with the calculated values.
    Args:
        temperature_calculated (float): The calculated temperature.
        salinity_calculated (float): The calculated salinity.
        liquid_fraction_calculated (float): The calculated liquid fraction.
    Returns:
        tuple: A tuple containing the following state variables:
            - temperature_new (float): The current temperature for iteration.
            - salinity_new (float): The current salinity for iteration.
            - liquid_fraction_new (float): The current liquid fraction for iteration.
    """

    temperature_new = (
        temperature_calculated  # Set current T as previous T for next iteration
    )
    salinity_new = salinity_calculated  # Set current S as previous S for next iteration
    liquid_fraction_new = (
        liquid_fraction_calculated  # Set current phi as previous phi for next iteration
    )
    if a_p_temperature is not None:
        a_p_temperature = a_p_temperature
        temp_factor3 = temp_factor3
    else:
        a_p_temperature = np.ones(len(temperature_new)) * 0.0
        temp_factor3 = np.ones(len(temperature_new))
    return (
        temperature_new,
        salinity_new,
        liquid_fraction_new,
        a_p_temperature,
        temp_factor3,
    )


# Function to compute errors for convergence
def compute_error_for_convergence(
    temperature_calculated,
    temperature_previous,
    salinity_calculated,
    salinity_previous,
    liquid_fraction_calculated,
    liquid_fraction_previous,
    voller=False,
    **kwargs,
):
    """Computes the errors for convergence between the calculated and previous values of temperature, salinity, and liquid fraction.

    Args:
        temperature_calculated (numpy.ndarray): Array of calculated temperature values.
        temperature_previous (numpy.ndarray): Array of previous temperature values.
        salinity_calculated (numpy.ndarray): Array of calculated salinity values.
        salinity_previous (numpy.ndarray): Array of previous salinity values.
        liquid_fraction_calculated (numpy.ndarray): Array of calculated liquid fraction values.
        liquid_fraction_previous (numpy.ndarray): Array of previous liquid fraction values.
    Returns:
        tuple: A tuple containing the following error values:
            - temperature_error_max (float): Maximum temperature error for convergence check.
            - temperature_error_all (numpy.ndarray): Full temperature error for convergence check.
            - salinity_error_max (float): Maximum salinity error for convergence check.
            - salinity_err_all (numpy.ndarray): Full salinity error for convergence check.
            - liquid_fraction_error_max (float): Maximum liquid fraction error for convergence check.
            - liquid_fraction_error_all (numpy.ndarray): Full liquid fraction error for convergence check.
    """
    # Function implementation goes here

    temperature_error_max = np.max(
        abs(temperature_calculated[1:-1] - temperature_previous[1:-1])
    )  # Compute maximum T error for convergence check
    temperature_error_all = abs(
        temperature_calculated - temperature_previous
    )  # Compute full T error for convergence check
    salinity_error_max = np.max(
        salinity_calculated[1:-1] - salinity_previous[1:-1]
    )  # Compute maximum S error for convergence check
    salinity_err_all = (
        salinity_calculated - salinity_previous
    )  # Compute full S error for convergence check
    liquid_fraction_error_max = np.max(
        liquid_fraction_calculated[1:-1] - liquid_fraction_previous[1:-1]
    )  # Compute maximum phi error for convergence check
    liquid_fraction_error_all = (
        liquid_fraction_calculated - liquid_fraction_previous
    )  # Compute full phi error for convergence check
    residual_sum = voller_residual_scheme(
        temperature_calculated,
        temperature_previous,
        liquid_fraction_calculated,
        liquid_fraction_previous,
        kwargs,
    )

    return (
        temperature_error_max,
        temperature_error_all,
        salinity_error_max,
        salinity_err_all,
        liquid_fraction_error_max,
        liquid_fraction_error_all,
        residual_sum,
    )


def voller_residual_scheme(
    temperature_calculated,
    temperature_previous,
    liquid_fraction_calculated,
    liquid_fraction_previous,
    kwargs,
):
    a_matrix = kwargs["A_matrix"]
    a_matrix[0], a_matrix[-1] = 0, 0
    a_res = a_matrix @ temperature_calculated
    mush_cond = (liquid_fraction_calculated >= 0.05) & (
        liquid_fraction_calculated <= 0.95
    )
    mush_indx = np.where(mush_cond)
    phi_diff = kwargs["phi_initial"] - liquid_fraction_calculated
    rhs_res = kwargs["t_initial"] + kwargs["temp_factor3"] * phi_diff
    residual = a_res[1:-1] - rhs_res[1:-1]

    residual_sum = np.sum(np.abs(residual))

    residual_sum = np.sum(np.abs(residual))
    return residual_sum


# Function to reset errors for while loop
def reset_error_for_while_loop(
    temperature_tolerance, salinity_tolerance, liquid_fraction_tolerance
):
    """Resets the error values for a while loop based on the given tolerances.

    Args:
        temperature_tolerance (float): The tolerance for temperature error.
        salinity_tolerance (float): The tolerance for salinity error.
        liquid_fraction_tolerance (float): The tolerance for liquid fraction error.
    Returns:
        tuple: A tuple containing the reset error values for temperature, salinity, and liquid fraction.
    """

    temperature_error = (
        1 + temperature_tolerance
    )  # Set initial T error to value greater than tolerance
    salinity_error = (
        1 + salinity_tolerance
    )  # Set initial S error to value greater than tolerance
    liquid_fraction_error = (
        1 + liquid_fraction_tolerance
    )  # Set initial phi error to value greater than tolerance
    return temperature_error, salinity_error, liquid_fraction_error


# Function to initialize state variables
def set_initial_statevariables(temperature, salinity, liquid_fraction):
    """Initializes the state variables for the given temperature, salinity, and liquid fraction.

    Args:
        temperature (float): The initial temperature.
        salinity (float): The initial salinity.
        liquid_fraction (float): The initial liquid fraction.
    Returns:
        tuple: A tuple containing the initialized state variables:
            - temperature_initial (float): The initial temperature.
            - temperature_new (float): The current temperature for iteration.
            - temperature_previous (float): The previous temperature for initialization.
            - salinity_initial (float): The initial salinity.
            - salinity_new (float): The current salinity for iteration.
            - salinity_previous (float): The previous salinity for initialization.
            - phi_initial (float): The initial liquid fraction.
            - liquid_fraction_new (float): The current liquid fraction for iteration.
            - phi_prev (float): The previous liquid fraction for initialization.
    """

    phi_initial = liquid_fraction  # Store initial phi for later use
    liquid_fraction_new = (
        liquid_fraction  # Set previous phi as current phi for iteration
    )
    phi_prev = liquid_fraction  # Set previous phi as current phi for initialization
    temperature_initial = temperature  # Store initial T for later use
    temperature_previous = temperature  # Set previous T as current T for initialization
    temperature_new = temperature  # Set current T as previous T for iteration
    salinity_initial = salinity  # Store initial S for later use
    salinity_new = salinity  # Set current S as previous S for initialization
    salinity_previous = salinity
    return (
        temperature_initial,
        temperature_new,
        temperature_previous,
        salinity_initial,
        salinity_new,
        salinity_previous,
        phi_initial,
        liquid_fraction_new,
        phi_prev,
    )
