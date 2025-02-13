import numpy as np

from src.spyice.parameters.user_input import UserInput

ui = UserInput()
(
    rho_br,
    rho_i,
    c_br,
    c_i,
    k_i,
    k_br,
    porosity_solid,
    L,
    diffusivity_solid,
    m,
    rho_w,
    c_w,
    k_w,
) = (
    ui.constants.rho_br,
    ui.constants.rho_i,
    ui.constants.c_br,
    ui.constants.c_i,
    ui.constants.k_i,
    ui.constants.k_br,
    ui.constants.P_s,
    ui.constants.L,
    ui.constants.D_s,
    ui.constants.m,
    ui.constants.rho_w,
    ui.constants.c_w,
    ui.constants.k_w,
)

def update_coefficients(argument, x_initial, w, phi, nz, salinity_initial):
    """Updates of coefficients required to solve the Advection Reaction Diffusion Equation for each time step for temperature or salinity

    Args:
        argument (str): Either 'temperature' or 'salinity'
        x_initial (float): Initial value for salinity or temperature [ppt] or [K]
        w (float): Brine velocity [ms-1]
        phi (float): Liquid fraction [-]
        nz (int): Number of computational nodes
    Returns:
        tuple: A tuple containing the following coefficients:
            - a (numpy.ndarray): 'temperature': heat capacity | 'salinity': liquid fraction
            - b (numpy.ndarray): Brine velocity [ms-1]
            - c (numpy.ndarray): 'temperature': thermal conductivity | 'salinity': salt diffusivity
            - d (numpy.ndarray): 'temperature': latent heat | 'salinity': factor to determine salinity increase due to liquid fraction decrease


    """
    assert argument in [
        "temperature",
        "salinity",
    ], "invalid input for argument in `update coefficients`"

    if argument == "temperature":
        temperature_initial = x_initial  # noqa: F841
        d = np.ones(nz, dtype=np.float64) * (rho_i * L)
        if salinity_initial in ["S0", "S1"]:  # no brine
            a = np.ones(nz, dtype=np.float64) * _effective_parameter(
                rho_w * c_w, rho_i * c_i, phi, nz
            )
            b = np.ones(nz, dtype=np.float64) * rho_w * c_w * w
            c = np.ones(nz, dtype=np.float64) * _effective_parameter(k_w, k_i, phi, nz)
        else:
            a = np.ones(nz, dtype=np.float64) * _effective_parameter(
                rho_br * c_br, rho_i * c_i, phi, nz
            )
            b = np.ones(nz, dtype=np.float64) * rho_br * c_br * w
            c = np.ones(nz, dtype=np.float64) * _effective_parameter(k_br, k_i, phi, nz)
    elif argument == "salinity":
        salinity_initial = x_initial
        a = np.ones(nz, dtype=np.float64) * phi
        b = np.ones(nz, dtype=np.float64) * w
        c = np.ones(nz, dtype=np.float64) * _salt_diffusivity(phi, nz)
        d = (
            np.ones(nz, dtype=np.float64) * salinity_initial * rho_i / rho_br
        )  # * (1 - P_s)

    return a, b, c, d


def _salt_diffusivity(phi, nz):
    """Computes effective diffusion coefficient of salt based on Archies law

    Args:
        phi (float): Ice volume fraction [-]
        nz (int): Number of computational nodes
    Returns:
        numpy.ndarray: Effective diffusion coefficient [s2m-1]
    """
    diffusion_effective = np.zeros(nz, dtype=np.float64)
    for i in range(nz):
        diffusion_effective[i] = diffusivity_solid * phi[i] ** m  # Archies law
    return diffusion_effective


def _effective_parameter(brine_coefficient, ice_coefficient, phi, nz):
    """Equation for effective parameters based on brine and ice physical properties based on liquid fraction

    Args:
        brine_coefficient (float): Brine property coefficient
        ice_coefficient (float): Ice property coefficient
        phi (float): Liquid fraction
        nz (int): Number of computational nodes
        float: Effective coefficient

    """
    eff = np.zeros(nz, dtype=np.float64)
    # phi = np.where(phi < 0.5, 0, np.where(phi > 0.5, 1, phi))
    eff = phi * brine_coefficient + (1 - phi) * ice_coefficient
    return eff

