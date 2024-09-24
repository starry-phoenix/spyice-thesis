from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

"""
This file lists all constants used for the computations. 
The names are (more or less) equivalent to the ones used in the Code by Buffo_et_al. (2018) 
available on GitHub (https://github.com/jbuffo/SlushFund-2.0-interactive).
They use several formulations and constant provided by Griewank and Notz (2013)
"""


@dataclass(frozen=True, slots=True)
class RealConstants:
    """Class representing real-valued constants used in the model.

    Args:
        phi_ini (int): Initial freezing value.
        phi_ini_Stefan (int): Initial freezing value according to Stefan condition.
        beta (float): Coefficient for density dependence on salinity.
        kappa (float): Thermal diffusivity.
        mu (float): Viscosity.
        Ra_c (float): Critical Rayleigh number.
        alpha (float): Linear coefficient for Rayleigh number driven advection.
        k_i (float): Thermal conductivity of ice.
        k_br (float): Thermal conductivity of brine.
        k_w (float): Thermal conductivity of water.
        D_s (float): Diffusivity for salt.
        c_i (int): Specific heat of ice.
        c_br (int): Specific heat of seawater.
        c_w (int): Specific heat of water.
        L (int): Latent heat of fusion between ice and water.
        rho_br (int): Density of ocean.
        rho_i (int): Density of ice.
        rho_w (int): Density of water.
        m (int): Cementation exponent for Archie's equation.
        g (float): Earth gravity.
        phi_c (float): Critical porosity.
        P_s (float): Partition coefficient.
        a_phi (float): Coefficient a for porosity calculation.
        b_phi (float): Coefficient b for porosity calculation.
        critical_depth (float): Critical depth value.
    """

    ...

    phi_ini: int = 1  # FREEZING
    phi_ini_Stefan: int = 1  # FREEZING
    beta: float = 0.0005836  # Coefficient for density dependence on salinity
    kappa: float = 1.37 * 10 ** (-7)  # Thermal Diffusivity
    mu: float = 1.88 * 10 ** (-3)  # Viscosity
    Ra_c: float = 1.01  # Critical Rayleigh
    alpha: float = 1.56 * 10 ** (
        -3
    )  # Linear coeff for Rayleigh number driven advection
    k_i: float = 2.0  # Thermal conductivity (ice) [W/m/K]
    k_br: float = 0.6  # Thermal conductivity (brine) [W/m/K]
    k_w: float = k_br  # 2.0
    D_s: float = 2 * 10 ** (-9)  # Diffusivity for Salt
    c_i: int = 2000  # Specific heat of ice
    c_br: int = 2000  # 3985  # Specific heat of seawater (J/kg/K)
    c_w: int = 2000  # 4200  # specific heat of water
    L: int = 334774  # Latent heat of fusion ice<->water (J/Kg)
    rho_br: int = 917  # 1028  # Density of Ocean (used in volume averaging - 1D grav. drainage uses delta S) 34ppt NaCl-1027
    rho_i: int = 917  # Density of Ice (Kg/m^3)                                     #  12.3ppt MgSO4-1012, 100pppt-1103, 282ppt-1323
    rho_w: int = 917  # 1000
    m: int = 2  # Cementation exponent for Archies equation
    g: float = 9.8  # Earth Gravity
    phi_c: float = 0.06823  # Critical Porosity
    P_s: float = 1 / 100  #    Partition coefficient
    a_phi: float = 0.0000059
    b_phi: float = 1 / a_phi
    critical_depth = 0.01
