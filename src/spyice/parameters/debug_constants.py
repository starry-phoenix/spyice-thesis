from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


@dataclass(frozen=True, slots=True)
class DebugConstants:
    """Class representing debug constants used in the model.

    Args:
        phi_ini (int): Freezing constant.
        phi_ini_Stefan (int): Freezing constant for Stefan condition.
        beta (int): Coefficient for density dependence on salinity.
        kappa (int): Thermal diffusivity.
        mu (int): Viscosity.
        Ra_c (int): Critical Rayleigh number.
        alpha (int): Linear coefficient for Rayleigh number driven advection.
        k_i (int): Thermal conductivity of ice [W/m/K].
        k_br (int): Thermal conductivity of brine [W/m/K].
        k_w (int): Thermal conductivity of water.
        D_s (int): Diffusivity for salt.
        c_br (int): Specific heat of seawater (J/kg/K).
        c_i (int): Specific heat of ice.
        c_w (int): Specific heat of water.
        L (int): Latent heat of fusion ice<->water (J/Kg).
        rho_i (int): Density of ice (Kg/m^3).
        rho_br (int): Density of ocean used in volume averaging.
        rho_w (int): Density of water.
        m (int): Cementation exponent for Archie's equation.
        g (int): Gravity constant.
        phi_c (int): Constant for phi.
        P_s (int): Constant for P_s.
        a_phi (int): Constant for a_phi.
        b_phi (int): Constant for b_phi.
    """

    # Class implementation goes here

    phi_ini: int = 1  # FREEZING
    phi_ini_Stefan: int = 1  # FREEZING
    beta: int = 1  # Coefficient for density dependence on salinity
    kappa: int = 1  # Thermal Diffusivity
    mu: int = 1  # Viscosity
    Ra_c: int = 0  # Critical Rayleigh
    alpha: int = 1  # Linear coeff for Rayleigh number driven advection
    k_i: int = 1  # Thermal conductivity (ice) [W/m/K]
    k_br: int = 1  # Thermal conductivity (brine) [W/m/K]
    k_w: int = 1
    D_s: int = 1e-4  # Diffusivity for Salt
    c_br: int = 1  # Specific heat of seawater (J/kg/K)
    c_i: int = 1  # Specific heat of ice
    c_w: int = 1  # specific heat of water
    L: int = 1  # Latent heat of fusion ice<->water (J/Kg)
    rho_i: int = 1  # Density of Ice (Kg/m^3)
    rho_br: int = 1  # Density of Ocean (used in volume averaging - 1D grav. drainage uses delta S) 34ppt NaCl-1027                                     #  12.3ppt MgSO4-1012, 100pppt-1103, 282ppt-1323
    rho_w: int = 1
    m: int = 1  # Cementation exponent for Archies equation
    g: int = 0
    phi_c: int = 1
    P_s: int = 1
    a_phi: int = 1
    b_phi: int = 1
