import numpy as np

from src.spyice.parameters.user_input import UserInput

ui = UserInput()

_temperature_melt, _boundary_top_temperature = (
    ui.temperature_melt,
    ui.boundary_top_temperature,
)

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
        if isinstance(_salinity, float):
            _melting_temperature_seawater = _temperature_melt * _salinity + (T_s - T_m)*_salinity/S_br
        else:
            _melting_temperature_seawater = _temperature_melt * np.ones(len(_salinity)) + (T_s - T_m)*_salinity/S_br
    elif _liquid_relation == "Frezchem":
        # _melting_temperature_seawater = _temperature_melt + (
        #     -(9.1969758 * (1e-05) * _salinity**2) - 0.03942059 * _salinity
        # )
        _melting_temperature_seawater = 272.63617665 + (
            -(9.1969758 * (1e-05) * _salinity**2) - 0.03942059 * _salinity
        )
    return _melting_temperature_seawater