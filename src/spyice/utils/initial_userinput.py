
def calculate_initial_melt_temperature_onephase(boundary_salinity):
    return 273.15 - 1.853 * boundary_salinity / 28.0

def calculate_initial_temperature_frezchem(boundary_salinity):
    return (
            -(9.1969758 * (1e-05) * boundary_salinity**2)
            - 0.03942059 * boundary_salinity
            + 272.63617665
        )

def calculate_initial_melt_temperature_twophase(boundary_salinity):
    T_m = 273.15  # melt temperature as solidus temperature
    T_s = 252.05  # eutectic temperature for Sbr = 233ppt
    S_br = 233.0   # brine salinity in ppt
    return 273.15 + (T_s - T_m)*boundary_salinity/S_br

def calculate_initial_melt_temperature(boundary_salinity, method):
    """
    Calculate the initial melt temperature based on the method specified.
    
    Parameters:
    - boundary_salinity: The salinity of the boundary.
    - method: The method to use for calculation. Options are:
        - 'onephase': One phase model
        - 'Frezchem': Frezchem model
        - 'twophase': Two phase model
    
    Returns:
    - Initial melt temperature in Kelvin.
    """
    if method == 'onephase':
        return calculate_initial_melt_temperature_onephase(boundary_salinity)
    elif method == 'Frezchem':
        return calculate_initial_temperature_frezchem(boundary_salinity)
    elif method == 'twophase':
        return calculate_initial_melt_temperature_twophase(boundary_salinity)
    else:
        raise ValueError("Invalid method specified.")