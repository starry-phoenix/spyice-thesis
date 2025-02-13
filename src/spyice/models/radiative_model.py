import numpy as np

def radiative_source_term(radiative_algae=0.0, radiative_ice=0.0, radiative_organicmatter=0.0):
    """Calculate the radiative source term for the model.

    Args:
        radiative_algae (float): The radiative source term for algae.
        radiative_ice (float): The radiative source term for ice.
        radiative_organicmatter (float): The radiative source term for organic matter.

    Returns:
        float: The radiative source term for the model.
    """

    return radiative_algae + radiative_ice + radiative_organicmatter

def radiative_ice(depth):
    albedo_ice = 0.58
    F_sw_solarradiance = 5 # W/m^2
    i0_surface_transmission_parameter = 0.17 
    kappa_attenuation_coefficient = 1.5 # of ice in m-1 + organic matter (TODO: Set actual value of organic matter in m-1)
    
    I_0 = i0_surface_transmission_parameter*(1-albedo_ice)*F_sw_solarradiance
    I_z = I_0*np.exp(-kappa_attenuation_coefficient*depth)
    radiative_ice = I_z*kappa_attenuation_coefficient

    return radiative_ice

def radiative_organicmatter(depth):
    albedo_organicmatter = 0.58
    F_sw_solarradiance = 5 # W/m^2
    i0_surface_transmission_parameter = 0.17 
    a_det_organic_matter = 0.01 # TODO: Set actual value of organic matter in m-1
    kappa_attenuation_coefficient = 1.5 + a_det_organic_matter # of organic matter in m-1 + ice in m-1

    I_0 = i0_surface_transmission_parameter*(1-albedo_organicmatter)*F_sw_solarradiance
    I_z = I_0*np.exp(-kappa_attenuation_coefficient*depth)
    radiative_organicmatter = I_z*a_det_organic_matter

    return radiative_organicmatter

def calculate_radiative_terms(depth, thickness_index, radiative_algae=0.0):
    radiative_ice_array = radiative_ice(depth)
    radiative_organicmatter = 0.0 # radiative_organicmatter(depth)
    # calculate for single layer BAL
    radiative_algae_all_depth = np.zeros(len(depth))
    radiative_algae_all_depth[thickness_index] = radiative_algae

    return radiative_source_term(radiative_algae=radiative_algae_all_depth, radiative_ice=radiative_ice_array, radiative_organicmatter=radiative_organicmatter)

def calculate_local_rayleigh_number(thickness_index, thickness, salinity, phi, grid_size,):
    gravity = 9.8 # m/s^2
    beta = 0.0005836
    kappa = 1.37*1e-7  # thermal diffusivity of ice in m^2/s
    mu = 1.8*1e-3 # dynamic viscosity of ice in Pa*s
    rho_sw = 1028 # kg/m^3
    phi_i = phi[:thickness_index]
    pi_i = calculate_permeability(phi_i)
    S_sw = salinity[thickness_index]
    z = thickness[thickness_index]
    h_i = np.linspace(0, z, thickness_index)  # height at each layer until the interface within the sea ice
    Ra_i = np.zeros(thickness_index)
    
    for i in range(thickness_index):
        if len(h_i) == 0:
            min_pi_i = min(pi_i[i:])
            Ra_i[i] = (gravity*rho_sw*beta*(salinity[i]-S_sw)*min_pi_i*(z-0))/(kappa*mu)
        else:
            min_pi_i = min(pi_i[i:])
            Ra_i[i] = (gravity*rho_sw*beta*(salinity[i]-S_sw)*min_pi_i*(z-h_i[i]))/(kappa*mu)

    return np.absolute(Ra_i)

def calculate_salinity_flux(dz,dt, Ra_c, thickness_index, thickness, salinity, phi):
    
    alpha = 1.56*1e-01  # Linear coeff for Rayleigh number driven advection
    Ra_i = calculate_local_rayleigh_number(thickness_index, thickness, salinity, phi, grid_size=dz)
    
    return alpha*(Ra_i - Ra_c)*(dz**3)*dt

def calculate_permeability(phi_i):
    
    return (1e-17)*((1e3)*phi_i)**3.1

def calculate_brine_flux(thickness_index, thickness, interface_depth, salinity, phi, dt, dz, Ra_c=10):
    brine_flux = np.zeros(thickness_index)
    depth_until_interface = np.arange(0, interface_depth, dz)
    Ra_i = calculate_local_rayleigh_number(thickness_index, thickness, salinity, phi, grid_size=dz)
    brine_flux = calculate_salinity_flux(dz,dt, Ra_c, thickness_index, thickness, salinity, phi)

    for i,z in enumerate(depth_until_interface):
        if z > interface_depth:
            brine_flux[i] = 0.0
        elif Ra_i[i]< Ra_c:
            brine_flux[i] = 0.0
        elif phi[i] < 0.05:
            brine_flux[i] = 0.0
        else:
            pass
    
    return brine_flux

def calculate_in_out_brine_flux(thickness_index, flux, phi):
    flux_in = np.zeros(thickness_index)
    flux_out = np.zeros(thickness_index)

    for i in range(thickness_index):
        if i == 0:
            flux_out[i] = 0.0
            if flux[i] == 0.0:
                flux_in[i] = 0.0
            else:
                flux_in[i] = flux[i]
        elif phi[i] >=0.95:
            flux_in[i] = 0.0
            flux_out[i] = 0.0
        else:
            flux_out[i] = flux_in[i-1]
            flux_in[i] = flux[i] + flux_out[i]

    return flux_in, flux_out

def calculate_salinity_source_term_from_brineflux(salinity, flux_in, flux_out, flux_channel):
    source_term = np.zeros(len(flux_channel))

    for i in range(len(source_term)):
        if i == 0 and len(source_term) > 1:
            source_term[i] = salinity[i+1]*flux_in[i+1] - salinity[i]*flux_channel[i] 
        elif len(source_term) > 0:
            source_term[i] = salinity[i+1]*flux_in[i] - salinity[i]*flux_channel[i] - salinity[i]*flux_out[i] 
        else:
            source_term[i] = 0.0

    return source_term

def get_salinity_source_term(thickness_index, thickness, salinity, phi, dt, dz):
    """
    Calculate the salinity source term from the brine flux.

    Args:
    - thickness_index (int): The index of the thickness.
    - thickness (np.array): The thickness of the sea ice.
    - salinity (np.array): The salinity of the sea ice.
    - phi (np.array): The porosity of the sea ice.
    - dt (float): The time step.
    - dz (float): The spatial step size.
    - alpha (float): The linear coefficient for Rayleigh number driven advection.
    """
    interface_depth = thickness[thickness_index]
    get_brine_flux = calculate_brine_flux(thickness_index, thickness, interface_depth, salinity, phi, dt, dz, Ra_c=10)
    get_influx, get_outflux = calculate_in_out_brine_flux(thickness_index, get_brine_flux, phi)
    get_source_term = calculate_salinity_source_term_from_brineflux(salinity, get_influx, get_outflux, get_brine_flux)
    
    if len(get_source_term) < len(salinity):
        get_source_term = np.append(get_source_term, np.zeros(len(salinity)-len(get_source_term)))
    
    return get_source_term

# TODO: calculate temperature source term from brine flux