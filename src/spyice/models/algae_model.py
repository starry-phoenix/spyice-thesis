# import packages

import warnings

import numpy as np

# Suppress runtime warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# salinity function
def fs_salinity(s_br_array):
    """
    Calculate the salinity function f_s as a function of salinity.
    The equation used is f_s = exp(-((2.16 - a - b) ** 2)), where:
    - a = 8.3 * 10 ** (-5) * s ** 2.11
    - b = 0.55 * log(s)
    - s is the salinity in g kg^-1

    Parameters:
    s_min (float): Minimum salinity in g kg^-1.
    s_max (float): Maximum salinity in g kg^-1.

    Returns:
    numpy.ndarray: Salinity function f
    """
    a_br = 8.3 * 10 ** (-5) * np.power(s_br_array, 2.11)
    b_br = 0.55 * np.log(s_br_array)
    f_s = np.exp(-((2.16 - a_br - b_br) ** 2))

    return f_s


# temperature function
def ft_temperature(t_c_array):
    """
    Model temperature function f_t as an exponential function of temperature.

    Parameters:
    t_min (float): Minimum temperature in degrees Celsius.
    t_max (float): Maximum temperature in degrees Celsius.

    Returns:
    numpy.ndarray: Temperature function f_t

    """
    r_g = 0.0633
    f_t = np.exp(r_g * (t_c_array - 273.15))

    return f_t


# nutrient function
def ln_nutrient(c_n_array, k: float):
    """A function to calculate the nutrient function as a function of nutrient concentration.
    The equation used is l_n = c_n / (k_n + c_n), where:
    - c_n is the nutrient concentration in mmol m^-3
    - k_n = 1.6 micro M
    - k_si = 3.9 micro M
    - k_p = 0.24 micro M
    - l_n = c_n / (k_n + c_n)

    Parameters:
    c_min (float): Minimum nutrient concentration in mmol m^-3.
    c_max (float): Maximum nutrient concentration in mmol m^-3.
    k (float): Nutrient concentration in mmol m^-3.

    Returns:
    numpy.ndarray: Nutrient function l_n
    """
    # c in mmol m-3
    k_n = 1.6  # in micrometers
    k_si = 3.9  # in micrometers
    k_p = 0.24  # in micrometers
    l_n = c_n_array / (k_n + c_n_array)
    l_si = c_n_array / (k_si + c_n_array)
    l_p = c_n_array / (k_p + c_n_array)
    # nutrient for k
    l_random = c_n_array / (k + c_n_array)

    return l_random


# light function


# PAR function
def photosynthetic_active_radiation(z, kappa=1.5, i_0=0.17, albedo=0.58, F_s_w=5):
    """
    Calculate the Photosynthetically Active Radiation (PAR) at a given depth in ice.
    The equation used is PAR = 4.91 * I(z), where I(z) = I_0 * exp(-kappa * z).
    The constants used are:
    - kappa = 1.5 m^-1
    - i_0 = 0.17 for snow free ice and =0 for snow covered ice
    - alpha = 1e-4 g C (g Chla h micro E m-2 s-1)^-1
    - F_s_w = 5 W m^-2 incoming solar irradiance 1361 https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://en.wikipedia.org/wiki/Solar_irradiance&ved=2ahUKEwiL5Z2Qw_6JAxV30wIHHWQXD4AQFnoECCEQAw&usg=AOvVaw3cpr_oSPFpcgj7Ny50mW76
    - I_0 = i_0 * (1 - alpha) * F_s_w
    - PAR = 4.91 * I(z)

    Parameters:
    - z (float): Depth in ice (m).
    - kappa (float): Attenuation coefficient (m^-1).
    - i_0 (float): Incident solar radiation (W m^-2).
    - alpha (float): Absorption coefficient (g C (g Chla h micro E m-2 s-1)^-1).
    - F_s_w (float): Incoming solar irradiance (W m^-2).

    Returns:
    - float: PAR at the given depth in ice (W m^-2).
    """
    I_array = irradiance(z, kappa, i_0, albedo, F_s_w)
    PAR = 4.91 * I_array  # in W m-2

    return PAR


def irradiance(z, kappa=1.5, i_0=0.17, albedo=0.58, F_s_w=5):
    I_0 = i_0 * (1 - albedo) * F_s_w
    I_array = I_0 * np.exp(-kappa * z)  # z is the depth in ice

    return I_array


# chlorophyll to carbon ratio function
def cholorophyl_to_C_ratio_par(PAR, ln, E=0.5, r_chl_c_max=0.05, r_chl_c_min=0.01):
    """
    Calculate the chlorophyll to carbon ratio r_chl_c_par as a function of depth in ice.
    The equation used is r_chl_c_par = r_chl_c_max - (r_chl_c_max - r_chl_c_min) * min(PAR / E, 1) * ln, where:
    - r_chl_c_max = 0.05
    - r_chl_c_min = 0.01
    - PAR = 4.91 * I(z)
    - E = 0.5 in micro E m-2 s-1
    - ln = c_n / (k_n + c_n)

    Parameters:
    r_chl_c_max (float): Maximum chlorophyll to carbon ratio.
    r_chl_c_min (float): Minimum chlorophyll to carbon ratio.
    PAR (float): Photosynthetically Active Radiation at the given depth in ice (W m^-2).
    E (float): Energy (W m^-2).
    ln (float): Nutrient function l_n.

    Returns:
    float: Chlorophyll to carbon ratio r_chl_c_par at the given depth in ice.
    """

    return r_chl_c_max - (r_chl_c_max - r_chl_c_min) * np.min([np.min(PAR / E), 1]) * ln


# maximum photosynthesis rate function
def Photosynthesis_rate_maximum_light(mu_m, fs, ft, ln, r_chl_c_par):
    """
    Calculate the maximum photosynthesis rate P_m as a function of depth in ice.
    The equation used is P_m = mu_m * fs * ft * ln / r_chl_c_par, where:
    - mu_m = 2e-5 in s^-1
    - fs = exp(-((2.16 - a - b) ** 2))
    - ft = exp(r_g * T_c)
    - ln = c_n / (k_n + c_n)
    - r_chl_c_par = r_chl_c_max - (r_chl_c_max - r_chl_c_min) * min(PAR / E, 1) * ln

    Parameters:
    mu_m (float): Maximum photosynthesis rate in s^-1.
    fs (numpy.ndarray): Salinity function f_s.
    ft (numpy.ndarray): Temperature function f_t.
    ln (numpy.ndarray): Nutrient function l_n.
    r_chl_c_par (numpy.ndarray): Chlorophyll to carbon ratio.

    Returns:
    numpy.ndarray: Maximum photosynthesis rate P_m at the given depth in ice.
    """
    return mu_m * fs * ft * ln / r_chl_c_par


def Ek_light(Photosynthesis_rate, alpha):
    """
    Calculate the light limitation Ek as a function of depth in ice.
    The equation used is Ek = P_m / alpha, where:
    - P_m = mu_m * fs * ft * ln / r_chl_c_par
    - mu_m = 2e-5 in s^-1
    - E = 0.5 in micro E m-2 s-1
    - alpha = 1e-4 in g C (g Chla h micro E m-2 s-1)-1
    - r_chl_c_max = 0.05
    - r_chl_c_min = 0.01
    - ln = c_n / (k_n + c_n)
    - fs = exp(-((2.16 - a - b) ** 2))
    - ft = exp(r_g * T_c)
    - r_chl_c_par = r_chl_c_max - (r_chl_c_max - r_chl_c_min) * min(PAR / E, 1) * ln
    - P_m = mu_m * fs * ft * ln / r_chl_c_par
    - E_k = P_m / alpha

    Parameters:
    - Photosynthesis_rate (float): Maximum photosynthesis rate in s^-1.
    - alpha (float): in g C (g Chla h micro E m-2 s-1)-1

    Returns:
    - float: Light limitation Ek at the given depth in ice.
    """

    return Photosynthesis_rate / alpha


def Ek_light_tanh(PAR_arr: np.array, Ek_arr: np.array):
    """

    Calculate the light limitation Ek as a function of depth in ice.
    The equation used is Ek = P_m / alpha, where:
    - P_m = mu_m * fs * ft * ln / r_chl_c_par
    - mu_m = 2e-5 in s^-1
    - E = 0.5 in micro E m-2 s-1
    - alpha = 1e-4 in g C (g Chla h micro E m-2 s-1)-1
    - r_chl_c_max = 0.05
    - r_chl_c_min = 0.01
    - ln = c_n / (k_n + c_n)
    - fs = exp(-((2.16 - a - b) ** 2))
    - ft = exp(r_g * T_c)
    - r_chl_c_par = r_chl_c_max - (r_chl_c_max - r_chl_c_min) * min(PAR / E, 1) * ln
    - P_m = mu_m * fs * ft * ln / r_chl_c_par
    - E_k = P_m / alpha

    Parameters:
    - PAR_arr (np.array): PAR at the given depth in ice.
    - Ek_arr (np.array): Light limitation Ek at the given depth in ice.

    Returns:
    - np.array: Light limitation Ek at the given depth in ice.

    """

    return np.tanh(PAR_arr / Ek_arr)


# photosynthetic rate function
def photosynthetic_rate(max_mu, fs, ft, ln, lpar):
    """
    Calculate the photosynthetic rate as a function of depth in ice.
    The equation used is mu = max_mu * fs * ft * ln * lpar, where:
    - max_mu = 2e-5 in s^-1
    - fs = exp(-((2.16 - a - b) ** 2)), where:
        - a = 8.3 * 10 ** (-5) * s ** 2.11
        - b = 0.55 * log(s)
        - s is the salinity in g kg^-1
    - ft = exp(r_g * T_c)
    - r_g = 0.0633
    - T_c is the temperature in degrees Celsius
    - ln = c_n / (k_n + c_n)
    - lpar = tanh(PAR / Ek)
    - Ek = P_m / alpha
    - P_m = mu_m * fs * ft * ln / r_chl_c_par
    - r_chl_c_par = r_chl_c_max - (r_chl_c_max - r_chl_c_min) * min(PAR / E, 1) * ln
    - E = 0.5 in micro E m-2 s^-1
    - alpha = 1e-4 in g C (g Chla h micro E m-2 s^-1)-1
    - mu = max_mu * fs * ft * ln * lpar

    Parameters:
    - s_min (float): Minimum salinity in g kg^-1.
    - s_max (float): Maximum salinity in g kg^-1.
    - t_min (float): Minimum temperature in degrees Celsius.
    - t_max (float): Maximum temperature in degrees Celsius.
    - c_min (float): Minimum nutrient concentration in mmol m^-3.
    - c_max (float): Maximum nutrient concentration in mmol m^-3.
    - k (float): Nutrient concentration in mmol m^-3.
    - z_min (float): Minimum depth in ice in m.
    - z_max (float): Maximum depth in ice in m.
    - alpha (float): in g C (g Chla h micro E m^-2 s^-1)-1
    - index (int): index of the depth in ice.

    Returns:
    - float: Photosynthetic rate at the given depth in ice.
    """
    return max_mu * fs * ft * ln * lpar


# f_s * f_t function
def fs_ft(f_s, f_t):
    """
    Calculate the temperature x salinity function f_s * f_t.
    The equation used is f_s * f_t = f_t * f_s, where:
    - f_t = exp(r_g * T_c)
    - f_s = exp(-((2.16 - a - b) ** 2)), where:
        - a = 8.3 * 10 ** (-5) * s ** 2.11
        - b = 0.55 * log(s)
        - s is the salinity in g kg^-1
    - r_g = 0.0633

    Parameters:
    t_min (float): Minimum temperature in degrees Celsius.

    Returns:
    numpy.ndarray: Temperature x salinity function f_s * f_t
    """
    return f_t * f_s


def fs_ft_ln(f_s, f_t, l_n):
    """
    Calculate the salinity x temperature x nutrient function f_s * f_t * l_n.
    The equation used is f_s * f_t * l_n = f_s * f_t * l_n, where:
    - f_s = exp(-((2.16 - a - b) ** 2)), where:
        - a = 8.3 * 10 ** (-5) * s ** 2.11
        - b = 0.55 * log(s)
        - s is the salinity in g kg^-1

    - f_t = exp(r_g * T_c)
    - r_g = 0.0633
    - T_c is the temperature in degrees Celsius

    - l_n = c_n / (k_n + c_n), where:
        - c_n is the nutrient concentration in mmol m^-3
        - k_n = 1.6 micrometers

    Parameters:
    s_min (float): Minimum salinity in g kg^-1.
    s_max (float): Maximum salinity in g kg^-1.
    t_min (float): Minimum temperature in degrees Celsius.
    t_max (float): Maximum temperature in degrees Celsius.
    c_min (float): Minimum nutrient concentration in mmol m^-3.
    c_max (float): Maximum nutrient concentration in mmol m^-3.
    k (float): Nutrient concentration in mmol m^-3.

    Returns:
    numpy.ndarray: Salinity x temperature x nutrient function f_s * f_t * l_n
    """

    return f_s * f_t * l_n


def model_algae_processes(
    s_br_array,
    t_c_array,
    c_n_array,
    z_array,
    k=1.6,
    mu_m=0.86 / (3600 * 24),
    alpha=1e-4,
    r_chl_c_max=0.05,
    r_chl_c_min=0.01,
    r_const=True,
    r_n_c=0.12,
    i_0=0.17,
    E=0.5,
    kappa=1.5,
    F_s_w=5,
    albedo=0.58,
):
    """
    Model the algae processes in sea ice.

    Parameters:
    - s_br_array (np.array): Salinity in g kg^-1.
    - t_c_array (np.array): Temperature in degrees Celsius.
    - c_n_array (np.array): Nutrient concentration in mmol m^-3.
    - z_array (np.array): Depth in ice in m.
    - k (float): half saturation constant in micro M = mmol m^-3.
    - mu_m (float): Maximum photosynthesis rate in s^-1.
    - alpha (float): Absorption coefficient in g C (g Chla h micro E m-2 s-1)^-1.
    - r_chl_c_max (float): Maximum chlorophyll to carbon ratio.
    - r_chl_c_min (float): Minimum chlorophyll to carbon ratio.
    - r_const (bool): Chlorophyll to carbon ratio constant.
    - r_n_c (float): Nutrient to carbon ratio in mmol m^-3. Here nutrient is dissolved silica
    - i_0 (float): Incident solar radiation in W m^-2.
    - E (float): Energy in W m^-2.
    - kappa (float): Attenuation coefficient in m^-1.
    - F_s_w (float): Incoming solar irradiance in W m^-2.

    Returns:
    - np.array: Photosynthetic rate at the given depth in ice in s^-1.
    """
    # field functions
    # TODO: change alpha value
    salinity_function = fs_salinity(s_br_array)
    temperature_function = ft_temperature(t_c_array)
    nutrient_function = ln_nutrient(c_n_array, k)

    # light function dependent on z_array
    PAR = photosynthetic_active_radiation(z_array, kappa, i_0, albedo, F_s_w)

    if r_const:
        r_chl_c_par = r_chl_c_max
    else:
        r_chl_c_par = cholorophyl_to_C_ratio_par(PAR, nutrient_function)
    P_m = Photosynthesis_rate_maximum_light(
        mu_m, salinity_function, temperature_function, nutrient_function, r_chl_c_par
    )
    Ek = Ek_light(P_m, alpha)
    l_PAR = Ek_light_tanh(PAR, Ek)

    # photosynthetic rate
    mu = photosynthetic_rate(
        mu_m, salinity_function, temperature_function, nutrient_function, l_PAR
    )

    return mu, PAR, nutrient_function


def ode_update_carbon_nutrient_uptake(
    dt,
    salinity_list,
    temperature_list,
    nutrient_list,
    depth_list,
    cc_old,
    cn_old,
    r_n_c=0.12,
):
    """
    Perform an ODE update for carbon and nutrient uptake in algae.

    Parameters:
    - dt (float): Time step in days.
    - salinity_list (np.array): Salinity in g kg^-1.
    - temperature_list (np.array): Temperature in degrees Celsius.
    - nutrient_list (np.array): Nutrient concentration in mmol m^-3.
    - depth_list (np.array): Depth in ice in node-wise and needs to be greater than 3.
    - cc_old (float): Old carbon concentration in mmol m^-3.
    - cn_old (float): Old nutrient concentration in mmol m^-3.
    - r_n_c (float): Nutrient to carbon ratio in mmol m^-3. Here nutrient is dissolved silica. [Vancoppennolle et al. 2010 case of dissolved silica]

    """

    # cc = 1  # TODO: find value - one algal group concentration in mmol m-3
    # cn = 15  # nutrient concentration in mmol m-3 as per Vancoppenoelle et al. 2010
    f = 1  # nutrient conservation
    lambda_ = 0.15 / (3600 * 24)  # carbon loss rate in s-1

    mu, PAR, _ = model_algae_processes(
        salinity_list, temperature_list, nutrient_list, depth_list
    )

    if dt > 0:
        cc_tp1 = cc_old + dt * (mu - lambda_) * cc_old
        cn_tp1 = cn_old + dt * (-mu + f * lambda_) * cc_old * r_n_c
        mu, PAR, _ = model_algae_processes(
            salinity_list, temperature_list, cn_tp1, depth_list
        )

    return cc_tp1, cn_tp1, mu, PAR


def chla_algae(PAR, nutrient_function, c_bulk_tracer):
    carbon_molar_mass = 12  # g mol-1
    r_chl_c_par = cholorophyl_to_C_ratio_par(PAR, nutrient_function)
    chla_bulk = c_bulk_tracer * r_chl_c_par * carbon_molar_mass

    return chla_bulk


def radiation_algae(chla_bulk_z, I_array):
    absorption_coefficient = (
        0.008  # m-1 (mg Chla m-3)-1 (Arrigo et. al 1991) of ice algae
    )
    # absorption_coefficient = 0.02   # m-1 (mg Chla m-3)-1 (Lavoie et. al 2005) of ice algae
    return chla_bulk_z * I_array * absorption_coefficient


def get_bulk_tracer_concentration(liquid_fraction, brine_concentration):
    return liquid_fraction * brine_concentration


def biogeochemical_model(
    temperature,
    salinity,
    liquid_fraction,
    nutrient_concentration,
    carbon_concentration,
    dt,
    thickness_index,
    thickness,
):
    # biogeochemical model for single BAL at interface
    # calculate photosynthetic rate
    # biologically active layer at a fixed region in ice which is at the interface of ice and ocean
    biologically_active_layer = int(thickness_index)
    # TODO: model for a dynamic BAL
    algae_salinity, algae_temperature, algae_liquid_fraction = (
        salinity[biologically_active_layer],
        temperature[biologically_active_layer],
        liquid_fraction[biologically_active_layer],
    )
    nutrient_concentration_at_interface = nutrient_concentration[biologically_active_layer]
    carbon_concentration_at_interface = carbon_concentration[biologically_active_layer]

    photosynthetic_rate_mu, PAR, nutrient_function = model_algae_processes(
        algae_salinity,
        algae_temperature,
        nutrient_concentration_at_interface,
        thickness,
    )
    # update carbon and nutrient uptake
    (
        carbon_concentration_at_interface,
        nutrient_concentration_at_interface,
        photosynthetic_rate_mu,
        PAR,
    ) = ode_update_carbon_nutrient_uptake(
        dt,
        algae_salinity,
        algae_temperature,
        nutrient_concentration_at_interface,
        thickness,
        carbon_concentration_at_interface,
        nutrient_concentration_at_interface,
    )
    # calculate bulk concentration in ice
    bulk_tracer_concentration = get_bulk_tracer_concentration(
        algae_liquid_fraction, nutrient_concentration_at_interface
    )
    # calculate chla concentration
    # TODO: is chla bulk dependent on bulk tracer or carbon??
    chla_bulk = chla_algae(PAR, nutrient_function, bulk_tracer_concentration)
    I_array = irradiance(thickness)
    radiation = radiation_algae(chla_bulk, I_array)
    # TODO: radiation/(rho_i*c_i)

    # reset the values and map to interface
    photosynthetic_rate_mu_all, radiation_all, chla_bulk_all = np.zeros(len(nutrient_concentration)),np.zeros(len(nutrient_concentration)),np.zeros(len(nutrient_concentration))
    nutrient_concentration[thickness_index] = nutrient_concentration_at_interface
    carbon_concentration[thickness_index] = carbon_concentration_at_interface
    photosynthetic_rate_mu_all[thickness_index] = photosynthetic_rate_mu
    radiation_all[thickness_index] = radiation
    chla_bulk_all[thickness_index] = chla_bulk

    return (
        carbon_concentration,
        nutrient_concentration,
        photosynthetic_rate_mu_all,
        radiation_all,
        chla_bulk_all,
    )


def biogeochemical_model_at_alldepths(
    temperature,
    salinity,
    liquid_fraction,
    nutrient_concentration,
    carbon_concentration,
    dt,
    thickness,
    thickness_index,
):
    # biogeochemical model for single BAL at interface
    # calculate photosynthetic rate
    # TODO: model for a dynamic BAL
    (
        algae_salinity,
        algae_temperature,
        algae_liquid_fraction,
        photosynthetic_rate_mu,
        chla_bulk,
        radiation,
    ) = (
        salinity[:thickness_index],
        temperature[:thickness_index],
        liquid_fraction[:thickness_index],
        np.zeros(len(thickness)),
        np.zeros(len(thickness)),
        np.zeros(len(thickness)),
    )
    nutrient_concentration_at_interface = nutrient_concentration[:thickness_index]
    carbon_concentration_at_interface = carbon_concentration[:thickness_index]

    photosynthetic_rate_mu_at_interface, PAR, nutrient_function = model_algae_processes(
        algae_salinity,
        algae_temperature,
        nutrient_concentration_at_interface,
        thickness[:thickness_index],
    )
    # update carbon and nutrient uptake
    (
        carbon_concentration_at_interface,
        nutrient_concentration_at_interface,
        photosynthetic_rate_mu_at_interface,
        PAR,
    ) = ode_update_carbon_nutrient_uptake(
        dt,
        algae_salinity,
        algae_temperature,
        nutrient_concentration_at_interface,
        thickness[:thickness_index],
        carbon_concentration_at_interface,
        nutrient_concentration_at_interface,
    )
    # calculate bulk concentration in ice
    bulk_tracer_concentration = get_bulk_tracer_concentration(
        algae_liquid_fraction, nutrient_concentration_at_interface
    )
    # calculate chla concentration
    # TODO: is chla bulk dependent on bulk tracer or carbon??
    chla_bulk_at_interface = chla_algae(
        PAR, nutrient_function, bulk_tracer_concentration
    )
    I_array = irradiance(thickness[:thickness_index])
    radiation_at_interface = radiation_algae(chla_bulk_at_interface, I_array)
    # TODO: radiation/(rho_i*c_i)

    # reset the values of the field arrays
    nutrient_concentration[:thickness_index] = nutrient_concentration_at_interface
    carbon_concentration[:thickness_index] = carbon_concentration_at_interface
    photosynthetic_rate_mu[:thickness_index] = photosynthetic_rate_mu_at_interface
    chla_bulk[:thickness_index] = chla_bulk_at_interface
    radiation[:thickness_index] = radiation_at_interface

    return (
        carbon_concentration,
        nutrient_concentration,
        photosynthetic_rate_mu,
        radiation,
        chla_bulk,
    )
