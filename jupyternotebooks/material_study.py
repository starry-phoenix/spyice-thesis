import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# .style.use("spyice.utils.custom")
plt.style.use("jupyternotebooks/custom.mplstyle")
# plt.rcParams.update(
#     {
#         "text.usetex": True,
#     }
# )
#plt.rcParams["pgf.texsystem"] = "pdflatex"
plt.rcParams["text.latex.preamble"].join(
    [
        r"\usepackage{dashbox}",
        r"\setmainfont{xcolor}",
    ]
)
plt.rcParams["animation.convert_path"] = Path(
    "C:/Program Files/ImageMagick-7.1.1-Q16-HDRI/magick.exe"
)
plt.rcParams.update({'font.size': 18})
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18


# # -----------------------------------------------------------------------------
# # Material Properties

S_boundary = 34.0

c_s = lambda T: 2112.2 + 7.6973*T
c_l = lambda T: 4208.8 + 111.71*T + 3.5611*T**2 + 0.052168*T**3
k_s = lambda T: 2.21 - 1e-02*T + 3.44*1e-05*T**2
k_l = lambda T: 0.52325*(1- S_boundary/1000) + 0.01256*T + 5.8604*1e-05*T**2 
rho = lambda S: 1000.3 + 0.7882327*S + 2.8008*1e-04*S**2

temperature = np.arange(-20, 0, 1)
salinity = np.arange(1, 230, 1)

plt.figure(figsize=(6.4, 5))
plt.plot((temperature+ 273.15), c_s(temperature), label=r'$c_{solid}$', color='red', linewidth=2)
plt.plot((temperature+273.15), c_l(temperature), label=r'$c_{liquid}$', color='blue', linewidth=2)
plt.legend()
plt.ylabel(r'Specific heat capacity [J/kg/K]')
plt.xlabel(r'Temperature [K]')
plt.grid(True)
plt.savefig('jupyternotebooks/specific_heat_capacity.pdf', backend='pgf')
plt.close()

plt.figure(figsize=(6.4, 5))
plt.plot((temperature+ 273.15), k_s(temperature), label=r'$k_{solid}$', color='red', linewidth=2)
plt.plot((temperature+ 273.15), k_l(temperature), label=r'$k_{liquid}$', color='blue', linewidth=2)
plt.legend()
plt.ylabel(r'Heat conductivity [W/m/K]')
plt.xlabel(r'Temperature [K]')
plt.grid(True)
plt.savefig('jupyternotebooks/heat_conductivity.pdf', backend='pgf')
plt.close()

plt.figure(figsize=(6.4, 5))
plt.plot(salinity, rho(salinity), label=r'$\rho_{liquid}$', color='blue', linewidth=2)
plt.legend()
plt.ylabel(r'Density [kg/m$^3$]')
plt.xlabel(r'Salinity [ppt]')
plt.grid(True)
plt.savefig('jupyternotebooks/density.pdf', backend='pgf')
plt.close()

# -----------------------------------------------------------------------------
# Enthalpy Equations 

# constants 
rho_l = 1028.0
rho_s =  917.0
c_l = 3985.0
c_s = 2000.0
S_br = 233.0
T_m = 273.15 # melt temperature 
T_s = 252.05 # eutectic temperature for Sbr = 233ppt 
T_l = lambda S: (T_m + (T_s - T_m)/S_br*S)
T_l_34ppt = T_l(34.0)
L = 334774

T_range = np.arange(240,300, 1)

def H_phi(phi, T):
    H_s = rho_s*c_s*(T-T_s)    
    H_l = rho_l*L + rho_l*c_l*(T-T_l_34ppt)
    H_phi = (1-phi)*H_s + phi*H_l 
    return H_phi 

phi_k = np.zeros(len(T_range))
nz = abs(T_s - T_l_34ppt)
for i, t in enumerate(T_range):
    if t <= T_s:
        phi_k[i] = 0.0
    elif t> T_s and t<T_l_34ppt:
        phi_k[i] = phi_k[i-1] + 1/nz
    elif t>= T_l_34ppt:
        phi_k[i] = 1.0
    else:
        print('t out of range')


solid_index = np.where((phi_k > 0) & (phi_k < 1.0))[0][0] - 1
liquid_index = np.where(phi_k == 1.0)[0][0] - 1

H_k = H_phi(phi_k, T_range)

H_solid = H_k[solid_index]
H_liquid = H_k[liquid_index]

fig1,(ax1) = plt.subplots(figsize=(6.4, 5))
ax1.plot(T_range,H_k/1e06, color='blue', linewidth=2)
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Enthalpy [MJ/kg]')
# ax1.axhline(
#     y=H_solid/1e06,
#     color='r',
#     linestyle='dashed',
#     label=r'$H_{solid}$',
# )
# ax1.axhline(
#     y=H_liquid/1e06,
#     color='black',
#     linestyle='dashed',
#     label=r'$H_{liquid}$',
# )

ax1.axvline(
    x=T_range[solid_index],
    color='r',
    linestyle='dashed',
    label=r'$T_s$',
)
# ax1.annotate(r'$T_s$', xy=(T_range[solid_index], 0), xytext=(T_range[solid_index], -0.5),
#              textcoords='offset points', ha='center', color='r')

ax1.axvline(
    x=T_range[liquid_index],
    color='black',
    linestyle='dashed',
    label=r'$T_l$',
)
# ax1.annotate(r'$T_l$', xy=(T_range[liquid_index], 0), xytext=(T_range[liquid_index], -0.5),
#              textcoords='offset points', ha='center', color='black')

plt.grid(True)
plt.legend()
fig1.savefig('jupyternotebooks/enthalpy.pdf', backend='pgf')
plt.close()

#     ax3 = ax1.twinx()
#     ax3.axvline(
#         T_melt,
#         color='b',
#         linestyle='dashed',
#         label=f'T_melt:{str(round(T_melt, 2))}',
#     )

# def f_s_salinity(s_min, s_max):
#     """
#     Calculate the salinity function f_s as a function of salinity.
#     The equation used is f_s = exp(-((2.16 - a - b) ** 2)), where:
#     - a = 8.3 * 10 ** (-5) * s ** 2.11
#     - b = 0.55 * log(s)
#     - s is the salinity in g kg^-1
    
#     Parameters:
#     s_min (float): Minimum salinity in g kg^-1.
#     s_max (float): Maximum salinity in g kg^-1.

#     Returns:
#     numpy.ndarray: Salinity function f
#     """
#     s_br = np.linspace(s_min, s_max, 150)
#     a_br = 8.3 * 10 ** (-5) * np.power(s_br, 2.11)
#     b_br = 0.55 * np.log(s_br)
#     f_s = np.exp(-((2.16 - a_br - b_br) ** 2))

#     return f_s

# f_s = f_s_salinity(20, 250)

# # plt.plot(np.linspace(1, 230, 150),f_s, color='black', linewidth=2)
# # plt.ylabel(r'$f_S$')
# # plt.xlabel(r'Salinity [ppt]')
# # plt.grid(True)
# # plt.savefig('jupyternotebooks/salinity_fs.pdf', backend='pgf')
# # plt.close()

# def ft_temperature(t_min, t_max):
#     """
#     Model temperature function f_t as an exponential function of temperature.

#     Parameters:
#     t_min (float): Minimum temperature in degrees Celsius.
#     t_max (float): Maximum temperature in degrees Celsius.

#     Returns:
#     numpy.ndarray: Temperature function f_t
    
#     """
#     T_c = np.linspace(t_min, t_max, 150)
#     r_g = 0.0633
#     f_t = np.exp(r_g * T_c)

#     return f_t

# f_T = ft_temperature(-20, 0)

# # plt.plot((np.linspace(-20, 0, 150)+273.15), f_T, color='black', linewidth=2)
# # plt.ylabel(r'$f_T$')
# # plt.xlabel(r'Temperature [K]')
# # plt.grid(True)
# # plt.savefig('jupyternotebooks/temperature_fT.pdf', backend='pgf')
# # plt.close()

# def nutrient_function(c_min, c_max, k):
#     """A function to calculate the nutrient function as a function of nutrient concentration.
#     The equation used is l_n = c_n / (k_n + c_n), where:
#     - c_n is the nutrient concentration in mmol m^-3
#     - k_n = 1.6 micrometers
#     - k_si = 3.9 micrometers
#     - k_p = 0.24 micrometers
#     - l_n = c_n / (k_n + c_n)

#     Parameters:
#     c_min (float): Minimum nutrient concentration in mmol m^-3.
#     c_max (float): Maximum nutrient concentration in mmol m^-3.
#     k (float): Nutrient concentration in mmol m^-3.

#     Returns:
#     numpy.ndarray: Nutrient function l_n
#     """
#     # c in mmol m-3
#     c_n = np.linspace(c_min, c_max, 150)
#     k_n = 1.6  # in micro M = mmol m-3
#     k_si = 3.9  # in micro M = mmol m-3
#     k_p = 0.24  # in micro M = mmol m-3
#     l_n = c_n / (k_n + c_n)
#     l_si = c_n / (k_si + c_n)
#     l_p = c_n / (k_p + c_n)
#     # nutrient for k
#     # l_random = c_n / (k + c_n)

#     return l_n, l_si, l_p

# l_n, l_si, l_p = nutrient_function(0, 100, 0)

# plt.plot(np.linspace(0, 100, 150),l_n, color='black', linewidth=2, label=r'$N$')
# plt.plot(np.linspace(0, 100, 150),l_si, color='black', linestyle='--', linewidth=2, label=r'$Si$')
# plt.plot(np.linspace(0, 100, 150),l_p, color='black', linestyle=':', linewidth=2, label=r'$P$')
# plt.legend()
# plt.ylabel(r'$L_N$')
# plt.xlabel(r'Nutrient concentration [mmol m$^{-3}$]')
# plt.xscale('log')
# plt.grid(True)
# plt.savefig('jupyternotebooks/nutrient_limitation_ln.pdf', backend='pgf')
# plt.close()

# def PAR(z):
#     """
#     Calculate the Photosynthetically Active Radiation (PAR) at a given depth in ice.
#     The equation used is PAR = 4.91 * I(z), where I(z) = I_0 * exp(-kappa * z).
#     The constants used are:
#     - kappa = 1.5 m^-1
#     - i_0 = 0.17 for snow free ice and =0 for snow covered ice
#     - alpha = 1e-4 g C (g Chla h micro E m-2 s-1)^-1
#     - F_s_w = 5 W m^-2 incoming solar irradiance 1361
#     - I_0 = i_0 * (1 - alpha) * F_s_w
#     - PAR = 4.91 * I(z)

#     Parameters:
#     z (float): Depth in ice (m).

#     Returns:
#     float: PAR at the given depth in ice (W m^-2).
#     """
#     kappa = 1.5  # in m^-1
#     i_0 = 0.17  # for snow free ice and =0 for snow covered ice
#     alpha = 0.58  # in g C (g Chla h micro E m-2 s-1)-1
#     F_s_w = 5  # in W m-2 incoming solar irradiance 1361 https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://en.wikipedia.org/wiki/Solar_irradiance&ved=2ahUKEwiL5Z2Qw_6JAxV30wIHHWQXD4AQFnoECCEQAw&usg=AOvVaw3cpr_oSPFpcgj7Ny50mW76
#     I_0 = i_0 * (1 - alpha) * F_s_w
#     I = lambda z: I_0 * np.exp(-kappa * z)  # z is the depth in ice
#     PAR = lambda z: 4.91 * I(z)  # in W m-2

#     return PAR(z)

# z = np.linspace(0, 1, 150)
# PAR_arr = PAR(z)
# plt.plot(PAR_arr,z, color='black', linewidth=2)
# plt.xlabel(r'PAR [$\mu E$ $m^{-2}$ $s^{-1}$]')
# plt.ylabel(r'Depth [m]')
# plt.gca().invert_yaxis()
# plt.grid(True)
# plt.savefig('jupyternotebooks/PAR_along_depth.pdf', backend='pgf')
# plt.close()

# def Ek_light_z(s_min, s_max, t_min, t_max, c_min, c_max, z_min, z_max, k, alpha, index):
#     mu_m = 2e-5  # in s^-1
#     E = 0.5  # in micro E m-2 s-1
#     # alpha = 1e-4  # in g C (g Chla h micro E m-2 s-1)-1
#     r_chl_c_max = 0.05
#     r_chl_c_min = 0.01

#     ln, _, _ = nutrient_function(c_min, c_max, k)
#     fs = f_s_salinity(s_min, s_max)
#     ft = ft_temperature(t_min, t_max)
#     fs_ft_ln_ = fs * ft * ln
    
#     # r_chl_c_par = (  # noqa: E731
#     #     lambda _z: r_chl_c_max
#     #     - (r_chl_c_max - r_chl_c_min) * np.min([np.min(PAR_div_E(_z)), 1]) * ln
#     # )

#     r_chl_c_par = 0.05

#     P_m_par = mu_m * fs_ft_ln_ / r_chl_c_par  # noqa: E731

#     E_k =  P_m_par / alpha  # noqa: E731

#     return np.array(E_k)

# def lpar(index, s_min=20, s_max=250, t_min=-20, t_max=0, c_min=0, c_max=100, k=0, z_min=0, z_max=1, alpha=1e-4):
#     """
#     Calculate the photosynthetic rate as a function of depth in ice.
#     The equation used is mu = max_mu * fs * ft * ln * lpar, where:
#     - max_mu = 2e-5 in s^-1
#     - fs = exp(-((2.16 - a - b) ** 2)), where:
#         - a = 8.3 * 10 ** (-5) * s ** 2.11
#         - b = 0.55 * log(s)
#         - s is the salinity in g kg^-1
#     - ft = exp(r_g * T_c)
#     - r_g = 0.0633
#     - T_c is the temperature in degrees Celsius
#     - ln = c_n / (k_n + c_n)
#     - lpar = tanh(PAR / Ek)
#     - Ek = P_m / alpha
#     - P_m = mu_m * fs * ft * ln / r_chl_c_par
#     - r_chl_c_par = r_chl_c_max - (r_chl_c_max - r_chl_c_min) * min(PAR / E, 1) * ln
#     - E = 0.5 in micro E m-2 s-1
#     - alpha = 1e-4 in g C (g Chla h micro E m-2 s-1)-1
#     - mu = max_mu * fs * ft * ln * lpar

#     Parameters:
#     s_min (float): Minimum salinity in g kg^-1.
#     s_max (float): Maximum salinity in g kg^-1.
#     t_min (float): Minimum temperature in degrees Celsius.
#     t_max (float): Maximum temperature in degrees Celsius.
#     c_min (float): Minimum nutrient concentration in mmol m^-3.
#     c_max (float): Maximum nutrient concentration in mmol m^-3.
#     k (float): Nutrient concentration in mmol m^-3.
#     z_min (float): Minimum depth in ice in m.
#     z_max (float): Maximum depth in ice in m.
#     alpha (float): in g C (g Chla h micro E m-2 s-1)-1
#     index (int): index of the depth in ice.

#     Returns:
#     float: Photosynthetic rate at the given depth in ice.
#     """
#     PAR_arr = np.linspace(0, 30, 150)
#     Ek_arr = Ek_light_z(s_min, s_max, t_min, t_max, c_min, c_max, z_min, z_max, k, alpha, index)
#     lpar_arr = np.tanh(PAR_arr/ Ek_arr)

#     return lpar_arr

# lpar_bottom = lpar(index=149)
# plt.plot(np.linspace(0, 30, 150), lpar_bottom, color='black', linewidth=2, label='z=0m')
# # plt.plot(PAR_arr, lpar_mid, color='black', linestyle=':', linewidth=2, label='z=0.5m')
# # plt.plot(PAR_arr, lpar_bottom, color='black', linewidth=2, linestyle='--', label='z=1m')
# plt.ylabel(r'$L_{PAR}$')
# plt.xlabel(r'PAR [$\mu E$ $m^{-2}$ $s^{-1}$]')
# plt.grid(True)
# plt.savefig('jupyternotebooks/lpar_tanh.pdf', backend='pgf')
# plt.close()


# def photosynthetic_rate(lpar, fs, ft, ln):
#     max_mu = 2e-5  # in s^-1
#     mu = max_mu * lpar * fs * ft * ln

#     return mu

# mu = photosynthetic_rate(lpar_bottom, f_s, f_T, l_n)
# plt.plot(np.linspace(0, 30, 150), mu, color='black', linewidth=2)
# plt.ylabel(r'photosynthetic rate $\mu$')
# plt.xlabel(r'Salinity [ppt]')
# plt.grid(True)
# plt.savefig('jupyternotebooks/photosynthetic_rate.pdf', backend='pgf')
# plt.close()