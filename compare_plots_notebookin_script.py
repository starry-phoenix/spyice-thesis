# %% 

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as pd

np.seterr(divide="ignore", invalid="ignore")
# .style.use("spyice.utils.custom")
plt.style.use("jupyternotebooks/custom.mplstyle")
# plt.rcParams.update(
#     {
#         "text.usetex": True,
#     }
# )
# plt.rcParams["pgf.texsystem"] = "pdflatex"
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
plt.rcParams.update({'lines.linewidth': 2})


# %%

tempmushPT1_0p01_1000iter_voller1 = np.load(
    r"C:\Users\sneha\Documents\MBDHiWi\MBDHiwi\spyicedir\spyicedir\outputs\2024-09-27\093344_real_47.0_1000_0.01_buffo_pt1\temperature_mushy.npy",
    allow_pickle=True,
)
tempmushPT2_0p01_1000iter_voller1 = np.load(
    r"C:\Users\sneha\Documents\MBDHiWi\MBDHiwi\spyicedir\spyicedir\outputs\2024-09-27\094942_real_47.0_1000_0.01_buffo_pt2\temperature_mushy.npy",
    allow_pickle=True,
)

iter_arr = np.arange(0, len(tempmushPT2_0p01_1000iter_voller1[1]), 1)
tempmushPT1_0p01_1000iter_voller1_arr = np.array(tempmushPT1_0p01_1000iter_voller1[1])[
    :, 1
]
tempmushPT2_0p01_1000iter_voller1_arr = np.array(tempmushPT2_0p01_1000iter_voller1[1])[
    :, 1
]
step_arr = (
    np.ones(len(tempmushPT2_0p01_1000iter_voller1_arr))
    * tempmushPT1_0p01_1000iter_voller1_arr[-1]
)
pt1_array_length = len(tempmushPT1_0p01_1000iter_voller1_arr)
tempmushPT1_0p01_1000iter_voller1_arr = np.append(
    tempmushPT1_0p01_1000iter_voller1_arr, step_arr[pt1_array_length:]
)

plt.figure(figsize=(8, 5))
plt.plot(step_arr, label=r"$T_{melt}$", color="black")
plt.scatter(iter_arr, step_arr, color="black", linewidth=1)
plt.plot(
    tempmushPT1_0p01_1000iter_voller1_arr,
    label=r"$T_{mush}$ for $\phi_k = f(T_{k})$",
    color="black",
)
plt.scatter(iter_arr, tempmushPT1_0p01_1000iter_voller1_arr, color="black", linewidth=1)
plt.plot(
    tempmushPT2_0p01_1000iter_voller1_arr,
    label=r"$T_{mush}$ for $\phi_k = f(T_{k-1})$",
    color="black",
)
plt.scatter(iter_arr, tempmushPT2_0p01_1000iter_voller1_arr, color="black", linewidth=1)
plt.xlabel(r"Iterations")
plt.ylabel(r"Temperature [K]")
# plt.title(r"Interface temperature at t=0.5H before convergence")
plt.legend(loc="lower right")
plt.grid(True)
plt.savefig("Buffo_interfacetrack_response_dz0p01_time0p5h_test.pdf", backend="pgf",)
plt.show()


# # %%

# t_voller_1p0 = np.load(
#     "outputs/2024-09-27/124638_real_47.0_1000_0.01_voller1_pt2/temperature_mushy.npy",
#     allow_pickle=True,
# )[1]
# t_voller_0p6 = np.load(
#     "outputs/2024-09-27/131506_real_47.0_1000_0.01_voller0p6_pt2/temperature_mushy.npy",
#     allow_pickle=True,
# )[1]
# t_voller_0p4 = np.load(
#     "outputs/2024-09-27/132047_real_47.0_1000_0.01_voller0p4_pt2/temperature_mushy.npy",
#     allow_pickle=True,
# )[1]
# t_voller_0p8 = np.load(
#     "outputs/2024-09-27/132248_real_47.0_1000_0.01_voller0p8_pt2/temperature_mushy.npy",
#     allow_pickle=True,
# )[1]

# iter_arr = np.arange(0, len(t_voller_1p0), 1)
# t_voller_1p0_arr = np.array(t_voller_1p0)[:, 1]
# t_voller_0p8_arr = np.array(t_voller_0p8)[:, 1]
# t_voller_0p6_arr = np.array(t_voller_0p6)[:, 1]
# t_voller_0p4_arr = np.array(t_voller_0p4)[:, 1]
# max_array_length = 45 #len(t_voller_1p0_arr)
# step_arr = np.ones(max_array_length) * t_voller_1p0_arr[-1]
# t_voller_0p8_arr = np.append(t_voller_0p8_arr, step_arr[len(t_voller_0p8_arr) :])
# t_voller_0p6_arr = np.append(t_voller_0p6_arr, step_arr[len(t_voller_0p6_arr) :])
# t_voller_0p4_arr = np.append(t_voller_0p4_arr, step_arr[len(t_voller_0p4_arr) :])


# # %%

# plt.figure(figsize=(8,5))
# plt.plot(step_arr, label=r"$T_{melt}$", color="black", linewidth=1.5)
# # plt.scatter(iter_arr, step_arr, color="black", alpha=0.6)
# plt.plot(
#     t_voller_1p0_arr[:max_array_length],
#     label=r"$\omega = 1.0 $",
#     color="black",
# )
# # plt.scatter(
# #     iter_arr,
# #     t_voller_1p0_arr,
# #     color="red",
# #     alpha=0.6,
# # )
# plt.plot(
#     t_voller_0p8_arr,
#     label=r"$\omega = 0.8 $",
#     color="black",
# )
# # plt.scatter(
# #     iter_arr,
# #     t_voller_0p8_arr,
# #     color="blue",
# #     alpha=0.6,
# # )
# # plt.plot(
# #     t_voller_0p6_arr,
# #     label=r"$\omega = 0.6 $",
# #     color="black",
# #     linewidth=1.5
# # )
# # plt.scatter(
# #     iter_arr,
# #     t_voller_0p6_arr,
# #     color="green",
# #     alpha=0.6,
# # )
# plt.plot(
#     t_voller_0p4_arr,
#     label=r"$\omega = 0.4 $",
#     color="black",
# )
# # plt.scatter(
# #     iter_arr,
# #     t_voller_0p4_arr,
# #     color="orange",
# #     alpha=0.6,
# # )
# plt.xlabel(r"Iterations")
# plt.ylabel(r"Temperature [K]")
# # plt.title(r"$T_{mush}$ for $\phi_k = f(T_{k-1})$ at t=0.5H before convergence ")
# plt.legend()
# plt.grid()
# plt.savefig("Omega_value_find_Voller_interfacetrack_responsePT2_dz0p01_time0p5h.pdf", backend="pgf")
# plt.show()

# %%

# buffo_res_pt1 = np.load(
#     "outputs/2024-09-27/093344_real_47.0_1000_0.01_buffo_pt1/residuals.npy",
#     allow_pickle=True,
# )[1]
# buffo_res_pt2 = np.load(
#     "outputs/2024-09-27/094942_real_47.0_1000_0.01_buffo_pt2/residuals.npy",
#     allow_pickle=True,
# )[1]
# # voller_res_1p0_pt1 = np.load(
# #     "outputs/2024-09-27/122035_real_47.0_25000_0.01_voller1_pt1/residuals.npy",
# #     allow_pickle=True,
# # )
# # voller_res_1p0_pt2 = np.load(
# #     "outputs/2024-09-27/124638_real_47.0_1000_0.01_voller1_pt2/residuals.npy",
# #     allow_pickle=True,
# # )
# voller_res_1p4_pt1 = np.load(
#     "outputs/2024-09-27/130550_real_47.0_1000_0.01_voller1p4_pt1/residuals.npy",
#     allow_pickle=True,
# )[1]
# voller_res_0p4_pt2 = np.load(
#     "outputs/2024-09-27/132047_real_47.0_1000_0.01_voller0p4_pt2/residuals.npy",
#     allow_pickle=True,
# )[1]

# iter_arr = np.arange(0, len(buffo_res_pt2), 1)
# max_array_length = len(buffo_res_pt2)
# step_arr = np.ones(max_array_length) * 0.01
# buffo_res_pt1 = np.append(buffo_res_pt1, step_arr[len(buffo_res_pt1) :])
# voller_res_1p4_pt1 = np.append(voller_res_1p4_pt1, step_arr[len(voller_res_1p4_pt1) :])
# voller_res_0p4_pt2 = np.append(voller_res_0p4_pt2, step_arr[len(voller_res_0p4_pt2) :])

# plt.figure(figsize=(8, 5))
# plt.plot(
#     buffo_res_pt1,
#     label=r"Buffo; $\phi_k = f(T_{k})$",
#     color="black",
# )
# # plt.scatter(
# #     iter_arr,
# #     buffo_res_pt1,
# #     color="black",
# # )
# plt.plot(
#     buffo_res_pt2,
#     label=r"Buffo; $\phi_k = f(T_{k-1})$",
#     color="black",
# )
# # plt.scatter(
# #     iter_arr,
# #     buffo_res_pt2,
# #     color="black",
# # )
# plt.plot(
#     voller_res_1p4_pt1,
#     label=r"Voller; $\phi_k = f(T_{k}); \omega=0.4$",
#     color="black",
# )
# # plt.scatter(
# #     iter_arr,
# #     voller_res_1p4_pt1,
# #     color="black",
# # )
# plt.plot(
#     voller_res_0p4_pt2,
#     label=r"Voller; $\phi_k = f(T_{k-1}); \omega=0.4$",
#     color="black",
# )
# # plt.scatter(
# #     iter_arr,
# #     voller_res_0p4_pt2,
# #     color="black",
# # )

# plt.xlabel(r"Iterations")
# plt.ylabel(r"$\Sigma Abs(Res_p)$")
# # plt.title(r"$T_{mush}$ for $\phi_k = f(T_{k-1})$ at t=0.5H before convergence ")
# plt.grid(True)
# plt.legend()
# plt.savefig("Residuals_dz0p01_time0p5h.pdf", backend="pgf")
# plt.show()
# plt.close()


# %%

# path_to_dir = r"C:\Users\sneha\Documents\MBDHiWi\MBDHiwi\spyicedir\spyicedir\outputs\2025-08-24\liquidus_curves_comparison"
# diff_salinty_constdens_lin = r"C:\Users\sneha\Documents\MBDHiWi\MBDHiwi\spyicedir\spyicedir\outputs\2025-08-24\liquidus_curves_comparison\080418_real_47.0_1000_diff_salinitz_constdens_linliquidus\Temperature_S34_Dirichlet_0.01_47.0_1000_nonconst_dens-mushfix\Temperature0.01.csv"
# diff_salinity_constdens_frezchem = r"C:\Users\sneha\Documents\MBDHiWi\MBDHiwi\spyicedir\spyicedir\outputs\2025-08-24\liquidus_curves_comparison\094158_real_47.0_1000_diff_salinity_consdens_frezchemliqduidus\Temperature_S34_Dirichlet_0.01_47.0_1000_nonconst_dens-mushfix\Temperature0.01.csv"
# diff_salinity_constdens_faden = r"C:\Users\sneha\Documents\MBDHiWi\MBDHiwi\spyicedir\spyicedir\outputs\2025-08-24\liquidus_curves_comparison\185600_real_47.0_1000_0.01_diffusion_only_const_dens_fadenliquidus\Temperature_S34_Dirichlet_0.01_47.0_1000_nonconst_dens-mushfix\Temperature.csv"

# # Read CSVs
# df_lin = pd.read_csv(diff_salinty_constdens_lin, header=None)
# df_frezchem = pd.read_csv(diff_salinity_constdens_frezchem, header=None)
# df_faden = pd.read_csv(diff_salinity_constdens_faden, header=None)

# depth_ = np.append(np.arange(0, 1, 0.01), 1.0)
# # Plot at time index 5
# plt.figure(figsize=(6.5, 5))
# plt.plot(df_lin.iloc[5, 1:], depth_, label="Linear", color="black", linewidth=2)
# plt.plot(df_frezchem.iloc[5, 1:], depth_, label="Frezchem", color="black",  linewidth=2)
# plt.plot(df_faden.iloc[5, :], depth_, label="Faden", color="black",  linewidth=2)
# plt.gca().invert_yaxis()
# plt.xlabel(r'Temperature [K]')
# plt.ylabel(r"Depth [$m$]")
# plt.legend()
# plt.grid(True)
# plt.savefig("LiquidusCurvesComparison.pdf", backend="pgf")
# plt.show()