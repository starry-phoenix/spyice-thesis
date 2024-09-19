import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pandas as pd

from src.spyice.parameters.results_params import ResultsParams
from src.spyice.parameters.user_input import UserInput
from src.spyice.postprocess.analysis import Analysis

ui = UserInput()

Tm_w, S_sw, iter_max, dz, cap_dens = (
    ui.temperature_melt,
    ui.boundary_salinity,
    ui.max_iterations,
    ui.grid_resolution_dz,
    ui.output_suffix,
)

np.seterr(divide="ignore", invalid="ignore")
# .style.use("spyice.utils.custom")
# plt.rcParams.update(
#     {
#         "text.usetex": True,
#     }
# )


class VisualiseModel:
    """A class for visualizing sea ice model results."""

    def __init__(
        self,
        user_input_dataclass: UserInput,
        results_dataclass: ResultsParams,
        error_analysis_dataclass: Analysis,
    ) -> None:
        """
        Args:
            user_input_dataclass (UserInput): An instance of the UserInput class containing user input data.
            results_dataclass (ResultsParams): An instance of the ResultsParams class containing results data.
            error_analysis_dataclass (Analysis): An instance of the Analysis class containing error analysis data.
        """
        self.results_object = results_dataclass
        self.error_analysis_object = error_analysis_dataclass
        self.ui_object = user_input_dataclass
        print("Visualisation object created...")

    def phi_slope(self, iteration):
        """Calculates the indices of the mushy regions based on the phi values.

        Args:
            iteration (int): The iteration number.
        Returns:
            numpy.ndarray: The indices of the mushy regions.
        """

        phi = self.results_object.phi_k_list[iteration]
        mush_idx = np.intersect1d(np.where(phi < 1.0), np.where(phi > 0.0))
        if mush_idx.size == 0 and phi.any() == 1.0:
            mush_idx = np.where(phi == 1.0)[0]
        elif mush_idx.size == 0 and phi.all() == 0.0:
            mush_idx = np.where(phi == 0.0)[0]
        return mush_idx

    def plot_error_temp_diff(self, zoom_x, savefig="True"):
        """Plots the temperature differences between two consecutive iterations.

        Args:
            zoom_x (int): The maximum value for the x-axis.
            savefig (str, optional): Indicates whether to save the figure or not. Defaults to "True".
        """

        print("Plotting Temperature differences between two consecutive iterations...")
        x_axis_iter = np.arange(0, zoom_x, 1)
        # fig, ax1 = plt.subplots(figsize=(10, 6))
        fig, ax1 = plt.subplots()
        # plt.grid()
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.error_analysis_object.t_k_stefan_diff_l2norm[x_axis_iter],
            "k",
            label="$L_2$",
        )
        ax1.set_yscale("log")
        self._extracted_from_plot_error_temp_8(
            ax1, r"$\lvert T_\text{Stefan} - T_k \rvert$"
        )
        ax1.tick_params(axis="y")
        ax1.legend()
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        color = "teal"
        ax2.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.results_object.thickness_list[: len(x_axis_iter)],
            label="Thickness in m",
            color=color,
            alpha=0.7,
            linestyle="dashed",
        )
        ax2.set_ylabel("Depth in m", color=color)
        ax2.tick_params(axis="y", labelcolor=color)
        ax1.autoscale()
        ax2.autoscale()
        # fig.tight_layout()

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/Temperature_error_diff_num_ana.pdf",
                backend="pgf",
            )

    def plot_error_temp(self, zoom_x: int, norm: str = "inf", savefig: bool = True):
        """Plots the temperature errors using the specified norm.

        Args:
            zoom_x (int): The maximum value for the x-axis.
            norm (str, optional): The norm to be used for plotting. Defaults to "inf".
            savefig (bool, optional): Whether to save the figure. Defaults to True.
        """

        print(f"Plotting Temperature errors using {norm} norm...")
        # fig1, (ax1) = plt.subplots(figsize=(10, 6))
        fig1, (ax1) = plt.subplots()
        # plt.grid()
        x_axis_iter = np.arange(0, zoom_x, 1)
        # plt.figure(figsize=(10,10))
        if norm == "inf":
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                self.error_analysis_object.t_stefan_diff_infnorm[x_axis_iter],
                "r--",
                label=r"$L_\infty$ Analytical",
            )
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                self.error_analysis_object.t_k_diff_infnorm[x_axis_iter],
                "k",
                label=r"$L_\infty$ Numerical",
            )
        elif norm == "2":
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                self.error_analysis_object.t_stefan_diff_l2norm[x_axis_iter],
                "r--",
                label="$L_2$ Analytical",
            )
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                self.error_analysis_object.t_k_diff_l2norm[x_axis_iter],
                "k",
                label="$L_2$ Numerical",
            )
        else:
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                self.error_analysis_object.t_stefan_diff_l1norm[x_axis_iter],
                "r--",
                label="$L_1$ Analytical",
            )
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                self.error_analysis_object.t_k_diff_l1norm[x_axis_iter],
                "k",
                label="$L_1$ Numerical",
            )
        self._extracted_from_plot_error_temp_8(
            ax1, "Temperature Error of Analytical and Numerical Results"
        )

        ax1.legend()
        color = "blue"
        ax3 = ax1.twinx()
        ax3.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.results_object.thickness_list[: len(x_axis_iter)],
            color=color,
            alpha=0.7,
        )
        # ax3.set_ylabel("Thickness in m")
        ax3.tick_params(axis="y", labelcolor=color)

        # fig1.tight_layout()

        if savefig:
            fig1.savefig(
                self.ui_object.dir_output_name + "/Temperature_error_num_and_ana",
                backend="pgf",
            )
        plt.close(fig1)

    # TODO Rename this here and in `plot_error_temp_diff` and `plot_error_temp`
    def _extracted_from_plot_error_temp_8(self, ax1, arg1):
        """Sets the x-label, y-label, and title of the given axis.

        Args:
            ax1 (matplotlib.axes.Axes): The axis to set the labels and title for.
            arg1 (str): The title of the plot.
        Returns:
            None
        """

        ax1.set_xlabel(r"$t$ [hours]")
        ax1.set_ylabel(r"$\Delta T$")
        ax1.set_title(arg1)

    def plot_depth_over_time(self, savefig: bool = False):
        """Plots the depth over time.

        Args:
            savefig (bool, optional): Whether to save the figure. Defaults to False.
        """

        print("Plotting Depth over time...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        depth_mush = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        mush_list_y1 = np.array(
            [
                [depth_mush[self.phi_slope(i)[0]], depth_mush[self.phi_slope(i)[-1]]]
                for i in x_axis_iter
            ]
        )
        # mush_list_y2 = [depth_mush[self.phi_slope(i)[-1]] for i in x_axis_iter]
        depth = np.array(self.results_object.thickness_list_buffo[: len(x_axis_iter)])
        # fig, ax = plt.subplots(figsize=(10, 6))
        fig, ax = plt.subplots()
        # ax.grid()
        ax.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.results_object.thickness_list[: len(x_axis_iter)],
            "r--",
            label="Numerical Depth",
        )
        ax.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.results_object.depth_stefan_all[: len(x_axis_iter)],
            "k",
            label="Analytical Depth",
        )
        if self.ui_object.is_buffo is True:
            ax.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                self.results_object.thickness_list_buffo[: len(x_axis_iter)],
                "b-.",
                alpha=0.5,
                label="Buffo Depth",
            )
        ax.fill_between(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            mush_list_y1[:, 0],
            mush_list_y1[:, 1],
            color="gray",
            alpha=0.2,
            label="Mushy Layer",
        )
        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.legend()
        ax.set_title(r"Numerical Depth Vs Analytical Depth")
        ax.set_yscale("log")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/Numerical_Analytical_Depth.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_depth_over_time_heatmap(self, savefig: bool = False):
        """Plots the depth over time.

        Args:
            savefig (bool, optional): Whether to save the figure. Defaults to False.
        """

        print("Plotting Heatmap Depth over time...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        depth_mush = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        mush_list_y1 = np.array(
            [
                [depth_mush[self.phi_slope(i)[0]], depth_mush[self.phi_slope(i)[-1]]]
                for i in x_axis_iter
            ]
        )
        heatmap_data = self.results_object.phi_k_list

        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="viridis",
            aspect="auto",
            interpolation="gaussian",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                self.results_object.depth_stefan_all[len(x_axis_iter) - 1],
                0,
            ],
            norm=colors.LogNorm(),
        )
        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        fig.colorbar(cax, ax=ax, label="Thickness [m]")
        # mush_list_y2 = [depth_mush[self.phi_slope(i)[-1]] for i in x_axis_iter]
        ax1 = ax.twinx()
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.results_object.depth_stefan_all[: len(x_axis_iter)],
            "k",
            label="Analytical Depth",
        )
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.results_object.thickness_list[: len(x_axis_iter)],
            "k--",
            label="Numerical Depth",
        )
        ax1.legend()
        ax1.invert_yaxis()
        ax1.set_yscale("log")
        ax1.set_ylim(
            self.results_object.depth_stefan_all[len(x_axis_iter) - 1],
            self.results_object.depth_stefan_all[0],
        )
        ax1.set_title(r"Numerical Depth Vs Analytical Depth")
        # ax1.set_yscale("log")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name
                + "/Numerical_Analytical_Depth_heatmap.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_temperature(
        self, z_depth: float, savefig: bool = True, Buffo_matlab: bool = False
    ):
        """Plots the temperature evolution at a given depth.

        Args:
            z_depth (float): The depth at which to plot the temperature evolution.
            savefig (bool, optional): Whether to save the figure. Defaults to True.
            Buffo_matlab (bool, optional): Whether to include Buffo-matlab data in the plot. Defaults to False.
        """

        print(f"Plotting Temperature evolution at {z_depth}m...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)[1:]
        # x_axis_iter = np.arange(0,22970,1)
        T_k_ = self.results_object.t_k_list[
            :, int(z_depth * self.ui_object.grid_resolution_dz)
        ]
        T_stefan_ = self.results_object.t_stefan_list[
            :, int(z_depth * self.ui_object.grid_resolution_dz)
        ]
        T_err = np.abs(T_k_ - T_stefan_)
        if self.ui_object.is_buffo is True:
            T_k_buffo_ = self.results_object.t_k_buffo_list[
                :, int(z_depth * self.ui_object.grid_resolution_dz)
            ]
            T_err_buffo = np.abs(T_k_ - T_k_buffo_)
        index = z_depth * self.ui_object.grid_resolution_dz
        if Buffo_matlab is True:
            df = pd.read_csv("MatlabData/temp_dirichletSalinity01.csv", sep=",")
            df_depth = np.array(
                df[
                    int(z_depth * self.ui_object.grid_resolution_dz) - 1 : int(
                        z_depth * self.ui_object.grid_resolution_dz
                    )
                ]
            ).reshape(25000, 1)

        # fig1, (ax1) = plt.subplots(figsize=(10, 6))
        fig1, (ax1) = plt.subplots()
        # plt.grid()
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            T_k_[x_axis_iter],
            "r--",
            label=r"Numerical Temperature",
        )
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            T_stefan_[x_axis_iter],
            "k",
            label=r"Analytical Temperature",
        )
        if self.ui_object.is_buffo is True:
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                T_k_buffo_[x_axis_iter],
                "b--",
                alpha=0.5,
                label=r"Buffo",
            )
        if Buffo_matlab is True:
            ax1.plot(
                x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
                df_depth[:24999],
                "b",
                alpha=0.2,
                label=r"Buffo-matlab",
            )
        ax1.set_xlabel(r"$t$ [hours]")
        ax1.set_ylabel(r"Temperature [K]")
        ax1.set_yscale("log")
        # ax1.legend()
        ax1.set_title(rf"Temperature evolution at {z_depth}m")
        color = "gray"
        ax3 = ax1.twinx()
        ax3.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            T_err[x_axis_iter],
            color=color,
            label=r"$\lvertT_\text{Numerical}-T_\text{Analytical}\rvert$",
            linestyle="dashed",
            linewidth=1,
        )
        ax3.tick_params(axis="y", labelcolor=color)
        ax3.set_ylabel(r"$T_k-T_Stefan$", color=color)
        ax3.set_yscale("log")
        ax1.legend()
        # fig1.tight_layout()
        if savefig:
            fig1.savefig(
                self.ui_object.dir_output_name
                + "/Temperature evolution at"
                + str(z_depth)
                + "m.pdf",
                backend="pgf",
            )
        plt.close(fig1)

        df_temp = pd.DataFrame(self.results_object.t_k_list)
        df_temp.to_csv(
            self.ui_object.dir_output_name
            + "/Temperature"
            + str(self.ui_object.grid_resolution_dz)
            + ".csv"
        )

    def plot_H_iter(self, h, t, param_name="Temperature", unit="K", savefig=False):
        iters = h.shape[0]
        iters_arr = np.linspace(0, iters - 1, iters)
        plt.grid()
        plt.plot(h[:, 0], label=r"cell Solid", color="turquoise")
        plt.scatter(iters_arr, h[:, 0], color="turquoise")
        plt.plot(h[:, 1], "--", label=r"cell Mushy", color="red", alpha=0.6)
        plt.scatter(iters_arr, h[:, 1], color="red", alpha=0.6)
        plt.plot(h[:, 2], ":", label=r"cell Liquid", color="teal")
        plt.scatter(iters_arr, h[:, 2], color="teal")
        plt.xlabel(r"iteration before convergence")
        plt.ylabel(rf"{param_name} [{unit}]")
        plt.title(rf"{param_name} at t={t}H")
        plt.legend()

        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name
                + "/"
                + param_name
                + "_iter"
                + self.ui_object.output_suffix
                + "_"
                + str(t)
                + "m.pdf",
                backend="pgf",
            )
        plt.close()

    def plot_H_iter_heatmap(
        self, h, p, t, param_name="Temperature", unit="K", savefig=False
    ):
        iters = h.shape[0]
        iters_arr = np.linspace(0, iters - 1, iters)
        plt.grid()
        Z = h.T
        fig1, (ax1) = plt.subplots()
        heatmap = ax1.imshow(
            Z,
            cmap="viridis",
            aspect="auto",
            interpolation="spline16",
            extent=[0, iters - 1, h[:, 0][-1], h[:, 2][-1]],
        )
        ax1.set_xlabel(r"iteration before convergence")
        ax1.set_ylabel(r"Liquid Fraction $\phi$")
        ax1.set_title(rf"{param_name} at t={t}H")
        ax1.set_yticks(
            [h[:, 2][-1], h[:, 1][-1], h[:, 0][-1]], ["Solid", "Mushy", "Liquid"]
        )
        fig1.colorbar(heatmap, ax=ax1, label=rf"{param_name} [{unit}]")
        ax1.grid(None)
        ax2 = ax1.twinx()
        ax2.plot(h[:, 0])
        ax2.scatter(iters_arr, h[:, 0])
        ax2.plot(h[:, 1], "--", label=r"cell Mushy", color="black", alpha=0.6)
        ax2.scatter(iters_arr, h[:, 1], color="black", alpha=0.6)
        ax2.plot(h[:, 2], ":")
        ax2.scatter(iters_arr, h[:, 2])
        ax2.set_ylim(h[:, 2][-1], h[:, 0][-1])
        ax2.invert_yaxis()
        ax2.grid(None)

        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name
                + "/"
                + param_name
                + "_iter"
                + self.ui_object.output_suffix
                + "_"
                + str(t)
                + "m_heatmap.pdf",
                backend="pgf",
            )
        plt.close(fig1)

    def plot_all_phi_mush(self, phi_mush, t, savefig=False):
        # plot all mush for len(phi_mush) iterations
        iters = phi_mush.shape[0]
        iters_arr = np.linspace(0, iters - 1, iters)
        plt.grid()
        plt.plot(phi_mush)
        plt.scatter(iters_arr, phi_mush)
        plt.xlabel(r"iteration before convergence")
        plt.ylabel(r"No. of mushy cells")
        plt.title(rf"Mushy Cells at t={t}H")

        if savefig:
            plt.savefig(
                f"{self.ui_object.dir_output_name}/LiquidFractionMush_iter_{cap_dens}_{str(t)}m.pdf",
                backend="pgf",
            )
        plt.close()

    def plot_H_iter_all(self, savefig=False):
        temperature_mushy_before_convergence = self.results_object.t_k_iter_all
        liquidfraction_mushy_before_convergence = self.results_object.phi_k_iter_all
        liquidfraction_before_convergence = self.results_object.all_phi_iter_all

        for h, p, t in zip(
            temperature_mushy_before_convergence,
            liquidfraction_mushy_before_convergence,
            [0.1, 0.5, 10, 100, 200, 300],
            strict=False,
        ):
            self.plot_H_iter(np.array(h), t, savefig=savefig)
            self.plot_H_iter_heatmap(np.array(h), np.array(p), t, savefig=savefig)

        for phi_mush, t in zip(
            liquidfraction_mushy_before_convergence,
            [0.1, 0.5, 10, 100, 200, 300],
            strict=False,
        ):
            self.plot_H_iter(
                np.array(phi_mush),
                t,
                param_name="Liquid-Fraction",
                unit="phi",
                savefig=savefig,
            )

        for phi, t in zip(
            liquidfraction_before_convergence,
            [0.1, 0.5, 10, 100, 200, 300],
            strict=False,
        ):
            self.plot_all_phi_mush(np.array(phi), t, savefig=savefig)

    # TODO: Add the following methods
    # def plot_phi(self, timestep, savefig=True):
    #     print(f"Plotting Liquid Fraction evolution at {timestep}H...")
    #     x_axis_iter = np.arange(0,self.nz,1)
    #     #x_axis_iter = np.arange(0,22970,1)
    #     index = timestep*3600/self.ui_object.grid_timestep_dt
    #     phi_k = self.phi_k_list[int(index)]
    #     if self.is_buffo is True:
    #         phi_buffo = self.phi_buffo_list[int(index)]

    #     fig1,(ax1) = plt.subplots(figsize=(10,6))
    #     plt.grid()
    #     ax1.plot(x_axis_iter*self.ui_object.grid_timestep_dt/3600, phi_k, 'r--',label='Numerical Temperature')
    #     if self.is_buffo is True:
    #         ax1.plot(x_axis_iter*self.ui_object.grid_timestep_dt/3600, phi_buffo, 'b--',alpha=0.5,label='Buffo')
    #     ax1.set_xlabel("depth in nodes")
    #     ax1.set_ylabel("Phi")
    #     #ax1.legend()
    #     ax1.set_title(f"Temperature evolution at {timestep}H")
    #     color = 'gray'
    #     ax1.legend()
    #     fig1.tight_layout()
    #     if savefig:
    #         fig1.savefig(self.dir_output_name + "/Liquid Fraction evolution at"+ str(timestep) +"s.png")

    # def plot_salinity(self, z_depth, savefig=True):
    #     x_axis_iter = np.arange(0, iter_max-1,1)
    #     index = z_depth*self.nz
    #     S_k = self.S_k_list[:,int(index)]
    #     if self.is_buffo is True:
    #         S_buffo = self.S_buffo_list[:,int(index)]
    #     fig1,(ax1) = plt.subplots(figsize=(10,6))
    #     plt.grid()
    #     ax1.plot(x_axis_iter*self.ui_object.grid_timestep_dt/3600, S_k[x_axis_iter], 'r--',label='Numerical Temperature')
    #     if self.is_buffo is True:
    #         ax1.plot(x_axis_iter*self.ui_object.grid_timestep_dt/3600, S_buffo[x_axis_iter], 'b--',alpha=0.5,label='Buffo')
    #     ax1.set_xlabel("t in hours")
    #     ax1.set_ylabel("Salinity in ppt")
    #     #ax1.legend()
    #     ax1.set_title(f"Salinity evolution at {z_depth}m")
    #     color = 'gray'
    #     ax1.legend()
    #     fig1.tight_layout()
    #     if savefig:
    #         fig1.savefig(self.dir_output_name + "/Salinity evolution at"+ str(z_depth) +"m.png")

    # def plot_enthalpy(self, timestep, savefig=False):
    #     print(f"Plotting Liquid Fraction Vs Temperature at {timestep}H...")
    #     x_axis_iter = np.arange(0, iter_max-1,1)
    #     index = timestep*3600/self.ui_object.grid_timestep_dt
    #     mush = self.phi_slope(int(index))
    #     phi = self.phi_k_list[int(index)]
    #     H_k = self.H_k_list[int(index)]
    #     H_solid = self.H_solid_list[int(index)]
    #     T_k = self.T_k_list[int(index)]
    #     H = (H_k - H_solid)/334774
    #     T_melt = Tm_w
    #     T_interface = self.Temp_interface[int(index)]
    #     fig1,(ax1) = plt.subplots(figsize=(10,6))
    #     plt.grid()
    #     ax1.plot(T_k,phi, 'r--',label='Phi')
    #     ax1.fill_betweenx(H, T_k[mush][0], T_k[mush][-1], color='gray', alpha=0.2, label='Mushy Layer')
    #     ax1.set_xlabel("Temperature in K")
    #     ax1.set_ylabel(r"Liquid Fraction $\phi$")
    #     ax1.set_title(f"Liquid Fraction Vs Temperature at {timestep}h")
    #     color = 'gray'
    #     ax3 = ax1.twinx()
    #     ax3.axvline(
    #         T_melt,
    #         color='b',
    #         linestyle='dashed',
    #         label=f'T_melt:{str(round(T_melt, 2))}',
    #     )
    #     ax3.axvline(
    #         T_interface,
    #         color='g',
    #         linestyle='dashed',
    #         label=f'T_interface:{str(round(T_interface, 2))}',
    #     )
    #     fig1.tight_layout()
    #     ax1.legend(loc=0)
    #     ax3.legend(loc=1)
    #     if savefig:
    #         fig1.savefig(self.dir_output_name + "/LiqVsTemp evolution at"+ str(timestep) +"h.png")
