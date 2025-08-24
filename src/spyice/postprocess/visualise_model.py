import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # Add this import at the top of your file
from functools import partial
from pathlib import Path

from src.spyice.parameters.results_params import ResultsParams
from src.spyice.parameters.user_input import UserInput
from src.spyice.postprocess.analysis import Analysis

np.seterr(divide="ignore", invalid="ignore")
# .style.use("spyice.utils.custom")
plt.style.use("src.spyice.utils.custom")
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
        depth = self.results_object.depth_stefan_all[len(x_axis_iter) - 1]
        index = int(depth / self.ui_object.grid_resolution_dz)
        heatmap_data = self.results_object.phi_k_list[:, :index]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Blues",
            aspect="auto",
            interpolation="gaussian",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                index * self.ui_object.grid_resolution_dz,
                0.0,
            ],
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
            label="Voller Numerical Depth",
        )
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            self.results_object.thickness_list_buffo[: len(x_axis_iter)],
            "k:",
            label="Buffo Numerical Depth",
        )
        ax1.legend()
        ax1.set_ylim(
            index * self.ui_object.grid_resolution_dz,
            0.0,
        )
        ax1.set_yticks([])
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

    def plot_temperature_heatmap(self, savefig: bool = True):
        """Plots the temperature heatmap."""

        print("Plotting Temperature heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = self.results_object.t_k_list[1:]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Blues",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                1.0,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label="Temperature [K]")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/Temperature_heatmap.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_salinity_heatmap(self, savefig: bool = True):
        """Plots the temperature heatmap."""

        print("Plotting Temperature heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = self.results_object.s_k_list[1:]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Blues",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                1.0,
                0,
            ],
        )

                # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label="Salinity in ppt")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/Salinity_heatmap.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_liquidfraction_heatmap(self, savefig: bool = True):
        """Plots the temperature heatmap."""

        print("Plotting Temperature heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = self.results_object.phi_k_list[1:]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Blues",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                1.0,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"Liquid Fraction")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/Liquidfraction.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_temperature_heatmap_as_gif(self):
        fps = 10
        nseconds = 10
        depth = self.results_object.depth_stefan_all[self.ui_object.max_iterations - 1]
        index = int(depth / self.ui_object.grid_resolution_dz) + 1
        fig, (ax) = plt.subplots()
        data = self.results_object.t_k_list[:, :index]
        frames_plot = int(self.ui_object.max_iterations / 1)
        data = data.reshape(frames_plot, 1, index)
        hmap = ax.imshow(
            data[1].T,
            cmap="Blues",
            aspect="auto",
            extent=[
                0,
                5 * 1 * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [m]")

        def animation_function(i):
            hmap.set_data(data[10 * i + 1].T)
            hmap.set_extent(
                [
                    0,
                    (10 * i + 1) * self.ui_object.grid_timestep_dt / 3600,
                    depth,
                    0,
                ]
            )

            return [hmap]

        fig.colorbar(hmap, ax=ax, label="Temperature [K]")

        ani = animation.FuncAnimation(
            fig,
            animation_function,
            repeat=True,
            frames=(fps * nseconds) - 2,
            interval=1000 / fps,
        )

        # To save the animation using Pillow as a gif
        ani.save(self.ui_object.dir_output_name + "/test_new.gif", writer="imagemagick")
        # ani.save("test_new_html.html", writer="html")
        plt.close(fig)

    def plot_H_iter(
        self, h, t, s=None, param_name="Temperature", unit="K", savefig=False
    ):
        iters = h.shape[0]
        iters_arr = np.linspace(0, iters - 1, iters)
        if s is not None:
            s_melt = s[:, 1]
            t_melt = (
                -(9.1969758 * (1e-05) * s_melt**2) - 0.03942059 * s_melt + 272.63617665
            )
        plt.grid()
        plt.plot(h[:, 0], label=r"cell Solid", color="turquoise")
        plt.scatter(iters_arr, h[:, 0], color="turquoise")
        plt.plot(h[:, 1], "--", label=r"cell Mushy", color="red", alpha=0.6)
        plt.scatter(iters_arr, h[:, 1], color="red", alpha=0.6)
        plt.plot(h[:, 2], ":", label=r"cell Liquid", color="teal")
        plt.scatter(iters_arr, h[:, 2], color="teal")
        if s is not None:
            plt.plot(t_melt, label=r"Temperature Melt", color="black", alpha=0.6)
            plt.scatter(iters_arr, t_melt, color="black", alpha=0.6)
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

    def plot_H_iter_heatmap_mushy(
        self, h, p, temp_all, t, param_name="Temperature", unit="K", savefig=False
    ):
        iters = h.shape[0]
        iters_arr = np.linspace(0, iters - 1, iters)
        plt.grid()
        Z = h.T
        fig1, (ax1) = plt.subplots()
        heatmap = ax1.imshow(
            Z,
            cmap="Blues",
            aspect="auto",
            extent=[0, iters - 1, h[:, 2][-1], h[:, 0][-1]],
        )
        ax1.set_xlabel(r"iteration before convergence")
        ax1.set_ylabel(r"Depth [m]")
        ax1.set_title(rf"{param_name} at t={t}H")

        fig1.colorbar(heatmap, ax=ax1, label=rf"{param_name} [{unit}]")
        ax2 = ax1.twinx()
        ax2.plot(h[:, 0])
        ax2.scatter(iters_arr, h[:, 0])
        ax2.plot(h[:, 1], "--", label=r"cell Mushy", color="black", alpha=0.6)
        ax2.scatter(iters_arr, h[:, 1], color="black", alpha=0.6)
        ax2.plot(h[:, 2], ":")
        ax2.scatter(iters_arr, h[:, 2])
        ax2.set_ylim(h[:, 2][-1], h[:, 0][-1])
        ax2.set_yticks(
            [h[:, 2][-1], h[:, 1][0], h[:, 0][-1]], ["Solid", "Mushy", "Liquid"]
        )
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
                + "m_heatmap_mushonly.pdf",
                backend="pgf",
            )
        plt.close(fig1)

    def plot_H_iter_heatmap(
        self, h, p, temp_all, d, t, param_name="Temperature", unit="K", savefig=False
    ):
        iters = h.shape[0]
        iters_arr = np.linspace(0, iters - 1, iters)
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        depth = self.results_object.depth_stefan_all[len(x_axis_iter) - 1]
        depth_mush = (d[0][0] - 1) * self.ui_object.grid_resolution_dz
        index = int(depth / self.ui_object.grid_resolution_dz) + 1
        temp_mush_array = np.array(h[:, 1])
        temp_mush_normalised = (
            temp_mush_array
            / (np.sqrt(temp_mush_array.T * temp_mush_array))
            * depth_mush
        )
        plt.figure(figsize=(10, 8))
        plt.grid()
        fig1, (ax1) = plt.subplots()
        heatmap = ax1.imshow(
            temp_all[:, :index].T,
            cmap="Blues",
            aspect="auto",
            extent=[
                0,
                iters - 1,
                self.results_object.depth_stefan_all[len(x_axis_iter) - 1],
                0,
            ],
        )
        ax1.set_xlabel(r"iteration before convergence")
        ax1.set_ylabel(r"Depth")
        ax1.set_title(rf"{param_name} at t={t}H")
        ax1.set_yticks(
            [
                0.0,
                depth_mush,
                self.results_object.depth_stefan_all[len(x_axis_iter) - 1],
            ],
            ["Solid", "Mushy", "Liquid"],
        )
        fig1.colorbar(heatmap, ax=ax1, label=rf"{param_name} [{unit}]")
        ax2 = ax1.twinx()
        # ax2.plot(h[:, 0])
        # ax2.scatter(iters_arr, h[:, 0])
        ax2.plot(
            temp_mush_normalised, "--", label=r"cell Mushy", color="black", alpha=1.0
        )
        ax2.scatter(iters_arr, temp_mush_normalised, color="red", alpha=0.6)
        # ax2.plot(h[:, 2], ":")
        # ax2.scatter(iters_arr, h[:, 2])
        ax2.set_ylim(0, self.results_object.depth_stefan_all[len(x_axis_iter) - 1])
        # ax2.set_yticks([h[:, 2][-1], h[:, 0][-1]], ["Solid", "Liquid"])
        ax2.set_yticks([])
        # ax2.legend()
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
                f"{self.ui_object.dir_output_name}/LiquidFractionMush_iter_{self.ui_object.output_suffix}_{str(t)}m.pdf",
                backend="pgf",
            )
        plt.close()

    def plot_H_iter_all(self, savefig=False):
        temperature_mushy_before_convergence = self.results_object.t_k_iter_all
        liquidfraction_mushy_before_convergence = self.results_object.phi_k_iter_all
        liquidfraction_before_convergence = self.results_object.all_phi_iter_all
        temperature_all_before_convergence = (
            self.results_object.t_k_before_convergence_all
        )
        salinity_mushy_before_convergence = self.results_object.s_k_iter_all
        depth_all = self.results_object.mush_indx_list_all

        for h, p, temp_all, s, d, t in zip(
            temperature_mushy_before_convergence,
            liquidfraction_mushy_before_convergence,
            temperature_all_before_convergence,
            salinity_mushy_before_convergence,
            depth_all,
            [0.1, 0.5, 10, 100, 200, 300],
            strict=False,
        ):
            self.plot_H_iter(np.array(h), t, np.array(s), savefig=savefig)
            self.plot_H_iter_heatmap(
                np.array(h),
                np.array(p),
                np.array(temp_all),
                d,
                t,
                savefig=savefig,
            )

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

    def plot_response_pt1_pt2(self, tempmushPT1, tempmushPT2, savefig=False):
        iter_arr = np.arange(0, len(tempmushPT2[1]), 1)
        tempmushPT1_arr = np.array(tempmushPT1[1])[:, 1]
        tempmushPT2_arr = np.array(tempmushPT2[1])[:, 1]
        step_arr = np.ones(len(tempmushPT2_arr)) * tempmushPT1_arr[-1]
        pt1_array_length = len(tempmushPT1_arr)
        tempmushPT1_arr = np.append(tempmushPT1_arr, step_arr[pt1_array_length:])

        plt.figure(figsize=(12, 10))
        plt.grid()
        plt.plot(step_arr, "--", label=r"$T_{melt}$", color="black", alpha=0.6)
        plt.scatter(iter_arr, step_arr, color="black", alpha=0.6)
        plt.plot(
            tempmushPT1_arr,
            "--",
            label=r"$T_{mush} \text{ for } \phi_k = f(T_{k})$",
            color="blue",
            alpha=0.6,
        )
        plt.scatter(iter_arr, tempmushPT1_arr, color="blue", alpha=0.6)
        plt.plot(
            tempmushPT2_arr,
            "--",
            label=r"$T_{mush} \text{ for } \phi_k = f(T_{k-1})$",
            color="red",
            alpha=0.6,
        )
        plt.scatter(iter_arr, tempmushPT2_arr, color="red", alpha=0.6)
        plt.xlabel(r"No. of Iterations")
        plt.ylabel(r"Temperature [K]")
        plt.title(r"Interface temperature at t=0.5H before convergence")
        plt.legend()

        if savefig:
            plt.savefig(
                f"{self.ui_object.dir_output_name}/Temperature_mush_response_{self.ui_object.output_suffix}.pdf",
                backend="pgf",
            )

        plt.close()

    def plot_temperature_liquid_solid_evolution(self, z_depth: float, savefig=False):
        """Plots the temperature evolution at a given depth.

        Args:
            z_depth (float): The depth at which to plot the temperature evolution.
            savefig (bool, optional): Whether to save the figure. Defaults to True.
            Buffo_matlab (bool, optional): Whether to include Buffo-matlab data in the plot. Defaults to False.
        """

        print(f"Plotting Liqiudus-Solidus Temperature evolution at {z_depth}m...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)[1:]
        # x_axis_iter = np.arange(0,22970,1)
        T_k_liquid = self.results_object.temperature_liquid[
            :, int(z_depth * self.ui_object.grid_resolution_dz)
        ]
        T_k_solid = self.results_object.temperature_solid[
            :, int(z_depth * self.ui_object.grid_resolution_dz)
        ]
        index = z_depth * self.ui_object.grid_resolution_dz

        # fig1, (ax1) = plt.subplots(figsize=(10, 6))
        fig1, (ax1) = plt.subplots()
        # plt.grid()
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            T_k_solid[x_axis_iter],
            "k--",
            label=r"Solidus Temperature",
        )
        ax1.plot(
            x_axis_iter * self.ui_object.grid_timestep_dt / 3600,
            T_k_liquid[x_axis_iter],
            "k:",
            label=r"Liquidus Temperature",
        )

        ax1.set_xlabel(r"$t$ [hours]")
        ax1.set_ylabel(r"Temperature [K]")
        ax1.set_yscale("log")
        # ax1.legend()
        ax1.set_title(rf"Liqduius-Solidus Temperature evolution at {z_depth}m")
        color = "gray"
        ax1.legend()
        # fig1.tight_layout()
        if savefig:
            fig1.savefig(
                self.ui_object.dir_output_name
                + "/Liquidus-Solidus Temperature evolution at"
                + str(z_depth)
                + "m.pdf",
                backend="pgf",
            )
        plt.close(fig1)
    
    # TODO: place create_2d_array in utils script
    def create_2d_array(self, data):
        data_2d = np.zeros([self.ui_object.max_iterations, 101])

        for t in range(len(data)):
            index = int(self.results_object.thickness_index_total[t])
            data_2d[t, index] = data[t]

        return data_2d

    def plot_carbon_concentration(self, savefig: bool = True):
        carbon_concentration = self.create_2d_array(self.results_object.carbon_concentration)
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting carbon concentration heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = carbon_concentration[1:, :(depth_index)]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Greens",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"carbon concentration $[mmol/m^3]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/carbon_concentration.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_carbon_concentration_multiplelayers(self, savefig: bool = True):
        carbon_concentration = self.results_object.carbon_concentration_multiplelayers
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting carbon concentration multiple layers...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = carbon_concentration[1:, :(depth_index)]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Greens",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )

        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"carbon concentration $[mmol/m^3]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/carbon_concentration_multiplelayers.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_chla_bulk_concentration(self, savefig: bool = True):
        chla_bulk_concentration = self.create_2d_array(self.results_object.chla_bulk)
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting carbon concentration heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = chla_bulk_concentration[1:, :(depth_index)]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Greens",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"Chla bulk concentration $[mg/m^3]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/chla_bulk_concentration.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_chla_bulk_concentration_multiplelayers(self, savefig: bool = True):
        chla_bulk_concentration = self.results_object.chla_bulk_multiplelayers
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting chla concentration multiple layers...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = chla_bulk_concentration[1:, :(depth_index)]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Greens",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"Chla bulk concentration $[mg/m^3]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/chla_bulk_concentration_multiplelayers.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_nutrient_concentration(self, savefig: bool = True):
        nutrient_concentration = self.create_2d_array(self.results_object.nutrient_concentration)
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting nutrient concentration heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = nutrient_concentration[1:, :depth_index]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Greens",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )

        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"Nutrient concentration $[mmol/m^3]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/nutrient_concentration.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_nutrient_concentration_multiplelayers(self, savefig: bool = True):
        nutrient_concentration = self.results_object.nutrient_concentration_multiplelayers

        print("Plotting nutrient concentration of all layers: heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = nutrient_concentration[1:]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Greens",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                1.0,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"Nutrient concentration $[mmol/m^3]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/nutrient_concentration_alllayers.pdf",
                backend="pgf",
            )
        plt.close(fig)


    def plot_photosynthetic_rate(self, savefig: bool = True):
        photosynthetic_rate_T_S_PAR = self.create_2d_array(self.results_object.photosynthetic_rate)
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting photosynthetic rate heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = photosynthetic_rate_T_S_PAR[1:, :depth_index]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Reds",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"photosynthetic rate $[mmol/m^3/s]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/photosynthetic_rate.pdf",
                backend="pgf",
            )
        plt.close(fig)
    
    def plot_photosynthetic_rate_multiplelayers(self, savefig: bool = True):
        photosynthetic_rate_T_S_PAR = self.results_object.photosynthetic_rate_multiplelayers
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting photosynthetic rate multiple layers...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = photosynthetic_rate_T_S_PAR[1:, :depth_index]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Reds",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"photosynthetic rate $[mmol/m^3/s]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/photosynthetic_rate_multiplelayers.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_radiation_algae(self, savefig: bool = True):
        radiation_algae = self.create_2d_array(self.results_object.radiation_algae)

        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting algae radiation heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = radiation_algae[1:, :depth_index]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Reds",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"algae radiation $[W/m^2]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/algae_radiation.pdf",
                backend="pgf",
            )
        plt.close(fig)
    
    def plot_radiation_algae_multiplelayers(self, savefig: bool = True):
        radiation_algae = self.results_object.radiation_algae_multiplelayers

        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        print("Plotting algae radiation multiple layers...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = radiation_algae[1:, :depth_index]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Reds",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"algae radiation $[W/m^2]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/algae_radiation_multiplelayers.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_radiation_algae_dt_by_rho_c(self, savefig: bool = True):
        depth = self.results_object.thickness_list[-1]
        depth_index = int(self.results_object.thickness_index_total[-1]) + 1

        phi = self.results_object.phi_k_list[:, depth_index]
        rho_c_eff = self.ui_object.constants.rho_br*phi + self.ui_object.constants.rho_i*(1-phi)
        radiation_algae_by_rho_c_eff = self.results_object.radiation_algae / rho_c_eff
        radiation_algae_by_rho_c_eff_2d = self.create_2d_array(radiation_algae_by_rho_c_eff)

        print("Plotting algae radiation heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = radiation_algae_by_rho_c_eff_2d[1:, :depth_index]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Reds",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                depth,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"algal-induced radiative heating $[K s^{-1}]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/algae_radiation_heating.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_radiation_all(self, savefig: bool = True):
        radiation_algae = self.results_object.radiation_multiplelayers

        print("Plotting radiation ice and algae heatmap...")
        x_axis_iter = np.arange(0, self.ui_object.max_iterations - 1, 1)
        heatmap_data = radiation_algae[1:]
        fig, ax = plt.subplots()
        cax = ax.imshow(
            heatmap_data.T,
            cmap="Reds",
            aspect="auto",
            extent=[
                0,
                len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600,
                1.0,
                0,
            ],
        )
        # Add contour lines
        X = np.linspace(0, len(x_axis_iter) * self.ui_object.grid_timestep_dt / 3600, heatmap_data.shape[0])
        Y = np.linspace(0, 1.0, heatmap_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        contour = ax.contour(
            X,
            Y,
            heatmap_data.T,
            colors='k',
            linewidths=0.5,
        )
        ax.clabel(contour, fmt="%.1f", inline=True, fontsize=5, use_clabeltext=True)

        ax.set_xlabel(r"$t$ [hours]")
        ax.set_ylabel(r"Depth [$m$]")
        ax.set_title(r"Freezing Overtime")
        fig.colorbar(cax, ax=ax, label=r"radiation ice and algae $[W/m^2]$")

        if savefig:
            fig.savefig(
                self.ui_object.dir_output_name + "/radiation_algae_heatmap_multiplelayers.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_nutrient_cn_profile(self, savefig: bool = True):
        time = len(self.results_object.nutrient_concentration_multiplelayers)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)

        plt.figure(figsize=(4, 4))
        
        plt.grid()
        plt.plot(self.results_object.nutrient_concentration_multiplelayers[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.nutrient_concentration_multiplelayers[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.nutrient_concentration_multiplelayers[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'nutrient concentration [mmol m^{-3}]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/nutrient_concentration_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()            

    def plot_salinity_profile(self, savefig: bool = True):
        time = len(self.results_object.s_k_list)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        
        plt.figure()
        plt.grid()
        plt.plot(self.results_object.s_k_list[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.s_k_list[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.s_k_list[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'salinity concentration [ppt]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/salinity_concentration_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()   

    def plot_liquid_fraction_profile(self, savefig: bool = True):
        time = len(self.results_object.phi_k_list)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        
        plt.figure(figsize=(4, 4))
        plt.grid()
        plt.plot(self.results_object.phi_k_list[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.phi_k_list[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.phi_k_list[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'liquid fraction')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/liquidfraction_concentration_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()        

    def plot_temperature_profile(self, savefig: bool = True):
        time = len(self.results_object.t_k_list)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        
        plt.figure(figsize=(4, 4))
        plt.grid()
        plt.plot(self.results_object.t_k_list[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.t_k_list[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.t_k_list[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'Temperature [K]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/temperature_3_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()     

    def plot_carbon_concentration_profile(self, savefig: bool = True):

        time = len(self.results_object.carbon_concentration_multiplelayers)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        
        plt.figure(figsize=(4, 4))
        plt.grid()
        plt.plot(self.results_object.carbon_concentration_multiplelayers[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.carbon_concentration_multiplelayers[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.carbon_concentration_multiplelayers[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'carbon concentration [mmol m^{-3}]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/carbon_concentration_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()              

    
    def plot_salinity_sourceterm_profile(self, savefig: bool = True):

        time = len(self.results_object.salinity_source_term)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        
        plt.figure(figsize=(4, 4))
        plt.grid()
        plt.plot(self.results_object.salinity_source_term[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.salinity_source_term[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.salinity_source_term[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'Salinity source term [kg]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/salinitysourceterm_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()      

    def plot_radiation_profile(self, savefig: bool = True):

        time = len(self.results_object.radiation_multiplelayers)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        
        plt.figure(figsize=(4, 4))
        plt.grid()
        plt.plot(self.results_object.radiation_multiplelayers[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.radiation_multiplelayers[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.radiation_multiplelayers[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'radiation [W/m^2]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/radiation_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()    
    
    def plot_liquid_salinity_profile(self, savefig: bool = True):
        time = len(self.results_object.s_k_list)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        liquid_salinity = self.results_object.s_k_list/self.results_object.phi_k_list

        plt.figure()
        plt.grid()
        plt.plot(liquid_salinity[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(liquid_salinity[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(liquid_salinity[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'salinity concentration [ppt]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/liquidsalinity_concentration_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()   
    
    def plot_brinevelocity_profile(self, savefig: bool = True):

        time = len(self.results_object.brine_velocity_list)
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        
        plt.figure(figsize=(4, 4))
        plt.grid()
        plt.plot(self.results_object.brine_velocity_list[5, :], depth_, ':',label=r't=5s', color='black', linewidth=1)
        plt.plot(self.results_object.brine_velocity_list[int(time/2), :], depth_, label=rf't={int(time/2)}s', color='black', linewidth=1)
        plt.plot(self.results_object.brine_velocity_list[int(time -5), :], depth_, '--',label=rf't={time-5}s', color='black', linewidth=1)
        plt.gca().invert_yaxis()
        plt.xlabel(r'Brine velocity [m/s]')
        plt.ylabel(r'depth in [m]')
        plt.legend()
        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/brinevelocity_vertical_profile.pdf",
                backend="pgf",
            )
        plt.close()     

    def plot_temperature_3D(self, savefig:bool = True):
        time = len(self.results_object.t_k_list[1:])
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        t_k_arr = np.array(self.results_object.t_k_list[1:])

        # Take every 100th time step
        t_k_arr = t_k_arr[::100, :]
        dt = self.ui_object.grid_timestep_dt  # time step in seconds
        time_indices = np.arange(0, time, 100)
        time_hours = time_indices * dt / 3600  # convert to hours
        D, T = np.meshgrid(depth_, time_hours)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(
            D, T, t_k_arr, cmap='Blues', edgecolor='black', alpha=0.9
        )
        ax.set_ylabel(r'$t$ [hours]')
        ax.set_xlabel('Depth $[m]$')
        ax.set_zlabel('Temperature [K]')
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Temperature [K]')

        ax.set_xlim(ax.get_xlim()[::-1])
        ax.set_ylim(ax.get_ylim()[::-1])

        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/temperature_3d_surface_profile.pdf",
                backend="pgf",
            )
        plt.close(fig)

    def plot_temperature_3d_contours(self, savefig=True):
        time = len(self.results_object.t_k_list[1:])
        depth_ = np.append(np.arange(0, 1, self.ui_object.grid_resolution_dz), 1.0)
        t_k_arr = np.array(self.results_object.t_k_list[1:])

        # Take every 100th time step
        t_k_arr = t_k_arr[::100, :]
        time_indices = np.arange(0, time, 100)
        T, D = np.meshgrid(time_indices, depth_)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        # Choose contour levels (e.g., 10 evenly spaced between min and max)
        levels = np.linspace(np.nanmin(t_k_arr), np.nanmax(t_k_arr), 10)

        # Plot 3D contours
        contours = ax.contour3D(
            T, D, t_k_arr, levels=levels, cmap='viridis', linewidths=2
        )

        ax.set_ylabel('Depth [m]')
        ax.set_xlabel('Time step')
        ax.set_zlabel('Temperature [K]')
        ax.set_title('Temperature 3D contours (all time steps)')

        fig.colorbar(contours, ax=ax, shrink=0.5, aspect=10, label='Temperature [K]')

        ax.set_xlim(ax.get_xlim()[::-1])
        ax.set_ylim(ax.get_ylim()[::-1])

        if savefig:
            plt.savefig(
                self.ui_object.dir_output_name + "/temperature_3d_contour_profile.pdf",
                backend="pgf",
            )
        plt.close(fig)


    # def residual_plot(self):
    #     res_0p01_1000iter_voller1 = np.load("outputs/2024-09-26/123923_real_47.0_1000_0.01/residuals.npy", allow_pickle=True)

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
