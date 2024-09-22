from __future__ import annotations

import os
from pathlib import Path

from omegaconf import DictConfig

from src.spyice.models.sea_ice_model import SeaIceModel
from src.spyice.postprocess.analysis import Analysis
from src.spyice.postprocess.visualise_model import VisualiseModel
from src.spyice.preprocess.pre_process import PreProcess


# TODO: Type hinting for all functions
# TODO: Return type hinting for all functions
class MainProcess:
    """Main class to run the model."""

    def __init__(
        self,
        config,
        hyd_output_dir: Path | str = Path(os.path.join(os.getcwd(), "outputs")),
        project_path=os.getcwd(),
    ):
        """
        Args:
            config: The configuration object.
            hyd_output_dir (Path | str): The directory path for the hydraulic output. Defaults to the 'outputs' directory in the current working directory.
            project_path (Path | str): The project path. Defaults to the current working directory.
        """
        self.config: DictConfig = config
        # self.config = ConfigSort.getconfig_dataclass(
        #     self.config_raw, config_type="default"
        # )
        self.project_path: Path | str = project_path
        self.out_dir_final = hyd_output_dir

    def run_model(self) -> None:
        """Runs the model using the provided configuration and output directory.

        Args:
            None
        Returns:
            None
        """
        # apply boundary and initial conditions during the pre-processing stage and get the pre-processed dataclass
        preprocess_data, userinput_data = PreProcess.get_variables(
            self.config, self.out_dir_final
        )
        # run the sea ice model and get the results dataclass

        results_data = SeaIceModel.get_results(preprocess_data, userinput_data)
        # error analysis of results and get the analysis dataclass

        analysis_data = Analysis.get_error_results(
            t_k_diff=results_data.t_k_diff, t_stefan_diff=results_data.t_stefan_diff
        )
        # plot the sea ice model using the user input, results, and analysis dataclasses
        self.plot_model(userinput_data, results_data, analysis_data)

    def plot_model(self, userinput_data, results_data, analysis_data):
        """Plots various visualizations of the model.

        Args:
            userinput_data (UserInputData): The user input data.
            results_data (ResultsData): The results data.
            analysis_data (AnalysisData): The error analysis data.
        Returns:
            None
        Raises:
            None
        """

        print("Postprocessing...")
        model_visualization_object = VisualiseModel(
            user_input_dataclass=userinput_data,
            results_dataclass=results_data,
            error_analysis_dataclass=analysis_data,
        )
        # model_visualization_object.plot_error_temp(100, norm="inf", savefig=False)
        # model_visualization_object.plot_depth_over_time(savefig=True)
        model_visualization_object.plot_depth_over_time_heatmap(savefig=True)
        # model_visualization_object.plot_temperature(
        #     z_depth=0.1, savefig=True, Buffo_matlab=False
        # )
        model_visualization_object.plot_H_iter_all(savefig=True)
        model_visualization_object.plot_temperature_heatmap(savefig=True)
        model_visualization_object.plot_temperature_heatmap_as_gif()
        print("Postprocessing done.")
