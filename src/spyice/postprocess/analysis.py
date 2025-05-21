from dataclasses import dataclass

import numpy as np
from pathlib import Path
from src.spyice.utils.error_norms import ErrorNorms


@dataclass
class AnalysisData:
    """Represents the analysis data.

    Attributes:
        all_variables (dict): A dictionary containing all the variables.
    """

    all_variables: dict


class Analysis:
    """Represents an analysis object that performs error analysis on temperature differences."""

    def __init__(self, t_k_diff, t_stefan_diff):
        """
        Args:
            t_k_diff (float): The temperature difference in Kelvin.
            t_stefan_diff (float): The temperature difference in Stefan-Boltzmann units.
        """

        self.t_k_diff = t_k_diff
        self.t_stefan_diff = t_stefan_diff

    def set_analysis(self):
        """Sets up the analysis for the current object.

        This method initializes the `error_norms_object` attribute with an instance of the `ErrorNorms` class,
        passing the `t_k_diff` and `t_stefan_diff` attributes as arguments. It then calculates the numerical
        analytical difference using the `numerical_analytical_diff` method of the `error_norms_object`.

        Args:
            None
        Returns:
            None
        """

        self.error_norms_object = ErrorNorms(self.t_k_diff, self.t_stefan_diff)
        self.num_ana_temperature_diff = (
            self.error_norms_object.numerical_analytical_diff()
        )

    def error_analytical_numerical(self):
        """Calculates the errors of Numerical and Analytical using error norms one, two and infinity.

        Args:
            None

        Returns:
            tuple: A tuple containing the errors calculated using different error norms:
                - T_k_Stefan_diff_L1norm (float): L1 norm of the difference between T_k_list and T_Stefan_list.
                - T_k_Stefan_diff_infnorm (float): Infinity norm of the difference between T_k_list and T_Stefan_list.
                - T_k_Stefan_diff_L2norm (float): L2 norm of the difference between T_k_list and T_Stefan_list.
                - T_Stefan_diff_L1norm (float): L1 norm of the difference between consecutive T_Stefan values.
                - T_Stefan_diff_infnorm (float): Infinity norm of the difference between consecutive T_Stefan values.
                - T_Stefan_diff_L2norm (float): L2 norm of the difference between consecutive T_Stefan values.
                - T_k_diff_infnorm (float): Infinity norm of the difference between consecutive T_k values.
                - T_k_diff_L2norm (float): L2 norm of the difference between consecutive T_k values.
                - T_k_diff_L1norm (float): L1 norm of the difference between consecutive T_k values.
        """
        # Calculates the errors of Numerical and Analytical using error norms one, two and infinity
        # norm(|T_k_list - T_Stefan_list|)
        print("Calculating errors...")
        self.num_ana_temperature_diff = self.num_ana_temperature_diff
        T_k_Stefan_diff_L1norm, T_k_Stefan_diff_infnorm, T_k_Stefan_diff_L2norm = (
            self.calculate_errors(
                self.num_ana_temperature_diff, self.error_norms_object
            )
        )
        # Calculates the errors of Analytical temperature differences between two consecutive iterations
        # norm(|T_Stefan[i+1] - T_Stefan[i]|)
        T_Stefan_diff_L1norm, T_Stefan_diff_infnorm, T_Stefan_diff_L2norm = (
            self.calculate_errors(self.t_stefan_diff, self.error_norms_object)
        )
        # Calculates the errors of Numerical temperature differences between two consecutive iterations
        # norm(|T_k[i+1] - T_k[i]|)
        T_k_diff_infnorm, T_k_diff_L2norm, T_k_diff_L1norm = self.calculate_errors(
            self.t_k_diff, self.error_norms_object
        )

        self.store_field_errors(
            T_k_Stefan_diff_L1norm,
            T_k_Stefan_diff_infnorm,
            T_k_Stefan_diff_L2norm,
            T_Stefan_diff_infnorm,
            T_Stefan_diff_L2norm,
            T_Stefan_diff_L1norm,
            T_k_diff_infnorm,
            T_k_diff_L2norm,
            T_k_diff_L1norm,
        )

    def calculate_errors(self, field_array, error_norms_object):
        """Calculates the errors of a given field array using the provided error norms object.

        Args:
            field_array (numpy.ndarray): The field array to calculate the errors for.
            error_norms_object (ErrorNorms): The error norms object used to calculate the errors.
        Returns:
            tuple: A tuple containing the one-norm error, infinity-norm error, and two-norm error.
        """

        one_norm_error = error_norms_object.one_norm(field_array)
        infinity_norm_error = error_norms_object.infinity_norm(field_array)
        two_norm_error = error_norms_object.two_norm(field_array)
        return one_norm_error, infinity_norm_error, two_norm_error

    def store_field_errors(
        self,
        T_k_Stefan_diff_L1norm,
        T_k_Stefan_diff_infnorm,
        T_k_Stefan_diff_L2norm,
        T_Stefan_diff_infnorm,
        T_Stefan_diff_L2norm,
        T_Stefan_diff_L1norm,
        T_k_diff_infnorm,
        T_k_diff_L2norm,
        T_k_diff_L1norm,
    ):
        """Stores the field errors.

        Args:
            T_k_Stefan_diff_L1norm (float): The L1 norm of the temperature and concentration difference.
            T_k_Stefan_diff_infnorm (float): The infinity norm of the temperature and concentration difference.
            T_k_Stefan_diff_L2norm (float): The L2 norm of the temperature and concentration difference.
            T_Stefan_diff_infnorm (float): The infinity norm of the temperature difference.
            T_Stefan_diff_L2norm (float): The L2 norm of the temperature difference.
            T_Stefan_diff_L1norm (float): The L1 norm of the temperature difference.
            T_k_diff_infnorm (float): The infinity norm of the temperature and concentration difference.
            T_k_diff_L2norm (float): The L2 norm of the temperature and concentration difference.
            T_k_diff_L1norm (float): The L1 norm of the temperature and concentration difference.
        """
        self.T_k_Stefan_diff_L1norm = T_k_Stefan_diff_L1norm
        self.T_k_Stefan_diff_infnorm = T_k_Stefan_diff_infnorm
        self.t_k_stefan_diff_l2norm = T_k_Stefan_diff_L2norm
        self.t_stefan_diff_infnorm = T_Stefan_diff_infnorm
        self.t_stefan_diff_l2norm = T_Stefan_diff_L2norm
        self.t_stefan_diff_l1norm = T_Stefan_diff_L1norm
        self.t_k_diff_infnorm = T_k_diff_infnorm
        self.t_k_diff_l2norm = T_k_diff_L2norm
        self.t_k_diff_l1norm = T_k_diff_L1norm

    @classmethod
    def get_error_results(
        cls,
        t_k_diff: np.ndarray,
        t_stefan_diff: np.ndarray,
        residual: np.ndarray,
        temperature_mushy: np.ndarray,
        phi_mushy: np.ndarray,
        salinity_mushy: np.ndarray,
        output_dir: Path | str,
    ) -> AnalysisData:
        """Runs error analysis on the given temperature differences.

        Args:
            cls: The class object.
            t_k_diff: The temperature difference for k.
            t_stefan_diff: The temperature difference for Stefan.
        Returns:
            AnalysisData: An instance of the AnalysisData class containing the error analysis results.
        """

        print("Running error analysis...")
        error_analysis_object = cls(t_k_diff, t_stefan_diff)
        error_analysis_object.set_analysis()
        error_analysis_object.error_analytical_numerical()
        error_analysis_object.export_residuals(
            residual, temperature_mushy, phi_mushy, salinity_mushy, output_dir
        ) 
        all_vars = dict(vars(error_analysis_object))
        return cls.set_dataclass(all_vars, AnalysisData)

    @staticmethod
    def set_dataclass(data_to_be_converted: dict, dataclass: dataclass):
        """Sets the values of a dataclass object using a dictionary.

        Args:
            data_to_be_converted (dict): A dictionary containing the values to be set.
            dataclass (dataclass): The dataclass object to be updated.
        Returns:
            dataclass: The updated dataclass object.
        """

        for key, value in data_to_be_converted.items():
            setattr(dataclass, key, value)

        return dataclass

    def export_residuals(
        self, residuals: list, temperature_mushy, phi_mushy, salinity_mushy, output_dir
    ):
        """Exports the residuals to a file.

        Args:
            residuals (np.array): The residuals to export.
            output_dir (str): The output directory to save the residuals to.
        Returns:
            None
        """

        residuals = np.array(residuals, dtype=object)
        temperature_mushy = np.array(temperature_mushy, dtype=object)
        phi_mushy = np.array(phi_mushy, dtype=object)
        salinity_mushy = np.array(salinity_mushy, dtype=object)
        np.save(output_dir + "/residuals.npy", residuals)
        np.save(output_dir + "/temperature_mushy.npy", temperature_mushy)
        np.save(output_dir + "/phi_mushy.npy", phi_mushy)
        np.save(output_dir + "/salinity_mushy.npy", salinity_mushy)
        print("Residuals exported successfully.")
