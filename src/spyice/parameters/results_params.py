from __future__ import annotations

import numpy as np


class ResultsParams:
    """Class to store the results of the simulation."""

    def __init__(self, iter_max, nz):
        """
        Args:
            iter_max (int): The maximum number of iterations.
            nz (int): The number of depth levels.
        Attributes:
            error_temperature (ndarray): An array of size nz to store the temperature errors.
            error_temperature_sum (ndarray): An array of size iter_max to store the sum of temperature errors.
            error_temperature_sum_weighted (ndarray): An array of size iter_max to store the weighted sum of temperature errors.
            temp_grad (None): A placeholder for the temperature gradient.
            thickness_index_total (ndarray): An array of size iter_max to store the total thickness index.
            depth_stefan_all (ndarray): An array of size iter_max to store the depth Stefan values.
            t_stefan_prev (ndarray): An array of size nz to store the previous Stefan temperature values.
            t_k_prev (ndarray): An array of size nz to store the previous temperature values.
            t_stefan_diff (ndarray): A 2D array of size (iter_max, nz) to store the differences in Stefan temperature values.
            t_k_diff (ndarray): A 2D array of size (iter_max, nz) to store the differences in temperature values.
            t_stefan_list (ndarray): A 2D array of size (iter_max, nz) to store the Stefan temperature values.
            t_k_list (ndarray): A 2D array of size (iter_max, nz) to store the temperature values.
            t_k_buffo_list (ndarray): A 2D array of size (iter_max, nz) to store the buffered temperature values.
            thickness_list (ndarray): An array of size iter_max to store the thickness values.
            thickness_list_buffo (ndarray): An array of size iter_max to store the buffered thickness values.
            phi_k_list (ndarray): A 2D array of size (iter_max, nz) to store the phi values.
            phi_buffo_list (ndarray): A 2D array of size (iter_max, nz) to store the buffered phi values.
            s_k_list (ndarray): A 2D array of size (iter_max, nz) to store the s values.
            s_buffo_list (ndarray): A 2D array of size (iter_max, nz) to store the buffered s values.
            h_k_list (ndarray): A 2D array of size (iter_max, nz) to store the h values.
            h_solid_list (ndarray): A 2D array of size (iter_max, nz) to store the solid h values.
            temp_interface (ndarray): An array of size iter_max to store the interface temperatures.
            t_k_iter (list): A list to store the temperature values at each iteration.
            phi_k_iter (list): A list to store the phi values at each iteration.
            all_phi_iter (ndarray): A 2D array of size (iter_max, nz) to store all phi values at each iteration.
            t_k_iter_all (ndarray): A 2D array of size (iter_max, nz) to store all temperature values at each iteration.
            phi_k_iter_all (ndarray): A 2D array of size (iter_max, nz) to store all phi values at each iteration.
            all_phi_iter_all (ndarray): A 2D array of size (iter_max, nz) to store all phi values at each iteration.
            all_t (ndarray): A 2D array of size (iter_max, nz) to store all temperature values.
            all_s (ndarray): A 2D array of size (iter_max, nz) to store all s values.
            all_phi (ndarray): A 2D array of size (iter_max, nz) to store all phi values.
            all_h (ndarray): A 2D array of size (iter_max, nz) to store all h values.
            all_h_solid (ndarray): A 2D array of size (iter_max, nz) to store all solid h values.
            all_w (ndarray): A 2D array of size (iter_max, nz) to store all w values.
            all_thick (ndarray): An array of size iter_max to store all thickness values.
            all_t_passed (ndarray): An array of size iter_max to store all passed temperature values.
        """
        # code implementation goes here

        self.error_temperature = np.zeros(nz)
        self.error_temperature_sum = np.zeros(iter_max)
        self.error_temperature_sum_weighted = np.zeros(iter_max)
        self.temp_grad = None
        self.thickness_index_total = np.zeros(iter_max)
        self.depth_stefan_all = np.zeros(iter_max)
        self.t_stefan_prev = np.zeros(nz)
        self.t_k_prev = np.zeros(nz)
        self.t_stefan_diff = np.zeros((iter_max, nz))
        self.t_k_diff = np.zeros((iter_max, nz))
        self.t_stefan_list = np.zeros((iter_max, nz))
        self.t_k_list = np.zeros((iter_max, nz))
        self.t_k_buffo_list = np.zeros((iter_max, nz))
        self.thickness_list = np.zeros(iter_max)
        self.thickness_list_buffo = np.zeros(iter_max)
        self.phi_k_list = np.zeros((iter_max, nz))
        self.phi_buffo_list = np.zeros((iter_max, nz))
        self.s_k_list = np.zeros((iter_max, nz))
        self.s_buffo_list = np.zeros((iter_max, nz))
        self.h_k_list = np.zeros((iter_max, nz))
        self.h_solid_list = np.zeros((iter_max, nz))
        self.temp_interface = np.zeros(iter_max)
        self.t_k_iter = []
        self.phi_k_iter = []
        self.all_phi_iter = np.zeros((iter_max, nz))
        self.t_k_iter_all = np.zeros((iter_max, nz))
        self.phi_k_iter_all = np.zeros((iter_max, nz))
        self.all_phi_iter_all = np.zeros((iter_max, nz))
        self.all_t = np.zeros([iter_max, nz])
        self.all_s = np.zeros([iter_max, nz])
        self.all_phi = np.zeros([iter_max, nz])
        self.all_h = np.zeros([iter_max, nz])
        self.all_h_solid = np.zeros([iter_max, nz])
        self.all_w = np.zeros([iter_max, nz])
        self.all_thick = np.zeros([iter_max])
        self.all_t_passed = np.zeros(iter_max)

    @staticmethod
    def store_results(
        results_dataclass, temp, s_sw, phi, h, h_solid, w, thickness, t_passed, t
    ):
        """Stores the results of the simulation in the given results_dataclass.

        Args:
            results_dataclass: An instance of the ResultsDataClass where the results will be stored.
            temp: The temperature values to be stored.
            s_sw: The s_sw values to be stored.
            phi: The phi values to be stored.
            h: The h values to be stored.
            h_solid: The h_solid values to be stored.
            w: The w values to be stored.
            thickness: The thickness value to be stored.
            t_passed: The t_passed value to be stored.
            t: The t value to be stored.
        Returns:
            The updated results_dataclass with the stored results.
        """

        results_dataclass.all_t[t, :] = temp
        results_dataclass.all_s[t, :] = s_sw
        results_dataclass.all_phi[t, :] = phi
        results_dataclass.all_h[t, :] = h
        results_dataclass.all_h_solid[t, :] = h_solid
        results_dataclass.all_w[t, :] = w
        results_dataclass.all_thick[t] = thickness
        results_dataclass.all_t_passed[t] = t_passed

        return results_dataclass

    @staticmethod
    def store_results_for_iter_t(
        results_dataclass,
        t,
        thickness_index,
        t_k,
        t_stefan,
        s_k,
        s_k_buffo,
        phi_k,
        phi_k_buffo,
        h_k,
        h_solid,
        thickness,
        thickness_buffo,
        thickness_stefan,
        t_k_buffo,
        buffo=False,
    ):
        """Stores the results for a given iteration 't' in the 'results_dataclass'.

        Args:
            results_dataclass: An instance of the results dataclass.
            t: The iteration number.
            thickness_index: The index of the thickness.
            t_k: The temperature values.
            t_stefan: The Stefan temperature values.
            s_k: The entropy values.
            s_k_buffo: The entropy values for buffo.
            phi_k: The phase fraction values.
            phi_k_buffo: The phase fraction values for buffo.
            h_k: The enthalpy values.
            h_solid: The solid enthalpy values.
            thickness: The thickness values.
            thickness_buffo: The thickness values for buffo.
            thickness_stefan: The Stefan thickness values.
            t_k_buffo: The temperature values for buffo.
            buffo: A boolean indicating if buffo is enabled (default is False).

        Returns:
            The updated 'results_dataclass' with the stored results.
        """
        results_dataclass.error_temperature[:thickness_index] = np.abs(
            np.abs(t_k[:thickness_index]) - np.abs(t_stefan[:thickness_index])
        )
        results_dataclass.error_temperature_sum[t] = np.sum(
            results_dataclass.error_temperature
        )
        results_dataclass.error_temperature_sum_weighted[t] = (
            results_dataclass.error_temperature_sum[t] / thickness_index
        )
        results_dataclass._update_diff_lists(t, t_k, t_stefan)
        results_dataclass._append_to_lists(
            t,
            t_k,
            t_stefan,
            s_k,
            s_k_buffo,
            phi_k,
            h_k,
            h_solid,
            thickness,
            thickness_buffo,
            thickness_stefan,
            t_k_buffo,
            phi_k_buffo,
            buffo,
        )
        results_dataclass.temp_interface[t] = t_k[thickness_index]

        return results_dataclass

    def _update_diff_lists(self, t, t_k, t_stefan):
        """Updates the difference lists for temperature values.

        Args:
            t (int): The current time step.
            t_k (ndarray): The temperature values at time step t.
            t_stefan (ndarray): The Stefan temperature values at time step t.

        Returns:
            None
        """
        self.t_stefan_diff[t, :] = t_stefan - self.t_stefan_prev
        self.t_stefan_prev = t_stefan
        self.t_k_diff[t, :] = t_k - self.t_k_prev
        self.t_k_prev = t_k

    def _append_to_lists(
        self,
        t,
        t_k,
        t_stefan,
        s_k,
        s_k_buffo,
        phi_k,
        h_k,
        h_solid,
        thickness,
        thickness_buffo,
        thickness_stefan,
        t_k_buffo,
        phi_k_buffo,
        buffo,
    ):
        """Appends the given parameters to the respective lists.

        Args:
            t (int): The index of the time step.
            t_k (numpy.ndarray): The temperature values.
            t_stefan (numpy.ndarray): The Stefan temperature values.
            s_k (numpy.ndarray): The solid fraction values.
            s_k_buffo (numpy.ndarray): The solid fraction values for buffo.
            phi_k (numpy.ndarray): The phase fraction values.
            h_k (numpy.ndarray): The enthalpy values.
            h_solid (numpy.ndarray): The solid enthalpy values.
            thickness (float): The thickness value.
            thickness_buffo (float): The thickness value for buffo.
            thickness_stefan (float): The Stefan thickness value.
            t_k_buffo (numpy.ndarray): The temperature values for buffo.
            phi_k_buffo (numpy.ndarray): The phase fraction values for buffo.
            buffo (bool): Indicates whether buffo is present or not.
        """
        self.t_k_list[t, :] = t_k
        self.t_stefan_list[t, :] = t_stefan
        if buffo:
            self.t_k_buffo_list[t, :] = t_k_buffo
            self.thickness_list_buffo[t] = thickness_buffo
            self.phi_buffo_list[t, :] = phi_k_buffo
        self.thickness_list[t] = thickness
        self.depth_stefan_all[t] = thickness_stefan
        self.phi_k_list[t, :] = phi_k
        self.s_k_list[t, :] = s_k
        self.s_buffo_list[t, :] = s_k_buffo
        self.h_k_list[t, :] = h_k
        self.h_solid_list[t, :] = h_solid
