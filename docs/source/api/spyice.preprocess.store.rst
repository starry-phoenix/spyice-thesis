src.spyice.preprocess.store
===========================

This module contains helper functions for initializing matrices,
storing simulation results, creating output directories,
and saving results to disk.

---

set_up_matrices
---------------

.. py:function:: set_up_matrices(iter_max, nz)

   Initialize matrices with one column for each time step.

   :param int iter_max: The maximum number of iterations (time steps).
   :param int nz: The number of vertical grid cells.
   :returns: A tuple containing initialized NumPy arrays:

      - **all_T** (*ndarray*) -- 2D array of shape ``(iter_max, nz)`` filled with zeros.
      - **all_S_sw** (*ndarray*) -- 2D array of shape ``(iter_max, nz)`` filled with zeros.
      - **all_phi** (*ndarray*) -- 2D array of shape ``(iter_max, nz)`` filled with zeros.
      - **all_H** (*ndarray*) -- 2D array of shape ``(iter_max, nz)`` filled with zeros.
      - **all_H_solid** (*ndarray*) -- 2D array of shape ``(iter_max, nz)`` filled with zeros.
      - **all_w** (*ndarray*) -- 2D array of shape ``(iter_max, nz)`` filled with zeros.
      - **all_thick** (*ndarray*) -- 1D array of length ``iter_max`` filled with zeros.
      - **all_t_passed** (*ndarray*) -- 1D array of length ``iter_max`` filled with zeros.

---

store_results
-------------

.. py:function:: store_results(T, S_sw, phi, H, H_solid, w, thickness, t_passed, all_T, all_S_sw, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed, t)

   Store simulation results for a given time step.

   :param numpy.ndarray T: Temperature values.
   :param numpy.ndarray S_sw: Salinity values.
   :param numpy.ndarray phi: Liquid fraction values.
   :param numpy.ndarray H: Enthalpy values.
   :param numpy.ndarray H_solid: Solid-state enthalpy values.
   :param numpy.ndarray w: Velocity values.
   :param float thickness: Ice thickness.
   :param float t_passed: Elapsed simulation time.
   :param numpy.ndarray all_T: Array storing temperature history.
   :param numpy.ndarray all_S_sw: Array storing salinity history.
   :param numpy.ndarray all_phi: Array storing liquid fraction history.
   :param numpy.ndarray all_H: Array storing enthalpy history.
   :param numpy.ndarray all_H_solid: Array storing solid enthalpy history.
   :param numpy.ndarray all_w: Array storing velocity history.
   :param numpy.ndarray all_thick: Array storing thickness history.
   :param numpy.ndarray all_t_passed: Array storing elapsed time history.
   :param int t: Current time step index.
   :returns: Tuple containing all updated arrays.

---

create_directory
----------------

.. py:function:: create_directory()

   Create a new directory using the current timestamp as its name.

   The directory is created inside a predefined parent directory.
   The generated path is appended to ``sys.path``.

   :returns: None

---

save_results
------------

.. py:function:: save_results(all_T, all_S, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed)

   Save simulation results to disk using compressed HKL files.

   :param list all_T: Temperature history.
   :param list all_S: Salinity history.
   :param list all_phi: Liquid fraction history.
   :param list all_H: Enthalpy history.
   :param list all_H_solid: Solid-state enthalpy history.
   :param list all_w: Velocity history.
   :param list all_thick: Thickness history.
   :param list all_t_passed: Time history.
   :returns: None
