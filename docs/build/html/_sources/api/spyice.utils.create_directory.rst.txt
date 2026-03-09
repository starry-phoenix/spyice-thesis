spyice.utils.create_directory
=================================

This module provides a utility function to **create an output directory** for storing simulation results such as temperature, salinity, or ice properties.  

Functions
=========

create_output_directory
-----------------------

.. py:function:: create_output_directory(hyd_dir, initial_salinity, boundary_condition_type, grid_resolution_dz, grid_timestep_dt, max_iterations, output_suffix)

   Creates an output directory based on simulation parameters.

   :param str hyd_dir: The base directory where the output directory will be created.
   :param float initial_salinity: Initial salinity of the system (used in folder naming).
   :param str boundary_condition_type: Type of boundary condition applied.
   :param float grid_resolution_dz: Spatial resolution of the simulation grid.
   :param float grid_timestep_dt: Time step of the simulation.
   :param int max_iterations: Maximum number of iterations in the simulation.
   :param str output_suffix: Suffix to append to the directory name (useful for versioning).
   :returns: Path of the created output directory.
   :rtype: str

   **Description:**  
   The function generates a directory name combining the input parameters.  
   If the directory does not exist, it is automatically created.

   **Example:**

   .. code-block:: python

       output_dir = create_output_directory(
           hyd_dir="results",
           initial_salinity=35.0,
           boundary_condition_type="top",
           grid_resolution_dz=0.01,
           grid_timestep_dt=60,
           max_iterations=1000,
           output_suffix="v1"
       )
       print(output_dir)
       # Output: results/Temperature_35.0_top_0.01_60_1000_v1