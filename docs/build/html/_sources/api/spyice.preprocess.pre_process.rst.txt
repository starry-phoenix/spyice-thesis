src.spyice.preprocess.pre_process
==========================================

This module handles **preprocessing of input data** for the sea ice model. It sets up:

- User configuration and input data.
- Geometry and grid settings.
- Initial and boundary conditions.
- Enthalpy calculations (solid and bulk).
- Ice thickness initialization.

Classes
=======

PreprocessData
--------------

.. py:class:: PreprocessData

   Data class representing preprocessing status.

   :ivar bool is_preprocessing: Flag indicating whether preprocessing is enabled. Defaults to True.


PreProcess
----------

.. py:class:: PreProcess(UserInput, GeometrySettings, ResultsParams)

   Main class for preprocessing input data before simulation.

   **Initialization**:

   .. py:method:: __init__(self, constants_type, config_data, output_dir)

      Initialize preprocessing with constants, geometry, and results parameters.

      :param str constants_type: Type of constants ("real" or other).
      :param ConfigData config_data: Configuration data object.
      :param str output_dir: Output directory.
      :returns: None

Methods
-------

preprocess
~~~~~~~~~~

.. py:method:: preprocess(self)

   Perform preprocessing steps:

   - Sets up iteration time tracking.
   - Initializes temperature, salinity, liquid fraction, and velocity arrays.
   - Computes solid and total enthalpy.
   - Initializes ice thickness.

   :returns: None


get_variables
~~~~~~~~~~~~~

.. py:classmethod:: get_variables(cls, config, out_dir_final: Path | str) -> tuple[PreprocessData, UserInput]

   Retrieve variables and user input data after preprocessing.

   :param cls: Class object.
   :param config: Configuration object.
   :param out_dir_final: Output directory.
   :returns: Tuple containing the dataclass with filtered variables and a UserInput object.


get_userinput
~~~~~~~~~~~~~

.. py:method:: get_userinput(self) -> UserInput

   Returns a `UserInput` object containing core input attributes such as constants, grid timestep, initial salinity, output directory, and maximum iterations.

   :returns: UserInput object


set_dataclass
~~~~~~~~~~~~~

.. py:staticmethod:: set_dataclass(data_to_be_converted: dict, dataclass: dataclass) -> dataclass

   Sets attributes of a dataclass object using a dictionary of values.

   :param dict data_to_be_converted: Dictionary of attributes to set.
   :param dataclass dataclass: Dataclass object to update.
   :returns: Updated dataclass object


Functions
=========

set_up_iter
-----------

.. py:function:: set_up_iter(iter_max, grid_timestep_dt)

   Compute and print iteration setup for the simulation.

   :param int iter_max: Maximum number of iterations.
   :param float grid_timestep_dt: Time step for the grid.
   :returns: 0 (placeholder return value)