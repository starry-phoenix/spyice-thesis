========================
Command Line Execution
========================

This section provides a step-by-step guide on how to run the SPyIce package
using the command line interface. It covers the necessary commands and options
to execute the main script and manage the simulation runs effectively.

Running SpyIce
==============

To run SpyIce from the command line:

.. code-block:: bash
   :caption: Terminal command

   python main_config.py

You can use the ``--help`` flag to see available options and arguments:

.. code-block:: bash
   :caption: Help command

   python main_config.py --help


Hydra Command-Line Overrides
============================

You can automate multiple runs with different configurations using Hydra's
command-line overrides for fast-run parameters.

Example: run the model with a specific time step and maximum iterations:

.. code-block:: bash
   :caption: Terminal command with overrides

   python main_config.py iter_max=iter_max1500 dt=dt47 dz=dz0p1 S_IC=S34


Fast-run parameters include:

- **iter_max**: Maximum number of iterations (e.g., ``iter_max1500``)
- **dt**: Time step size (e.g., ``dt47``)
- **dz**: Spatial step size (e.g., ``dz0p1``)
- **constants**: Debug mode or realistic mode (``real``)


Hydra Multirun
==============

You can run multiple jobs for all combinations of parameters using Hydra's
multirun feature (``-m`` flag).

Example 1: Keep ``dt`` fixed at 47 s

.. code-block:: bash
   :caption: Multirun example

   python main_config.py -m dt=dt47


Example 2: Keep ``iter_max`` fixed at 1500

.. code-block:: bash
   :caption: Multirun example

   python main_config.py -m iter_max=iter_max1500


Example 3: Keep ``dt`` and ``S_IC`` fixed

.. code-block:: bash
   :caption: Multirun example

   python main_config.py -m dt=dt47 S_IC=S34


Output Directory Structure
==========================

The outputs for each run will be located in the ``outputs`` directory:

.. code-block:: text

   YOUR_ROOT
   ├── ...
   ├── docs
   ├── outputs
   ├── src/spyice
   ├── ...


Adding New Parameters for Automation
=====================================

1. Create configuration files for the new parameter:

.. code-block:: text

   example/
   └── conf/
       └── parameter/
           ├── param1.yaml
           └── param2.yaml


2. Define values inside the YAML files.

Float example:

.. code-block:: yaml

   # param1.yaml
   parameter: 1.0


String argument example:

.. code-block:: yaml

   # param2.yaml
   parameter: "Argument"


3. Modify ``config_sort.py``

Add the new parameter inside the ``get_config_params()`` function
(similar to ``get_ownconfig_params()``):

.. code-block:: python

   def get_config_params():
       self.new_parameter = config.get(
           "parameter_directory", {}
       ).get("parameter_from_yaml", "default")
       return config_params


4. Modify ``pre_process.py``

Inside the ``PreprocessData`` class:

.. code-block:: python

   class PreprocessData:
       def __init__(self):
           ...

       if config_data.constants_type == "real":
           UserInput.__init__(
               self,
               new_parameter=config_data.new_parameter,
               ...
           )


Main Script Modifications
=========================

To modify the main execution workflow, refer to the implementation
inside ``main_config.py``.

.. note::

   Adding new parameters and modifying ``main_config.py`` is recommended
   only for developers.