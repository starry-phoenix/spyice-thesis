src.spyice.models.sea_ice_model
==========================================

Overview
--------

The ``SeaIceModel`` module implements the full numerical framework
for simulating thermodynamic and biogeochemical evolution of sea ice.

It includes:

- Ice–ocean interface tracking
- Convergence iteration loops
- Phase change handling (Stefan problem)
- Radiation and salinity source terms
- Algae and biogeochemical coupling
- Residual tracking and mushy layer diagnostics
- Results storage and export


Functions
---------

locate_ice_ocean_interface
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: spyice.models.sea_ice_model.locate_ice_ocean_interface


Classes
-------

SeaIceModel
^^^^^^^^^^^

.. autoclass:: spyice.models.sea_ice_model.SeaIceModel
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource


Main Methods
------------

The most important high-level methods are:

``run_sea_ice_model()``
   Executes the full time-stepping simulation.

``convergence_loop_iteration()``
   Performs a full nonlinear convergence iteration for a given timestep.

``run_while_convergence_iteration()``
   Core iterative solver loop until tolerance is met.

``calculate_source_terms()``
   Computes radiation, salinity, and algae-related source terms.

``choose_phase_type_iteration()``
   Selects between one-phase and two-phase Stefan problem.

``get_results()``
   Class method to run the model and return a ``ResultsParams`` object.


Dependencies
------------

This module relies on:

- ``spyice.models.algae_model``
- ``spyice.models.radiative_model``
- ``spyice.models.stefan_problem``
- ``spyice.statevariables``
- ``spyice.update_physical_values``
- ``spyice.utils.helpers``


Notes
-----

- Convergence is controlled by user-defined tolerances.
- A ``ConvergenceError`` is raised if iteration exceeds the counter limit.
- An ``InvalidPhaseError`` is raised if an unsupported Stefan phase type is selected.
- Mushy layer diagnostics are recorded during Stefan validation runs.
- Results are stored incrementally at each timestep.

``SeaIceModel.get_results(...)`` is the recommended entry point for external use.