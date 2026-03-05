src.spyice.models.stefan_problem
==========================================

This module provides analytical and semi-analytical solutions for the
Stefan phase-change problem in one-dimensional sea ice systems.

Both:

- One-phase Stefan formulation
- Two-phase (solid–liquid) formulation

are implemented.

The solutions are primarily used for:

- Model validation
- Benchmarking numerical solvers
- Comparing analytical and numerical ice growth


Mathematical Background
-----------------------

The classical one-phase Stefan problem describes phase boundary motion
governed by heat diffusion:

.. math::

   \frac{\partial T}{\partial t}
   = \alpha \frac{\partial^2 T}{\partial z^2}

with a moving boundary determined by the Stefan condition:

.. math::

   \rho L \frac{ds}{dt}
   = -k \frac{\partial T}{\partial z}

The analytical solution involves the similarity variable:

.. math::

   s(t) = 2 \lambda \sqrt{\alpha t}

where :math:`\lambda` is obtained from a transcendental equation.


Class: StefanProblem
====================

.. py:class:: StefanProblem()

   A utility class providing static methods to compute
   one-phase and two-phase Stefan solutions.


Static Methods
==============


stefan_problem
--------------

.. py:staticmethod:: stefan_problem(t, ui)

   Compute the one-phase Stefan solution depth.

   :param float t: Time.
   :param UserInput ui: User input configuration.
   :returns: Ice thickness depth.

   The solution uses a Newton method to determine
   the similarity parameter :math:`\lambda`.


calculate_temperature_profile
-----------------------------

.. py:staticmethod:: calculate_temperature_profile(depth_stefan, t, dz, nz, ui)

   Compute the temperature profile for the one-phase solution.

   :param float depth_stefan: Ice thickness.
   :param float t: Time.
   :param float dz: Grid spacing.
   :param int nz: Number of grid cells.
   :param UserInput ui: User input configuration.
   :returns: Temperature array.

   The profile is based on the error function solution:

   .. math::

      T(z,t)
      =
      T_b
      - (T_b - T_m)
      \frac{\operatorname{erf}(z / 2\sqrt{\alpha t})}
           {\operatorname{erf}(s / 2\sqrt{\alpha t})}


calculate_temperature_twophase_profiles
---------------------------------------

.. py:staticmethod:: calculate_temperature_twophase_profiles(depth_stefan, t, dz, nz, ui)

   Compute temperature and salinity profiles
   for the two-phase Stefan problem.

   :param float depth_stefan: Interface position.
   :param float t: Time.
   :param float dz: Grid spacing.
   :param int nz: Number of grid cells.
   :param UserInput ui: User input configuration.
   :returns: Tuple of:

      - Temperature array
      - Salinity array

   Uses complementary error functions for the liquid region.


stefan_problem_twophase
-----------------------

.. py:staticmethod:: stefan_problem_twophase(t, ui)

   Compute the moving interface for a two-phase
   solid–liquid system.

   :param float t: Time.
   :param UserInput ui: User configuration.
   :returns: Interface depth.

   The solution requires solving a nonlinear root
   problem via Newton iteration.


_plot_stefan_temp_twophase
--------------------------

.. py:staticmethod:: _plot_stefan_temp_twophase(z_depth=0.5)

   Utility method for plotting temperature evolution
   at a specified depth.

   :param float z_depth: Relative depth location.
   :returns: Temperature and salinity history arrays.

   Intended for validation and debugging.


Dependencies
============

- ``numpy``
- ``scipy.optimize``
- ``scipy.special.erfc``
- ``matplotlib``
- ``spyice.parameters.user_input``


Notes
=====

- Newton iteration tolerance is set to moderate precision.
- Analytical solutions assume constant thermophysical properties.
- Two-phase formulation accounts for salinity depression of melting temperature.
- Intended primarily for verification against the numerical solver.