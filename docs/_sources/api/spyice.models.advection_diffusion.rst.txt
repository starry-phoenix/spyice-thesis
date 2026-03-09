src.spyice.models.advection_diffusion
==========================================

This module implements a 1D implicit finite-difference solver for the
advection–diffusion equation:

.. math::

   a \frac{\partial U}{\partial t}
   + b \frac{\partial U}{\partial z}
   + \frac{\partial}{\partial z}
     \left(c \frac{\partial U}{\partial z}\right)
   + d \frac{\partial W}{\partial t}
   = 0

The model supports both **temperature** and **salinity** evolution and
includes optional boundary-condition treatments (Stefan, Buffo, Voller).

---

Class: AdvectionDiffusion
==========================

.. py:class:: AdvectionDiffusion(argument, X, source, X_initial, W, W_initial, thickness_index, w, dt, dz, nz, t_passed, S_IC, Stefan=False, Buffo=False, Voller=False, bc_neumann=None)

   Class representing a 1D advection–diffusion model.

   :param str argument: Either ``"temperature"`` or ``"salinity"``.
   :param ndarray X: Current field values.
   :param float source: Source term.
   :param ndarray X_initial: Field values from previous timestep.
   :param ndarray W: Current liquid fraction.
   :param ndarray W_initial: Previous liquid fraction.
   :param int thickness_index: Index of ice thickness boundary.
   :param ndarray w: Convective brine velocity.
   :param float dt: Time step.
   :param float dz: Spatial step.
   :param int nz: Number of spatial grid cells.
   :param float t_passed: Simulation time.
   :param float S_IC: Initial salinity condition.
   :param bool Stefan: Enable Stefan boundary formulation.
   :param bool Buffo: Enable Buffo boundary formulation.
   :param bool Voller: Enable Voller enthalpy formulation.
   :param float bc_neumann: Optional Neumann boundary gradient.

   :raises AssertionError: If invalid argument is provided.
   :raises AssertionError: If both Stefan and Buffo are enabled.

---

Methods
=======

TDMAsolver
----------

.. py:method:: TDMAsolver(a, b, c, d)

   Solve a tridiagonal system using the Thomas algorithm.

   :param array a: Sub-diagonal.
   :param array b: Main diagonal.
   :param array c: Super-diagonal.
   :param array d: RHS vector.
   :returns: Solution vector.

Buffosolver
-----------

.. py:method:: Buffosolver(a, b, c, f)

   Modified tridiagonal solver enforcing diagonal dominance.

   :param ndarray a: Lower diagonal.
   :param ndarray b: Upper diagonal.
   :param ndarray c: Main diagonal.
   :param ndarray f: RHS vector.
   :returns: Solution vector.
   :raises AssertionError: If matrix is not diagonally dominant.

set_up_tridiagonal
------------------

.. py:method:: set_up_tridiagonal()

   Construct diagonal vectors for the implicit scheme matrix.
   Handles separate logic for temperature and salinity,
   including upwinding and boundary adjustments.

assemble_tridiagonal
--------------------

.. py:method:: assemble_tridiagonal()

   Assemble full matrix from diagonal vectors.

   :returns: Assembled system matrix.

modify_tridiagonal_voller_scheme
--------------------------------

.. py:method:: modify_tridiagonal_voller_scheme()

   Modify diagonal dominance for mushy-layer cells when
   using the Voller enthalpy scheme.

voller_X_array_set_zero_to_melt
-------------------------------

.. py:method:: voller_X_array_set_zero_to_melt(X, X_boundary)

   Force solution values above melting threshold to boundary value.

set_source_term
---------------

.. py:method:: set_source_term()

   Compute discretized source term depending on field type.

   :returns: Source contribution.

unknowns_matrix
---------------

.. py:method:: unknowns_matrix(temperature_melt, non_constant_physical_properties=False)

   Solve the implicit linear system for the updated field.

   :param float temperature_melt: Melting temperature.
   :returns:
      - **X_new** (ndarray) — Updated solution
      - **X_wind** (ndarray) — Wind-corrected solution
      - **A_before_correction** (ndarray) — System matrix

factor_1
--------

.. py:method:: factor_1(argument, a, c, dt, dz, nz)

   Compute diffusion discretization factor.

   :returns:
      - ndarray for temperature
      - tuple of ndarrays for salinity (central, plus, minus)

factor_2
--------

.. py:method:: factor_2(a, b, dt, dz, nz)

   Compute advection discretization factor.

   :returns: ndarray

factor_3
--------

.. py:method:: factor_3(a, d, nz)

   Compute phase-coupling factor.

   :returns: ndarray

---

Numerical Characteristics
==========================

- Implicit time discretization
- Finite difference spatial scheme
- Tridiagonal linear system
- Upwind advection (salinity)
- Optional enthalpy formulation (Voller)
- Stefan and Buffo boundary options

---

Dependencies
=============

- ``spyice.parameters.user_input``
- ``spyice.coefficients.update_coefficients``
- ``spyice.rhs.apply_boundary_condition``
- ``spyice.rhs.correct_for_brine_movement``

---

Notes
=====

Only one of ``Stefan`` or ``Buffo`` may be active at a time.

Diagonal dominance is required for the Buffo solver.

Zero-division protection is implemented in factor calculations.