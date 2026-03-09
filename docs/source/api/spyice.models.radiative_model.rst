src.spyice.models.radiative_model
=============================================

Module: ``spyice.models.brine_flux_radiation``

This module contains functions for modeling **radiative fluxes, brine fluxes, and salinity transport** in sea ice. It is designed for process-based simulations of ice-ocean interactions, focusing on:

- Radiative source terms from algae, ice, and organic matter.
- Calculation of Photosynthetically Active Radiation (PAR) in ice.
- Brine flux driven by Rayleigh number instabilities.
- Salinity source terms from brine advection.
- Permeability and local Rayleigh number calculations.

Functions
=========

radiative_source_term
----------------------

.. py:function:: radiative_source_term(radiative_algae=0.0, radiative_ice=0.0, radiative_organicmatter=0.0)

   Compute the total radiative source term for the ice-algae system.

   :param float radiative_algae: Radiative contribution from algae.
   :param float radiative_ice: Radiative contribution from ice.
   :param float radiative_organicmatter: Radiative contribution from organic matter.
   :returns: Total radiative source term.


radiative_ice
-------------

.. py:function:: radiative_ice(depth)

   Compute radiative attenuation through ice.

   :param float or np.ndarray depth: Depth in ice (m).
   :returns: Radiative flux contribution from ice.

   Formula:

   .. math::

      I_0 = i_0 (1 - \alpha) F_{sw}, \quad
      I(z) = I_0 \exp(-\kappa z), \quad
      R_\text{ice} = I(z) \kappa


radiative_organicmatter
------------------------

.. py:function:: radiative_organicmatter(depth)

   Compute radiative flux attenuation due to organic matter in ice.

   :param float or np.ndarray depth: Depth in ice (m).
   :returns: Radiative flux contribution from organic matter.


calculate_radiative_terms
--------------------------

.. py:function:: calculate_radiative_terms(depth, thickness_index, radiative_algae=0.0, algae_model_depth_type="single")

   Combine radiative contributions from algae, ice, and organic matter for the given depth.

   :param np.ndarray depth: Depth array in ice.
   :param int thickness_index: Index corresponding to interface layer.
   :param float radiative_algae: Radiative flux from algae.
   :param str algae_model_depth_type: 'single' or 'all' for biologically active layer modeling.
   :returns: Total radiative source term.


calculate_local_rayleigh_number
-------------------------------

.. py:function:: calculate_local_rayleigh_number(thickness_index, thickness, salinity, phi, grid_size)

   Calculate the **local Rayleigh number** for brine convection in sea ice.

   :param int thickness_index: Index at the ice-ocean interface.
   :param np.ndarray thickness: Ice thickness array.
   :param np.ndarray salinity: Salinity profile in sea ice.
   :param np.ndarray phi: Porosity profile.
   :param float grid_size: Grid spacing (dz).
   :returns: Array of local Rayleigh numbers.


calculate_salinity_flux
-------------------------

.. py:function:: calculate_salinity_flux(dz, dt, Ra_c, thickness_index, thickness, salinity, phi)

   Compute brine-driven salinity flux based on Rayleigh number exceeding critical value.

   :param float dz: Grid spacing in ice.
   :param float dt: Time step.
   :param float Ra_c: Critical Rayleigh number.
   :param int thickness_index: Interface index.
   :param np.ndarray thickness: Ice thickness array.
   :param np.ndarray salinity: Salinity profile.
   :param np.ndarray phi: Porosity profile.
   :returns: Salinity flux array.


calculate_permeability
------------------------

.. py:function:: calculate_permeability(phi_i)

   Compute brine permeability as a function of porosity.

   :param np.ndarray phi_i: Porosity array.
   :returns: Permeability array.

   Formula:

   .. math::

      \Pi = 10^{-17} (1000 \phi)^{3.1}


calculate_brine_flux
---------------------

.. py:function:: calculate_brine_flux(thickness_index, thickness, interface_depth, salinity, phi, dt, dz, Ra_c=10)

   Calculate brine flux in the ice column.

   :param int thickness_index: Interface layer index.
   :param np.ndarray thickness: Ice thickness array.
   :param float interface_depth: Depth at ice-ocean interface.
   :param np.ndarray salinity: Salinity profile.
   :param np.ndarray phi: Porosity profile.
   :param float dt: Time step.
   :param float dz: Grid spacing.
   :param float Ra_c: Critical Rayleigh number.
   :returns: Brine flux array.


calculate_in_out_brine_flux
-----------------------------

.. py:function:: calculate_in_out_brine_flux(thickness_index, flux, phi)

   Compute **incoming and outgoing brine fluxes** for each layer.

   :param int thickness_index: Interface layer index.
   :param np.ndarray flux: Brine flux array.
   :param np.ndarray phi: Porosity profile.
   :returns: Tuple (flux_in, flux_out)


calculate_salinity_source_term_from_brineflux
-----------------------------------------------

.. py:function:: calculate_salinity_source_term_from_brineflux(salinity, flux_in, flux_out, flux_channel)

   Compute the **salinity source term** from brine fluxes.

   :param np.ndarray salinity: Salinity profile.
   :param np.ndarray flux_in: Incoming brine flux.
   :param np.ndarray flux_out: Outgoing brine flux.
   :param np.ndarray flux_channel: Total brine flux channel.
   :returns: Salinity source term array.


get_salinity_source_term
--------------------------

.. py:function:: get_salinity_source_term(thickness_index, thickness, salinity, phi, dt, dz)

   High-level function to calculate **salinity source term** from brine flux for a given ice column.

   :param int thickness_index: Index at ice-ocean interface.
   :param np.ndarray thickness: Ice thickness array.
   :param np.ndarray salinity: Salinity profile.
   :param np.ndarray phi: Porosity profile.
   :param float dt: Time step.
   :param float dz: Grid spacing.
   :returns: Salinity source term array.