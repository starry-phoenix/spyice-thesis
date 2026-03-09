spyice.utils.initial_userinput
===================================

This module provides functions to **calculate the initial melt temperature of seawater or ice based on boundary salinity**.  
It supports multiple calculation methods: one-phase, Frezchem, and two-phase models.

Functions
=========

calculate_initial_melt_temperature_onephase
-------------------------------------------

.. py:function:: calculate_initial_melt_temperature_onephase(boundary_salinity)

   Calculates the initial melt temperature using a **one-phase model**.

   :param float boundary_salinity: Salinity at the boundary in ppt.
   :returns: Initial melt temperature in Kelvin.

   **Formula:**

   .. math::
      T_{melt} = 273.15 - 1.853 \frac{S}{28.0}

calculate_initial_temperature_frezchem
--------------------------------------

.. py:function:: calculate_initial_temperature_frezchem(boundary_salinity)

   Calculates the initial melt temperature using the **Frezchem model**.

   :param float boundary_salinity: Salinity at the boundary in ppt.
   :returns: Initial melt temperature in Kelvin.

   **Formula:**

   .. math::
      T_{melt} = 272.63617665 - 9.1969758 \times 10^{-5} S^2 - 0.03942059 S

calculate_initial_melt_temperature_twophase
-------------------------------------------

.. py:function:: calculate_initial_melt_temperature_twophase(boundary_salinity)

   Calculates the initial melt temperature using a **two-phase model**.

   :param float boundary_salinity: Salinity at the boundary in ppt.
   :returns: Initial melt temperature in Kelvin.

   **Parameters used:**

   - Solidus temperature: ``T_m = 273.15 K``
   - Eutectic temperature: ``T_s = 252.05 K``
   - Brine salinity: ``S_br = 233 ppt``

   **Formula:**

   .. math::
      T_{melt} = 273.15 + \frac{(T_s - T_m) \cdot S}{S_{br}}

calculate_initial_melt_temperature
----------------------------------

.. py:function:: calculate_initial_melt_temperature(boundary_salinity, method)

   Calculates the initial melt temperature using the specified method.

   :param float boundary_salinity: Salinity at the boundary in ppt.
   :param str method: Calculation method. Options are:

      - ``onephase``: One phase model
      - ``Frezchem``: Frezchem model
      - ``twophase``: Two-phase model

   :returns: Initial melt temperature in Kelvin.
   :raises ValueError: If an invalid method is specified.

   **Description:**  
   This function wraps the three specific calculation methods above and selects the appropriate formula based on the ``method`` argument.