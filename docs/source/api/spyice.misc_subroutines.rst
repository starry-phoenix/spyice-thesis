spyice.misc_subroutines
=====================================

This module provides functions to **calculate the melting temperature of seawater based on salinity**.  
It supports different liquidus relations such as "Normal" and "Frezchem".

Functions
=========

calculate_melting_temperature_from_salinity
-------------------------------------------

.. py:function:: calculate_melting_temperature_from_salinity(_salinity, _temperature_melt=_temperature_melt, _liquid_relation="Normal")

   Calculates the melting temperature of seawater based on salinity.

   :param numpy.ndarray _salinity: Array of salinity values (in ppt).
   :param float _temperature_melt: Reference melting temperature (solidus) in Kelvin. Defaults to the UserInput value.
   :param str _liquid_relation: Type of liquidus relation. Must be either "Normal" or "Frezchem". Defaults to "Normal".
   :returns: numpy.ndarray of melting temperatures corresponding to the input salinity values.
   :raises TypeError: If ``_liquid_relation`` is not "Normal" or "Frezchem".

   **Description:**

   - **Normal relation:**  
     Uses a linear interpolation between the solidus temperature and the eutectic temperature.  
     Eutectic temperature `T_s = 252.05 K`, brine salinity `S_br = 233 ppt`.

     For scalar salinity:
     
     .. math::
        T_{melt} = T_m + (T_s - T_m) * \frac{S}{S_{br}}

     For array salinity, the same formula is applied element-wise.

   - **Frezchem relation:**  
     Uses a quadratic fit from the Frezchem database:

     .. math::
        T_{melt} = 272.63617665 - 9.1969758 \times 10^{-5} S^2 - 0.03942059 S