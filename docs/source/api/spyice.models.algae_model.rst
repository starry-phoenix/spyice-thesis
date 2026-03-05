src.spyice.models.algae_model
==========================================

This module contains functions for modeling **sea-ice algae biogeochemistry**, including:

- Salinity, temperature, and nutrient-dependent growth functions.
- Photosynthetically active radiation (PAR) modeling through ice.
- Chlorophyll-to-carbon ratio calculations.
- Carbon and nutrient uptake dynamics.
- Radiation absorption by ice algae.
- Bulk tracer and chlorophyll calculations.
- ODE-based updates for biologically active layers (BAL).

The module is designed for process-based modeling in ice-ocean systems and supports depth-resolved calculations.


Functions
=========

fs_salinity
-------------

.. py:function:: fs_salinity(s_br_array)

   Compute the salinity limitation function :math:`f_s` for algae growth.

   :param numpy.ndarray s_br_array: Salinity in g kg^-1.
   :returns: Salinity function :math:`f_s`.

   Formula:

   .. math::

      f_s = \exp\Big(-(2.16 - a - b)^2\Big), 
      \quad
      a = 8.3 \cdot 10^{-5} s^{2.11}, 
      \quad
      b = 0.55 \log s


ft_temperature
----------------

.. py:function:: ft_temperature(t_c_array)

   Compute the temperature limitation function :math:`f_t`.

   :param numpy.ndarray t_c_array: Temperature in Celsius.
   :returns: Temperature function :math:`f_t`.

   Formula:

   .. math::

      f_t = \exp(r_g (T - 273.15)), \quad r_g = 0.0633


ln_nutrient
-------------

.. py:function:: ln_nutrient(c_n_array, k)

   Compute nutrient limitation function :math:`l_n`.

   :param numpy.ndarray c_n_array: Nutrient concentration in mmol m^-3.
   :param float k: Half-saturation constant (mmol m^-3).
   :returns: Nutrient limitation function :math:`l_n`.

   Formula:

   .. math::

      l_n = \frac{c_n}{k + c_n}


photosynthetic_active_radiation
---------------------------------

.. py:function:: photosynthetic_active_radiation(z, kappa=1.5, i_0=0.17, albedo=0.58, F_s_w=5)

   Calculate Photosynthetically Active Radiation (PAR) at depth.

   :param float or np.ndarray z: Depth in ice (m).
   :param float kappa: Attenuation coefficient (m^-1).
   :param float i_0: Incident solar radiation factor.
   :param float albedo: Ice/snow albedo.
   :param float F_s_w: Incoming solar irradiance (W m^-2).
   :returns: PAR in W m^-2.

   Formula:

   .. math::

      I(z) = I_0 \exp(-\kappa z), \quad
      PAR = 4.91 \cdot I(z)


cholorophyl_to_C_ratio_par
----------------------------

.. py:function:: cholorophyl_to_C_ratio_par(PAR, ln, E=0.5, r_chl_c_max=0.05, r_chl_c_min=0.01)

   Chlorophyll-to-carbon ratio as a function of PAR and nutrients.

   :param np.ndarray PAR: Photosynthetically active radiation.
   :param np.ndarray ln: Nutrient limitation function.
   :param float E: Half-saturation energy.
   :param float r_chl_c_max: Maximum chlorophyll-to-carbon ratio.
   :param float r_chl_c_min: Minimum chlorophyll-to-carbon ratio.
   :returns: Chlorophyll-to-carbon ratio.


Photosynthesis_rate_maximum_light
-----------------------------------

.. py:function:: Photosynthesis_rate_maximum_light(mu_m, fs, ft, ln, r_chl_c_par)

   Compute maximum potential photosynthesis rate :math:`P_m`.

   :param float mu_m: Maximum photosynthesis rate (s^-1).
   :param np.ndarray fs: Salinity function.
   :param np.ndarray ft: Temperature function.
   :param np.ndarray ln: Nutrient function.
   :param np.ndarray r_chl_c_par: Chlorophyll-to-carbon ratio.
   :returns: Maximum photosynthesis rate.


Ek_light
----------

.. py:function:: Ek_light(Photosynthesis_rate, alpha)

   Light limitation parameter :math:`E_k`.

   :param np.ndarray Photosynthesis_rate: Maximum photosynthesis rate.
   :param float alpha: Photosynthetic efficiency.
   :returns: Light limitation parameter :math:`E_k = P_m / \alpha`


Ek_light_tanh
---------------

.. py:function:: Ek_light_tanh(PAR_arr, Ek_arr)

   Tanh-based light limitation factor.

   :param np.ndarray PAR_arr: Photosynthetically active radiation.
   :param np.ndarray Ek_arr: Light limitation parameter.
   :returns: Light limitation as tanh(PAR / Ek)


photosynthetic_rate
---------------------

.. py:function:: photosynthetic_rate(max_mu, fs, ft, ln, lpar)

   Compute actual photosynthetic rate considering light, temperature, salinity, and nutrient limitation.

   :param float max_mu: Maximum photosynthesis rate.
   :param np.ndarray fs: Salinity function.
   :param np.ndarray ft: Temperature function.
   :param np.ndarray ln: Nutrient function.
   :param np.ndarray lpar: Light limitation factor.
   :returns: Photosynthetic rate.


fs_ft / fs_ft_ln
------------------

.. py:function:: fs_ft(f_s, f_t)
.. py:function:: fs_ft_ln(f_s, f_t, l_n)

   Combine limitation functions multiplicatively:

   - :math:`f_s \cdot f_t`
   - :math:`f_s \cdot f_t \cdot l_n`


model_algae_processes
------------------------

.. py:function:: model_algae_processes(s_br_array, t_c_array, c_n_array, z_array, ...)

   High-level function to compute depth-resolved algae growth:

   - Salinity, temperature, nutrient limitation
   - PAR calculation
   - Chlorophyll-to-carbon ratio
   - Maximum photosynthesis
   - Light limitation
   - Actual photosynthetic rate

   :returns: Tuple (mu, PAR, nutrient_function)


ode_update_carbon_nutrient_uptake
-----------------------------------

.. py:function:: ode_update_carbon_nutrient_uptake(dt, salinity_list, temperature_list, nutrient_list, depth_list, cc_old, cn_old, r_n_c=0.12)

   Update carbon and nutrient concentrations for algae over a timestep `dt`.


chla_algae
------------

.. py:function:: chla_algae(PAR, nutrient_function, c_bulk_tracer)

   Compute chlorophyll-a bulk concentration from carbon tracer.


radiation_algae
-----------------

.. py:function:: radiation_algae(chla_bulk_z, I_array)

   Compute absorbed radiation by algae in ice.


get_bulk_tracer_concentration
-------------------------------

.. py:function:: get_bulk_tracer_concentration(liquid_fraction, brine_concentration)

   Compute bulk tracer concentration considering liquid fraction.


biogeochemical_model / biogeochemical_model_at_alldepths
---------------------------------------------------------

.. py:function:: biogeochemical_model(temperature, salinity, liquid_fraction, nutrient_concentration, carbon_concentration, dt, thickness_index, thickness)
.. py:function:: biogeochemical_model_at_alldepths(temperature, salinity, liquid_fraction, nutrient_concentration, carbon_concentration, dt, thickness, thickness_index)

   Full depth-resolved or interface-level biogeochemical updates for algae:

   - Photosynthesis
   - Carbon and nutrient uptake
   - Bulk tracer and chlorophyll
   - Radiation absorption

   :returns: Updated field arrays:

      - Carbon concentration
      - Nutrient concentration
      - Photosynthetic rate
      - Radiation
      - Chlorophyll-a bulk