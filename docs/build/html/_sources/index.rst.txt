##############################
Sea-Ice Model Package (SPyIce)
##############################
   
The SPyIce package is a software tool that enables 1D finite difference simulation for vertical transport equations. It specifically focuses on thermal diffusion with the influence of salinity and physical properties. The package utilizes the Thomas tridiagonal solver as the solver algorithm. With SPyIce, users can model and analyze the behavior of temperature, salinity, and other relevant variables in a vertical system. It provides a comprehensive framework for studying the thermal diffusion process and its interaction with salinity in various scenarios. Hydra is used to automate the simulation runs of the Sea-Ice Model. It is used to manage and run sea ice simulations, making it easier for users to explore different scenarios and optimize their models.
You can also model microalgae bloom on sea ice and how nutrient transport, sunlight and radiation affects its growth with this package.


.. list-table::
   :align: center
   :widths: 33 33 33

   * - .. figure:: images/seaicemodel.png
         :width: 270px
         :align: center
         :height: 200px

         Sketch of sea ice model with biogeochemical processes.

     - .. figure:: images/seaicemodel_algae_layers.png
         :width: 200px
         :align: center
         :height: 200px

         Sketch of sea ice model with multiple microalgae layers near the ocean interface.

     - .. figure:: images/sunlight_brinedrainage_seaice.png
         :width: 200px
         :align: center
         :height: 200px

         Model with gravity drainage and sunlight penetration through sea ice.

#########
Contents
#########

.. toctree::
  :maxdepth: 1

  About <report/_index>
  Installation <quick_start/installation.rst>
  API <api/_index>
  User Guide <user_guide/_index>
  Examples <quick_start/example/_index>
