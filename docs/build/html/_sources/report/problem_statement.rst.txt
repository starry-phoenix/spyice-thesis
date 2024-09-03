Problem Statement
=================

.. figure:: images/SeaIceModelMush.png
    :align: center
    :width: 300
    
    caption

Modelling of Sea ice freezing involves ice-brine interface tracking which is modelling using the generalised Stefan problem. In 
this project, a 1D advective-diffusion model predominantly diffusive in nature is implemented using finite differences. The 
numerical equations models thermal and brine effects on the sea ice system with the following set of equations :eq:`eq_mushylayer`:

.. math::
    :label: eq_mushylayer

    (\rho c)_{eff} \frac{\partial T}{\partial t} &= \frac{\partial}{\partial z}\left(k_{eff} \frac{\partial T}{\partial z} \right) - \rho_{i}L \frac{\partial \phi}{\partial t} 

    \phi \frac{\partial S_{br}}{\partial t} &= \frac{\partial}{\partial z} \left( D_{eff} \frac{\partial S_{br}}{\partial z}\right) - \frac{\rho_{i}}{\rho_{br}} P_S S_{br} \frac{\partial \phi} {\partial t} 


Brine in a sea ice system propagates along a vertical column and the flux across boundaries is considered to be constant (CHECK!!) 
as shown in Fig(). The Dirichlet boundary conditions (essential boundary condition) is implemented at the top layer of the 
vertical column and the system is assumed to be at melting temperature of sea ice as :eq:`eq_liquidus_temp`:. 

.. math::
    :label: eq_liquidus_temp

    T_m(S) = T_m - \Gamma S

The linear numerical model is solved using the *Thomas tri-diagonal solver* **(REF)** which is for a positive definite diagonally dominant matrix **(REF)**. The 
conditional stability of the implicit numerical difference system is verified using the *Fourier stability criteria* **(REF)**. The model 
parameters used in this project are number of iterations, sea ice freezing duration, and initial salinity.

**Stefan Problem:** So :math:`$\lambda$` giving the minimum absolute result is determined and used to compute the location of the phase change. 

.. math::
    :label: eq_stefanProblem
    
    -\dfrac{2sl\mathrm{e}^{-l^2}}{\sqrt{{\pi}}\operatorname{erf}\left(l\right)}-\dfrac{2s\mathrm{e}^{-2l^2}}{{\pi}\operatorname{erf}^2\left(l\right)}-1
