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

        {\rho c}_{\text{eff}} \frac{\partial T}{\partial t} = -\rho_{l} c_l \textbf{w} \cdot   \nabla T + \nabla \cdot \left( \kappa_{\text{eff}} \nabla T \right) - \rho_{s} L \frac{\partial \phi}{\partial t} + Q

        \phi \frac{\partial S}{\partial t} = - \textbf{w} \cdot \nabla S + \nabla \cdot \left( D_{\text{eff}} \nabla S \right) - \frac{\rho_s}{\rho_l} P_S S \frac{\partial \phi}{\partial t} + Q'

        \nabla \cdot \textbf{w} = \left( \frac{\rho_s}{\rho_l} - 1 \right) \frac{\partial \phi}{\partial t}

        \mu \textbf{w} = \Pi (-\nabla p + \rho g)


Brine in a sea ice system propagates along a vertical column and the flux across boundaries is considered to be constant (CHECK!!) 
as shown in Fig(). The Dirichlet boundary conditions (essential boundary condition) is implemented at the top layer of the 
vertical column and the system is assumed to be at melting temperature of sea ice as :eq:`eq_liquidus_temp`:. 

.. math::
    :label: eq_liquidus_temp

    T_m(S) = T_m - \Gamma S

The linear numerical model is solved using the *Thomas tri-diagonal solver* **(REF)** which is for a positive definite diagonally dominant matrix **(REF)**. The 
conditional stability of the implicit numerical difference system is verified using the *Fourier stability criteria* **(REF)**. The model 
parameters used in this project are number of iterations, sea ice freezing duration, and initial salinity.

The instabilities due to temperature and salinity gradients causing expulsion of salt in sea ice causing gravity drainage can be described using a parametrised equation as a function of Rayleigh number.

.. math::
    :label: Gravity drainage

     Ra = \frac{g(h_i - z) \rho_w \beta_w (S(z) - S_w) \Pi(\phi_{min})}{\kappa_s \eta}

The transport of nutrients follows salt transport; therefore, the salt mass conservation equation can be used:

.. math::
    :label: Nutrient transport

    \frac{\partial C}{\partial t} = -w \nabla C + \nabla \cdot \left( D_{\text{eff}} \nabla C \right).

The micro-algal growth dynamics follows the Nicholson-Bailey model, which is a Malthusian model where the growth and depletion of algae and nutrients follows the law of exponential growth. Algal activity is associated with its byproduct, carbon concentration and the presence of algae with chlorophyll concentration,

.. math::
    :label: Nutrient transport

    \frac{d C_C}{d t} &= (\mu - \lambda)C_C,

    \frac{d C_N}{d t} &= r^N_C(-\mu + f\lambda)C_C,


**Stefan Problem:** So :math:`$\lambda$` giving the minimum absolute result is determined and used to compute the location of the phase change. 

.. math::
    :label: eq_stefanProblem
    
    -\dfrac{2sl\mathrm{e}^{-l^2}}{\sqrt{{\pi}}\operatorname{erf}\left(l\right)}-\dfrac{2s\mathrm{e}^{-2l^2}}{{\pi}\operatorname{erf}^2\left(l\right)}-1
