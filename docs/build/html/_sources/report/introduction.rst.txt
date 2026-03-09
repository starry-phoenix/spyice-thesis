************
Introduction
************

Multiphysics simulation has been growing in the past decades due to computational capabilities. It involves translating the 
physical partial differential equations to an efficient numerical code. Implementing good Scientific software development 
practices help in reducing computational power. Testing applicability of multiphysical computational models can prove to be 
challenging due to thermodynamic and numerical constraints. Therefore, improving the code reproducibility through 
automation and object oriented concepts is beneficial. 

.. figure:: images/SoftwareCycle2p0.png
    :align: center
    :height: 300
    :width: 300
    
    Caption

This project aims to improve an existing 1D advection-diffusion sea ice 
simulation model using python to improve the benchmark testing for various model parameters. The automation of the model 
is performed using hydra python package. Hatch python project manager is used for testing, static analysis checks, and for 
creating reproducible build ecosystem. OOP concepts are leveraged and the mediator design pattern is implemented 
to improve sustainability of code. 

