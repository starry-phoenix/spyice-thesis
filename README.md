# Sea-Ice Model Package (SPyIce)

The SPyIce package is a software tool that enables 1D finite difference simulation for vertical transport equations. It specifically focuses on thermal diffusion with the influence of salinity and physical properties. The package utilizes the Thomas tridiagonal solver as the solver algorithm. With SPyIce, users can model and analyze the behavior of temperature, salinity, and other relevant variables in a vertical system. It provides a comprehensive framework for studying the thermal diffusion process and its interaction with salinity in various scenarios. Hydra is used to automate the simulation runs of the Sea-Ice Model. It is used to manage and run sea ice simulations, making it easier for users to explore different scenarios and optimize their models.
Here is the link to the full documentation of SPyIce package: [SPyIce documentation](https://starry-phoenix.github.io/spyice-thesis/build/html/index.html). 

## Installation Guide

### Prerequisites

Before installing Spyice, make sure you have the following prerequisites installed:

- **Python 3.11 or above**: Spyice requires Python 3.11 or above to run. You can download the latest version of Python from the official website.

- **Hatch**: Hatch is a Python package manager that we'll use to create a virtual environment for Spyice. You can install Hatch by following the instructions on the [official website](https://hatch.pypa.io/latest/install/).

- **Sphinx**: Sphinx is a documentation generator that we'll use to build the Spyice documentation. You can install Sphinx using the appropriate package manager for your system by following the instructions on the [official website](https://www.sphinx-doc.org/en/master/usage/installation.html).

### Installation Steps

Follow these steps to install Spyice and set up the environment using hatch:

1. Clone this project repository to your local machine.

2. (Optional) The required wheels can be built with the help of `hatch` without the worries of cross-compilation and native architecture support.

    ```bash
    hatch build
    ```

    This command will build the wheels for the project and store them under the name `dist/spyice-1.0.0.dev0-py3-none-any.whl`

3. Create a the new default python virtual environment

    ```bash
    hatch env create
    hatch shell
    ```

    The project will be automatically installed in editable mode by `hatch` when the environment is created. Confirm the installation by running `pip show spyice` in the shell.

4. You're all set! You can now start using Spyice.

5. You can get started with SPyIce by running the following command in the terminal:

    ```bash
    python main_config.py
    ```

    You can use the `--help` flag to see the available options and arguments.

### To run Jupyter Notebooks


Enter into the spyice virtual environment using hatch (mentioned in the previous steps) and enter the following commands on the terminal. You can now find the spyice package renamed as spyice-jupyter as a kernel option. 

   ```console
   $pip install ipykernel
   $python -m ipykernel install --user --name spyice --display-name "Python (spyice-jupyter)"
   ```