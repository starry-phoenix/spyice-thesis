##############################
Sea-Ice Model Package (SPyIce)
##############################

The SPyIce package is a software tool that enables 1D finite difference simulation for vertical transport equations. It specifically focuses on thermal diffusion with the influence of salinity and physical properties. The package utilizes the Thomas tridiagonal solver as the solver algorithm. With SPyIce, users can model and analyze the behavior of temperature, salinity, and other relevant variables in a vertical system. It provides a comprehensive framework for studying the thermal diffusion process and its interaction with salinity in various scenarios. Hydra is used to automate the simulation runs of the Sea-Ice Model. It is used to manage and run sea ice simulations, making it easier for users to explore different scenarios and optimize their models.

Installation Guide
==================

Prerequisites
-------------

Before installing Spyice, make sure you have the following prerequisites installed:

- **Python 3.11 or above**: Spyice requires Python 3.11 or above to run. You can download the latest version of Python from the official website.

- **Hatch**: Hatch is a Python package manager that we'll use to create a virtual environment for Spyice. You can install Hatch by following the instructions on the `official website <https://hatch.pypa.io/latest/install/>`.

- **Sphinx**: Sphinx is a documentation generator that we'll use to build the Spyice documentation. You can install Sphinx using the appropriate package manager for your system by following the instructions on the `official website <https://www.sphinx-doc.org/en/master/usage/installation.html>`.

Installation Steps
------------------

Follow these steps to install Spyice and set up the environment using hatch:

1. Clone this project repository to your local machine.

2. (Optional) The required wheels can be built with the help of `hatch` without the worries of cross-compilation and native architecture support.

    .. code-block:: bash

        hatch build

    This command will build the wheels for the project and store them under the name `dist/spyice-1.0.0.dev0-py3-none-any.whl`

3. Create a the new default python virtual environment

    .. code-block:: bash

        hatch env create
        hatch shell

    The project will be automatically installed in editable mode by `hatch` when the environment is created. Confirm the installation by running `pip show spyice` in the shell.

4. You're all set! You can now start using Spyice.
5. You can get started with SPyIce by running the following command in the terminal:
   
       .. code-block:: bash

        python main_config.py 

    You can use the `--help` flag to see the available options and arguments.

Setting up a new parameter to automate
--------------------------

>!Note: The new parameter must belong to the already existing parameter directory which can be found in the `src/spyice/parameters/user_input.py`

    .. code-block:: bash
    src/spyice
    ├── ...
    ├── ...
    ├── main_process.py
    └── utils
        ├── config_sort.py
        ├── ...
        └── ...
    ├── preprocess/
        ├── pre_process.py
        ├── ...
        └── ...
    ├── parameters/
        ├── user_input.py
        ├── ...
        └── ...


1. You can setup configuration for a new parameter 
    .. code-block:: bash
    example/
    └── conf/
        └── parameter/
            └── param1.yaml
            └── param2.yaml   
    
2. Inside the  .yaml file you can set the parameters for values as a float/integer variable and arguments as strings, for example:
    .. code-block:: yaml

        # param1.yaml
        parameter:  1.0 # parameter value
    
    .. code-block:: yaml

        # param2.yaml
        parameter:  'Argument' # parameter value
3.  Make the following changes inside config_sort.py:
    add the new parameter in get_config_params() function and follow a similar format to get_ownconfig_params():
    .. code-block:: python
        ... 
        def get_config_params():
            self.new_parameter = config.get("parameter_directory", {}).get("parameter_from_yaml", "default")
            return config_params
    
4. Inside pre_process.py, make the following changes:
    .. code-block:: python
    class PreprocessData:
    ...
        def __init__():
            ...
            ...
        if config_data.constants_type == "real":
            UserInput.__init__(
                self,
                new_parameter=config_data.new_parameter,
                ...
            )

Viewing the documentation
--------------------------

Enter the project directory and run the following command in the terminal to view the pre-built static version of the sphinx documentation.

.. code-block:: bash

    python -m http.server --directory=docs/build/html/

This will start a local server at `http://localhost:8000/` where you can view the documentation.

1. (Optional) To test and run coverage and use hydra on command line, run the following command in the terminal:

## Running Tests

```console
hatch test tests\test_userinput.py
```

Checks for code coverage Coverage.py

```console
hatch test --cover
```

## Static Analysis

Hatch uses Ruff for linting

```console
hatch fmt
```

## Type Hinting

Hatch uses mypy for type hinting

```console
hatch run types:check
```