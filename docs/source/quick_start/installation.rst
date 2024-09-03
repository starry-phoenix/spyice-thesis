Installation Guide
==================

Prerequisites
-------------

Before installing Spyice, make sure you have the following prerequisites installed:

- **Python 3.11 or above**: Spyice requires Python 3.11 or above to run. You can download the latest version of Python from the official website.

- **Hatch**: Hatch is a Python package manager that we'll use to create a virtual environment for Spyice. You can install Hatch by following the instructions on the `official website <https://hatch.pypa.io/latest/install/>`_.

- **Sphinx**: Sphinx is a documentation generator that we'll use to build the Spyice documentation. You can install Sphinx using the appropriate package manager for your system by following the instructions on the `official website <https://www.sphinx-doc.org/en/master/usage/installation.html>`_.

Installation Steps
------------------

Follow these steps to install Spyice and set up the environment using hatch:

1. Clone this project repository to your local machine.

2. (Optional) The required wheels can be built with the help of `hatch` without the worries of cross-compilation and native architecture support.
    
    ```
    hatch build
    ```

This command will build the wheels for the project and store them under the name `dist/spyice-1.0.0.dev0-py3-none-any.whl`

3. Create a the new default python virtual environment

    ```
    hatch env create
    ```
    
The project will be automatically installed in editable mode by `hatch` when the environment is created. Confirm the installation by running `pip show spyice` in the shell.

4. You're all set! You can now start using Spyice.

Viewing the documentation
--------------------------

Enter the project directory and run the following command in the terminal to view the pre-built static version of the sphinx documentation.

    ```
    python -m http.server --directory=docs/build/html/
    ```

This will start a local server at `http://localhost:8000/` where you can view the documentation by opening the link in your browser.