import contextlib
import logging
import sys
import os
import hydra
from omegaconf import DictConfig

from spyice.main_process import MainProcess
from spyice.utils import SpyiceLogger

"""
This script is used to configure and run a main process for SeaIce application.
It imports necessary modules and defines a main function `my_app` that takes a configuration
object as input and performs the following steps:
1. Sets the output directory for the Hydra configuration.
2. Initializes an instance of the `CMain` class with the given configuration and output directory.
3. Sets user input for the process object.
4. Runs the preprocess step of the process object.
5. Plots the outputs of the process object and saves them in the output directory.

The script also includes a conditional block that checks if the script is being run as the main
module, and if so, calls the `my_app` function.

To run the script, use the following command in the terminal:
python main_config.py -m consts=real dt=dt47 S_IC=S34

Note: This script requires the `omegaconf`, `hydra`, and `subprocess` modules to be installed.
You can generate a `requirements.txt` file by running `pip freeze > requirements.txt`.
"""


log = logging.getLogger(__name__)


@contextlib.contextmanager
def custom_stdout_logger(logger):
    original_stdout = sys.stdout
    try:
        sys.stdout = SpyiceLogger(logger)
        yield
    finally:
        sys.stdout = original_stdout


@hydra.main(version_base=None, config_path="example/conf", config_name="config")
def my_app(cfg: DictConfig) -> None:
    # Set the output directory for the Hydra configuration
    with custom_stdout_logger(log):
        out_hydra_dir = hydra.core.hydra_config.HydraConfig.get().runtime.output_dir
        sys.stdout = SpyiceLogger(log)
        log.debug("Debug level message")
        # print(f"cfg is :{out_hydra_dir}")
        # log.info("Info level message")
        log.warning("Warning level message")
        main = MainProcess(cfg, out_hydra_dir)
        main.run_model()
        # main.postprocess()


if __name__ == "__main__":
    my_app()
