*********************************************************
Implementation with Hydra hierarchical configuration
*********************************************************

Hydra configuration tool is used to feed multiple input combinations to the Sea ice model and generate a structured output directory based on the given input combination. For example: all combinations of user defined input values for salinity, 
maximum iterations and time step size are sequentially run using hydra which improves the traceability of simulation along with 
logging files in their respective output directory. 

=============================================
Hydra configuration for main_config.py script
==============================================


Check out the `main_config.py` script for the implementation.

.. code-block:: python
    :caption: Python psuedo-example for main_config.py script


    @hydra.main(version_base=None, config_path="example/conf", config_name="config")
    def my_app(cfg: DictConfig) -> None:
        """Runs the model using the provided configuration and output directory."""
        with custom_stdout_logger(log):
            out_hydra_dir = hydra.core.hydra_config.HydraConfig.get().runtime.output_dir
            sys.stdout = SpyiceLogger(log)
            log.debug("Debug level message")
            log.warning("Warning level message")
            main = MainProcess(cfg, out_hydra_dir)
            main.run_model()

You can run for this main_config.py script using the following command in the terminal:
   
    .. code-block:: bash
        :caption: Terminal command
        
        python main_config.py 


===================================================
Hydra configuration for jupter notebook
===================================================

You can also use hydra configuration for jupyter notebook as well. Check out the `main_original.ipynb` notebook for the implementation.

.. code-block:: python
    :caption: Python psuedo-example for main_original.ipynb notebook
    :linenos:

    import os
    from pathlib import Path
    from hydra import (
        compose,
        initialize,
    )
    from omegaconf import OmegaConf

    # import the main process class
    from spyice.main_process import MainProcess

    with initialize(version_base=None, config_path="conf"):
        cfg = compose(
            config_name="config.yaml",
            overrides=["iter_max=iter_max1500", "dt=dt47", "S_IC=S34"],
        )
        out_hydra_dir = Path(output_base_dir, "with_hydra")
        main = MainProcess(cfg, out_hydra_dir)
        main.run_model()