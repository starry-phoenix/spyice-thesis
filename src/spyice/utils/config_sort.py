from omegaconf import DictConfig
from dataclasses import dataclass
from .helpers import set_dataclass


@dataclass
class ConfigData:
    """Class representing configuration data.

    Args:
        setup (bool): Indicates if the setup is enabled or not.
    """

    setup = True


def read_omegaconfig(config, parameter: str):
    """Reads the parameter from the OmegaConf configuration.

    Args:
        parameter (str): The parameter to read.
    Returns:
        str: The value of the parameter.
    Raises:
        None
    """

    return config.get(parameter, {}).get(parameter, "default")


class ConfigSort:
    def __init__(self, config: DictConfig):
        """A class that provides methods for retrieving configuration parameters.

        Args:
            config (DictConfig): The configuration dictionary.
        Attributes:
            config (DictConfig): The configuration dictionary.
        Methods:
            get_config_params(config: DictConfig): Retrieves configuration parameters using the 'consts', 'dt', 'S_IC', and 'iter_max' keys.
            get_ownconfig_params(config): Retrieves configuration parameters using the 'constants', 'dt', 'S_IC', and 'iter_max' keys.
            getconfig_dataclass(config, config_type="default"): Retrieves configuration parameters based on the specified config_type.
        """
        self.config = config

    def get_config_params(self, config: DictConfig):
        """Get the configuration parameters from the given `config` dictionary.

        Args:
            config (DictConfig): The configuration dictionary.
        Returns:
            None
        Raises:
            None
        """

        self.constants_type = read_omegaconfig("constants")
        # self.time_step = self.read_omegaconfig("dt")
        # self.initial_salinity = self.read_omegaconfig("S_IC")
        # self.max_iterations = self.read_omegaconfig("iter_max")
        # self.is_salinity_equation = self.read_omegaconfig("salinity")

    def get_ownconfig_params(self, config):
        """Retrieves the parameters from the given configuration dictionary.

        Args:
            config (dict): The configuration dictionary.
        Returns:
            None
        """

        self.constants_type = config.get("constants")
        self.time_step = config.get("dt")
        self.initial_salinity = config.get("S_IC")
        self.max_iterations = config.get("iter_max")

    @classmethod
    def getconfig_dataclass(
        cls, config: dataclass, config_type: str = "default"
    ) -> dataclass:
        """Retrieves configuration parameters based on the specified config_type.

        Args:
            config: The configuration dictionary.
            config_type (str): The type of configuration "default" or "jupyter".
                               "jupyter" is used for Jupyter notebook configurations.
        Returns:
            ConfigData: An instance of the ConfigData class.
        """

        config_method = cls(config)
        if config_type == "default":
            config_method.get_config_params(config)
        elif config_type == "jupyter":
            config_method.get_ownconfig_params(config)
        config_vars = dict(vars(config_method))
        return set_dataclass(config_vars, ConfigData)
