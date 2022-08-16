# ====================================================
# Configuration base class.
# ====================================================

# Library imports.
import yaml
from pathlib import Path


class Config:
    def __init__(self, cfg_path:str):
        cfg_path = Path(cfg_path)
        self.cfg = Config._parse_cfg(cfg_path)


    @staticmethod
    def _parse_cfg(cfg_path:Path) -> dict:
        """
        Parses the configuration yaml file into a configuration dictionary.

        :param cfg_path: Path to configuration .yml file.
        :return: Configuration dictionary.
        """
        with open(cfg_path, "r") as ymlfile:
            cfg = yaml.safe_load(ymlfile)

        return cfg


    def _check_keys(self):
        """
        Checks the keys to ensure that all configurations passed are valid.

        :return: Raises ValueError if a non-valid configuration is found.
        """
        for key, value in self.cfg.items():
            if key not in dir(self):
                raise ValueError(f'"{key}" is not a valid configuration.')


    def as_dict(self) -> dict:
        """
        Converts the configuration yml file to a dictionary.

        :return: The configurations as a dictionary.
        """
        return self.cfg

