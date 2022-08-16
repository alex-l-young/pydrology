# ==========================
# NEXRAD configuration.
# ==========================

# Library imports.
from pathlib import Path

# Local imports.
from pydrology.config import Config


class NexradConfig(Config):
    """
    Configuration class for NEXRAD processing.

    Parameters
    ---------------
    cfg_path
        - The path to the NEXRAD configuration file as a string.
    """

    def __init__(self, cfg_path:str):
        super().__init__(cfg_path)
        self._check_keys()


    @property
    def nexrad_directory(self):
        return self.cfg.get("nexrad_directory")

    @nexrad_directory.getter
    def nexrad_directory(self):
        return Path(self.cfg.get("nexrad_directory")) / str(self.cfg.get("date"))

    @property
    def output_files(self):
        return self.cfg.get("output_files")

    @output_files.getter
    def output_files(self):
        ld = self.nexrad_directory.glob('**/*')
        output_files = [self.nexrad_directory / f.name for f in ld if str(f)[-3:] == 'V06']
        return output_files

    @property
    def npool(self) -> int:
        return self.cfg.get("npool", 1)

    @property
    def date(self) -> str:
        return self.cfg.get("date")

    @property
    def site(self):
        return self.cfg.get("site")

    @property
    def variable(self):
        return self.cfg.get("variable")

    @property
    def grid_interp(self):
        return self.cfg.get("grid_interp")

    @property
    def bounds(self):
        return self.cfg.get("bounds")

    @property
    def num_points(self):
        return self.cfg.get("num_points")

    @property
    def spacing(self):
        return self.cfg.get("spacing")

    @property
    def sweep(self):
        return self.cfg.get("sweep")