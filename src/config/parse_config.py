
from typing import List, Dict
import configparser


class FaseInternalConfig:
    allowed_species: List[str]
    biomart_name_for_species: Dict[str, str]

    def __init__(self, internal_config_path):

        config = configparser.ConfigParser()
        config.read(internal_config_path)

        self.allowed_species = config["DEFAULT"]["AllowedSpecies"].split(" ")
        self.biomart_name_for_species = {k: v for k, v in [
            s.split(":") for s in config["DEFAULT"]["BiomartNameForSpecies"].split(" ")
        ]}


class FaseUserConfig:

    pass


class FaseRunConfig:

    pass
