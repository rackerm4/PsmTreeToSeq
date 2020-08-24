import os
import yaml


def cfg_load(config):
    path = os.path.join(os.getcwd(), 'src', config + '.yaml')
    if not os.path.isfile(path):
        print("Config file does not exist.")
        # self._exit_handler()
    with open(path) as f:
        config = yaml.safe_load(f)
    return config


def load_headers(config):
    return ['id'] + list(config['ProtractedSpeciationProcess']) + list(config['generate_sample']) + list(
        config['seq-gen'])
