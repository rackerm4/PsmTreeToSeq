import os
import sys
import numpy as np
import yaml


class Loader:
    def __init__(self, config, main=None):
        self.config = config
        self.main = main
    #     self.check_version()
    #
    # def check_version(self):
    #     if not dendropy.__version__.startswith('4'):
    #         sys.stderr.write('''ERROR: Your DendroPy version is not up to date, use: 'pip install dendropy --upgrade'
    #         ''')
    #         self._exit_handler()
    #
    # def _exit_handler(self):
    #     if self.main:
    #         sys.exit()
    #     else:
    #         pass

    def get_path(self):
        return os.path.join(os.getcwd(), 'src', self.config + '.yaml')

    def cfg_load(self):
        path = self.get_path()
        if not os.path.isfile(path):
            print("Config file does not exist.")
            # self._exit_handler()
        with open(path) as f:
            config = yaml.safe_load(f)
        return config

    def get_specific_config_values(self, getting):
        return self.cfg_load()[str(getting)]

    @staticmethod
    def generate_protracted_speciation_process_values():
        vital = ['incipient_species_extinction_rate', 'speciation_initiation_from_orthospecies_rate',
                 'speciation_initiation_from_incipient_species_rate', 'speciation_completion_rate',
                 'orthospecies_extinction_rate', 'aincipient_species_extinction_rate']
        rng_values = []
        n = 2  # number of digits after the decimal point
        for i in range(len(vital)):
            rng_values.append(round(np.random.uniform(0.001, 0.25), n))
        # test values below
        #return {'incipient_species_extinction_rate': 0.2, 'speciation_initiation_from_orthospecies_rate': 0.2, 'speciation_initiation_from_incipient_species_rate': 0.2, 'speciation_completion_rate': 0.2, 'orthospecies_extinction_rate': 0.2, 'aincipient_species_extinction_rate': 0.2}
        return dict(zip(vital, rng_values))

    # old function, reads in parameters from config file 'default.yaml'
    # def generate_protracted_speciation_process_values(self):
    #     psp = self.get_specific_config_values('ProtractedSpeciationProcess')
    #     vital = ['incipient_species_extinction_rate', 'speciation_initiation_from_orthospecies_rate',
    #              'speciation_initiation_from_incipient_species_rate', 'speciation_completion_rate',
    #              'orthospecies_extinction_rate', 'aincipient_species_extinction_rate']
    #     if 'ProtractedSpeciationProcess' in self.cfg_load():
    #         for arg in vital:
    #             if arg not in psp:
    #                 print("Missing ProtractedSpeciationProcess argument: " + arg + " in file " + self.get_path())
    #                 sys.exit()
    #         values_protractedspeciationprocess = psp
    #         return values_protractedspeciationprocess
    #     else:
    #         print("Missing ProtractedSpeciationProcess in config file.")
    #         sys.exit()

    def get_generate_sample_values(self):
        return self.get_specific_config_values('generate_sample')

    def generate_seq_gen_general_rates(self):
        n = 2 # number of digits after the decimal point
        return [round(np.random.uniform(0.001, 1), n) for _ in range(6)]

    def generate_seq_gen_state_freqs(self):
        """Generates 4 random numbers using dirichlet distribution. Return when sum = 1"""
        sum = 0
        vals = []
        while True:
            for i in np.random.dirichlet(np.ones(4)) * 1:
                rounded = np.round(i, 3)
                if rounded == 0:
                    sum = 0
                    vals = []
                    break
                vals.append(rounded)
                sum = np.sum(vals)
            if sum == 1:
                break
            sum = 0
            vals = []
        # vals = [int(i) for i in vals]
        return vals

    def get_seq_gen_values(self):
        """Grabs Seq-Gen parameters of config and generates random values."""
        config = self.get_specific_config_values('seq-gen')
        random_args = {k: v for k, v in config.items() if v == 1}
        for k in random_args.keys():
            if k == 'state_freqs':
                random_args['state_freqs'] = self.generate_seq_gen_state_freqs()
            elif k == 'general_rates':
                random_args['general_rates'] = self.generate_seq_gen_general_rates()
            else:
                random_args[k] = round(np.random.uniform(0.001, 1), 2)
        return {**config, **random_args}

    def load_headers(self):
        return list(self.cfg_load()['ProtractedSpeciationProcess']) + list(self.cfg_load()['generate_sample']) + list(
            self.cfg_load()['seq-gen'])