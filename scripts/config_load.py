#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import yaml
import sys


class Profile:
    def __init__(self, profile):
        self.profile = profile

    def get_path(self):# 'scripts',  <= ADD FOR DOCKER
        cfg_path = os.path.join(os.getcwd(), 'profiles', self.profile + '.yaml')
        return cfg_path

    def cfg_load(self):
        path = self.get_path()
        if not os.path.isfile(path):
            print("Config file does not exist.")
            sys.exit(1)
        with open(path) as f:
            config = yaml.safe_load(f)
        return config

    def get_specific_config_values(self, getting):
        cfg = self.cfg_load()
        return cfg[str(getting)]

    def get_generate_sample_values(self):
        values_generate_sample = self.get_specific_config_values('generate_sample')
        return values_generate_sample

    def get_protractedspeciationprocess_values(self):
        psp = self.get_specific_config_values('ProtractedSpeciationProcess')
        vital = ['incipient_species_extinction_rate', 'speciation_initiation_from_orthospecies_rate', 'speciation_initiation_from_incipient_species_rate', 'speciation_completion_rate', 'orthospecies_extinction_rate', 'aincipient_species_extinction_rate']
        if 'ProtractedSpeciationProcess' in self.cfg_load():
            for arg in vital:
                if arg not in psp:
                    print("Missing ProtractedSpeciationProcess argument: " + arg + " in file " + self.get_path())
                    sys.exit(1)
            values_protractedspeciationprocess = psp
            return values_protractedspeciationprocess
        else:
            print("Missing ProtractedSpeciationProcess in config file.")
            sys.exit(1)

    def get_seq_gen_values(self):
        values_seq_gen = self.get_specific_config_values('seq-gen')
        return values_seq_gen
