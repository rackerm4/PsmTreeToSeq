import numpy as np


def generate_protracted_speciation_process_values():
    vital = ['incipient_species_extinction_rate', 'speciation_initiation_from_orthospecies_rate',
             'speciation_initiation_from_incipient_species_rate', 'speciation_completion_rate',
             'orthospecies_extinction_rate', 'aincipient_species_extinction_rate']
    rng_values = []
    n = 2  # number of digits after the decimal point
    for i in range(len(vital)):
        rng_values.append(round(np.random.uniform(0.001, 0.3), n))
    # test values below
    #return {'incipient_species_extinction_rate': 0.2, 'speciation_initiation_from_orthospecies_rate': 0.2, 'speciation_initiation_from_incipient_species_rate': 0.2, 'speciation_completion_rate': 0.2, 'orthospecies_extinction_rate': 0.2, 'aincipient_species_extinction_rate': 0.2}
    return dict(zip(vital, rng_values))


def generate_seq_gen_general_rates():
    n = 2 # number of digits after the decimal point
    return [round(np.random.uniform(0.001, 1), n) for _ in range(6)]


def generate_seq_gen_state_freqs():
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
    return vals


def get_seq_gen_values(config):
    """Grabs Seq-Gen parameters of config and generates random values."""
    random_args = {k: v for k, v in config.items() if v == 1}
    for k in random_args.keys():
        if k == 'state_freqs':
            random_args['state_freqs'] = generate_seq_gen_state_freqs()
        elif k == 'general_rates':
            random_args['general_rates'] = generate_seq_gen_general_rates()
        else:
            random_args[k] = round(np.random.uniform(0.001, 1), 2)
    return {**config, **random_args}


def gen_sample_values(values):
    """
    Returns variables for psp_ini.generate_sample in call_sample_tree function. Joins args to dict and filters empty args
    :param :
    :return : args with parameters only
    """
    return {k: v for k, v in values.items() if v}