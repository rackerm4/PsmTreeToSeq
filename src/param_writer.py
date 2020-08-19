import os
import csv


def parameters_to_txt(config, args, headers):
    """Parameter to file handler."""
    z = {**config.get_generate_sample_values(), **config.generate_protracted_speciation_process_values(),
         **config.get_seq_gen_values()}
    write_data_to_txt(z, args.output, headers)
    # todo:
    # call params first, then write
    # when n runs 4 => written 4, but 8 generated !!


def write_data_to_txt(data_dict, output, headers):
    """Parameter to file writer."""
    csv_columns = headers
    csv_file = 'used_parameters.txt'
    path = os.path.join(output, csv_file)
    if os.path.isfile(path):
        with open(path, 'a') as csvfile:
            writer = csv.DictWriter(csvfile, lineterminator='\n', delimiter=' ', fieldnames=csv_columns)
            writer.writerow(data_dict)
    else:
        with open(path, 'a') as csvfile:
            writer = csv.DictWriter(csvfile, lineterminator='\n', delimiter=' ', fieldnames=csv_columns)
            writer.writeheader()
            writer.writerow(data_dict)

    # csv_columns = ['incipient_species_extinction_rate', 'speciation_initiation_from_orthospecies_rate',
    # 'speciation_initiation_from_incipient_species_rate', 'speciation_completion_rate',
    # 'orthospecies_extinction_rate', 'aincipient_species_extinction_rate'], ['max_time', 'num_extant_orthospecies',
    # 'num_extant_lineages', 'is_retry_on_total_extinction', 'max_retries'], ['state_freqs', 'general_rates',
    # 'MODEL', 'SEQUENCE_LENGTH', 'NUMBER_OF_DATASETS', 'NUMBER_OF_PARTITIONS', '-sSCALE', '-dSCALE',
    # 'CODON_POSITION_RATES', 'ALPHA', 'NUM_CATEGORIES', 'PROPORTION_INVARIABLE', 'STATE_FREQUENCIES',
    # 'TRANSITION_TRANSVERSION_RATIO', 'RATE_MATRIX_VALUES', 'ANCESTRAL_SEQUENCE_NUMBER', 'RANDOM_NUMBER_SEED', 'op',
    # 'or', True, 'of', 'TEXT_FILE_NAME', 'wa', 'wr', 'q', 'h']
