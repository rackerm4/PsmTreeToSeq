import csv
import os


def parameters_prep(generated_sample_parameters, generated_protracted_speciation_process_parameters,
                      seqgen_parameters, args, headers, file_name):
    """Combines all parameters dicts to one dict and pass dict + headers to write_data_to_txt"""
    id_file = {'id': file_name.split('.')[0]}
    data_dict = {**id_file, **generated_sample_parameters, **generated_protracted_speciation_process_parameters,
                 **seqgen_parameters}
    write_data_to_txt(data_dict, args.output, headers)


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
