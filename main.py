#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import time
import tempfile
import dendropy
import src.data_to_csv as dtc
import src.loader as cl
import src.db as data
from dendropy.model import protractedspeciation
from dendropy.interop import seqgen

"""
    author: Robin Ackermann
    version: -
"""


def main():
    parser = argparse.ArgumentParser(description='Generate sample tree data under the protracted speciation model')
    parser.add_argument('--output', '-o', required=False, help='Output dir')
    parser.add_argument('--schema', '-s', choices=['newick', 'nexus'], required=True,
                        help='Tree schema: Newick, Nexus')
    parser.add_argument('--profile', '-p', default="default", required=True, help='')
    parser.add_argument('--num_runs', '-n', default=1, type=int, required=False, help='')

    args = parser.parse_args()
    args.parser = parser

    config = cl.Loader(args.profile)
    db = data.DB(args.output)
    headers = config.load_headers()

    # start
    start = time.perf_counter()
    c = 0
    for _ in range(args.num_runs):
        # try:
        # getting trees
        get_trees = call_sample_tree(args, config)
        # generating Sequences & saving trees
        for names in file_output(get_trees, args, ["lineage", "orthospecies"]):
            print(names)
            call_seq_gen(args, names, config)

        # saving parameters
        parameters_to_txt(config, args, headers)
        # except:
        #     print("Run dismissed.")
        #     c += 1
        #     pass
    print('Runs dismissed:',c)
    print('\nProcess took %.2f seconds to complete.' % (time.perf_counter() - start))


def temp_file_name():
    """
    Generates random string
    :return : random string
    """
    temp_name = next(tempfile._get_candidate_names())
    return temp_name


def gen_sample_values(values):
    """
    Returns variables for psp_ini.generate_sample in call_sample_tree function. Joins args to dict and filters empty args
    :param :
    :return : args with parameters only
    """
    return {k: v for k, v in values.items() if v}


def call_sample_tree(args, config):
    """
    Calls ProtractedSpeciationProcess and generates sample trees.
    :param args, config:
    :return trees,
    generate_tree[0] lineage_tree (|Tree| instance) – A tree from the protracted speciation process, with all lineages
    (good species as well as incipient species).
    generate_tree[1] orthospecies_tree (|Tree| instance) – A tree from the protracted speciation process with only
    “good” species.:
    """
    try:
        # calling args
        values = gen_sample_values(config.get_generate_sample_values())
        # generate trees
        generated_trees = protractedspeciation.ProtractedSpeciationProcess(
            **config.generate_protracted_speciation_process_values()).generate_sample(**values)
        # generated_trees = protractedspeciation.ProtractedSpeciationProcess(
        #     **config.generate_protracted_speciation_process_values()).generate_sample(**values)
        return generated_trees
    except BaseException as e:
        print("Maximum number of runs to execute in the event of prematurely-terminated simulations due to all "
              "lineages going extinct. Once this number or re-runs is exceed, then TreeSimTotalExtinctionException "
              "is raised. Defaults to 1000. Set to None to never quit trying.\n" + str(e))


def file_output(trees, args, tree_names):
    """Stores output files."""
    output_dir = args.output
    # sanity check
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(output_dir):
        sys.exit()
    for i in range(len(trees)):
        try:
            file_name = tree_names[i] + '_' + temp_file_name() + "." + str(args.schema)
            tmp_path = os.path.join(output_dir, file_name)
            trees[i].write_to_path(tmp_path, suppress_edge_lengths=True,
                                   schema=args.schema)
            yield file_name
        except BaseException as e:
            return "Unexpected error while saving tree data:\n" + str(e)


def call_seq_gen(args, name, config):
    full_path = os.path.join(args.output, "seq_align_{}.txt".format(name.split('.')[0]))
    #try:
    path_to_tree = os.path.join(args.output, name)
    trees = dendropy.TreeList.get(path=path_to_tree, schema=args.schema)
    s = seqgen.SeqGen()

    # generate one alignment per tree
    # as substitution model is not specified, defaults to a JC model
    # will result in a DataSet object with one DnaCharacterMatrix per input tree
    # d0 = s.generate(trees)
    # print(len(d0.char_matrices))
    # print(d0.char_matrices[0].as_string("nexus"))
    # with open("seq_align.txt", "w") as f:
    #     f.write(d0.char_matrices[0].as_string("nexus"))

    # instruct Seq-Gen to scale branch lengths by factor of 0.1
    # note that this does not modify the input trees
    s.scale_branch_lens = 0.1

    # more complex model
    s.char_model = seqgen.SeqGen.GTR
    s.state_freqs = [0.4, 0.4, 0.1, 0.1]
    s.general_rates = [0.8, 0.4, 0.4, 0.2, 0.2, 0.1]
    d1 = s.generate(trees)
    #print(len(d1.char_matrices))
    print(d1.char_matrices[0].as_string(args.schema))
    with open(full_path, "w") as f:
        f.write(d1.char_matrices[0].as_string(args.schema))

    print("Seq-Gen file stored in {}".format(os.path.join(os.getcwd(), full_path)))
    # except BaseException as e:
    #     print("Unexpected error while using Seq-Gen:\n" + str(e))

#     #a = config.get_seq_gen_values()['state_freqs']
#     """Generates Sequences based on generated trees from ProtractedSpeciationProcess"""
#     file_path = os.path.join(args.output, "seq_align_{}.txt".format(name))
#     #try:
#     trees = dendropy.TreeList.get(path=tree, schema=args.schema)
#     s = seqgen.SeqGen()
#
#     # generate one alignment per tree
#     # as substitution model is not specified, defaults to a JC model
#     # will result in a DataSet object with one DnaCharacterMatrix per input tree
#     # d0 = s.generate(trees)
#     # print(len(d0.char_matrices))
#     # print(d0.char_matrices[0].as_string("nexus"))
#     # with open("seq_align.txt", "w") as f:
#     #     f.write(d0.char_matrices[0].as_string("nexus"))
#
#     # instruct Seq-Gen to scale branch lengths by factor of 0.1
#     # note that this does not modify the input trees
#     s.scale_branch_lens = 0.1
#
# # more complex model
#     s.char_model = seqgen.SeqGen.GTR
#     s.state_freqs = [0.4, 0.4, 0.1, 0.1]
#     s.general_rates = [0.8, 0.4, 0.4, 0.2, 0.2, 0.1]
#     d1 = s.generate(trees)
#     # print(len(d1.char_matrices))
#     # print(d1.char_matrices[0].as_string(args.schema))
#     with open(file_path, "w") as f:
#         f.write(d1.char_matrices[0].as_string(args.schema))
#

def parameters_to_txt(config, args, headers):
    z = {**config.get_generate_sample_values(), **config.generate_protracted_speciation_process_values(),
         **config.get_seq_gen_values()}
    dtc.write_data_to_txt(z, args.output, headers)


if __name__ == '__main__':
    main()
