#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import time
import dendropy
import tempfile

# import _thread as thread
from dendropy.model import protractedspeciation
from dendropy.interop import seqgen

"""
    author: Robin Ackermann
    version: -
    Script to generate sample tree & sequence data with protracted speciation model
"""
#
# def main():
#     start = time.perf_counter()
#     # thread.start_new_thread(waiting_std_output, (20, 0.5))
#
#     print("Creating trees and sequences.")
#     # getting trees
#     call_sample_tree(args)
#     print("\nTree files finished.")
#
#     # saving trees
#     # lineage_path = file_output(get_trees[0], args, "lineage")
#     # orthospecies_path = file_output(get_trees[1], args, "orthospecies")
#
#     # # generating Sequences
#     # print("\nStarting Seq-Gen")
#     # call_seq_gen(lineage_path, args, "lineage")
#     # call_seq_gen(orthospecies_path, args, "orthospecies")
#     #
#     # print("Finished.")
#     # print('\nProcess took %.2f seconds to complete.' % (time.perf_counter() - start))
#


def temp_file_name():
    temp_name = next(tempfile._get_candidate_names())
    return temp_name


def gen_sample_values(config):
    """
    Returns variable for psp_ini.generate_sample in call_sample_tree function. Joins args to dict and filters empty args
    :param :
    :return : args with values only
    """

    k = ["max_retries", "num_extant_lineages", "is_retry_on_total_extinction", "num_extant_orthospecies", "max_time"]
    v = [args.max_retries, args.num_extant_lineages, args.is_retry_on_total_extinction, args.num_extant_orthospecies,
         args.max_time]
    d = dict(zip(k, v))

    if 'generate_sample' in config:
        return {k: v for k, v in d.items() if v}


def call_sample_tree(args, parameters):
    """
    Calls ProtractedSpeciationProcess and generates sample trees.
    :param args:
    :return trees,
    generate_tree[0] lineage_tree (|Tree| instance) – A tree from the protracted speciation process, with all lineages
    (good species as well as incipient species).
    generate_tree[1] orthospecies_tree (|Tree| instance) – A tree from the protracted speciation process with only
    “good” species.:
    """

    psp_ini = protractedspeciation.ProtractedSpeciationProcess(**parameters)
    try:
        # calling for args with values
        values = gen_sample_values(args)
        # generate trees
        generated_trees = psp_ini.generate_sample(**values)
        return generated_trees
    except BaseException as e:
        print("Maximum number of runs to execute in the event of prematurely-terminated simulations due to all "
              "lineages going extinct. Once this number or re-runs is exceed, then TreeSimTotalExtinctionException "
              "is raised. Defaults to 1000. Set to None to never quit trying.\n" + str(e))


def file_output(trees, args, tree_names):
    """Stores output files."""

    output_dir = args.output + 'trees'
    # sanity check
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(output_dir):
        exit()
    for i in range(len(trees)):
        try:
            filename = temp_file_name() + tree_names[i] + "." + str(args.schema)
            tmp_path = os.path.join(os.getcwd(), output_dir, filename)
            trees[i].write_to_path(tmp_path, suppress_edge_lengths=True,
                               schema=args.schema)
        except BaseException as e:
            return "Unexpected error while saving tree data:\n" + str(e)


def call_seq_gen(tree, args, name):
    """Generates Sequences based on generated trees from ProtractedSpeciationProcess"""

    filename = "seq_align_{}.txt".format(name)
    output_dir = os.path.join(args.output, 'seqs')
    full_path = os.path.join(os.getcwd(), output_dir, filename)
    # sanity check
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(output_dir):
        exit()
    try:
        trees = dendropy.TreeList.get(path=tree, schema=args.schema)
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
        #print(d1.char_matrices[0].as_string(args.schema))
        with open(full_path, "w") as f:
            f.write(d1.char_matrices[0].as_string(args.schema))
        print("Seq-Gen file stored in {}".format(os.path.join(os.getcwd(), output_dir)))
    except BaseException as e:
        print("Unexpected error while using Seq-Gen:\n" + str(e))

# waiting animation
# def waiting_std_output(lenstr=1, zzz=0.5, dispstr='Creating trees and sequences'):
#     dots = '.' * lenstr
#     spaces = ' ' * lenstr
#     print(dispstr)
#     while True:
#         for i in range(lenstr):
#             time.sleep(zzz)
#             outstr = dots[:i] + spaces[i:]
#             sys.stdout.write('\b' * lenstr + outstr)
#             sys.stdout.flush()


#if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description='Generate sample tree data under the protracted speciation model')
    #
    # parser.add_argument('--output', '-o', required=False, help='Output dir')
    # parser.add_argument('--schema', '-s', choices=['newick', 'nexus'], required=True, help='Tree schema: Newick, Nexus')
    # # parser.add_argument('--override', '-r', action='store_true', required=False, help='Overrides existing files.')
    #
    # parser.add_argument('--incipient_species_extinction_rate', type=float, required=True,
    #                     help='')
    # parser.add_argument('--speciation_initiation_from_orthospecies_rate', type=float, required=True,
    #                     help='')
    # parser.add_argument('--speciation_initiation_from_incipient_species_rate', type=float, required=True,
    #                     help='')
    # parser.add_argument('--speciation_completion_rate', type=float, required=True,
    #                     help='')
    # parser.add_argument('--orthospecies_extinction_rate', type=float, required=True,
    #                     help='')
    # parser.add_argument('--aincipient_species_extinction_rate', type=float, required=True,
    #                     help='')
    # parser.add_argument('--max_time', type=float, required=False,
    #                     help='')
    # parser.add_argument('--num_extant_orthospecies', type=int, required=False,
    #                     help='')
    # parser.add_argument('--num_extant_lineages', type=int, required=False,
    #                     help='')
    # parser.add_argument('--is_retry_on_total_extinction', type=bool, required=False,
    #                     help='')
    # parser.add_argument('--max_retries', required=False, help='')

    # Seq-Gen args, default.yaml?
    #         s.state_freqs = [0.4, 0.4, 0.1, 0.1]
    #         s.general_rates = [0.8, 0.4, 0.4, 0.2, 0.2, 0.1]

    # args = parser.parse_args()
    # args.parser = parser
    #main()


# --schema nexus --incipient_species_extinction_rate 0.1 --speciation_initiation_from_orthospecies_rate 0.2 --speciation_initiation_from_incipient_species_rate 0.3 --speciation_completion_rate 0.4 --orthospecies_extinction_rate 0.3 --aincipient_species_extinction_rate 0.2 --max_time 140 -o data