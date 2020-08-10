#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import time
import dendropy
import tempfile
# import _thread as thread

import config_load as cl
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
    parser.add_argument('--override', '-r', action='store_true', required=False, help='Overrides existing files.')
    parser.add_argument('--profile', '-p', default="default", required=True, help='')
    parser.add_argument('--num_runs', '-n', default=1, type=int, required=False, help='')

    args = parser.parse_args()
    args.parser = parser
    config = cl.Profile(args.profile)

    start = time.perf_counter()
    # thread.start_new_thread(waiting_std_output, (20, 0.5))

    try:
        print("Creating trees and sequences.")
        # get_trees = []
        for _ in range(args.num_runs):
            # getting trees
            get_trees = call_sample_tree(args, config)
            #print("\nTree files finished.")
            # saving trees
            file_output(get_trees, args, ["lineage", "orthospecies"])

            # generating Sequences
            # print("\nStarting Seq-Gen")
            # call_seq_gen(lineage_path, args, "lineage")
            # call_seq_gen(orthospecies_path, args, "orthospecies")
            # print("Finished.")
            # print('\nProcess took %.2f seconds to complete.' % (time.perf_counter() - start))
    except BaseException as e:
        print(e)


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
    :return : args with values only
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
    psp_ini = protractedspeciation.ProtractedSpeciationProcess(**config.get_protractedspeciationprocess_values())
    try:
        # calling for args with values
        values = gen_sample_values(config.get_generate_sample_values())
        # generate trees
        generated_trees = psp_ini.generate_sample(**values)
        return generated_trees
    except BaseException as e:
        print("Maximum number of runs to execute in the event of prematurely-terminated simulations due to all "
              "lineages going extinct. Once this number or re-runs is exceed, then TreeSimTotalExtinctionException "
              "is raised. Defaults to 1000. Set to None to never quit trying.\n" + str(e))


def file_output(trees, args, tree_names):
    """Stores output files."""

    output_dir =  os.path.join(os.getcwd(), 'scripts', args.output + '/trees')
    #output_dir =  + args.output + '/trees'
    print(output_dir)
    # sanity check
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(output_dir):
        sys.exit()
    for i in range(len(trees)):
        try:
            filename = tree_names[i] + '_' + temp_file_name() + "." + str(args.schema)
            tmp_path = os.path.join(output_dir, filename)
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
        # print(len(d1.char_matrices))
        # print(d1.char_matrices[0].as_string(args.schema))
        with open(full_path, "w") as f:
            f.write(d1.char_matrices[0].as_string(args.schema))
        print("Seq-Gen file stored in {}".format(os.path.join(os.getcwd(), output_dir)))
    except BaseException as e:
        print("Unexpected error while using Seq-Gen:\n" + str(e))


if __name__ == '__main__':
    main()
