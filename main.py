#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import tempfile
import time
import glob
import dendropy

from dendropy.interop import seqgen
from dendropy.model import protractedspeciation

from src import param_writer
from src import loader

VERSION = "0.0.1"


def main():
    parser = argparse.ArgumentParser(description='Generate sample tree data under the protracted speciation model')
    parser.add_argument('--output', '-o', required=False, help='Output dir')
    parser.add_argument('--schema', '-s', choices=['newick', 'nexus'], required=True,
                        help='Tree schema: Newick, Nexus')
    parser.add_argument('--config', '-c', default="default", required=True, help='')
    parser.add_argument('--num_runs', '-n', default=1, type=int, required=False, help='')

    args = parser.parse_args()
    args.parser = parser
    config = loader.Loader(args.config, main=1)
    headers = config.load_headers()

    # start
    start = time.perf_counter()
    for _ in range(args.num_runs):
        # try:
        # getting trees
        get_trees = call_sample_tree(args, config)
        # generating Sequences & saving trees
        for names in file_output(get_trees, args, ["lineage", "orthospecies"]):
            # generating sequences
            call_seq_gen(args, names, config)
        # saving parameters
        param_writer.parameters_to_txt(config, args, headers)
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
    while True:
        try:
            # calling args
            values = gen_sample_values(config.get_generate_sample_values())
            # generate trees
            generated_trees = protractedspeciation.ProtractedSpeciationProcess(
                **config.generate_protracted_speciation_process_values()).generate_sample(**values)
            # generated_trees = protractedspeciation.ProtractedSpeciationProcess(
            #     **config.generate_protracted_speciation_process_values()).generate_sample(**values)
            return generated_trees
        except:
            continue


def file_output(trees, args, tree_names):
    """Stores output files."""
    output_dir = args.output
    # sanity check
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if not os.path.isdir(args.output):
        sys.exit()
    for i in range(len(trees)):
        if args.schema != 'nexus':
            file_name = tree_names[i][:3] + '_' + temp_file_name() + "." + str(args.schema)
            tmp_path = os.path.join(output_dir, file_name)
            trees[i].write_to_path(tmp_path, suppress_rooting=True, suppress_edge_lengths=True,
                                   schema=args.schema)
            new_fname = convert_newick_to_nexus(args, file_name)
            yield new_fname
        else:
            file_name = tree_names[i][:3]+ '_' + temp_file_name() + "." + str(args.schema)
            tmp_path = os.path.join(output_dir, file_name)
            trees[i].write_to_path(tmp_path, suppress_rooting=True, suppress_edge_lengths=True,
                                   schema=args.schema)
            yield file_name
        for f in glob.glob(os.path.join(output_dir, "*.newick")):
            os.remove(f)


def convert_newick_to_nexus(args, fname):
    tree = dendropy.Tree.get(path=os.path.join(args.output, fname), schema='newick')
    new_file_name = fname.split('.')[0] + '.nexus'
    path = os.path.join(args.output, new_file_name)
    tree.write_to_path(path, suppress_rooting=True, suppress_edge_lengths=True,
                                                 schema="nexus")
    return new_file_name

def call_seq_gen(args, name, config):
    """
    Passing tree names to Seq-Gen.
    Calls randomized parameter from Loader.
    """
    full_path = os.path.join(args.output, "seq_{}.txt".format(name.split('.')[0]))
    # try:
    path_to_tree = os.path.join(args.output, name)
    trees = dendropy.TreeList.get(path=path_to_tree, schema='nexus')
    s = seqgen.SeqGen()

    # generate one alignment per tree
    # as substitution model is not specified, defaults to a JC model
    # will result in a DataSet object with one DnaCharacterMatrix per input tree
    d0 = s.generate(trees)
    # instruct Seq-Gen to scale branch lengths by factor of 0.1
    # note that this does not modify the input trees
    s.scale_branch_lens = 0.1
    seqgen_vals = config.get_seq_gen_values()
    for k, v in seqgen_vals.items():
        seqgen_vals[k] = seqgen.SeqGen(v)
    d1 = s.generate(trees)
    with open(full_path, "w") as f:
        f.write(d1.char_matrices[0].as_string('nexus'))


if __name__ == '__main__':
    main()
