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

import src.gen as g
import src.cfg as c
from src import param_writer


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
    config = c.cfg_load(args.config)
    headers = c.load_headers(config)

    # start
    start = time.perf_counter()
    for _ in range(args.num_runs):
        # load run parameters
        generated_sample_parameters = g.gen_sample_values(config[('generate_sample')])
        generated_protracted_speciation_process_parameters = g.generate_protracted_speciation_process_values()
        seqgen_parameters = g.get_seq_gen_values(config[('seq-gen')])
        # getting trees
        get_trees = call_sample_tree(generated_sample_parameters, generated_protracted_speciation_process_parameters)
        # saving trees
        names = file_output(get_trees, args, ["lineage", "orthospecies"])
        # saving parameters
        for n in names:
            param_writer.parameters_prep(generated_sample_parameters, generated_protracted_speciation_process_parameters,
                              seqgen_parameters, args, headers, n)
        # generating sequences
        call_seq_gen(args, names, seqgen_parameters)
    print('\nProcess took %.2f seconds to complete.' % (time.perf_counter() - start))


def temp_file_name():
    """
    Generates random string
    :return : random string
    """
    temp_name = next(tempfile._get_candidate_names())
    return temp_name


def call_sample_tree(generated_sample_parameters, generated_protracted_speciation_process_parameters):
    """
    Calls ProtractedSpeciationProcess and generates sample trees.
    :return trees,
    generated_tree[0] lineage_tree (|Tree| instance) – A tree from the protracted speciation process, with all lineages
    (good species as well as incipient species).
    generated_tree[1] orthospecies_tree (|Tree| instance) – A tree from the protracted speciation process with only
    “good” species.:
    """
    while True:
        try:
            generated_trees = protractedspeciation.ProtractedSpeciationProcess(
                **generated_protracted_speciation_process_parameters).generate_sample(**generated_sample_parameters)
            return generated_trees
        except:
            continue


def file_output(trees, args, tree_names):
    """Stores output files."""
    output_dir = args.output
    re_names = []
    # sanity check
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if not os.path.isdir(args.output):
        sys.exit()
    a = temp_file_name()
    for i in range(len(trees)):
        if args.schema != 'nexus':
            file_name = tree_names[i][:3] + '_' + a + "." + str(args.schema)
            tmp_path = os.path.join(output_dir, file_name)
            trees[i].write_to_path(tmp_path, suppress_rooting=True, suppress_edge_lengths=True,
                                   schema=args.schema)
            new_fname = convert_newick_to_nexus(args, file_name)
            re_names.append(new_fname)
        else:
            file_name = tree_names[i][:3] + '_' + a + "." + str(args.schema)
            tmp_path = os.path.join(output_dir, file_name)
            trees[i].write_to_path(tmp_path, suppress_rooting=True, suppress_edge_lengths=True,
                                   schema=args.schema)
            re_names.append(file_name)
    for f in glob.glob(os.path.join(output_dir, "*.newick")):
        os.remove(f)
    return re_names


def convert_newick_to_nexus(args, fname):
    tree = dendropy.Tree.get(path=os.path.join(args.output, fname), schema='newick')
    new_file_name = fname.split('.')[0] + '.nexus'
    path = os.path.join(args.output, new_file_name)
    tree.write_to_path(path, suppress_rooting=True, suppress_edge_lengths=True,
                       schema="nexus")
    return new_file_name


def call_seq_gen(args, name, seqgen_vals):
    """
    Passing tree names to Seq-Gen.
    """
    for n in name:
        full_path = os.path.join(args.output, "seq_{}.txt".format(n.split('.')[0]))
        path_to_tree = os.path.join(args.output, n)
        trees = dendropy.TreeList.get(path=path_to_tree, schema='nexus')
        s = seqgen.SeqGen()
        # instruct Seq-Gen to scale branch lengths by factor of 0.1
        # note that this does not modify the input trees
        s.scale_branch_lens = 0.1
        for k, v in seqgen_vals.items():
            seqgen_vals[k] = seqgen.SeqGen(v)
        d1 = s.generate(trees)
        with open(full_path, "w") as f:
            f.write(d1.char_matrices[0].as_string('nexus'))


if __name__ == '__main__':
    main()
