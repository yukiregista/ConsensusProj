from Consensus import *
from Consensus._consensus import _compatible
import numpy as np
from importlib.resources import files
import sys

from bitstring import Bits
import matplotlib.pyplot as plt

import argparse


def check_n_comptaible_elements(input_trees_path: str):
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="newick")
    
    print(input_trees.n_taxa)
    
    edges = []
    for tree in input_trees:
        for internal_edge in tree.internal_edges(exclude_seed_edge=True):
            edges.append(internal_edge.bipartition.split_as_int())
    edges = np.unique(edges)
    edge_supports = np.array([input_trees.edge_dict[item] for item in edges])
    order = np.argsort(edge_supports)[::-1]
    edges = [Bits(uint=item, length=input_trees.n_taxa) for item in edges[order]]
    edge_supports = edge_supports[order]
    pick = np.sum(edge_supports > 500)
    
    n_bipars = [len(edges)]
    support_threshold = [100]
    edges_copy = edges.copy()
    
    for i in range(pick):
        bipar = edges[i]
        keep_indices = []
        for j in range(len(edges_copy)):
            if _compatible(bipar, edges_copy[j]):
                keep_indices.append(j)
        edges_copy = [edges_copy[item] for item in keep_indices]
        n_bipars.append(len(edges_copy))
        support_threshold.append(float(edge_supports[i])/10)
        
    print(n_bipars)
    print(support_threshold)
    
    fig, ax = plt.subplots()
    ax.scatter(support_threshold, n_bipars)
    ax.set_xlabel("Threshold of support value")
    ax.set_ylabel("Number of compatible bipartitions")
    plt.show()
    
parser = argparse.ArgumentParser()

parser.add_argument("--input_path")


if __name__ == "__main__":
    args = parser.parse_args()
    if args.input_path:
        check_n_comptaible_elements(args.input_path)
    else:
        input_trees_path = files('Consensus.example_data').joinpath('GTR_edit.tre')
        check_n_comptaible_elements(input_trees_path)