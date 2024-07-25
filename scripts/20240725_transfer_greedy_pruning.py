
import time
from importlib.resources import files
from Consensus import *
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Specify the number after 'high_astral'")
parser.add_argument('number', type=int, help="sample number")
parser.add_argument('--signal', type=str,choices=['high', 'low'], default="high", help="high or low")
args = parser.parse_args()

# Construct the file path using the provided number
ASTRAL_TREE_PATH = files("Consensus.sample100").joinpath(f"{args.signal}_astral{args.number}.tre")
INPUT_TREE_PATH = files("Consensus.sample100").joinpath(f"{args.signal}sample{args.number}.tre")

input_trees = TreeList_with_support.get(path = INPUT_TREE_PATH, schema = "newick")
consensus_tree = Tree_with_support.get(path = ASTRAL_TREE_PATH, schema = "newick",taxon_namespace = input_trees.taxon_namespace )

astral_std = consensus_tree.STD_greedy_pruning(input_trees, normalized=True,time_flag =True)
