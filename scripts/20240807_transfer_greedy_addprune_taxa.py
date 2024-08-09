import time
from importlib.resources import files
from Consensus import *
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="number and signal")
parser.add_argument('taxa',type=int,help="taxa number")
parser.add_argument('number', type=int, help="sample number")
parser.add_argument('--signal', type=str,choices=['high', 'low'], default="high", help="high or low")
args = parser.parse_args()

# Construct the file path using the provided number
INPUT_TREE_PATH = files(f"Consensus.sample{args.taxa}").joinpath(f"{args.signal}sample{args.number}.tre")

input_trees = TreeList_with_support.get(path = INPUT_TREE_PATH, schema = "newick")
majority = input_trees.majority_rule_consensus()
#t,ts = plot_example_func(majority)
#t.render(file_name="test.png",tree_style=ts)
start = time.process_time()
stdg = STDGreedyConsensus(input_trees)
stdg.specify_initial_tree(majority) # starting from majority tree
stdg.greedy(method="first", order="BS")
end = time.process_time()
print("TIME_all:", end = " ", flush=True)
print('{:.2f}'.format((end-start)), flush=True)
#t,ts = plot_example_func(stdg.return_current_tree())
#t.render(file_name="test2.png",tree_style=ts)
