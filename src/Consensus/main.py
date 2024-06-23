import argparse
import sys
import os
from ._consensus import *
from ._greedy import *


def _get_format_from_extension(filename):
    # Dictionary of extensions and their corresponding formats
    extension_format = {
        ".tre": "newick",
        ".nwk": "newick",
        ".nw": "newick",
        ".nex": "nexus",
        ".phy": "phylip",
        ".phylip": "phylip",
    }

    # Iterate over the dictionary and check if the filename ends with the extension
    for ext, fmt in extension_format.items():
        if filename.endswith(ext):
            return fmt

    # Return None if no matching extension is found
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_trees', type=str, help="Path to input trees")
    parser.add_argument('-is', '--input_schema', type=str, help="Schema of input trees", default="newick")
    parser.add_argument('-m', '--method', type=str, help="Greedy consensus method")
    parser.add_argument('-s', '--starting_tree', type=str, help="Starting tree")
    parser.add_argument('-ss', '--starting_tree_schema', type=str, help="Schema of starting tree", default="newick")
    parser.add_argument('-o', '--output_prefix')
    
    args = parser.parse_args()
    
    if args.input_trees is None:
        sys.exit("Please specify the path to input trees with -i option.")
    if args.input_schema is None:
        input_format = _get_format_from_extension(args.input_trees)
        if input_format is None:
            sys.exit("Please specify the schema of input trees by -is argument, or use the standard suffix.")
    else:
        input_format = args.input_schema
        if input_format not in ["newick", "nexus", "phylip"]:
            sys.exit("Please specify a correct schema for input trees (-is).")
    
    if args.method is None:
        sys.exit("Please specify the consensus method with -m option.")
    elif args.method not in ["STDG", "SQDG", "USTDG"]:
        sys.exit("Please choose the method from 'STDG")
     
    try:   
        input_trees = TreeList_with_support.get(path =args.input_trees, schema=input_format)
    except Exception as e:
        sys.exit(f"Error in reading input trees: \n{e}")
    
    if args.starting_tree is None:
        if args.method == "SQDG":
            sys.exit(f"You need to specify the path to (nearly) fully-resolved startnig tree with SQDG approach.")
        else:
            starting_tree = input_trees.majority_rule_consensus()
    else:
        if args.starting_tree_schema is None:
            ss_format = _get_format_from_extension(args.starting_tree)
            if ss_format is None:
                sys.exit("Please specify the schema of the starting tree by -ss argument, or use the standard suffix.")
        else:
            ss_format = args.starting_tree_schema
            if ss_format not in ["newick", "nexus", "phylip"]:
                sys.exit("Please specify a correct schema for the startinh tree (-ss).")
    try:
        starting_tree = Tree_with_support.get(path =args.starting_tree, schema=ss_format, taxon_namespace=input_trees.taxon_namespace)
    except Exception as e:
        sys.exit(f"Error in reading input trees: \n{e}")

    if args.output_prefix is None:
        output_file = args.input_trees + f".{args.method}.nex"
    else:
        output_file = args.output_prefix + ".nex"
        dir_name = os.path.dirname(output_file)
        if not os.path.exists(dir_name):
            sys.exit(f"Directory '{dir_name}' does not exist.")
    
    
    ## main computation
    
    method_dict = {
        "STDG": STDGreedyConsensus,
        "SUTDG": SUTDGreedyConsensus,
        "SQDG": SQDGreedyConsensus
    }
    
    gre = method_dict[args.method](input_trees)
    gre.specify_initial_tree(starting_tree)
    gre.greedy(method="most")
    
    gre.return_current_tree().write(path = output_file, schema="nexus")
    
    