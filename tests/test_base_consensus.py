import pytest

from Consensus import *


from importlib.resources import files
import numpy as np

def test_base_consensus():
    input_trees_path = files('Consensus.example_data').joinpath('GTR_edit.nex')
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="nexus")
    
    
    # majority
    maj = input_trees.majority_consensus()
    maj.compute_branch_support(input_trees)
    maj.compute_transfer_support(input_trees)
    print("\nmajority support sum:", np.sum(list(maj.branch_support.values())))
    
    # MCC
    MCC = input_trees.MCC_tree()
    MCC.compute_branch_support(input_trees)
    print("MCC support sum:", np.sum(list(MCC.branch_support.values())))
    
    # MAP
    MAP = input_trees.MAP()[0][0]
    MAP.compute_branch_support(input_trees)
    print("MAP support sum:", np.sum(list(MAP.branch_support.values())))
    
    
    