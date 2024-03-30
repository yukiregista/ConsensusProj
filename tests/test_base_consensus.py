import pytest

from Consensus import *


from importlib.resources import files
import numpy as np

def test_base_consensus():
    input_trees_path = files('Consensus.example_data').joinpath('GTR_edit.nex')
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="nexus")
    
    
    # majority
    maj = input_trees.majority_consensus()
    maj.compute_bootstrap(input_trees)
    print("\nmajority support sum:", np.sum(list(maj.bootstrap_support.values())))
    
    # MCC
    MCC = input_trees.MCC_tree()
    MCC.compute_bootstrap(input_trees)
    print("MCC support sum:", np.sum(list(MCC.bootstrap_support.values())))
    
    # MAP
    MAP = input_trees.MAP()[0][0]
    MAP.compute_bootstrap(input_trees)
    print("MAP support sum:", np.sum(list(MAP.bootstrap_support.values())))
    
    
    