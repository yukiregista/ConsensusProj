import pytest

from Consensus import *


from importlib.resources import files
import numpy as np


def test_STDGreedyConsensus():
    #t = Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.tre'), schema="newick")
    input_trees_path = files('Consensus.example_data').joinpath('GTRgamma_edit.nex')
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="nexus")
    majority = input_trees.majority_consensus()
    
    # from majority
    stdg_maj = STDGreedyConsensus(input_trees)
    stdg_maj.specify_initial_tree(majority)
    stdg_maj.greedy()
    
    current_tree2 = stdg_maj.return_current_tree()
    print(current_tree2)
    print(current_tree2.branch_resolution())
    print(stdg_maj.current_loss())

    
    # stdg_cls = STDGreedyConsensus(input_trees)
    # stdg_cls.greedy()

    # current_tree = stdg_cls.return_current_tree()
    
    
    # start from astral
    astral = Tree_with_support.get(path = files('Consensus.example_data').joinpath('astral_GTRgamma.tre'),schema="newick", taxon_namespace = input_trees.taxon_namespace)
    #stdg_astral = STDGreedyConsensus(input_trees)
    stdg_astral = stdg_maj
    stdg_astral.specify_initial_tree(astral)
    stdg_astral.greedy()
    current_tree = stdg_astral.return_current_tree()
    
    print(current_tree)
    print(current_tree.branch_resolution())
    print(stdg_astral.current_loss())
    
    stdg_astral.specify_initial_tree(majority)
    stdg_astral.greedy(method = "most")
    current_tree = stdg_astral.return_current_tree()
    
    print(current_tree)
    print(current_tree.branch_resolution())
    print(stdg_astral.current_loss())
    
    greedy = astral.std_greedy(input_trees)
    print(greedy)
    print(greedy.branch_resolution())
    print(np.sum(STD_fp_fn(greedy, input_trees)))
    
    
    