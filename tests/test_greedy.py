import pytest

from Consensus import *

import dendropy

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

def test_refine_majority():
    input_trees_path = files('Consensus.example_data').joinpath('GTRgamma_edit.nex')
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="nexus")
    majority = input_trees.majority_consensus()
    
    # from majority
    stdg_maj = STDGreedyConsensus(input_trees)
    stdg_maj.specify_initial_tree(majority)
    stdg_maj.greedy()
    
    tree1 = stdg_maj.return_current_tree()
    
    print("refine majority")
    
    stdg_maj.reset_initial_tree()
    stdg_maj.greedy(refine_majority=True)
    tree2 = stdg_maj.return_current_tree()
    
    print("DIFF should be 0", dendropy.calculate.treecompare.symmetric_difference(tree1, tree2))
    
    
    sqdg_maj = SQDGreedyConsensus(input_trees)
    sqdg_maj.specify_initial_tree(majority)
    sqdg_maj.greedy()
    
    tree1 = sqdg_maj.return_current_tree()
    
    print("refine majority")
    
    sqdg_maj.reset_initial_tree()
    sqdg_maj.greedy(refine_majority=True)
    tree2 = sqdg_maj.return_current_tree()
    
    print("DIFF should be 0", dendropy.calculate.treecompare.symmetric_difference(tree1, tree2))
    
    
    sutdg_maj = STDGreedyConsensus(input_trees)
    sutdg_maj.specify_initial_tree(majority)
    sutdg_maj.greedy()
    
    tree1 = sutdg_maj.return_current_tree()
    
    print("refine majority")
    
    sutdg_maj.reset_initial_tree()
    sutdg_maj.greedy(refine_majority=True)
    tree2 = sutdg_maj.return_current_tree()
    
    print("DIFF should be 0", dendropy.calculate.treecompare.symmetric_difference(tree1, tree2))

 
def test_SUTDGreedyConsensus():
    #t = Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.tre'), schema="newick")
    input_trees_path = files('Consensus.example_data').joinpath('GTRgamma_edit.nex')
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="nexus")
    majority = input_trees.majority_consensus()
    
    # from majority
    stdg_maj = SUTDGreedyConsensus(input_trees)
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
    
    greedy = astral.std_greedy(input_trees, normalized=False)
    stdg_astral.specify_initial_tree(greedy)
    print(stdg_astral.current_loss())
        
    
def test_SQDGreedyConsensus():
    input_trees_path = files('Consensus.example_data').joinpath('GTRgamma_edit.nex')
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="nexus")
    majority = input_trees.majority_consensus()
    
    sqd_maj = SQDGreedyConsensus(input_trees)
    sqd_maj.specify_initial_tree(majority)
    sqd_maj.greedy()