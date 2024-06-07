import pytest
from Consensus import *
from importlib.resources import files
import numpy as np



def test_metric_behaviour():
    input_trees_path = files('Consensus.example_data').joinpath('std_ex.tre')
    input_trees = TreeList_with_support.get(path = input_trees_path, schema="newick")
    tree1 = Tree_with_support.get(data = "((A,B),(C,D), ((E,F),G), H);", schema="newick", taxon_namespace=input_trees.taxon_namespace)
    
    # Check STD
    std_fp, std_fn = STD_fp_fn(tree1, input_trees)
    epsilon = 10**(-10)
    for i in range(3):
        assert np.abs(std_fp[i]) < epsilon
        assert np.abs(std_fn[i] - 1/3) < epsilon
    for i in range(3,5):
        assert np.abs(std_fp[i] - 1/2) < epsilon
        assert np.abs(std_fn[i] - 5/3) < epsilon
    std_fp, std_fn = STD_fp_fn(tree1, input_trees[0])
    assert np.abs(std_fp) < epsilon
    assert np.abs(std_fn - 1/3) < epsilon
    
    # Check SUTD
    sutd_fp, sutd_fn = SUTD_fp_fn(tree1, input_trees)
    for i in range(3):
        assert np.abs(sutd_fp[i]) < epsilon
        assert np.abs(sutd_fn[i] - 1) < epsilon
    for i in range(3,5):
        assert np.abs(sutd_fp[i] - 1) < epsilon
        assert np.abs(sutd_fn[i] - 3) < epsilon
    sutd_fp, sutd_fn = SUTD_fp_fn(tree1, input_trees[0])
    assert np.abs(sutd_fp) < epsilon
    assert np.abs(sutd_fn - 1) < epsilon

    # Check SQD
    sqd_trees = TreeList_with_support.get(data = "(((A,B),C), (D,E));\n((A,B), (D,(C,E)));", schema="newick")
    sqd_con1 = Tree_with_support.get(data = "((A,B),C,D,E);", schema="newick", taxon_namespace=sqd_trees.taxon_namespace)
    sqd_con2 = Tree_with_support.get(data = "(A,(B,C),D,E);", schema="newick", taxon_namespace=sqd_trees.taxon_namespace)
    sqd_fp, sqd_fn = SQD_fp_fn(sqd_con1, sqd_trees)
    assert np.abs(sqd_fp[0]) == 0
    assert np.abs(sqd_fp[1]) == 0
    assert np.abs(sqd_fn[0]) == 2
    assert np.abs(sqd_fn[1]) == 2
    sqd_fp, sqd_fn = SQD_fp_fn(sqd_con1, sqd_con2)
    assert sqd_fp == 3
    assert sqd_fn == 3

    # Check SBD
    sbd_fp, sbd_fn = SBD_fp_fn(tree1, input_trees)
    epsilon = 10**(-10)
    for i in range(3):
        assert np.abs(sbd_fp[i]) < epsilon
        assert np.abs(sbd_fn[i] - 1) < epsilon
    for i in range(3,5):
        assert np.abs(sbd_fp[i] - 1) < epsilon
        assert np.abs(sbd_fn[i] - 2) < epsilon
    sbd_fp, sbd_fn = SBD_fp_fn(tree1, input_trees[0])
    assert np.abs(sbd_fp) < epsilon
    assert np.abs(sbd_fn - 1) < epsilon