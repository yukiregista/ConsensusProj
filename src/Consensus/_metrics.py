import sys

import numpy as np
from bitstring import Bits
import dendropy

from ._consensus import Tree_with_support, TreeList_with_support, transfer_support, unnormalized_transfer_support, _tqdist_fp_fn, _tqdist_fp_fn_trees_str



def _SUTD1_fp_fn(true_tree, tree2: Tree_with_support):
    tree1_int_bipars = [node.bipartition for node in true_tree.postorder_internal_node_iter(exclude_seed_node=True)]
    tree2_int_bipars = [node.bipartition for node in tree2.postorder_internal_node_iter(exclude_seed_node=True)]
    fp = np.sum(unnormalized_transfer_support(tree2_int_bipars, [true_tree]))
    fn = np.sum(unnormalized_transfer_support(tree1_int_bipars, [tree2]))
    return fp, fn

def _STD1_fp_fn(true_tree, tree2:  Tree_with_support):
    tree1_int_bipars = [node.bipartition for node in true_tree.postorder_internal_node_iter(exclude_seed_node=True)]
    tree2_int_bipars = [node.bipartition for node in tree2.postorder_internal_node_iter(exclude_seed_node=True)]
    fp = np.sum((1-transfer_support(tree2_int_bipars, [true_tree])))
    fn = np.sum((1-transfer_support(tree1_int_bipars, [tree2])))
    return fp, fn

def SUTD_fp_fn(estimate, inputs):
    """Computes false positive and false negative of unnormalized STD.
        
    Parameters
    ----------
    estimate : Tree_with_support
        Estimates.
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate. 
        

    Returns
    -------
    (float, float) or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives.
    """
    
    if isinstance(inputs, Tree_with_support):
        return _SUTD1_fp_fn(inputs, estimate)
    elif isinstance(inputs, TreeList_with_support):
        fp_list = []
        fn_list = []
        for i in range(len(inputs)):
            fp_tmp, fn_tmp = _SUTD1_fp_fn(inputs[i], estimate)
            fp_list.append(fp_tmp)
            fn_list.append(fn_tmp)
        return np.array(fp_list), np.array(fn_list)
    else:
        print("input trees have to be either Tree_with_support or TreeList_with_support.")
        sys.exit(1)
        
def STD_fp_fn(estimate, inputs):
    """Computes false positive and false negative of STD.
        
    Parameters
    ----------
    estimate : Tree_with_support
        Estimates.
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate.      

    Returns
    -------
    (float, float)  or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives.
    """
    
    if isinstance(inputs, Tree_with_support):
        return _STD1_fp_fn(inputs, estimate)
    elif isinstance(inputs, TreeList_with_support):
        fp_list = []
        fn_list = []
        for i in range(len(inputs)):
            fp_tmp, fn_tmp = _STD1_fp_fn(inputs[i], estimate)
            fp_list.append(fp_tmp)
            fn_list.append(fn_tmp)
        return np.array(fp_list), np.array(fn_list)
    else:
        print("input trees have to be either Tree_with_support or TreeList_with_support.")
        sys.exit(1)
        
        
        
def SQD_fp_fn(estimate, inputs, parent_dir=None):
    """Compute false positives and false negatives of symmetric quartet distance.
    

    Parameters
    ----------
    estimate : Tree_with_support
        Estimates.
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate. 
    parent_dir : str, optional
        Path to place intermediate files, by default None (place files in the folder of execution).
        
    Returns
    -------
    (float, float) or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives respectively.
    """
    if isinstance(inputs, Tree_with_support):
        fp, fn = _tqdist_fp_fn(estimate, inputs, parent_dir=parent_dir)
        return fp, fn
    elif isinstance(inputs, TreeList_with_support):
        inputs_string = inputs.as_string("newick", suppress_rooting=True)
        fp, fn = _tqdist_fp_fn_trees_str(estimate, inputs_string, len(inputs), parent_dir = parent_dir)
        return fp, fn
    else:
        print("Please provide an instance of Tree_with_support or TreeList_with_support as inputs.")
        sys.exit(1)
   
   
def SBD_fp_fn(estimate, inputs):
    """Compute false positives and false negatives of symmetric bipartition distance.
    

    Parameters
    ----------
    estimate : Tree_with_support
        Estimates.
    inputs : Tree_with_support or TreeList_with_support
        Input trees to evaluate the estimate. 
    parent_dir : str, optional
        Path to place intermediate files, by default None (place files in the folder of execution).
        
    Returns
    -------
    (float, float) or (numpy.ndarray, numpy.ndarray)
        False positives and False negatives respectively.
    """
    if isinstance(inputs, Tree_with_support):
        fp, fn = dendropy.calculate.treecompare.false_positives_and_negatives(inputs, estimate)
        return fp, fn
    elif isinstance(inputs, TreeList_with_support):
        fps = []; fns=  []
        for tree in inputs:
            fp, fn = dendropy.calculate.treecompare.false_positives_and_negatives(tree, estimate)
            fps.append(fp); fns.append(fn)
        return np.array(fps), np.array(fns)