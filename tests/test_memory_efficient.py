import pytest
from importlib.resources import files
import ctypes
from bitstring import Bits
import Consensus
import numpy as np



class BipartitionMatchesHash(ctypes.Structure):
    _fields_ = [
        ("h1", ctypes.c_uint),
        ("h2", ctypes.c_uint),
        ("topo_depth", ctypes.c_int),
        ("matched_ids", ctypes.POINTER(ctypes.c_int)),
        ("td", ctypes.POINTER(ctypes.c_int)),
        ("input_count", ctypes.c_int)
    ]

class BipartitionDictHash(ctypes.Structure):
    _fields_ = [
        ("entries", ctypes.POINTER(ctypes.POINTER(BipartitionMatchesHash))),
        ("num_entries", ctypes.c_int),
        ("num_matches", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
    ]
    

class BipartitionSupportHash(ctypes.Structure):
    _fields_ = [
        ("h1", ctypes.c_uint),
        ("h2", ctypes.c_uint),
        ("is_external", ctypes.c_bool),
        ("support", ctypes.c_double),
        ("node_id", ctypes.c_int),
        ("bipartition", ctypes.POINTER(ctypes.c_uint))
    ]
    
class BipartitionSupportArrayHash(ctypes.Structure):
    _fields_ = [
        ("num_entries", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
        ("bipartition_supports", ctypes.POINTER(ctypes.POINTER(BipartitionSupportHash))),
        ("num_alt_trees", ctypes.c_int)
    ]
    
    
class TbeSupportMatchResult(ctypes.Structure):
    _fields_ = [
        ("bdict", ctypes.POINTER(BipartitionDictHash)),
        ("bsupp_arr", ctypes.POINTER(BipartitionSupportArrayHash))
    ]
    
    
class transferDistance(ctypes.Structure):
    _fields_ = [
        ("dist", ctypes.c_int),
        ("h1", ctypes.c_uint),
        ("h2", ctypes.c_uint),
        ("node_id", ctypes.c_int)
    ]

class transferDistanceAll(ctypes.Structure):
    _fields_ = [
        ("td", ctypes.POINTER(ctypes.POINTER(transferDistance))),
        ("n_elements", ctypes.c_int),
        ("node_h1", ctypes.c_uint),
        ("node_h2", ctypes.c_uint),
        ("topo_depth", ctypes.c_int)
    ]




def test_call_booster_tbe_support_and_match(K=20):
    # example trees
    print(files('Consensus.example_data').joinpath('sample1.nex'))
    input_trees = Consensus.TreeList_with_support.get(path = files('Consensus.example_data').joinpath('boot10.tre'), schema="newick")
    reference_tree = Consensus.Tree_with_support.get(path = files('Consensus.example_data').joinpath('astral_GTRgamma.tre'), schema="newick",
                                            taxon_namespace = input_trees.taxon_namespace)
    num_trees = len(input_trees)
    
    ## For using booster
    input_trees_str = input_trees.as_string("newick", suppress_rooting=True)
    input_trees_lines = input_trees_str.split("\n")
    reference_tree_str = reference_tree.as_string("newick", suppress_rooting=True)
    
    booster_lib = Consensus.load_booster()
    booster_lib.tbe_support_and_match.restype = ctypes.POINTER(TbeSupportMatchResult)

    # Specify the argument types of the function
    booster_lib.tbe_support_and_match.argtypes = [
        ctypes.c_char_p,                    # char* nh_initial_tree
        ctypes.POINTER(ctypes.c_char_p),    # char** nh_input_trees
        ctypes.c_int,                       # int num_trees
        ctypes.c_int                        # int k
    ]
    
    reference_tree_cstr = ctypes.c_char_p(reference_tree_str.encode('utf-8'))
    input_trees_array = (ctypes.c_char_p * len(input_trees_lines))(*[line.encode('utf-8') for line in input_trees_lines])
    result = booster_lib.tbe_support_and_match(reference_tree_cstr, input_trees_array, num_trees, K)
    
    print(hex(ctypes.addressof(result.contents)))
    
    # do the support computation in pure Python
    reference_tree.compute_transfer_support(input_trees)
    # list of support values
    support_list = []
    for edge in reference_tree.internal_edges(exclude_seed_edge = True):
        key = edge.bipartition.split_as_int()
        support_list.append(reference_tree.transfer_support[key])
    
    # list of support values from C implementation 
    bsupp_arr = result.contents.bsupp_arr.contents
    # Number of entries in the bipartition_supoprts array
    num_entries = bsupp_arr.num_entries
    C_support_list = []
    for i in range(num_entries):
        bipartition_support = bsupp_arr.bipartition_supports[i].contents
        support = bipartition_support.support
        if not bipartition_support.is_external:
            C_support_list.append(support)
    
    assert( np.all(np.sort(support_list) == np.sort(C_support_list)) ) # the sets of support values are the same
    
    
    ## Check if matches contain reasonable values
    bdict = result.contents.bdict.contents
    print("Number of bipartitions", bdict.num_entries)
    print("Number of matches", bdict.num_matches)
    print("Number of taxa", bdict.n_taxa)
    h1_list = []
    h2_list = []
    for i in range(bdict.num_entries):
        h1 = bdict.entries[i].contents.h1
        h1_list.append(h1)
        h2 = bdict.entries[i].contents.h2
        h2_list.append(h2)
        topo_depth = bdict.entries[i].contents.topo_depth
        for j in range(bdict.num_matches):
            td = bdict.entries[i].contents.td[j]
            matched_id = bdict.entries[i].contents.matched_ids[j]
            # print(h1,h2,topo_depth, td,matched_id)
    
    print(set(h1_list), len(set(h1_list)))
    print(set(h2_list), len(set(h1_list)))

    ## freeing
    booster_lib.free_tbe_support_and_match.argtypes = [ctypes.POINTER(TbeSupportMatchResult)]
    booster_lib.free_tbe_support_and_match.restype = None
    booster_lib.free_tbe_support_and_match(result)
    
    print("freed memory")


def test_recomputation(K=10): # example trees
    print(files('Consensus.example_data').joinpath('sample1.nex'))
    input_trees = Consensus.TreeList_with_support.get(path = files('Consensus.example_data').joinpath('boot10.tre'), schema="newick")
    reference_tree = Consensus.Tree_with_support.get(path = files('Consensus.example_data').joinpath('astral_GTRgamma.tre'), schema="newick",
                                            taxon_namespace = input_trees.taxon_namespace)
    
    reference_tree.STD_greedy_pruning(input_trees)
    
    num_trees = len(input_trees)
    
    ## For using booster
    input_trees_str = input_trees.as_string("newick", suppress_rooting=True)
    input_trees_lines = input_trees_str.split("\n")
    reference_tree_str = reference_tree.as_string("newick", suppress_rooting=True)
    
    booster_lib = Consensus.load_booster()
    booster_lib.tbe_support_and_match.restype = ctypes.POINTER(TbeSupportMatchResult)

    # Specify the argument types of the function
    booster_lib.tbe_support_and_match.argtypes = [
        ctypes.c_char_p,                    # char* nh_initial_tree
        ctypes.POINTER(ctypes.c_char_p),    # char** nh_input_trees
        ctypes.c_int,                       # int num_trees
        ctypes.c_int                        # int k
    ]
    
    reference_tree_cstr = ctypes.c_char_p(reference_tree_str.encode('utf-8'))
    input_trees_array = (ctypes.c_char_p * len(input_trees_lines))(*[line.encode('utf-8') for line in input_trees_lines])
    result = booster_lib.tbe_support_and_match(reference_tree_cstr, input_trees_array, num_trees, K)
    
    print(hex(ctypes.addressof(result.contents)))
    
    # do the support computation in pure Python
    reference_tree.compute_transfer_support(input_trees)
    # list of support values
    support_list = []
    for edge in reference_tree.internal_edges(exclude_seed_edge = True):
        key = edge.bipartition.split_as_int()
        support_list.append(reference_tree.transfer_support[key])
    
    # list of support values from C implementation 
    bsupp_arr = result.contents.bsupp_arr.contents
    # Number of entries in the bipartition_supoprts array
    num_entries = bsupp_arr.num_entries
    C_support_list = []
    for i in range(num_entries):
        bipartition_support = bsupp_arr.bipartition_supports[i].contents
        support = bipartition_support.support
        if not bipartition_support.is_external:
            C_support_list.append(support)
    
    print(np.sort(support_list))
    print(np.sort(C_support_list))
    
    assert( np.all(np.sort(support_list) == np.sort(C_support_list)) ) # the sets of support values are the same
    
    
    ## Check if matches contain reasonable values
    bdict = result.contents.bdict.contents
    print("Number of bipartitions", bdict.num_entries)
    print("Number of matches", bdict.num_matches)
    print("Number of taxa", bdict.n_taxa)
    h1_list = []
    h2_list = []
    input_count_list = []
    for i in range(bdict.num_entries):
        h1 = bdict.entries[i].contents.h1
        h1_list.append(h1)
        h2 = bdict.entries[i].contents.h2
        h2_list.append(h2)
        input_count_list.append(bdict.entries[i].contents.input_count)
        topo_depth = bdict.entries[i].contents.topo_depth
        for j in range(bdict.num_matches):
            td = bdict.entries[i].contents.td[j]
            matched_id = bdict.entries[i].contents.matched_ids[j]
            # print(h1,h2,topo_depth, td,matched_id)
    
    print(set(h1_list), len(set(h1_list)))
    print(set(h2_list), len(set(h1_list)))
    print(set(input_count_list))
    
    
    # recomputation for some branches
    booster_lib.recompute.restype = ctypes.POINTER(transferDistanceAll)

    # Specify the argument types of the function
    booster_lib.recompute.argtypes = [
        ctypes.c_uint,                       # unsigned int h1
        ctypes.c_uint                        # unsigned int h2
    ]
    
    recomp_res = booster_lib.recompute(h1_list[0], h2_list[0])
    # check if some reasonable values are returned.
    
    assert(recomp_res.contents.node_h1 == h1_list[0])
    assert(recomp_res.contents.node_h2 == h2_list[0])
    print(f"topological depth {recomp_res.contents.topo_depth}")
    for i in range(recomp_res.contents.n_elements):
        print(recomp_res.contents.td[i].contents.dist, recomp_res.contents.td[i].contents.h1, recomp_res.contents.td[i].contents.h2)
    
    

    ## freeing
    booster_lib.free_tbe_support_and_match.argtypes = [ctypes.POINTER(TbeSupportMatchResult)]
    booster_lib.free_tbe_support_and_match.restype = None
    booster_lib.free_tbe_support_and_match(result)
    
    booster_lib.free_TDA_LIST.argtypes = []
    booster_lib.free_TDA_LIST.restype=None
    booster_lib.free_TDA_LIST()
    
    print("freed memory")
    
    
def test_c_prune(K=20):
    #Consensus.c_prune(inittree_file=files('Consensus.example_data').joinpath('astral_GTRgamma.tre'), inputtrees_file=files('Consensus.example_data').joinpath('boot10.tre'), K=K)
    Consensus.c_prune(inittree_file=files('Consensus.example_data').joinpath('astral_GTRgamma.tre'), inputtrees_file=files('Consensus.example_data').joinpath('boot10.tre'), K=2)