from ._booster import load_booster, prepare_tbe_support_and_match, prepare_tbe_support_and_match_args, organize_support_and_match, OrganizedTransferResult, _generate_bitmask, prepare_recompute, transferDistanceAll, prepare_prune_and_return_newick, prepare_prune_and_return_newick_args
import dendropy
import numpy as np
import time
import ctypes
from tqdm import trange
from bitstring import Bits

import threading, psutil, os
from ._dev_utils import monitor_cpu

import gc


class PruningVars():
    """
        Class to hold all variables used in pruning.
    """
    def __init__(self, FP: np.ndarray, FN_diff: np.ndarray, E: set, R1: np.ndarray, R2: np.ndarray, W1: list[set], W2: list[set], OSR: OrganizedTransferResult, booster_lib: ctypes.CDLL):
        self.FP: np.ndarray = FP # False Positive Loss caused by each internal bipartitoin int the ref tree.
        self.FN_diff: np.ndarray = FN_diff # False Negaive Loss increase that will be induced by pruning each branch in the ref tree.
        self.E: set = E # Set of indices of ref_tree to keep
        self.pruned: set = set() # Set of indices that have been pruned.
        self.R1: np.ndarray = R1 # Rank of the first match (starts with index 0). 
        self.R2: np.ndarray = R2 # Rank of the second match (initialized with 1 at first).
        self.W1: list[set] = W1 # List of size maximum n-3, containing the set of branches matched to each internal bipar in ref tree.
        self.W2: list[set] = W2 # List of size maximum n-3, containing the set of branches second matched to each internal bipar in ref tree.
        self.L: np.ndarray = self.FP - self.FN_diff # Utility of pruning branches
        self.OSR: OrganizedTransferResult = OSR # OrganizedTransferResult
        self.use_K: list[bool] = [True for _ in range(len(self.R1))] # Whether to use the originally returned K matches, or the result of recomputation.
        self.booster_lib = booster_lib
        
        # some useful variables
        self.n_internal_ref = len(self.FP) # Number of internal edges in the initial tree
        self.m = len(self.R1) # Number of unique internal branches in the input trees
        self.K = self.OSR.M.shape[1] # K
        
        # For storing recomputation results
        self.RECOMP: dict[int, RecomputationResults] = dict()


class RecomputationResults():
    
    def __init__(self, M: np.ndarray, S: np.ndarray):
        self.M = M
        self.S = S

    # import tracemalloc

    # # we hope! that inittree_file and inputtrees_file have the same taxons
    # inittree = dendropy.Tree.get(path = inittree_file, schema="newick")
    # taxon_labels = np.sort([item.label for item in inittree.taxon_namespace])
    
    
def transfer_distance(a: tuple, b: tuple, bits: int, n_taxa: int):
    assert(len(a)==len(b))
    spl = len(a)
    mask = _generate_bitmask(bits)
    mask_last = _generate_bitmask(n_taxa % bits)
    hamdist1 = 0
    for i in range(spl):
        hamdist1 += bin(a[i] ^ b[i]).count('1')
    
    hamdist2 = 0
    for i in range(spl-1):
        hamdist2 += bin(a[i] ^ (~b[i] & mask)).count('1')
    hamdist2 += bin(a[spl-1] ^  (~b[spl-1] & mask_last)).count('1')
    
    return min(hamdist1, hamdist2)
        

def recompute_td(osr: OrganizedTransferResult, bipar_id: int, booster_lib: ctypes.CDLL):
    h1 = osr.H1[bipar_id]
    h2 = osr.H2[bipar_id]
    # if bipar_id==75: print("h1, h2:", h1, h2)
    tda: transferDistanceAll = booster_lib.recompute(ctypes.c_uint(h1), ctypes.c_uint(h2)).contents
    
    # We need to reconcile the order between the original OSR result and tda result.
    node_id_list = []
    tda_dict = dict()
    use_external=False
    for i in range(tda.n_elements):
        if (tda.td[i].contents.dist >= tda.topo_depth-1):
            continue
        else:  
            bipar_loc = osr.REF_BIPAR_ID_to_location[ tda.td[i].contents.node_id ]
            node_id_list.append( bipar_loc )
            tda_dict[bipar_loc] = (tda.td[i].contents.dist / (tda.topo_depth-1))


    if (len(node_id_list)) == 0:
        # Others are all filled with externals
        new_M = np.concatenate([osr.M[bipar_id], [-1 for _ in range(tda.n_elements - len(osr.M[bipar_id]))]]).astype(int)
        new_S = np.concatenate([osr.S[bipar_id], [1 for _ in range(tda.n_elements - len(osr.M[bipar_id]))]])
    else:
        first_K_IDS = set(osr.M[bipar_id])
        IDs_to_add = [item for item in node_id_list if item not in first_K_IDS]
        Values = [tda_dict[item] for item in IDs_to_add]
        # sort the values
        order = np.argsort(Values)
        M_mid = [IDs_to_add[x] for x in order]
        S_mid = [Values[x] for x in order]
        rem_elements = tda.n_elements - len(osr.M[bipar_id]) - len(M_mid)
        new_M = np.concatenate([ osr.M[bipar_id], M_mid, [-1 for _ in range(rem_elements)] ]).astype(int)
        new_S = np.concatenate([ osr.S[bipar_id], S_mid, [1 for _ in range(rem_elements)] ])
    
    return RecomputationResults(new_M, new_S)
            
    
     
def compute_all_normalized_td(bipar, matched, K, reverse_id_ref, bipar_val, topo_depth, bits, n_taxa):
    processed_set = set(matched)
    remaining_branches = [(id,val) for id, val in reverse_id_ref.items() if id not in processed_set]
    td = [transfer_distance(bipar_val, val[1], bits, n_taxa) for val in remaining_branches]
    order = np.argsort(td)
    ns = np.array([td[item]/(topo_depth - 1.0) for item in order])
    ids = [remaining_branches[item][0] for item in order]
    
    ext_index = 0
    for i in range(len(ns)):
        if ns[i] > topo_depth-2:
            ext_index=i
            break
    
    ids = [remaining_branches[item][0] for item in order[:ext_index]] + [-1 for _ in range(len(ns)-ext_index)]
    ns[ext_index:] = 1.0
    return ids, ns
    
    
def renew_match(pv: PruningVars, b: int):
    """
        renew match by pruning b.
    """
    pv.L[b] = -1 # Set to -1 so that it will not be picked again.
    ## update bipartitions in pv.W1[b]
    # print(f"W1 len: {len(pv.W1[b])}, W2 len: {len(pv.W2[b])}")
    for bipar_id in pv.W1[b]: # bipartitions in the first match
        # print("W1:", bipar_id)
        # Promote second match to the fist match
        second_match = pv.OSR.M[bipar_id, pv.R2[bipar_id]] if pv.use_K[bipar_id] else pv.RECOMP[bipar_id].M[pv.R2[bipar_id]]
        if (second_match!=-1): # If the original second match of bipar_id is not external one.
            pv.W1[second_match].add(bipar_id) # Add bipar_id to the original second match's W1
            pv.W2[second_match].remove(bipar_id) # Remove bipar_id from the original second match's W2
            pv.R1[bipar_id] = pv.R2[bipar_id]
        
            # Find new second match
            now_rank = pv.R1[bipar_id]+1
            while True:
                if pv.use_K[bipar_id]:
                    if now_rank > pv.K-1: # Then do recomputation
                        # recomputation
                        
                        pv.RECOMP[bipar_id] = recompute_td(pv.OSR, bipar_id, pv.booster_lib)
                        pv.use_K[bipar_id] = False
                    elif pv.OSR.M[bipar_id, now_rank] in pv.pruned: # Then proceed to next rank
                        # it is already pruned
                        now_rank += 1
                    else: # Then store the new second match and modify losses
                        if pv.OSR.M[bipar_id, now_rank]!=-1: # if matched to one of the internal branches
                            pv.W2[pv.OSR.M[bipar_id, now_rank]].add(bipar_id)
                        pv.R2[bipar_id] = now_rank
                        # Adjust the FN_diff and L of the first match
                        diff = (pv.OSR.S[bipar_id, pv.R2[bipar_id]] - pv.OSR.S[bipar_id, pv.R1[bipar_id]]) * pv.OSR.C[bipar_id]
                        pv.FN_diff[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] += diff
                        pv.L[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] -= diff
                        break
                        
                else:
                    if pv.RECOMP[bipar_id].M[now_rank] in pv.pruned: # Then proceed to next rank
                        now_rank += 1
                    else: # Then store the new second match and modify losses
                        if pv.RECOMP[bipar_id].M[now_rank] != -1:
                            pv.W2[pv.RECOMP[bipar_id].M[now_rank]].add(bipar_id)
                        pv.R2[bipar_id] = now_rank
                        # Adjust the FN_diff and L of the first match
                        diff = (pv.RECOMP[bipar_id].S[pv.R2[bipar_id]] - pv.RECOMP[bipar_id].S[pv.R1[bipar_id]]) * pv.OSR.C[bipar_id]
                        pv.FN_diff[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] += diff
                        pv.L[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] -= diff
                        break
     
    for bipar_id in pv.W2[b]: # Bipartitions that had b as the second match -> renew second match
        # print("W2:", bipar_id)
        now_rank = pv.R2[bipar_id] + 1
        while True:
            if pv.use_K[bipar_id]:
                if now_rank > pv.K-1: # Then do recomputation
                    # recomputation
                    pv.RECOMP[bipar_id] = recompute_td(pv.OSR, bipar_id, pv.booster_lib)
                    pv.use_K[bipar_id] = False
                elif pv.OSR.M[bipar_id, now_rank] in pv.pruned: # Then proceed to next rank
                    now_rank += 1
                else: # Then store the new second match and modify losses
                    if pv.OSR.M[bipar_id, now_rank]!=-1: # if matched to one of the internal branches
                        pv.W2[pv.OSR.M[bipar_id, now_rank]].add(bipar_id)
                    old_rank = pv.R2[bipar_id]
                    pv.R2[bipar_id] = now_rank
                    # Adjust the FN_diff and L of the first match
                    diff = (pv.OSR.S[bipar_id, pv.R2[bipar_id]] - pv.OSR.S[bipar_id, old_rank]) * pv.OSR.C[bipar_id]
                    # if (pv.OSR.M[bipar_id, pv.R1[bipar_id]]==77):
                    #     # print("diff: ", diff, "C:", pv.OSR.C[bipar_id], "new rank: ", now_rank, " new S: ", pv.OSR.S[bipar_id, pv.R2[bipar_id]])
                    pv.FN_diff[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] += diff
                    pv.L[pv.OSR.M[bipar_id, pv.R1[bipar_id]]] -= diff
                    break
            else:
                if pv.RECOMP[bipar_id].M[now_rank] in pv.pruned: # Then proceed to next rank
                    now_rank += 1
                else: # Then store the new second match and modify losses
                    if pv.RECOMP[bipar_id].M[now_rank] != -1:
                        pv.W2[pv.RECOMP[bipar_id].M[now_rank]].add(bipar_id)
                    old_rank = pv.R2[bipar_id]
                    pv.R2[bipar_id] = now_rank
                    # Adjust the FN_diff and L of the first match
                    diff = (pv.RECOMP[bipar_id].S[pv.R2[bipar_id]] - pv.RECOMP[bipar_id].S[old_rank]) * pv.OSR.C[bipar_id]
                    # if (pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]==77):
                    #     print("diff: ", diff, "C:", pv.OSR.C[bipar_id], "new rank: ", now_rank, " new S: ", pv.RECOMP[bipar_id].S[pv.R2[bipar_id]])
                    pv.FN_diff[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] += diff
                    pv.L[pv.RECOMP[bipar_id].M[pv.R1[bipar_id]]] -= diff
                    break      
                        
        
          
        
        
        
def greedy_pruning(pv: PruningVars):
    from copy import deepcopy
    old_L = deepcopy(pv.L)
    while True:
        old_L = deepcopy(pv.L)
        best_index = np.argmax(pv.L)
        # print(f"Pruning {best_index} with loss reduction: {pv.L[best_index]}. Loss of 77: {pv.L[77]}.")
        if (pv.L[best_index] <= 0): break # No edges induce loss reduction. Exiting the while loop
        # Otherwise, prune best_index
        pv.E.remove(best_index)
        pv.pruned.add(best_index)
        renew_match(pv, best_index)
    
    

# def greedy_pruning(FP_loss, FN_diff, count_arr, normalized_scores, match_arr, topo_depth_list, K, reverse_id_ref, reverse_id_alt, bits, n_taxa):
#     print("pruning process started")
#     benefit = FP_loss + FN_diff
#     nb_init_edges = len(FP_loss)
#     nb_bipars = len(count_arr)
    
#     which_match = np.zeros(nb_bipars, dtype=int) # 何番目のmatchと今matchしているか。
#     which_second_match = np.ones(nb_bipars, dtype=int) # 何番目に良いものがsecond matchになっているか。
    
#     reverse_match = [list() for _ in range(nb_init_edges)]
#     reverse_second_match = [list() for _ in range(nb_init_edges)]
#     for i in range(nb_bipars):
#         reverse_match[match_arr[i,0]].append(i)
#         reverse_second_match[match_arr[i,1]].append(i)
    
    
#     all_match_dict = dict()
#     all_normalized_score = dict()
#     use_K = [True for _ in range(nb_bipars)]
    
    
#     pruned_set = set(); current_set = set([i for i in range(nb_init_edges)])
#     for i in trange(nb_init_edges):
#         max_ind = np.argmax(benefit)
#         if benefit[max_ind] < 0:
#             break # break from loop, end pruning
        
#         # prune
#         benefit[max_ind] = -1 # setting benefit of pruned edge to be -1. Might be inefficient.
#         pruned_set.add(max_ind)
#         current_set.remove(max_ind)
        
#         # For bipartitions in reverse_match[max_ind]
#         for bipar in reverse_match[max_ind]:
#             which_match[bipar]  = which_second_match[bipar]
#             if use_K[bipar]:
#                 if match_arr[bipar, which_match[bipar]] == -1:
#                     continue
#                 # move second match to first match
#                 reverse_match[match_arr[bipar, which_match[bipar]]].append(bipar)
#             else:
#                 if all_match_dict[bipar][which_match[bipar]] == -1:
#                     continue
#                 # move second match to first match
#                 reverse_match[ all_match_dict[bipar][which_match[bipar]] ].append(bipar)
            
#             # search new second match  
#             now_matched = which_match[bipar] + 1
#             while True:
#                 if use_K[bipar]:
#                     if now_matched > K-1:
#                         # compute everything
#                         matches, scores = compute_all_normalized_td(bipar, match_arr[bipar], K, reverse_id_ref, reverse_id_alt[bipar], topo_depth_list[bipar], bits, n_taxa)
#                         all_match_dict[bipar] = np.append(match_arr[bipar], matches)
#                         all_normalized_score[bipar] = np.append(normalized_scores[bipar], scores)
#                         use_K[bipar] = False
#                     elif match_arr[bipar, now_matched] in pruned_set:
#                         now_matched += 1
#                     else:
#                         # match_arr[bipar, now_matched] is a valid second match.
#                         if match_arr[bipar, now_matched] != -1:
#                             reverse_second_match[match_arr[bipar, now_matched]].append(bipar)
#                         which_second_match[bipar] = now_matched
#                         benefit[which_match[bipar]] += normalized_scores[bipar, which_match[bipar]] - normalized_scores[bipar, which_second_match[bipar]]
#                         break
#                 else:
#                     if all_match_dict[bipar][now_matched] in pruned_set:
#                         now_matched += 1
#                     else:
#                         if all_match_dict[bipar][now_matched] != -1:
#                             reverse_second_match[match_arr[bipar, now_matched]].append(bipar)
#                         which_second_match[bipar] = now_matched
#                         benefit[which_match[bipar]] += all_normalized_score[bipar][which_match[bipar]] - all_normalized_score[bipar][which_second_match[bipar]]
#                         break
        
#         # For bipartitions in reverse_second_match[max_ind]     
#         for bipar in reverse_second_match[max_ind]:
#             if use_K[bipar]:
#                 old_second_score = normalized_scores[bipar, which_second_match[bipar]]
#             else:
#                 old_second_score = all_normalized_score[bipar][which_second_match[bipar]]
#             # search new second match  
#             now_matched = which_second_match[bipar] + 1
#             while True:
#                 if use_K[bipar]:
#                     if now_matched > K-1:
#                         # compute everything
#                         matches, scores = compute_all_normalized_td(bipar, match_arr[bipar], K, reverse_id_ref, reverse_id_alt[bipar], topo_depth_list[bipar], bits, n_taxa)
#                         all_match_dict[bipar] = np.append(match_arr[bipar], matches)
#                         all_normalized_score[bipar] = np.append(normalized_scores[bipar], scores)
#                         use_K[bipar] = False
#                     elif match_arr[bipar, now_matched] in pruned_set:
#                         now_matched += 1
#                     else:
#                         # match_arr[bipar, now_matched] is a valid second match.
#                         if match_arr[bipar, now_matched] != -1:
#                             reverse_second_match[match_arr[bipar, now_matched]].append(bipar)
#                         which_second_match[bipar] = now_matched
#                         benefit[match_arr[bipar, which_match[bipar]]] += old_second_score - normalized_scores[bipar, which_second_match[bipar]]
#                         break
#                 else:
#                     if all_match_dict[bipar][now_matched] in pruned_set:
#                         now_matched += 1
#                     else:
#                         which_second_match[bipar] = now_matched
#                         if all_match_dict[bipar][now_matched] != -1:
#                             reverse_second_match[all_match_dict[bipar][now_matched]].append(bipar)
#                         benefit[all_match_dict[bipar][which_match[bipar]]] += old_second_score - all_normalized_score[bipar][which_second_match[bipar]]
#                         break        
        
#     return current_set



def _prune_first_comp(res: OrganizedTransferResult, booster_lib: ctypes.CDLL):
    # False Positive Loss
    FP_loss = (1-np.array(res.TS)) * res.num_alt_trees
    
    # False Negative Losses and W1, W2
    W1 = [set() for _ in range(len(res.TS))]
    W2 = [set() for _ in range(len(res.TS))]
    FN_diff = np.zeros_like(FP_loss)
    for i in range(res.S.shape[0]):
        if res.M[i,0]!= -1:
            FN_diff[res.M[i,0]] += (res.S[i,1] - res.S[i,0]) * res.C[i]
            W1[res.M[i,0]].add(i)
            if res.M[i,1]!=-1:
                W2[res.M[i,1]].add(i)   
    # Initialize Other variables
    E = [i for i in range(len(res.TS))]
    R1 = np.zeros(res.S.shape[0], dtype=int)
    R2 = np.full(res.S.shape[0], 1, dtype=int)
    
    return PruningVars(FP_loss, FN_diff, E, R1, R2, W1, W2, res, booster_lib)
           
    
    
    # print(np.array(match_list))
    
    
def split_set_to_bipar_encodings(bipar_set, n_taxa, bits):
    split_length = (n_taxa-1) // bits + 1
    last_bits = n_taxa % bits
    
    # first take unique list of all unsigned integers to convert.
    all_integers_pre = [item for integer_list in bipar_set for item in integer_list[:split_length]]
    all_integers_pre = np.unique(all_integers_pre)
    
    all_integers_last = [integer_list[-1] for integer_list in bipar_set]
    all_integers_last = np.unique(all_integers_last)
    
    # split to bipar
    pre_bits = [''.join('1' if item & (1 << i) else '0' for i in range(bits)) for item in all_integers_pre]
    pre_dict = dict(zip(all_integers_pre, pre_bits))
    last_bits = [''.join('1' if item & (1 << i) else '0' for i in range(last_bits)) for item in all_integers_last]
    last_dict = dict(zip(all_integers_last, last_bits))
    
    bipartition_list = []
    for bipar in bipar_set:
        bipartition_bin = ''.join([pre_dict[item] for item in bipar[:split_length-1]] + [last_dict[bipar[split_length-1]]])
        bipartition_list.append(dendropy.Bipartition(leafset_bitmask = Bits(bin=bipartition_bin).uint, tree_leafset_bitmask=2**n_taxa-1))
    return bipartition_list
 
def report_memory(process, print_string):
    memory_info = process.memory_info() 
    rss_memory = memory_info.rss  # Resident Set Size: Memory currently in use by the process (in bytes)
    print(f"{print_string}: {rss_memory / (1024 ** 3):.2f} GB")

def c_prune(inittree_file: str, inputtrees_file: str, K=30):
    # tracemalloc.start()
    # event = threading.Event()
    # initial_time = time.time()
    # m = threading.ThreaTBEd(target=monitor_cpu,args=((initial_time,event)))
    # m.start() # 開始
    # Get the current process using psutil
    process = psutil.Process(os.getpid())
    report_memory(process, "memory used before anything")
    
    booster_lib = load_booster()
    booster_lib = prepare_tbe_support_and_match(booster_lib)
    booster_lib = prepare_recompute(booster_lib)
    booster_lib = prepare_prune_and_return_newick(booster_lib)
    
    
    with open(inputtrees_file, "r") as f:
        inputtrees = f.read()
    with open(inittree_file, "r") as f:
        inittree = f.read()
    args = prepare_tbe_support_and_match_args(inittree, inputtrees, K)
    first_done = 0
    final_tree = None
    try:
        # Step 1: tbe_support_and_match
        start = time.perf_counter()
        res_ptr = booster_lib.tbe_support_and_match(*args)
        middle = time.perf_counter()
        print(f"tbe support and match time with K={K}: {middle-start}")    
        first_done += 1

        # copying results: organize_support_and_match
        organized_results: OrganizedTransferResult = organize_support_and_match(res_ptr)
        end = time.perf_counter()
        print(f"organizing time with K={K}: {end - middle}")
        
        # Step2: initial gain and initializations
        pruning_vars: PruningVars = _prune_first_comp(organized_results, booster_lib)
        first_comp_end = time.perf_counter()
        print(f"First computation ended in {first_comp_end - end}")
        
        
        print(f"Initial number of branches {len(pruning_vars.E)}")
        
        # Step3: Greedy Pruning
        greedy_pruning(pruning_vars)
        pruning_end = time.perf_counter()
        print(f"Pruning ended in {pruning_end - first_comp_end}")
        
        # Show Final Result: print the tree to stderr
        loc_to_id = dict()
        for id, loc in pruning_vars.OSR.REF_BIPAR_ID_to_location.items():
            loc_to_id[loc] = id
        id_list = [loc_to_id[item] for item in pruning_vars.pruned]
        args = prepare_prune_and_return_newick_args(id_list)
        newick_string = booster_lib.prune_and_return_newick(*args)
        newick_string_decoded = ctypes.cast(newick_string, ctypes.c_char_p).value.decode('utf-8')
        print(f"Newick string at last: {newick_string_decoded}")
        
        booster_lib.free_buffer(newick_string)
        
        print(f"Number of remaining edges: {len(pruning_vars.E)}")
        
    finally:
        booster_lib.free_tbe_support_and_match(res_ptr)
        booster_lib.free_TDA_LIST()
        
        # event.set()

    
    
    return final_tree