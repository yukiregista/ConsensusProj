from ._booster import load_booster, prepare_tbe_match, prepare_tbe_match_args, match_score_list, prepare_tbe_support_args, prepare_tbe_support, tbe_support, _generate_bitmask
import dendropy
import numpy as np
import time
from tqdm import trange
from bitstring import Bits

import threading, psutil, os
from ._dev_utils import monitor_cpu

import gc
import tracemalloc

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
    
    

def greedy_pruning(FP_loss, FN_diff, count_arr, normalized_scores, match_arr, topo_depth_list, K, reverse_id_ref, reverse_id_alt, bits, n_taxa):
    print("pruning process started")
    benefit = FP_loss + FN_diff
    nb_init_edges = len(FP_loss)
    nb_bipars = len(count_arr)
    
    which_match = np.zeros(nb_bipars, dtype=int) # 何番目のmatchと今matchしているか。
    which_second_match = np.ones(nb_bipars, dtype=int) # 何番目に良いものがsecond matchになっているか。
    
    reverse_match = [list() for _ in range(nb_init_edges)]
    reverse_second_match = [list() for _ in range(nb_init_edges)]
    for i in range(nb_bipars):
        reverse_match[match_arr[i,0]].append(i)
        reverse_second_match[match_arr[i,1]].append(i)
    
    
    all_match_dict = dict()
    all_normalized_score = dict()
    use_K = [True for _ in range(nb_bipars)]
    
    
    pruned_set = set(); current_set = set([i for i in range(nb_init_edges)])
    for i in trange(nb_init_edges):
        max_ind = np.argmax(benefit)
        if benefit[max_ind] < 0:
            break # break from loop, end pruning
        
        # prune
        benefit[max_ind] = -1 # setting benefit of pruned edge to be -1. Might be inefficient.
        pruned_set.add(max_ind)
        current_set.remove(max_ind)
        
        # For bipartitions in reverse_match[max_ind]
        for bipar in reverse_match[max_ind]:
            which_match[bipar] += 1
            if use_K[bipar]:
                if match_arr[bipar, which_match[bipar]] == -1:
                    continue
                # move second match to first match
                reverse_match[match_arr[bipar, which_match[bipar]]].append(bipar)
            else:
                if all_match_dict[bipar][which_match[bipar]] == -1:
                    continue
                # move second match to first match
                reverse_match[ all_match_dict[bipar][which_match[bipar]] ].append(bipar)
            
            # search new second match  
            now_matched = which_match[bipar] + 1
            while True:
                if use_K[bipar]:
                    if now_matched > K-1:
                        # compute everything
                        matches, scores = compute_all_normalized_td(bipar, match_arr[bipar], K, reverse_id_ref, reverse_id_alt[bipar], topo_depth_list[bipar], bits, n_taxa)
                        all_match_dict[bipar] = np.append(match_arr[bipar], matches)
                        all_normalized_score[bipar] = np.append(normalized_scores[bipar], scores)
                        use_K[bipar] = False
                    elif match_arr[bipar, now_matched] in pruned_set:
                        now_matched += 1
                    else:
                        # match_arr[bipar, now_matched] is a valid second match.
                        if match_arr[bipar, now_matched] != -1:
                            reverse_second_match[match_arr[bipar, now_matched]].append(bipar)
                        which_second_match[bipar] = now_matched
                        benefit[which_match[bipar]] += normalized_scores[bipar, which_match[bipar]] - normalized_scores[bipar, which_second_match[bipar]]
                        break
                else:
                    if all_match_dict[bipar][now_matched] in pruned_set:
                        now_matched += 1
                    else:
                        if all_match_dict[bipar][now_matched] != -1:
                            reverse_second_match[match_arr[bipar, now_matched]].append(bipar)
                        which_second_match[bipar] = now_matched
                        benefit[which_match[bipar]] += all_normalized_score[bipar][which_match[bipar]] - all_normalized_score[bipar][which_second_match[bipar]]
                        break
        
        # For bipartitions in reverse_second_match[max_ind]     
        for bipar in reverse_second_match[max_ind]:
            if use_K[bipar]:
                old_second_score = normalized_scores[bipar, which_second_match[bipar]]
            else:
                old_second_score = all_normalized_score[bipar][which_second_match[bipar]]
            # search new second match  
            now_matched = which_second_match[bipar] + 1
            while True:
                if use_K[bipar]:
                    if now_matched > K-1:
                        # compute everything
                        matches, scores = compute_all_normalized_td(bipar, match_arr[bipar], K, reverse_id_ref, reverse_id_alt[bipar], topo_depth_list[bipar], bits, n_taxa)
                        all_match_dict[bipar] = np.append(match_arr[bipar], matches)
                        all_normalized_score[bipar] = np.append(normalized_scores[bipar], scores)
                        use_K[bipar] = False
                    elif match_arr[bipar, now_matched] in pruned_set:
                        now_matched += 1
                    else:
                        # match_arr[bipar, now_matched] is a valid second match.
                        if match_arr[bipar, now_matched] != -1:
                            reverse_second_match[match_arr[bipar, now_matched]].append(bipar)
                        which_second_match[bipar] = now_matched
                        benefit[match_arr[bipar, which_match[bipar]]] += old_second_score - normalized_scores[bipar, which_second_match[bipar]]
                        break
                else:
                    if all_match_dict[bipar][now_matched] in pruned_set:
                        now_matched += 1
                    else:
                        which_second_match[bipar] = now_matched
                        if all_match_dict[bipar][now_matched] != -1:
                            reverse_second_match[match_arr[bipar, now_matched]].append(bipar)
                        benefit[match_arr[bipar, which_match[bipar]]] += old_second_score - all_normalized_score[bipar][which_second_match[bipar]]
                        break        
        
    return current_set



def _prune_first_comp(support_list, count_list, match_list, score_list, topo_depth_list, n_trees: int):
    # False Positive Loss
    FP_loss = (1-np.array(support_list)) * n_trees
    
    # False Negative Losses
    FN_diff = np.zeros_like(FP_loss)
    count_arr = np.array(count_list)
    score_arr = np.array(score_list)
    topo_depth_arr = np.array(topo_depth_list)
    normalized_scores = score_arr / (topo_depth_arr.reshape(-1,1)-1)
    match_arr = np.array(match_list)
    
    for i in range(normalized_scores.shape[0]):
        if match_arr[i,0]!= -1:
            FN_diff[match_arr[i,0]] += (normalized_scores[i,0] - normalized_scores[i,1]) * count_list[i]

    return FP_loss, FN_diff, count_arr, normalized_scores, match_arr
           
    
    
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
    tracemalloc.start()
    # event = threading.Event()
    # initial_time = time.time()
    # m = threading.Thread(target=monitor_cpu,args=((initial_time,event)))
    # m.start() # 開始
    # Get the current process using psutil
    process = psutil.Process(os.getpid())
    report_memory(process, "memory used before anything")
    
    booster_lib = load_booster()
    booster_lib = prepare_tbe_match(booster_lib)
    booster_lib = prepare_tbe_support(booster_lib)
    with open(inputtrees_file, "r") as f:
        inputtrees = f.read()
    with open(inittree_file, "r") as f:
        inittree = f.read()
    args = prepare_tbe_match_args(inputtrees, inittree, K)
    args2 = prepare_tbe_support_args(inittree, inputtrees)
    first_done = 0
    final_tree = None
    try:
        # the other way: tbe_match
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)
        
        report_memory(process, "memory used before tbe_match")
        start = time.perf_counter()
        res_ptr = booster_lib.tbe_match(*args)
        middle = time.perf_counter()
        print(f"tbe match time with K={K}: {middle-start}")
        first_done += 1
        
        report_memory(process, "memory used after tbe_match")
        
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)
        
        # the normal way: tbe_support
        res_ptr2 = booster_lib.tbe_support(*args2)
        id_dict_ref, support_list, bits, n_taxa = tbe_support(res_ptr2)
        end = time.perf_counter()
        print(f"tbe support time: {end-middle}")
        
        report_memory(process, "memory used after tbe_support")
        
        # copying results: match_score_list
        reverse_id_ref = {val:key for key, val in id_dict_ref.items()}
        id_dict_alt, count_list, match_list, score_list, topo_depth_list, num_trees = match_score_list(res_ptr, id_dict_ref, K)
        reverse_id_alt = {val:key for key, val in id_dict_alt.items()}
        copy_end = time.perf_counter()
        print(f"copying time with K={K}: {copy_end-end}")
        
        report_memory(process, "memory used after copying into Python object")
        
        
        # compute the initial gain
        FP_loss, FN_diff, count_arr, normalized_scores, match_arr = _prune_first_comp(support_list, count_list, match_list, score_list, topo_depth_list, num_trees)
        initial_gain_end = time.perf_counter()
        print(f"initial gain computation time with K={K}: {initial_gain_end-copy_end}")
        
        # greedy pruning
        edge_id_set = greedy_pruning(FP_loss, FN_diff, count_arr, normalized_scores, match_arr, topo_depth_list, K, reverse_id_ref, reverse_id_alt, bits, n_taxa)
        pruning_end = time.perf_counter()
        print(f"pruning time with K={K}: {pruning_end - initial_gain_end}")
        
        split_set = [reverse_id_ref[item] for item in edge_id_set]
        bipartition_encodings = split_set_to_bipar_encodings(split_set, n_taxa, bits)
        bipar_en_end = time.perf_counter()
        print(f"Bipartition encoding time: {bipar_en_end - pruning_end}")
        
        # load initial tree
        init_tree = dendropy.Tree.get(path = inittree_file, schema="newick")
        init_load_end = time.perf_counter()
        print(f"Load init tree time: {init_load_end - bipar_en_end}")
        
        taxon_labels = [item for item in init_tree.taxon_namespace]
        taxon_namespace = dendropy.TaxonNamespace(np.sort(taxon_labels))
        
        final_tree = dendropy.Tree.from_bipartition_encoding(bipartition_encodings, taxon_namespace)
        final_tree_end = time.perf_counter()
        print(f"Final Tree Creation time: {final_tree_end - init_load_end}")
        
        report_memory(process, "memory used after everything")
        
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)
        
        
    finally:
        booster_lib.after_tbe_match(res_ptr)
        if first_done:
            booster_lib.after_tbe_support(res_ptr2)
        
        mem = psutil.virtual_memory() 
        print(mem)
        
        report_memory(process, "memory used after freeing C memory")
        
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)
        
        
        # event.set()

    
    
    return final_tree