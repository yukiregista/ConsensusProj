from ._consensus import TreeList_with_support, Tree_with_support, _compatible
from bitstring import Bits
import numpy as np
from itertools import combinations
import distance
import sys
from collections import OrderedDict
import time
import dendropy


def _MinHammingDist(a : Bits, b : Bits):
    return min(distance.hamming(a.bin,b.bin), distance.hamming((~a).bin,b.bin))

class STDGreedyConsensus():
    
    def _CreateBipartitionListAndCounts(self):
        # create BipartitionDict
        BipartitionDict = OrderedDict()
        for tree in self.input_trees:
            dict_keys = BipartitionDict.keys() # no duplicate keys in one tree
            tree.encode_bipartitions()
            for branch in tree.internal_edges(exclude_seed_edge=True):
                key = branch.bipartition.split_as_int()
                if key in dict_keys:
                    BipartitionDict[key] += 1
                else:
                    BipartitionDict[key] = 1
        # Create BipartitionList, BipartitionCounts, Index2Bipartition and BipartitionP
        BipartitionList = []; BipartitionCounts = []; BipartitionP = []; BipartitionBitsList = []
        
        for key, value in BipartitionDict.items():
            BipartitionList.append(key)
            BipartitionCounts.append(value)
            bitstr = Bits(uint = key, length=self.n_taxa)
            BipartitionBitsList.append(bitstr)
            BipartitionP.append(min(bitstr.count(0), bitstr.count(1))) # size of the light side
        
        Bipartition2Index = {BipartitionList[i]:i for i in range(len(BipartitionList))}
        return np.array(BipartitionList), np.array(BipartitionCounts), BipartitionBitsList, np.array(BipartitionP), Bipartition2Index
    
    
    # def _CreateCompatibilityGraph(self):
    #     # intiialize CompatibilityGraph
    #     CompatibilityGraph = dict(zip([i for i in range(self.n_bipartitions)], [set() for _ in range(self.n_bipartitions)]))
    #     for key in CompatibilityGraph.keys():
    #         CompatibilityGraph[key].add(key)
        
    #     for i, j in combinations([i for i in range(self.n_bipartitions)], 2):
    #         if not _compatible(self.BipartitionBitsList[i], self.BipartitionBitsList[j]):
    #             CompatibilityGraph[i].add(j)
    #             CompatibilityGraph[j].add(i)
    #     return CompatibilityGraph
           
    def _ComputeDIST(self, method="simple"):
        DIST = np.zeros((self.n_bipartitions+1, self.n_bipartitions+1))
        if method == "simple":
            for i, j in combinations([i for i in range(self.n_bipartitions)], 2):
                bit_i = self.BipartitionBitsList[i]
                bit_j = self.BipartitionBitsList[j]
                DIST[i,j] = DIST[j,i] = _MinHammingDist(bit_i, bit_j)
        else:
            sys.exit("Invalid method.")
        DIST[self.n_bipartitions, :-1] = self.BipartitionP - 1
        DIST[:-1, self.n_bipartitions] = self.BipartitionP - 1
        return DIST

    
    def _ComputeTransferDissimilarityCost(self):
        TD = np.zeros(self.n_bipartitions)
        for tree in self.input_trees:
            internal_bipartitions = [item.bipartition.split_as_int() for item in tree.internal_edges(exclude_seed_edge=True)]
            indices = [self.Bipartition2Index[item] for item in internal_bipartitions]
            for bipar in [i for i in range(self.n_bipartitions)]:
                minDist = np.min(self.DIST[bipar, indices])
                TD[bipar] += min (minDist / (self.BipartitionP[bipar]-1), 1)
        return TD
            
    
    def __init__(self, input_trees : TreeList_with_support):
        self.input_trees = input_trees
        self.taxon_namespace = input_trees.taxon_namespace
        self.n_taxa = len(input_trees.taxon_namespace)
        
        # Create BipartitionList etc
        print("Creating Bipartition List etc...", end=" ", flush=True)
        time1 = time.time()
        self.BipartitionList, self.BipartitionCounts, self.BipartitionBitsList, self.BipartitionP, self.Bipartition2Index \
            = self._CreateBipartitionListAndCounts()
        time2 = time.time()
        print(f"time passed : {time2-time1}", flush=True)
        
        
        self.n_bipartitions = len(self.BipartitionList)
        
        # # Create CompatibilityGraph
        # print("Creating CompatibilityGraph...", end = " ", flush=True)
        # self.CompatibilityGraph = self._CreateCompatibilityGraph()
        # time3 = time.time()
        # print(f"time passed : {time3-time2}", flush=True) -> probably very inefficient, no need for IncompatibilityGraph.
        
        # Create DIST
        print("Creating DIST...", end = " ", flush=True)
        self.DIST = self._ComputeDIST(method="simple")
        time4 = time.time()
        print(f"time passed : {time4-time2}", flush=True)
        
        # Compute Transfer Dissimilarity Cost
        print("Creating TD...", end = " ", flush=True)
        self.TD = self._ComputeTransferDissimilarityCost()
        time5 = time.time()
        print(f"time passed : {time5-time4}", flush=True)
        
        
        # Initialize self.MATCH and self.SECOND_MATCH
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        # self.n_bipartitions is used to represent the match with leaf bipartitions.
        
        # Initialize current tree
        self.current_tree_included = np.zeros(self.n_bipartitions)
        
    def specify_initial_tree(self, initial_tree : Tree_with_support):
        # initial tree should only have those bipartitions in the input trees
        assert (initial_tree.taxon_namespace == self.taxon_namespace)
        internal_branches = initial_tree.internal_edges(exclude_seed_edge = True)
        internal_bipar_keys = [item.bipartition.split_as_int() for item in internal_branches]
        if not np.all([item in self.BipartitionList for item in internal_bipar_keys]):
            sys.exit("Initial tree should only have those bipartitions in the input trees")
        
        # renew self.current_tree_included
        self.current_tree_included = np.zeros(self.n_bipartitions)
        internal_indices = np.array([self.Bipartition2Index[item] for item in internal_bipar_keys])
        for item in internal_indices:
            self.current_tree_included[item] = 1
            
        # renew MATCH and SECOND_MATCH
        self.MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        self.SECOND_MATCH = np.array([self.n_bipartitions for _ in range(self.n_bipartitions)])
        
        if len(internal_indices) == 0:
            return
        
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            min_ind = internal_indices[np.argmin(self.DIST[bipar_ind, internal_indices])]
            if (self.DIST[bipar_ind, min_ind]) / (self.BipartitionP[bipar_ind]-1) >= 1:
                continue
            self.MATCH[bipar_ind] = min_ind
            candidate_indices = internal_indices[internal_indices != min_ind]
            if len(candidate_indices) == 0:
                break
            min2_ind = candidate_indices[np.argmin(self.DIST[bipar_ind, candidate_indices])]
            if (self.DIST[bipar_ind, min2_ind]) / (self.BipartitionP[bipar_ind]-1) >= 1:
                continue
            self.SECOND_MATCH[bipar_ind] = min2_ind
        
        
        
    
    def is_compatible(self, BiparBit):
        
        included_indices = np.nonzero(self.current_tree_included)[0]
        
        if len(included_indices) == 0:
            return True
        elif len(included_indices) == self.n_taxa - 3:
            return False
        else:
            # check compatibility with each edge.
            compatible = True
            for index in included_indices:
                current = self.BipartitionBitsList[index]
                if not _compatible(BiparBit, current):
                    return False
        return compatible
                
    def _remove_benefit(self, branch_index):
        benefit_fp = self.TD[branch_index]
        # compute benefit_fn
        benefit_fn = 0
        branch = self.BipartitionList[branch_index]
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            # check if match
            if self.MATCH[bipar_ind] == branch_index:
                benefit_fn += ( (self.DIST[bipar_ind, branch_index] - self.DIST[bipar_ind,self.SECOND_MATCH[bipar_ind]]) / (self.BipartitionP[bipar_ind]-1) ) * self.BipartitionCounts[bipar_ind]
        
        return benefit_fp + benefit_fn 

    
    def _add_benefit(self, branch_index):
        benefit_fp = - self.TD[branch_index]
        # compute benefit_fn
        benefit_fn = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            dist_diff = self.DIST[self.MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
            if dist_diff > 0:
                # risk reduction occurs for this bipartition
                benefit_fn += dist_diff / (self.BipartitionP[bipar_ind] - 1) * self.BipartitionCounts[bipar_ind]
        
        return benefit_fn + benefit_fp

    
    def benefit(self, branch_index, method = "add"):
        if method == "add":
            return self._add_benefit(branch_index)
        elif method == "remove":
            return self._remove_benefit(branch_index)
        else:
            sys.exit("Invalid method specified.")
    
    
    def add(self, branch_index):
        self.current_tree_included[branch_index] = 1
        # renew match and second match
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            first_dist_diff = self.DIST[self.MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
            if first_dist_diff > 0:
                # renew match and second match
                self.SECOND_MATCH[bipar_ind] = self.MATCH[bipar_ind]
                self.MATCH[bipar_ind] = branch_index
            else:
                second_dist_diff = self.DIST[self.SECOND_MATCH[bipar_ind], bipar_ind] - self.DIST[branch_index, bipar_ind]
                if second_dist_diff > 0:
                    # renew second match
                    self.SECOND_MATCH[bipar_ind] = branch_index
        
    def remove(self, branch_index):
        self.current_tree_included[branch_index] = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            # check if match
            if self.MATCH[bipar_ind] == branch_index:
                # renew match and search for the new second match
                self.MATCH[bipar_ind] = self.SECOND_MATCH[bipar_ind]
                ## search for second match
                if self.MATCH[bipar_ind] == self.n_bipartitions:
                    self.SECOND_MATCH[bipar_ind] == self.n_bipartitions # equivalent to doing nothing
                else:
                    current_indices = np.nonzero(self.current_tree_included)[0]
                    candidate_indices = current_indices[current_indices!=self.MATCH[bipar_ind]]
                    cand = candidate_indices[np.argmin(self.DIST[bipar_ind,candidate_indices])]
                    if self.DIST[bipar_ind, cand]/(self.BipartitionP[bipar_ind] - 1) >= 1:
                        self.SECOND_MATCH[bipar_ind] = self.n_bipartitions
                    else:
                        self.SECOND_MATCH[bipar_ind] = cand
            if self.SECOND_MATCH[bipar_ind] == branch_index:
                # just renew the second match
                current_indices = np.nonzero(self.current_tree_included)[0]
                candidate_indices = current_indices[current_indices!=self.MATCH[bipar_ind]]
                cand = candidate_indices[np.argmin(self.DIST[bipar_ind,candidate_indices])]
                if self.DIST[bipar_ind, cand]/(self.BipartitionP[bipar_ind] - 1) >= 1:
                    self.SECOND_MATCH[bipar_ind] = self.n_bipartitions
                else:
                    self.SECOND_MATCH[bipar_ind] = cand

    def _first_greedy(self, sorted_index):
        improving = True
        while improving:
            improving = False
            for bipar_ind in sorted_index:
                if self.current_tree_included[bipar_ind] == 1:
                    if self.benefit(bipar_ind, method="remove") > 0:
                        self.remove(bipar_ind)
                        print("branch removed")
                        improving =True
                        break
                elif self.is_compatible(self.BipartitionBitsList[bipar_ind]):
                    if self.benefit(bipar_ind, method="add") > 0:
                        self.add(bipar_ind)
                        print("branch added")
                        improving=True
                        break
        
    def _risk_greedy(self, sorted_index):
        while True:
            benefit = np.zeros(self.n_bipartitions)
            for bipar_ind in sorted_index:
                if self.current_tree_included[bipar_ind] == 1:
                    benefit[bipar_ind] = self.benefit(bipar_ind, method="remove")
                elif self.is_compatible(self.BipartitionBitsList[bipar_ind]):
                    benefit[bipar_ind] = self.benefit(bipar_ind, method="add")
            most_benefit_arg = np.argmax(benefit)
            if benefit[most_benefit_arg] > 0:
                if self.current_tree_included[most_benefit_arg] == 1:
                    self.remove(most_benefit_arg)
                    print("branch removed with benefit of ", benefit[most_benefit_arg])
                else:
                    self.add(most_benefit_arg)
                    print("branch added with benefit of ", benefit[most_benefit_arg])
            else:
                break
        
    def greedy(self, method = "first", order="BS"):
        
        if order == "BS":
            # order branch by branch support
            sorted_index = np.argsort(self.BipartitionCounts)[::-1]
            if method == "first":
                self._first_greedy(sorted_index)
            if method == "most":
                self._risk_greedy(sorted_index)
            
        
    def return_current_tree(self):
        bipartition_keys = self.BipartitionList[np.nonzero(self.current_tree_included)[0]]
        bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.n_taxa - 1) for item in bipartition_keys]
        dendropy_tree = dendropy.Tree.from_bipartition_encoding(bipartitions, self.taxon_namespace)
        return Tree_with_support(dendropy_tree)  
    
    def current_loss(self):
        nonzero = np.nonzero(self.current_tree_included)[0]
        fp = np.sum(self.TD[nonzero])
        fn = 0
        for bipar_ind, bipar in enumerate(self.BipartitionList):
            fn += ((self.DIST[bipar_ind, self.MATCH[bipar_ind]]) / (self.BipartitionP[bipar_ind] - 1)) * self.BipartitionCounts[bipar_ind]
        return (fp + fn) / len(self.input_trees)
            
            
                
        
            
                   
            
        
        
        
        
        
        
        
        
    

