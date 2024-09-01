import pytest
import Consensus
from importlib.resources import files
import dendropy
import ctypes
from bitstring import Bits

import numpy as np

class Split(ctypes.Structure):
    _fields_ = [("data", ctypes.POINTER(ctypes.c_uint))]

class BipartitionMatches(ctypes.Structure):
    _fields_ = [
        ("bipartition", ctypes.POINTER(ctypes.c_uint)),
        ("matches", ctypes.POINTER(ctypes.POINTER(ctypes.c_uint))),
        ("td", ctypes.POINTER(ctypes.c_int)),
        ("is_external", ctypes.c_bool)
    ]

class BipartitionDict(ctypes.Structure):
    _fields_ = [
        ("entries", ctypes.POINTER(ctypes.POINTER(BipartitionMatches))),
        ("num_entries", ctypes.c_int),
        ("bipartition_size", ctypes.c_size_t),
        ("num_matches", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
        ("uint_bits", ctypes.c_int)
    ]


def int_to_reversed_bin(x, bit_length):
    '''
    This conversion is inefficient!: 
    It would be better if we also use the similar structure (array of uints)
    in Python as well.
    '''
    return ''.join('1' if x & (1 << i) else '0' for i in range(bit_length))


def to_bitvector(bipartition: ctypes.POINTER(ctypes.c_uint), n_tips: int=100, bit_length: int=32) -> Bits:
    n_uints = (n_tips-1) // bit_length + 1
    bipartition_array = [bipartition[i] for i in range(n_uints)]
    bit_size = [bit_length for _ in range(n_uints-1)]
    bit_size.append(n_tips - bit_length * (n_uints-1))
    bit_length = 32
    binary_strings = [int_to_reversed_bin(x, y) for x, y in zip(bipartition_array, bit_size)]
    concatenated_bits = ''.join(binary_strings)
    return Bits(bin=concatenated_bits[::-1])


def test_first_K_match(K=20):
    
    # example trees
    print(files('Consensus.example_data').joinpath('sample1.nex'))
    tree1 = Consensus.Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.nex'), schema="nexus")
    tree2 = Consensus.Tree_with_support.get(path = files('Consensus.example_data').joinpath('astral_GTRgamma.tre'), schema="newick",
                                            taxon_namespace = tree1.taxon_namespace)
    
    
    
    ## For using booster
    tree1_str = tree1.as_string("newick", suppress_rooting=True)
    tree2_str = tree2.as_string("newick", suppress_rooting=True)
    print(tree2_str)
    
    booster_lib = Consensus.load_booster()
    booster_lib.tbe_match.restype = ctypes.c_void_p # To get pointer accurately
    booster_lib.print_bdict.argtypes = ctypes.c_void_p,
    booster_lib.after_tbe_match.argtypes = ctypes.c_void_p,
    
    tbe_match_args = ['program_name', '-k', f'{K}']
    argc = len(tbe_match_args)
    argv_ctypes = (ctypes.POINTER(ctypes.c_char) * argc)()
    for i, str in enumerate(tbe_match_args):
        enc_str = str.encode('utf-8')
        argv_ctypes[i] = ctypes.create_string_buffer(enc_str)
    
    tree1_ctypes = ctypes.create_string_buffer(tree1_str.encode('utf-8'))
    tree2_ctypes = (ctypes.POINTER(ctypes.c_char) * 1)()
    tree2_ctypes[0] = ctypes.create_string_buffer(tree2_str.encode('utf-8'))
    
    try:
        res = booster_lib.tbe_match(argc, argv_ctypes, tree1_ctypes ,tree2_ctypes)
        print(booster_lib.tbe_match.restype)
        import time
        time.sleep(5)
        print(hex(res))
        booster_lib.after_tbe_match(res)
        print("here")
        print("loading")
        print(tbe_match_args)
        argc = len(tbe_match_args)
        argv_ctypes = (ctypes.POINTER(ctypes.c_char) * argc)()
        for i, str in enumerate(tbe_match_args):
            enc_str = str.encode('utf-8')
            argv_ctypes[i] = ctypes.create_string_buffer(enc_str)
        
        
        tree1_ctypes = ctypes.create_string_buffer(tree1_str.encode('utf-8'))
        tree2_ctypes = (ctypes.POINTER(ctypes.c_char) * 1)()
        tree2_ctypes[0] = ctypes.create_string_buffer(tree2_str.encode('utf-8'))
        
        print(argc,argv_ctypes, tree1_ctypes, tree2_ctypes)
        print(booster_lib.tbe_match.restype)
        res = booster_lib.tbe_match(argc, argv_ctypes, tree1_ctypes ,tree2_ctypes)
        if not res:
            raise RuntimeError("NULL pointer received from tbe_match")
        print(hex(res))
        import time
        time.sleep(5)
        # print(hex(res))
        
        print("here")
        #0x14eb7c0e0<- 1(C) 0x14eb7c0e0
        #0x1050e7f50 <- 2
        #EXC_BAD_ACCESS (code=1, address=0x50e7f50)
        
        bipartition_dict = ctypes.cast(res, ctypes.POINTER(BipartitionDict)).contents
        print("there")
        
        print(bipartition_dict.num_entries)
        
        ## Check if set of internal bipartitions are the same.
        tree1_internals = [bipar.split_as_bitstring() for bipar in tree1.bipartition_encoding if not bipar.is_trivial()]
        tree1_internals_C = []
        tree1_internals_C2 = []
        match_dict = dict()
        score_dict = dict()
        for i in range(bipartition_dict.num_entries): 
            entry = bipartition_dict.entries[i] # 各枝に対応: BipartitionMatches インスタンスへのポインタ
            if(entry.contents.is_external): 
                continue # 外部枝にはmatchの情報がありません。
            matches = entry.contents.matches # matchのポインタ配列
            bipartition = to_bitvector(entry.contents.bipartition) # 枝のbipartitionのBitsインスタンス
            tree1_internals_C2.append(bipartition.bin)
            if bipartition[-1]:
                bipartition = ~bipartition
            tree1_internals_C.append(bipartition.bin)
            #print(bipartition.bin) # Bipartitionのプリント
            match_dict[bipartition.bin] = []; score_dict[bipartition.bin]=[]
            for i in range(K):
                match_bipar = to_bitvector(entry.contents.matches[i])
                if match_bipar[-1]:
                    match_bipar = ~match_bipar
                match_dict[bipartition.bin].append(match_bipar.bin)
                score_dict[bipartition.bin].append(entry.contents.td[i])
        
        assert(len(set(tree1_internals).difference(tree1_internals_C))==0) # two set of bipartitions are the same.
        assert(len(set(tree1_internals_C).difference(tree1_internals))==0)
        
        # tree2.plot_Bio()
        
        # compute by Python
        tree2_bipars = [bipar.split_as_bitstring() for bipar in tree2.bipartition_encoding] #list of strings
        for internal, internal_orig in zip(tree1_internals_C, tree1_internals_C2):
            all_hamming_dist = np.array([Consensus._greedy._MinHammingDist(Bits(bin=internal), Bits(bin=item)) for item in tree2_bipars]) #bitstring.Bits
            order = np.argsort(all_hamming_dist)
            
            before_K = order[all_hamming_dist[order] < all_hamming_dist[order][K]]
            Kplus = order[all_hamming_dist[order] <= all_hamming_dist[order][K-1]]
            
            matches = [tree2_bipars[i] for i in before_K]
            matches2 = [tree2_bipars[i] for i in Kplus]
            scores = all_hamming_dist[order][:K]
            
            # for mat in match_dict[internal]:
            #     if not mat in tree2_bipars:
            #         print("this does not exist", [item.label for i, item in enumerate(tree1.taxon_namespace) if mat[99-i]=="1"])
            
            # check if before_K is included in the matches returned from C
            assert(len(set(matches).difference(set(match_dict[internal]))) == 0)
            
            # check if Kplus includes the matches from C
            assert(len(set(match_dict[internal]).difference(matches2))==0)
            
            assert((scores == score_dict[internal]).all()) # check if scores are the same
        
    finally:   
        print("freeing memory")
        booster_lib.after_tbe_match(res) # free memory
    
    print(type(bipartition.bin))
    
    ## Function to compute first K match using Python code.
    ### Consensus._greedy._MinHammingDist() を使って一つずつ計算 -> 最初のKこを取り出す。
    #print([bipar.split_as_bitstring().count("1") for bipar in tree1.bipartition_encoding if not bipar.is_trivial()])
    #print([bipar.split_as_bitstring().count("1") for bipar in tree1.bipartition_encoding if bipar.is_trivial()]) # is_trivial also contains root edge.
    
    
    
    ## Function to compute first K match using booster lib
    
    
    ## Compare two results.
    