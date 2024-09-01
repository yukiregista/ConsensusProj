import pytest
import Consensus
from importlib.resources import files
import dendropy
import ctypes
from bitstring import Bits

from copy import deepcopy

import numpy as np
from tqdm import tqdm

from bitarray import bitarray


import time

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
        ("bits_uint", ctypes.c_int)
    ]


def int_to_reversed_bin(x, bit_length):
    '''
    This conversion is inefficient!: 
    It would be better if we also use the similar structure (array of uints)
    in Python as well.
    '''
    return ''.join('1' if x & (1 << i) else '0' for i in range(bit_length))

def int_to_reversed_bin2(x, bit_length):
    """
    Convert an integer to a reversed binary string, but avoid string operations.
    Instead of generating a string and then reversing it, we use bitwise operations.
    """
    reversed_bin = 0
    for i in range(bit_length):
        if x & (1 << i):
            reversed_bin |= (1 << (bit_length - 1 - i))
    return reversed_bin


def to_bitvector(bipartition: ctypes.POINTER(ctypes.c_uint), n_tips: int=100, bit_length: int=32) -> Bits:
    n_uints = (n_tips-1) // bit_length + 1
    bipartition_array = [bipartition[i] for i in range(n_uints)]
    bit_size = [bit_length for _ in range(n_uints-1)]
    bit_size.append(n_tips - bit_length * (n_uints-1))
    bit_length = 32
    binary_strings = [int_to_reversed_bin(x, y) for x, y in zip(bipartition_array, bit_size)]
    concatenated_bits = ''.join(binary_strings)
    return Bits(bin=concatenated_bits[::-1])

def to_bitvector2(bipartition: ctypes.POINTER(ctypes.c_uint), n_tips: int=100, bit_length: int=32) -> Bits:
    start_time = time.perf_counter()
    
    n_uints = (n_tips - 1) // bit_length + 1
    bipartition_array = [bipartition[i] for i in range(n_uints)]
    mid_time1 = time.perf_counter()

    bit_size = [bit_length for _ in range(n_uints - 1)]
    bit_size.append(n_tips - bit_length * (n_uints - 1))
    mid_time2 = time.perf_counter()

    concatenated_bits = 0
    total_length = 0
    
    for x, y in zip(bipartition_array, bit_size):
        reversed_bin = int_to_reversed_bin2(x, y)
        concatenated_bits |= reversed_bin << total_length
        total_length += y
    mid_time3 = time.perf_counter()
    
    result = Bits(uint=concatenated_bits, length=n_tips)
    end_time = time.perf_counter()
    
    print(f"Array extraction time: {mid_time1 - start_time:.6f} seconds")
    print(f"Bit size setup time: {mid_time2 - mid_time1:.6f} seconds")
    print(f"Bit manipulation time: {mid_time3 - mid_time2:.6f} seconds")
    print(f"Bits object creation time: {end_time - mid_time3:.6f} seconds")
    print(f"Total time: {end_time - start_time:.6f} seconds")
    
    return result


def to_bitvector_optimized(bipartition: ctypes.POINTER(ctypes.c_uint), n_tips: int=100, bit_length: int=32) -> bitarray:
    # Start timing
    start_time = time.perf_counter()
    
    # Calculate the number of uints needed to represent the bits
    n_uints = (n_tips - 1) // bit_length + 1
    mid_time1 = time.perf_counter()

    # Prepare a bitarray of the required size
    bits = bitarray(n_tips)
    bits.setall(0)  # Initialize all bits to 0
    mid_time2 = time.perf_counter()

    # Process each integer in the bipartition array
    for i in range(n_uints):
        current_uint = bipartition[i]
        for j in range(bit_length):
            if i * bit_length + j < n_tips:
                bits[i * bit_length + j] = (current_uint >> j) & 1
    end_time = time.perf_counter()
    
    # Time monitoring
    print(f"Calculation of n_uints time: {mid_time1 - start_time:.6f} seconds")
    print(f"Bitarray initialization time: {mid_time2 - mid_time1:.6f} seconds")
    print(f"Bit manipulation time: {end_time - mid_time2:.6f} seconds")
    print(f"Total time: {end_time - start_time:.6f} seconds")
    
    return bits

def temp_bitvector(bipartition: ctypes.POINTER(ctypes.c_uint), n_tips: int=100, bit_length: int=32):
    n_uints = (n_tips-1) // bit_length + 1
    bipartition_array = [bipartition[i] for i in range(n_uints)]


def get_first_K_match(tree1: str, tree2: str, booster_lib, K=10):
    booster_lib.tbe_match.restype = ctypes.c_void_p # To get pointer accurately
    booster_lib.print_bdict.argtypes = ctypes.c_void_p,
    booster_lib.after_tbe_match.argtypes = ctypes.c_void_p,
    tbe_match_args = ['program_name', '-k', f'{K}']
    argc = len(tbe_match_args)
    argv_ctypes = (ctypes.POINTER(ctypes.c_char) * argc)()
    for i, str in enumerate(tbe_match_args):
        enc_str = str.encode('utf-8')
        argv_ctypes[i] = ctypes.create_string_buffer(enc_str)
    
    tree1_ctypes = ctypes.create_string_buffer(tree1.encode('utf-8'))
    tree2_ctypes = (ctypes.POINTER(ctypes.c_char) * 1)()
    tree2_ctypes[0] = ctypes.create_string_buffer(tree2.encode('utf-8'))
    
    try:
        start_match = time.perf_counter()
        res = booster_lib.tbe_match(argc, argv_ctypes, tree1_ctypes ,tree2_ctypes)
        end_match = time.perf_counter()
        bipartition_dict = ctypes.cast(res, ctypes.POINTER(BipartitionDict)).contents
        
        # convert to python object
        
        # method 1: Store as a tuple of integers
        match_dict = dict()
        score_dict = dict()
        for i in range(bipartition_dict.num_entries): 
            entry = bipartition_dict.entries[i] # 各枝に対応: BipartitionMatches インスタンスへのポインタ
            if(entry.contents.is_external): 
                continue # 外部枝にはmatchの情報がありません。
            # bipartition = to_bitvector(entry.contents.bipartition) # 枝のbipartitionのBitsインスタンス
            bipartition = entry.contents.bipartition
            if bipartition[0] %2:
                bipartition = (entry.contents.bipartition[0], entry.contents.bipartition[1], 
                            entry.contents.bipartition[2], entry.contents.bipartition[3])
            else:
                bipartition = (~entry.contents.bipartition[0], ~entry.contents.bipartition[1], 
                            ~entry.contents.bipartition[2], ~entry.contents.bipartition[3])
            
            # #print(bipartition.bin) # Bipartitionのプリント
            match_dict[bipartition] = []; score_dict[bipartition]=[]
            for i in range(K):
                match_bipar = entry.contents.matches[i]
                if match_bipar[0]%2:
                    match_bipar = (match_bipar[0], match_bipar[1], 
                            match_bipar[2], match_bipar[3])
                else:
                    match_bipar = (~match_bipar[0], ~match_bipar[1], 
                            ~match_bipar[2], ~match_bipar[3])
                match_dict[bipartition].append(match_bipar)
                score_dict[bipartition].append((entry.contents.td[i]))
        # method 1: end
            
        # # ## method 2 :naive -> very slow! 
        # match_dict = dict()
        # score_dict = dict()
        # for i in range(bipartition_dict.num_entries): 
        #     entry = bipartition_dict.entries[i] # 各枝に対応: BipartitionMatches インスタンスへのポインタ
        #     if(entry.contents.is_external): 
        #         continue # 外部枝にはmatchの情報がありません。
        #     bipartition = to_bitvector2(entry.contents.bipartition) # 枝のbipartitionのBitsインスタンス
        #     if bipartition[-1]:
        #         bipartition = ~bipartition
        #     # #print(bipartition.bin) # Bipartitionのプリント
        #     match_dict[bipartition.bin] = []; score_dict[bipartition.bin]=[]
        #     for i in range(K):
        #         match_bipar = to_bitvector2(entry.contents.matches[i])
        #         if match_bipar[-1]:
        #             match_bipar = ~match_bipar
        #         match_dict[bipartition.bin].append(match_bipar.bin)
        #         score_dict[bipartition.bin].append((entry.contents.td[i]))
        # convert_fin = time.perf_counter() 
        
        
        ### method 3 using bitarray -> still very slow!
                # ## method 2 :naive -> very slow! 
        # match_dict = dict()
        # score_dict = dict()
        # for i in range(bipartition_dict.num_entries): 
        #     entry = bipartition_dict.entries[i] # 各枝に対応: BipartitionMatches インスタンスへのポインタ
        #     if(entry.contents.is_external): 
        #         continue # 外部枝にはmatchの情報がありません。
        #     bipartition = to_bitvector_optimized(entry.contents.bipartition) # 枝のbipartitionのBitsインスタンス
        #     if bipartition[-1]:
        #         bipartition = ~bipartition
        #     # #print(bipartition.bin) # Bipartitionのプリント
        #     match_dict[bipartition.tobytes()] = []; score_dict[bipartition.tobytes()]=[]
        #     for i in range(K):
        #         match_bipar = to_bitvector_optimized(entry.contents.matches[i])
        #         if match_bipar[-1]:
        #             match_bipar = ~match_bipar
        #         match_dict[bipartition.tobytes()].append(match_bipar.tobytes())
        #         score_dict[bipartition.tobytes()].append((entry.contents.td[i]))
        # convert_fin = time.perf_counter() 
                
        # match_dict = deepcopy(match_dict)
        # score_dict = deepcopy(score_dict)
        print(bipartition_dict.num_entries)
        
    finally:
        booster_lib.after_tbe_match(res)
    
    
    # return match_dict, score_dict
    return end_match - start_match, convert_fin - end_match
    
    
    


def test_match_and_convert(K=10):
    # example trees
    print(files('Consensus.example_data').joinpath('sample1.nex'))
    tree1 = Consensus.Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.nex'), schema="nexus")
    tree2 = Consensus.TreeList_with_support.get(path = files('Consensus.example_data').joinpath('GTRgamma_edit.nex'), schema="nexus",
                                            taxon_namespace = tree1.taxon_namespace)
    
        ## For using booster
    tree1_str = tree1.as_string("newick", suppress_rooting=True)
    tree2_str = tree2.as_string("newick", suppress_rooting=True)
    
    booster_lib = Consensus.load_booster()
    time_match = 0; time_convert = 0
    for i, line in tqdm(enumerate(tree2_str.split("\n"))):
        booster_lib = Consensus.load_booster()
        if not len(line): continue
        a, b = get_first_K_match(tree1_str, line, booster_lib, K=K)
        time_match += a
        time_convert += b
    
    print(f"Time for match (K={K}): {time_match}\n Time for convert (K={K}): {time_convert}")
        
        
if __name__ == "__main__":
    test_match_and_convert()
    