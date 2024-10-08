import ctypes
import importlib.util
import os
import ctypes

from bitstring import Bits

class Split(ctypes.Structure):
    _fields_ = [("data", ctypes.POINTER(ctypes.c_uint))]

class BipartitionMatches(ctypes.Structure):
    _fields_ = [
        ("bipartition", ctypes.POINTER(ctypes.c_uint)),
        ("topo_depth", ctypes.c_int),
        ("matches", ctypes.POINTER(ctypes.POINTER(ctypes.c_uint))),
        ("td", ctypes.POINTER(ctypes.c_int)),
        ("is_external", ctypes.c_bool)
    ]

class BipartitionDict(ctypes.Structure):
    _fields_ = [
        ("entries", ctypes.POINTER(ctypes.POINTER(BipartitionMatches))),
        ("num_entries", ctypes.c_int),
        ("bipartition_size", ctypes.c_size_t), # number of unsigned integers used to represent the bipartition
        ("num_matches", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
        ("uint_bits", ctypes.c_int) # number of bits in one unsigned integer
    ]


class BipartitionDictArray(ctypes.Structure):
    _fields_ = [
        ("bdict", ctypes.POINTER(BipartitionDict)),
        ("num_trees", ctypes.c_int)
    ]



class BipartitionSupport(ctypes.Structure):
    _fields_ = [
        ("bipartition", ctypes.POINTER(ctypes.c_uint)),
        ("is_external", ctypes.c_bool),
        ("support", ctypes.c_double)
    ]
    
class BipartitionSupportArray(ctypes.Structure):
    _fields_ = [
        ("num_entries", ctypes.c_int),
        ("bipartition_size", ctypes.c_size_t),
        ("n_taxa", ctypes.c_int),
        ("uint_bits", ctypes.c_int),
        ("bipartition_supports", ctypes.POINTER(ctypes.POINTER(BipartitionSupport)))
    ]



def load_booster() -> ctypes.CDLL:
    spec = importlib.util.find_spec("booster")

    if spec is None or spec.origin is None:
        raise ImportError("Cannot find the 'booster' shared library")

    # Load the shared library
    booster_lib = ctypes.CDLL(spec.origin)
    
    return booster_lib

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


def prepare_tbe_match(booster_lib: ctypes.CDLL):
    """Prepare `tbe_match` and `after_tbe_match`

    Args:
        booster_lib (ctypes.CDLL): booster library.
    """
    booster_lib.tbe_match.restype = ctypes.c_void_p # To get pointer accurately
    booster_lib.after_tbe_match.argtypes = ctypes.c_void_p,
    return booster_lib

def prepare_tbe_support(booster_lib: ctypes.CDLL):
    """Prepare `tbe_match` and `after_tbe_match`

    Args:
        booster_lib (ctypes.CDLL): booster library.
    """
    booster_lib.tbe_support.restype = ctypes.c_void_p # To get pointer accurately
    booster_lib.after_tbe_support.argtypes = ctypes.c_void_p,
    return booster_lib


def prepare_tbe_match_args(reftrees_str:str, alttree_str: str, K: int):
    tbe_match_args = ['program_name', '-k', f'{K}']
    argc = len(tbe_match_args)
    argv_ctypes = (ctypes.POINTER(ctypes.c_char) * argc)()
    for i, str in enumerate(tbe_match_args):
        enc_str = str.encode('utf-8')
        argv_ctypes[i] = ctypes.create_string_buffer(enc_str)  
        
    tree1_lines = reftrees_str.split("\n")
    tree1_lines = [item for item in tree1_lines if len(item)>0]
    num_trees = len(tree1_lines)
    tree1_ctypes = (ctypes.POINTER(ctypes.c_char) * num_trees)()
    tree2_ctypes = ctypes.create_string_buffer(alttree_str.encode('utf-8'))
    for i in range(num_trees):
        buffer = ctypes.create_string_buffer(tree1_lines[i].encode('utf-8'))
        tree1_ctypes[i] = buffer
    
    num_trees_ctypes = ctypes.c_int(num_trees)    
    return argc, argv_ctypes, tree1_ctypes, tree2_ctypes, num_trees_ctypes


def prepare_tbe_support_args(reftree_str:str, alttrees_str: str):
    tree2_lines = alttrees_str.split("\n")
    tree2_lines = [item for item in tree2_lines if len(item)>0]
    num_trees = len(tree2_lines)
    tree2_ctypes = (ctypes.POINTER(ctypes.c_char) * num_trees)()
    tree1_ctypes = ctypes.create_string_buffer(reftree_str.encode('utf-8'))
    for i in range(num_trees):
        buffer = ctypes.create_string_buffer(tree2_lines[i].encode('utf-8'))
        tree2_ctypes[i] = buffer
    
    num_trees_ctypes = ctypes.c_int(num_trees)    
    return tree1_ctypes, tree2_ctypes, num_trees_ctypes


def _generate_bitmask(bits: int)-> int:
    return (1<<bits) - 1

def match_score_list(res_ptr: int, id_dict_ref: dict, K: int):
    bdict_arr: BipartitionDictArray = ctypes.cast(res_ptr, ctypes.POINTER(BipartitionDictArray)).contents  
    num_trees: int = bdict_arr.num_trees
    
    id = 0
    
    match_list = []
    score_list = []
    id_dict = dict()
    count_list = []
    topo_depth_list = []
    n_bits = bdict_arr.bdict[0].uint_bits
    mask = _generate_bitmask(n_bits)
    mask_last = _generate_bitmask(bdict_arr.bdict[0].n_taxa % n_bits)
    for i_tree in range(num_trees):
        bipartition_dict = bdict_arr.bdict[i_tree]
        for i in range(bipartition_dict.num_entries): 
            entry = bipartition_dict.entries[i] # 各枝に対応: BipartitionMatches インスタンスへのポインタ
            size = bipartition_dict.bipartition_size
            if(entry.contents.is_external): 
                continue # 外部枝にはmatchの情報がありません。
            # bipartition = to_bitvector(entry.contents.bipartition) # 枝のbipartitionのBitsインスタンス
            
            bipartition = entry.contents.bipartition
            if bipartition[0] %2:
                bipartition = tuple(entry.contents.bipartition[item] for item in range(size))
            else:
                bipartition = tuple([~entry.contents.bipartition[item] & mask for item in range(size-1)] + [~entry.contents.bipartition[size-1] & mask_last])
            
            # #print(bipartition.bin) # Bipartitionのプリント
            if bipartition in id_dict.keys():
                count_list[id_dict[bipartition]] += 1
                continue
            
            id_dict[bipartition] = id
            topo_depth = entry.contents.topo_depth
            better_than_external = True
            matches = []; scores =[]
            for i in range(K):
                # print(entry.contents.td[i])
                if not better_than_external:
                    matches.append(-1)
                    scores.append(topo_depth-1)
                elif topo_depth-2 < entry.contents.td[i]:
                    matches.append(-1)
                    scores.append(topo_depth-1)
                    better_than_external = False
                else:
                    match_bipar = entry.contents.matches[i]
                    if match_bipar[0]%2:
                        match_bipar = id_dict_ref[tuple(match_bipar[i] for i in range(size))]
                    else:
                        try:
                            match_bipar = id_dict_ref[tuple([~match_bipar[i] & mask for i in range(size-1)] + [~match_bipar[size-1] & mask_last])]
                        except:
                            import sys
                            print([entry.contents.bipartition[item] for item in range(size)])
                            print([entry.contents.matches[i][item] for item in range(size)])
                            sys.exit()
                            print(i)
                            print(topo_depth-2 , entry.contents.td[i:10])
                            print(to_bitvector(match_bipar, 10000, 32).bin.count("1"), to_bitvector(match_bipar, 10000, 32).bin.count("0"))
                            print(to_bitvector(entry.contents.bipartition, 10000, 32).bin.count("1"))
                            print((to_bitvector(match_bipar, 10000, 32) & to_bitvector(entry.contents.bipartition, 10000, 32)).bin.count("1"))
                            print(tuple(match_bipar[i] for i in range(size)) in id_dict_ref)
                            similar = []
                            for item in id_dict_ref:
                                bad=False
                                n = len(item)
                                for s in range(n-1):
                                    if item[s] != ~match_bipar[i] & mask:
                                        bad=True
                                        break
                                if not bad:
                                    similar.append(item[-1])
                            print(similar)
                            print(match_bipar[-1], ~match_bipar[size-1] & mask_last)
                            

                                    
                            continue
                    matches.append(match_bipar)
                    scores.append((entry.contents.td[i]))
            match_list.append(matches)
            score_list.append(scores)
            count_list.append(1)
            topo_depth_list.append(entry.contents.topo_depth)
            id += 1
    
    return id_dict, count_list, match_list, score_list, topo_depth_list, num_trees

def tbe_support(res_ptr: int):
    bdict_supp: BipartitionSupportArray = ctypes.cast(res_ptr, ctypes.POINTER(BipartitionSupportArray)).contents
    num_edges = bdict_supp.num_entries
    size = bdict_supp.bipartition_size
    support_list = []; id_dict = dict()
    mask = _generate_bitmask(bdict_supp.uint_bits)
    mask_last = _generate_bitmask(bdict_supp.n_taxa % bdict_supp.uint_bits)
    id = 0
    for i in range(num_edges):
        if bdict_supp.bipartition_supports[i].contents.is_external:
            continue
        bipartition = bdict_supp.bipartition_supports[i].contents.bipartition
        if bipartition[0] %2:
            bipartition = tuple(bipartition[item] for item in range(size))
        else:
            bipartition = tuple([~bipartition[item] & mask for item in range(size-1)] + [~bipartition[size-1] & mask_last])
        support_list.append(bdict_supp.bipartition_supports[i].contents.support)
        id_dict[bipartition] = id
        id += 1
    return id_dict, support_list, bdict_supp.uint_bits, bdict_supp.n_taxa
        
    
def get_bipartition_dict_array(res_ptr: int):
    bdict_arr: BipartitionDictArray = ctypes.cast(res_ptr, ctypes.POINTER(BipartitionDictArray)).contents
    print(bdict_arr)
    print(type(bdict_arr))
    print(bdict_arr.bdict[1])
    print(type(bdict_arr.bdict[1]))
    
    

