import pytest
from importlib.resources import files
import ctypes
from bitstring import Bits

import numpy as np



class BipartitionMatchesHash(ctypes.Structure):
    _fields_ = [
        ("h1", ctypes.c_uint),
        ("h2", ctypes.c_uint),
        ("topo_depth", ctypes.c_int),
        ("matched_ids", ctypes.POINTER(ctypes.c_int)),
        ("td", ctypes.POINTER(ctypes.c_int))
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
        ("node_id", ctypes.c_int)
    ]
    
class BipartitionSupportArrayHash(ctypes.Structure):
    _fields_ = [
        ("num_entries", ctypes.c_int),
        ("n_taxa", ctypes.c_int),
        ("bipartition_supoprts", ctypes.POINTER(ctypes.POINTER(BipartitionSupportHash))),
    ]
    
    
class TbeSupportMatchResult(ctypes.Structure):
    __fields__ = [
        ("bdict", ctypes.POINTER(BipartitionDictHash)),
        ("bsupp_arr", ctypes.POINTER(BipartitionSupportArrayHash))
    ]




