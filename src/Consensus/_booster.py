import ctypes
import importlib.util
import os
import ctypes

from bitstring import Bits

def load_booster():
    spec = importlib.util.find_spec("booster")

    if spec is None or spec.origin is None:
        raise ImportError("Cannot find the 'booster' shared library")

    # Load the shared library
    booster_lib = ctypes.CDLL(spec.origin)
    
    return booster_lib

