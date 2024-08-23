import pytest
import Consensus
from importlib.resources import files
import dendropy



def test_first_K_match(K=10):
    
    # example trees
    print(files('Consensus.example_data').joinpath('sample1.nex'))
    tree1 = Consensus.Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.nex'), schema="nexus")
    tree2 = Consensus.Tree_with_support.get(path = files('Consensus.example_data').joinpath('astral_GTRgamma.tre'), schema="newick",
                                            taxon_namespace = tree1.taxon_namespace)
    
    ## Function to compute first K match using Python code.
    ### Consensus._greedy._MinHammingDist() を使って一つずつ計算 -> 最初のKこを取り出す。
    
    ## Function to compute first K match using booster lib
    
    
    ## Compare two results.
    