import pytest 


from Consensus import *


def test_plot():
    t = Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.tre'), schema="newick")
    newick_str = t.as_string("newick", suppress_rooting=True)
    print(newick_str) # newickで出力
    
    print(plot_example_func(newick_str))