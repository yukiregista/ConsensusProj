import pytest

from Consensus import *


from importlib.resources import files


def test_read_example_data():
    t = Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.tre'), schema="newick")
    print(files('Consensus.example_data').joinpath('true.tre'))


def test_taxon_namespace_order():
    # class Tree
    t = Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.tre'), schema="newick")
    print(t.taxon_namespace) # make sure that it is sorted
    t2 = Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.tre'), schema="newick", taxon_namespace = t.taxon_namespace)
    assert(t.taxon_namespace == t2.taxon_namespace) # make sure that the two taxon_namespace is the same.
    
    # class TreeList
    t = TreeList_with_support.get(path = files('Consensus.example_data').joinpath('GTR_edit.nex'), schema="nexus")
    print(t.taxon_namespace) # make sure that it is sorted
    t2 = TreeList_with_support.get(path = files('Consensus.example_data').joinpath('GTR_edit.nex'), schema="nexus", taxon_namespace = t.taxon_namespace)
    assert(t.taxon_namespace == t2.taxon_namespace) # make sure that the two taxon_namespace is the same.
    
    
