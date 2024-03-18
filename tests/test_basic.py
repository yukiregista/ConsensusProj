import pytest

from Consensus import *


from importlib.resources import files


def test_read_example_data():
    t = Tree_with_support.get(path = files('Consensus.example_data').joinpath('true.tre'), schema="newick")
    print(files('Consensus.example_data').joinpath('true.tre'))


