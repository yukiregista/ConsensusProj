import dendropy
import numpy as np
from bitstring import Bits
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt
import warnings
from collections import OrderedDict
import sys
from copy import deepcopy
import itertools
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import to_tree
from scipy.cluster.hierarchy import linkage

def _get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
    
    Code copied from https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format.

    Parameters
    ----------
    node : _type_
        output of sciply.cluster.hierarchy.to_tree()
    parent_dist : _type_
        output of sciply.cluster.hierarchy.to_tree().dist
    leaf_names : List
        list of leaf names
    newick : str, optional
        leave empty, this variable is used in recursion., by default ''

    Returns
    -------
    str
        tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = _get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = _get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

def TBE(bipartitions, treelist):
    """Compute TBE of given set of bipartitions against input trees.

    Parameters
    ----------
    bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        Bipartitions to evaluate TBE support.
    treelist : TreeList_with_support or Iterable of either Tree_with_support or 
        Input trees to evaluate TBE support against.
    
    Notes
    -----
    Trees that produce input bipartitions and treelist need to have the same taxon_namespace.
    
    Otherwise, the produced TBE support will not be a valid one.
    
    Please check this manually, as we don't check it inside the function.

    Returns
    -------
    numpy.ndarray
        Array of TBE support. 
    """
    
    n_taxa = len(treelist[0].taxon_namespace)
    refinfo = _create_refinfo(bipartitions, n_taxa)
    n_bipartitions = len(bipartitions)
    totalSupport = np.zeros(n_bipartitions)

    for tree in treelist:
        if tree.bipartition_encoding is None:
            tree.encode_bipartitions()
        tree_bipartition_ints = [bipartition.split_as_int() for bipartition in tree.bipartition_encoding]
        for i, bipartition in enumerate(bipartitions):
            bipar_int = bipartition.split_as_int()
            if bipar_int in tree_bipartition_ints:
                support = 1
            elif refinfo[bipar_int][1] == 2:
                support  = 0
            else:
                support = 1 - _minDist(refinfo[bipar_int], tree) / (refinfo[bipar_int][1] - 1)
            totalSupport[i] += support
    
    totalSupport = totalSupport / len(treelist)
    return totalSupport

def unnormalized_TBE(bipartitions, treelist):
    """Unnormalized version of transfer bootstrap expectation.
    
    Due to the lack of normalization, contrary to the usual TBE, the lower value (close to zero) indicates well supported edges.

    Parameters
    ----------
    bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        Bipartitions to evaluate unnormalized TBE.
    treelist : TreeList_with_support or Iterable of either Tree_with_support or 
        Input trees to evaluate unnormalized TBE against.

    Returns
    -------
    numpy.ndarray
        Array of unnormalized TBE. 
    """
    n_taxa = len(treelist[0].taxon_namespace)
    refinfo = _create_refinfo(bipartitions, n_taxa)
    n_bipartitions = len(bipartitions)
    totalDistance = np.zeros(n_bipartitions)

    for tree in treelist:
        if tree.bipartition_encoding is None:
            tree.encode_bipartitions()
        tree_bipartition_ints = [bipartition.split_as_int() for bipartition in tree.bipartition_encoding]
        for i, bipartition in enumerate(bipartitions):
            bipar_int = bipartition.split_as_int()
            if bipar_int in tree_bipartition_ints:
                distance = 0
            elif refinfo[bipar_int][1] == 2:
                distance  = 1
            else:
                distance = _minDist(refinfo[bipar_int], tree)
            totalDistance[i] += distance
    
    totalDistance = totalDistance/ len(treelist)
    return totalDistance


def _create_refinfo(bipartitions, n_taxa):
    """Helper function for TBE

    Parameters
    ----------
    bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        _description_
    n_taxa : int
        _description_

    Returns
    -------
    dict
        bipar_int as key an tuple of (bitstr, p) as value.
    """
    
    # bipartitions: iterable of dendropy.Bipartition
    refinfo = dict()
    for bipartition in bipartitions:
        edge_bitstr = Bits(uint = bipartition.split_as_int(), length=n_taxa)
        counts = (edge_bitstr.count(0), edge_bitstr.count(1))
        if counts[0] < counts[1]:
            # zero lready assigned to light side
            bitstr = edge_bitstr
            p = counts[0]
        else:
            # reverse bitstring
            bitstr = (~edge_bitstr)
            p = counts[1]
        refinfo[bipartition.split_as_int()] = (bitstr, p)
    # O( len(bipartitions) *  n_taxa )
    return refinfo

def _minDist(refinfo_b, tree):
    """Helper minDist function for computing TBE.

    Parameters
    ----------
    refinfo_b : dict
        refinfo of bipartition b created by _create_refinfo function.
    tree : Tree_with_support or dendropy.datamodel.treemodel.Tree
        Tree to compute support of the bipartition
        
    Notes
    -----
    This function expects input bipartition (``refinfo_b``) to be NOT included in ``tree``.

    Returns
    -------
    int
        Hamming distance.
    """
    # bipartition: Bipartition object
    # tree: Tree_with_support or Tree object
    p = refinfo_b[1]
    b_bitstr = refinfo_b[0]

    d = p - 1
    #edge_bitstr = Bits(uint = bipartition.split_as_int(), length=self.n_taxa)
    #p = min(edge_bitstr.count(1), edge_bitstr.count(0))
    #m = len(tree.edge)

    # WE ASSUME THAT bipartition and tree has the exact same TAXON_NAMESPACE 
    taxon_namespace = tree.taxon_namespace
    taxon_labels = [taxon.label for taxon in taxon_namespace]
    n_taxa = len(taxon_labels)
    taxon_labels_to_bitstr_digit = dict(zip(taxon_labels, [i for i in range(n_taxa)]))
    node_biparint_to_postorder_index = dict()

    if tree.bipartition_encoding is None:
        tree.encode_bipartitions()
    
    m = len(tree.bipartition_encoding)
    CountOnesSubtree = np.zeros(m)
    CountZerosSubtree = np.zeros(m)

    for i, node in enumerate(tree.postorder_node_iter(exclude_seed_node=True)):
        node_biparint_to_postorder_index[node.bipartition.split_as_int()] = i
        if node.is_leaf():
            digit = taxon_labels_to_bitstr_digit[node.taxon.label]
            CountOnesSubtree[i] =  int(b_bitstr[- digit - 1])
            CountZerosSubtree[i] = 1 - CountOnesSubtree[i]
        else:
            for child in node.child_node_iter(exclude_seed_node=True):
                CountOnesSubtree[i] += CountOnesSubtree[ node_biparint_to_postorder_index[child.bipartition.split_as_int()] ]
                CountZerosSubtree[i] += CountZerosSubtree[ node_biparint_to_postorder_index[child.bipartition.split_as_int()] ]
            actDist = p - CountZerosSubtree[i] + CountOnesSubtree[i]
            if actDist > n_taxa / 2:
                actDist = n_taxa - actDist
            if actDist < d:
                d = actDist
                if d == 1:
                    return d
    
    return d

class Tree_with_support(dendropy.Tree):
    """Child class of `dendropy.datamodel.treemodel.Tree`, with supportfor support values.
    """
    def __init__(self, *args, **kwargs):
        bootstrap_support = kwargs.pop("bootstrap_support", None)
        TBE_support = kwargs.pop("TBE_support", None)
        super().__init__(*args, **kwargs)
        self.bootstrap_support = bootstrap_support
        self.TBE_support = TBE_support
        self.n_taxa = len(self.taxon_namespace)
        self.encode_bipartitions() #some methods (from_bipartitions etc) seem to produce a tree with bad bipartition_encoding, so renew it 

    def compute_bootstrap(self, treelist: TreeList_with_support):
        """Compute bootstrap values of self against input ``treelist`` and update `self.bootstrap_support`

        Parameters
        ----------
        treelist : TreeList_with_support
            List of trees to evaluate bootstrap support against.

        Returns
        -------
        None
        """
        assert self.taxon_namespace == treelist.taxon_namespace
        bootstrap_support = dict()
        if self.bipartition_encoding is None:
            self.encode_bipartitions()
        for bipar in self.bipartition_encoding:
            key = bipar.split_as_int()
            if key in treelist.edge_dict.keys():
                bootstrap_support[key] = treelist.edge_dict[key]/len(treelist)
            else:
                bootstrap_support[key] = 0
        self.bootstrap_support = bootstrap_support
        return None
    
    
    def compute_TBE(self, treelist: TreeList_with_support):
        """Compute TBE support of self against input ``treelist`` and update `self.TBE_support`

        Parameters
        ----------
        treelist : TreeList_with_support
            List of trees to evaluate TBE support against.

        Returns
        -------
        None
        """
        assert self.taxon_namespace == treelist.taxon_namespace
        if self.bipartition_encoding  is None:
            self.encode_bipartitions()
        bipartitions = self.bipartition_encoding
        bipartitions_ints = [bipar.split_as_int() for bipar in bipartitions]
        TBE_list = TBE(bipartitions, treelist)
        self.TBE_support = dict(zip(bipartitions_ints, TBE_list))
        return None

    def clone(self, depth=1):
        cloned_tree = super().clone(depth=1)
        return Tree_with_support(cloned_tree, bootstrap_support=self.bootstrap_support, TBE_support = self.TBE_support)
        

    def _make_Bio_tree(self, support = None):
        """Helper function to convert self to Bio tree.

        Parameters
        ----------
        support : dict or None, optional
            Support values to provide to Bio tree.   
            If dict, it should have integer representing each bipartition as key and support values as value.
            If None, no support values are provided to Bio tree.
            By default None

        Returns
        -------
        Bio.Phylo.Newick.Tree
            Newly created Bio tree.
        """
        tree_Bio = Phylo.read(file = StringIO(self.as_string(schema="newick")) , format = "newick" )
        if support is not None:
            ##### start -- add support information to each edge
            for clade in tree_Bio.find_clades(order="postorder"):
                terminals = [item.name for item in clade.get_terminals()]
                taxonnames_array = np.array([item.label for item in self.taxon_namespace])
                digits = [ np.where( item == taxonnames_array )[0][0] for item in terminals ]
                clade_bool = [False for i in range(self.n_taxa)]
                for digit in digits:
                    clade_bool[-(digit+1)] = True
                clade_bit = Bits(clade_bool) # bitstring corresponding to the clade
                if int(clade_bit.bin[-1]) == 1:
                    clade_bit = (~clade_bit) # convert to least significant 0 format
                if clade_bit.uint in support.keys():
                    clade_support = support[clade_bit.uint]
                    clade.confidence = clade_support
            #####  end  -- add support information to each edge
        return tree_Bio  

    def plot_Bio(self, ax = None, type = "plain", decimals = None):
        """Function to plot the tree using Bio.Phylo package.
        
        This function converts self to Bio.Phylo tree, and plot using Bio.Phylo's functionality.

        Parameters
        ----------
        ax : matplotlib.axes or None, optional
            Axes to plot the tree on. If None, new figure will be created.
        type : str, optional
            The type of support values to plot, by default "plain" (no support values are plotted).
            If not "plain", it should be one of the followings:
            
            - "bootstrap" : bootstrap support
            - "TBE" : TBE support
            
            If these support values are not computed before calling this function, it will use the option "plain".
            
            Please compute the support values beforehand.
        decimals : int or None, optional
            Decimal precision of support values, by default None (no rounding).
        """

        ## type: one of ["plain", "bootstrap", "TBE"]. If other name is provided, it will automatically plot "plain" one.
        ## decimals: int, if provided, round confidence number at that number

        if ax is None:
            fig, ax = plt.subplots(figsize = (20,20))
        tree_Bio = Phylo.read(file = StringIO(self.as_string(schema="newick")) , format = "newick" )

        if type == "bootstrap":
            if self.bootstrap_support is None:
                warnings.warn("Bootstrap values has not been provided to the instance. Plotting the plain tree.")
            if decimals is not None:
                boot = {k:np.round(v, decimals=decimals) for k, v in self.bootstrap_support.items()}
            else:
                boot = self.bootstrap_support
            tree_Bio = self._make_Bio_tree(boot)
            Phylo.draw( tree_Bio, axes = ax)
        elif type == "TBE":
            if self.TBE_support is None:
                warnings.warn("TBE values has not been provided to the instance. Plotting the plain tree.")
            if decimals is not None:
                tbe_supp = {k:np.round(v, decimals=decimals) for k, v in self.TBE_support.items()}
            else:
                tbe_supp = self.TBE_support
            tree_Bio = self._make_Bio_tree(tbe_supp)
            Phylo.draw( tree_Bio, axes = ax)
        else:
            tree_Bio = self._make_Bio_tree()
            Phylo.draw( tree_Bio, axes = ax) 
        return ax
    
    def branch_resolution(self):
        """Returns branch resolution of self.
        
        The branch resolution is defined by:
        
        .. math::

            \begin{frac}{\text{\# internal edges in self}}{\text{self.n_taxa}-3}

        Returns
        -------
        _type_
            _description_
        """
        # returns how resolved tree is w.r.t. #internal branches
        possible_num = self.n_taxa - 3 # dendropy Tree stores n external edges + 1 external edge connecting to seed_node. There are n-3 internal edges.
        internal_num = len(self.inernal_edges(exclude_seed_edge=True))
        return internal_num / possible_num


    @classmethod
    def get(cls, **kwargs):
        tree = dendropy.Tree.get(**kwargs)
        bootstrap_support = kwargs.pop("bootstrap_support", None)
        TBE_support = kwargs.pop("TBE_support", None)
        return cls(tree, bootstrap_support = bootstrap_support, TBE_support = TBE_support)


def _compatible(bits_a : Bits, bits_b : Bits):
    """Check if two bitsting represented clades are compatible.

    Parameters
    ----------
    bits_a : Bits
        Bitstring representation of first bipartition.
    bits_b : Bits
        Bitstring representation of second bipartition.

    Returns
    -------
    int
        1 indicates compatible, 0 indicates incompatible.
    """
    count1 = ( bits_a & bits_b ).count(1)
    count2 = ((~bits_a) & bits_b).count(1)
    count3 = (bits_a & (~bits_b)).count(1)
    count4 = ((~bits_a) & (~bits_b)).count(1)
    if np.min([count1, count2, count3, count4]) == 0:
        return 1 # compatible
    else:
        return 0 # incompatible

def _Rstar(trees,root_name : str):
    """Computes R* consensus tree using Bio.Phylo package.

    Parameters
    ----------
    trees : Iterable of Bio.Phylo.Newick.Tree
        Input trees to compute R* consensus.
    root_name : str
        Name of the root node picked.

    Returns
    -------
    Bio.Phylo.Newick.Tree
        R* consensus tree
    """
    # root: clade instance of the leaf we make a root.
    # currently very inefficient and messy
    rooted_trees = []
    for tree in trees:
        rooted_tree = deepcopy(tree)
        terminals = rooted_tree.get_terminals()
        root_ind = np.where(np.array([item.name for item in terminals]) == root_name)[0][0]
        root = terminals[root_ind]
        rooted_tree.root_with_outgroup(root)
        rooted_tree.prune(root)
        rooted_trees.append(rooted_tree)
    
    all_taxa_list = [item.name for item in rooted_trees[0].get_terminals()]
    all_taxa_set = set(all_taxa_list)

    
    # Initialization: dict_RT_count O(n^3)
    dict_RT_count = {}
    for combi in itertools.combinations(all_taxa_list, 3):
        dict_RT_count[ ( frozenset((combi[0], combi[1])) ,combi[2] ) ] = 0
        dict_RT_count[ ( frozenset((combi[1], combi[2])) ,combi[0] ) ] = 0
        dict_RT_count[ ( frozenset((combi[2], combi[0])) ,combi[1] ) ] = 0
    # For each tree, we traverse by DFS and add counts of rooted triples O(Kn^3)
    for tree in rooted_trees:
        for clade in tree.find_clades(order="postorder"):
            nonterminals = clade.get_nonterminals()
            if len(nonterminals) == 0:
                continue # leaf clade
            
            # obtain the list of subtrees
            subtree_list = clade.clades
            # obtain list of list of taxon in each subtree
            subtree_taxa_list = [ [item.name for item in  tr.get_terminals()] for tr in subtree_list]
            # set of taxon in this current clade and the other side of bipartition
            taxa_set = {item.name for item in clade.get_terminals()}
            taxa_other = all_taxa_set - taxa_set
            # if this corresponds to root, skip
            if len(taxa_other)==0:
                continue

            for st_a, st_b in itertools.combinations(subtree_taxa_list, 2):
                for taxon_a in st_a:
                    for taxon_b in st_b:
                        for taxon_c in taxa_other:
                            # add count of (a,b)|c
                            dict_RT_count[(frozenset({taxon_a, taxon_b}), taxon_c)] += 1
    # Initialization of Rmaj
    #print(np.max(list(dict_RT_count.values())))
    Rmaj = {}
    for combi in itertools.combinations(all_taxa_set, 2):
        Rmaj[frozenset({combi[0], combi[1]})] = list()
    # again iterate over all rooted triples
    for combi in itertools.combinations(all_taxa_set, 3):
        count_0 = dict_RT_count[(frozenset({combi[0], combi[1]}) , combi[2])]
        count_1 = dict_RT_count[(frozenset({combi[1], combi[2]}) , combi[0])]
        count_2 = dict_RT_count[(frozenset({combi[2], combi[0]}) , combi[1])]
        counts = [count_0, count_1, count_2]
        argmax_ind = np.where(counts == np.max(counts))[0]
        if len(argmax_ind) > 1:
            # multiple candidates, so add nothing to Rmaj
            continue
        # else, add to Rmaj
        Rmaj[ frozenset( { combi[argmax_ind[0]], combi[(argmax_ind[0] + 1)%3 ] } ) ].append(combi[(argmax_ind[0] + 2)%3])
    n_taxa = len(all_taxa_list)
    s_maj = np.diag([-1 for i in range(n_taxa)])# initialization
    for i in range (n_taxa):
        for j in range(i):
            s_maj[i,j] = s_maj[j,i] = len( Rmaj[ frozenset( { all_taxa_list[i], all_taxa_list[j] } ) ] )
    condensed_similarity = squareform(s_maj, checks=False) # condensed SIMILARITY matrix
    condensed_dissimilarity = - condensed_similarity + n_taxa # transforming into nonnegative value as well
    linkage_matrix = linkage(condensed_dissimilarity, method = "single", metric="IGNORED")
    slink_tree = to_tree(linkage_matrix)
    slink_newick = _get_newick(slink_tree, slink_tree.dist, leaf_names=all_taxa_list)
    slink = Phylo.read(StringIO(slink_newick), "newick")
    #unit branch length version
    slink_unit = deepcopy(slink)
    for clade in slink_unit.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = 1
    collapse_clades = []
    for clade in slink_unit.find_clades(order = "postorder"):
        if clade.is_terminal():
            # leaf node
            continue
        terminals = {item.name for item in clade.get_terminals()}
        outsiders = all_taxa_set - terminals
        if len(outsiders) ==0:
            # trivial clade
            continue
        end = False
        for a, b in itertools.combinations(terminals, 2):
            for c in outsiders:
                if c not in Rmaj[ frozenset( {a, b} ) ]:
                    collapse_clades.append(clade)
                    end = True
                    break
            if end:
                break
    slink_unit.collapse_all(lambda c: c in collapse_clades)
    # Phylo.write(slink_unit, f"simiiformes/rstar/multi_{root_name}.txt", format="newick")
    return slink_unit


class TreeList_with_support(dendropy.TreeList):
    """Child class of dendropy.datamodel.treecollectionmodel.TreeList, with support for several consensus constructions and support values.
    """
    # EXPECTS TreeList with the same leaves. 
    # EXPECTS TreeList of unrooted trees

    def _make_edge_dict(self):
        edge_dict = OrderedDict()
        for tree in self:
            dict_keys = edge_dict.keys() # no duplicate keys in one tree
            for bipartition in tree.encode_bipartitions():
                key = bipartition.split_as_int()
                if key in dict_keys:
                    edge_dict[key] += 1
                else:
                    edge_dict[key] = 1
        return edge_dict
    
    def _make_all_TBE_dict(self):
        keys = list(self.edge_dict.keys())
        bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.n_taxa - 1, compile_bipartition=True) for item in keys]
        tbe_list = self.TBE(bipartitions)
        bipartition_ints = [bipar.split_as_int() for bipar in bipartitions]
        self.all_TBE_dict =  dict(zip(bipartition_ints, tbe_list))
        return None

    def __init__(self, *args, **kwargs):
        edge_dict = kwargs.pop("edge_dict", None)
        self.all_TBE_dict =  kwargs.pop("all_TBE_dict", None)
        if len(args) ==  1 and isinstance(args[0], dendropy.TreeList):
            # just reference it
            self.__dict__.update(args[0].__dict__)
        else:
            super().__init__(*args, **kwargs)
        trees_init = True
        if len(args) == 1:
            try:
                if len(args[0]) > 0 and all(isinstance(x, Tree_with_support) for x in args[0]):
                    # the given trees are Tree_with_support, so initialize it with this (not to lose additional info)
                    self._trees = list(args[0])
                    trees_init = False
                    #self.taxon_namespace = args[0][0].taxon_namespace
            except:
                pass
        if trees_init:
            self._trees = [Tree_with_support(tree) for tree in self._trees]
        self.edge_dict = edge_dict
        self.n_trees = len(self)
        self.n_taxa = len(self[0].leaf_nodes())
        # self.bipartition_encodings_list = [tree.encode_bipartitions() for tree in self]
        if self.edge_dict is None:
            self.edge_dict = self._make_edge_dict()
    
    @classmethod
    def get(cls, **kwargs):
        edge_dict = kwargs.pop("edge_dict", None)
        tree = dendropy.TreeList.get(**kwargs)
        return cls(tree, edge_dict = edge_dict)



    def bootstrap_support(self, bipartition):
        """Compute bootstrap support of a bipartition against self.

        Parameters
        ----------
        bipartition : `dendropy.datamodel.treemodel.Bipartition`
            Bipartition to evaluate bootstrap support

        Returns
        -------
        float
            Bootstrap support.
        """
        # computes bootstrap support of an edge.
        # bipartition: Bipartition object

        key = bipartition.split_as_int()
        if key in self.edge_dict.keys():
            return self.edge_dict[key]/self.n_trees
        return 0
    
    def compute_all_TBE(self):
        """Compute all TBE of all bipartitions present in the treelist against self and save it in self.all_TBE_dict.

        Returns
        -------
        dict
            Dictionary that has integers representing bipartitions as keys and TBE support as values. 
        """
        self._make_all_TBE_dict()
        return self.all_TBE_dict
    
    def TBE(self, bipartitions):
        """Compute TBE support of bipartitions against self.

        Parameters
        ----------
        bipartitions : Iterable of `dendropy.datamodel.treemodel.Bipartition`
        Bipartitions to evaluate TBE support.

        Returns
        -------
        numpy.ndarray
            Array of TBE support
        """
        return TBE(bipartitions, self)
    
    def majority_consensus(self):
        """Computes majority rule consensus.

        Returns
        -------
        Tree_with_support
            Majority rule consensus tree.
        """
        tree_leafset_bitmask = 2**self.n_taxa - 1 # corresponds to "1111......1"
        bipartitions = []; bipar_ints = []; bipar_bs = []
        for k, v in self.edge_dict.items():
            if v > self.n_trees/2:
                bipartitions.append(dendropy.Bipartition(bitmask = k, leafset_bitmask = k, tree_leafset_bitmask = tree_leafset_bitmask, is_rooted=False))
                bipar_ints.append(k); bipar_bs.append(v/self.n_trees)
        majority_tree = dendropy.Tree.from_bipartition_encoding(bipartitions, self.taxon_namespace)
        bootstrap_support = dict(zip(bipar_ints, bipar_bs))
        # majority_tree.bootstrap_support = bootstrap_support
        return Tree_with_support(majority_tree, bootstrap_support = bootstrap_support, taxon_namespace = self.taxon_namespace)

    def MAP(self):
        """Returns the most appearing tree topology.

        Returns
        -------
        List
            List of the most appearing tree topologies (has length > 1 when more than two topologies have the same maximum count).
        
        int
            The number of occurrences of MAP trees.
        """
        # returns tree(s) with maximum number of appearances
        dict_topologies = dict()
        for tree_idx in range(len(self)):
            tree = self[tree_idx]
            bipar_set = frozenset({node.bipartition.split_as_int() for node in tree.postorder_internal_node_iter(exclude_seed_node=True)})
            if bipar_set in dict_topologies.keys():
                dict_topologies[bipar_set].append(tree_idx)
            else:
                dict_topologies[bipar_set] = [tree_idx]
        len_list = [len(item) for item in dict_topologies.values()]
        max_len = np.max(len_list)
        max_len_topology_indices = [item for item in dict_topologies.values() if len(item)==max_len]
        repr_trees = []
        #print(max_len_topology_indices)
        for item in max_len_topology_indices:
            repr_tree = self[item[0]].clone(depth=1)
            len_dict = dict()
            for edge in self[item[0]].postorder_internal_edge_iter(exclude_seed_edge=True):
                len_dict[edge._bipartition.split_as_int()] = edge.length/max_len
            for treeind in item[1:]:
                tree = self[treeind]
                for edge in tree.postorder_internal_edge_iter(exclude_seed_edge=True):
                    len_dict[edge._bipartition.split_as_int()] += edge.length/max_len
            for edge in repr_tree.postorder_internal_edge_iter(exclude_seed_edge=True):
                edge.length = len_dict[edge._bipartition.split_as_int()]
            repr_trees.append(repr_tree)
        return repr_trees, max_len
        

    
    def greedy_TBE_consensus(self, normalization : bool = False):
        """Greedily (in terms of TBE or its unnormalized version) add edges to construct consensus tree.
        
        The candidate bipartitions are the set of all bipartitions included in self.

        Parameters
        ----------
        normalization : bool, optional
            Whether to use the usual TBE or unnormalized TBE.

        Returns
        -------
        Tree_with_support
            Consensus tree.
        """
        all_edges = list(self.edge_dict.keys())
        all_bipartitions = [dendropy.Bipartition(leafset_bitmask = item, tree_leafset_bitmask = 2**self.n_taxa - 1) for item in all_edges]
        all_Bits = [Bits(uint = item.split_as_int(), length=self.n_taxa) for item in all_bipartitions]
        # compile all bipartitions
        for bipartition in all_bipartitions:
            bipartition.compile_leafset_bitmask()
        
        encoding_bipartitions = []; encoding_TS = []
        if normalization:
            if self.all_TBE_dict is None:
                tbe_list=  self.TBE(all_bipartitions)
            else:
                tbe_list = list(self.all_TBE_dict.values())
            sort_index = np.argsort(tbe_list)[::-1] # high to low
            all_Bits_sorted = [all_Bits[item] for item in sort_index]
            all_bipartitions_sorted =  [all_bipartitions[item] for item in sort_index]
            tbe_list_sorted = [tbe_list[item] for item in sort_index]
            while len(all_Bits_sorted) > 0:
                bipar_bits = all_Bits_sorted.pop(0)
                TS = tbe_list_sorted.pop(0)
                bipartition = all_bipartitions_sorted.pop(0)
                encoding_TS.append(TS)
                encoding_bipartitions.append(bipartition)
                # delete incompatible edges
                if TS == 1:
                    count_1 = bipar_bits.count(1)
                    if (count_1 <= 1) or (count_1 >= self.n_taxa - 1):
                        continue # compatible with all edges
                
                remain_index = [i for i in range(len(all_Bits_sorted)) if _compatible(all_Bits_sorted[i], bipar_bits) ]
                all_Bits_sorted = [all_Bits_sorted[item] for item in remain_index]
                all_bipartitions_sorted = [all_bipartitions_sorted[item] for item in remain_index]
                tbe_list_sorted = [tbe_list_sorted[item] for item in remain_index]
        else:
            if self.all_TBE_dict is None:
                tbe_list= unnormalized_TBE(all_bipartitions, self)
            else:
                normalized_tbe_list = np.array(list(self.all_TBE_dict.values()))
                tmp = 1 - normalized_tbe_list
                refinfos = _create_refinfo(all_bipartitions, self.n_taxa)
                normalizing_constants = np.array([refinfos[bipar_int][1] - 1 for bipar_int in all_edges])
                tbe_list = tmp * normalizing_constants

            sort_index = np.argsort(tbe_list) # low to high
            all_Bits_sorted = [all_Bits[item] for item in sort_index]
            all_bipartitions_sorted =  [all_bipartitions[item] for item in sort_index]
            tbe_list_sorted = list(tbe_list[sort_index])
            while len(all_Bits_sorted) > 0:
                bipar_bits = all_Bits_sorted.pop(0)
                TS = tbe_list_sorted.pop(0)
                bipartition = all_bipartitions_sorted.pop(0)
                encoding_TS.append(TS)
                encoding_bipartitions.append(bipartition)

                # delete incompatible edges
                if TS == 0:
                    count_1 = bipar_bits.count(1)
                    if (count_1 <= 1) or (count_1 >= self.n_taxa - 1):
                        continue # compatible with all edges
                remain_index = [i for i in range(len(all_Bits_sorted)) if _compatible(all_Bits_sorted[i], bipar_bits) ]
                all_Bits_sorted = [all_Bits_sorted[item] for item in remain_index]
                all_bipartitions_sorted = [all_bipartitions_sorted[item] for item in remain_index]
                tbe_list_sorted = [tbe_list_sorted[item] for item in remain_index]
        
        encoding_bipartition_ints = [bipar.split_as_int() for bipar in encoding_bipartitions]
        greedy_tree = dendropy.Tree.from_bipartition_encoding(encoding_bipartitions, self.taxon_namespace)
        if normalization:
            return Tree_with_support(greedy_tree, TBE_support = dict(zip( encoding_bipartition_ints, encoding_TS )))
        else:
            return Tree_with_support(greedy_tree)

    def Rstar_consensus(self, root_index = None, root_name = None):
        """Compute R* consensus.

        Parameters
        ----------
        root_index : int, optional
            The index of root in self.taxon_namespace, by default None.
            You need to specify either ``root_index`` or ``root_name``.
        root_name : _type_, optional
            The name of root taxa, by default None
            You need to specify either ``root_index`` or ``root_name``.

        Returns
        -------
        Tree_with_support
            R* consensus tree.
        """

        # specify either root_index or root_name
        # currently very inefficient and messy

        if root_index is not None:
            if root_name is not None:
                sys.exit("Error: Both root_index and root_name are specified. Specify only one of them.")
            root_name = self.taxon_namespace[root_index].label
        elif root_name is None:
            sys.exit("Error: You need to specify either root_index or root_name.")
        
        taxon_idx = np.where(np.array( [item.label for item in self.taxon_namespace] ) == root_name)[0][0]
        taxon = self.taxon_namespace[int(taxon_idx)]

        trees = [Phylo.read(file = StringIO(tree.as_string(schema="newick")) , format = "newick" ) for tree in self]
        
        tree = Tree_with_support.get(data = _Rstar(trees, root_name).__format__("newick"), schema="newick", taxon_namespace = self.taxon_namespace, rooting="force-rooted")
        node = dendropy.Node(taxon = taxon, edge_length = 1)
        tree.seed_node.add_child(node)
        tree.is_rooted = False
        return tree
    

    # def Qstar_consensus(self):
    #     dict_Q_count = count_quartets2(self)
    #     list_Qmaj = construct_Qmaj(dict_Q_count, self.taxon_namespace)
    #     write_Qmaj_quartfile(list_Qmaj, self.taxon_namespace)
    #     try:
    #         subprocess.run("./qstar")
    #     except:
    #         sys.exit("Error: Error running ./qstar. Make sure that qstar executable is in the same directory.")
    #     try:
    #         subprocess.run("./tree-pop")
    #     except:
    #         sys.exit("Error: Error running ./tree-pop. Make sure that qstar executable is in the same directory.")
    #     qstar_tree = Tree_with_support.get(path = "treefile", schema="newick", taxon_namespace = self.taxon_namespace)
    #     return qstar_tree
    
    # def Qstar_by_Rstar(self):
    #     # construct rstar consensus by picking first taxon as root.
    #     root_name = self.taxon_namespace[0].label
    #     taxon_idx = np.where(np.array( [item.label for item in self.taxon_namespace] ) == root_name)[0][0] # which should equal zero
    #     taxon = self.taxon_namespace[int(taxon_idx)]

    #     trees = [Phylo.read(file = StringIO(tree.as_string(schema="newick")) , format = "newick" ) for tree in self]
        
    #     tree = Tree_with_support.get(data = _Rstar(trees, root_name).__format__("newick"), schema="newick", taxon_namespace = self.taxon_namespace, rooting="force-rooted")
    #     node = dendropy.Node(taxon = taxon, edge_length = 1)
    #     tree.seed_node.add_child(node)
    #     print(len(tree.leaf_nodes()))
    #     tree.is_rooted = False
    #     tree.encode_bipartitions()

    #     # now construct majority consensus
    #     maj_tree = self.majority_consensus()
    #     maj_tree.encode_bipartitions()

    #     # store legit bipartitions
    #     ## first, mark legitimate edges with 1, unknown edges with -1
    #     bipartitions = []
    #     maj_bipar_ints = [edge.bipar.split_as_int() for edge in maj_tree.postorder_edge_iter()]
    #     for edge in tree.postorder_edge_iter():
    #         if edge.is_leaf():
    #             edge.annotations.add_new(name="in", value=1)
    #             bipartitions.append(edge.bipartition)
    #         elif edge.bipartition.split_as_int() in maj_bipar_ints:
    #             edge.annotations.add_new(name="in", value=1)
    #             bipartitions.append(edge.bipartition)
    #         else:
    #             edge.annotations.add_new(name="in", value=-1)
        
    #     ## Check if the unknown edges are legitimate
    #     for edge in tree.postorder_edge_iter():
    #         if edge.annotations['in'].value == 1:
    #             continue
    #         elif edge.annotations['in'].value == -1:
    #             # unknown, check if it's valid
    #             head = edge.head_node
    #             tail = edge.tail_node
    #             tree.reroot_at_node(tail)
    #             head_adjacents = head.child_nodes()
    #             tail_adjacents = tail.child_nodes() # this still includes head
    #             tail_adjacents = set(tail_adjacents) - {head}
    #             for s_a, s_b in itertools.combinations(head_adjacents, 2):
    #                 for s_c, s_d in itertools.combinations(tail_adjacents, 2):
    #                     for a in s_a.leaf_nodes():
    #                         for b in s_b.leaf_nodes():
    #                             for c in s_c.leaf_nodes():
    #                                 for d in s_d.leaf_nodes():
    #                                     pass
                                        

    # def Qstar_mode_consensus(self):
    #     dict_Q_count = count_quartets2(self)
    #     list_Qmaj = construct_Qmaj(dict_Q_count, self.taxon_namespace)
    #     write_Qmaj_quartfile(list_Qmaj, self.taxon_namespace)
    #     try:
    #         subprocess.run("./qstar")
    #     except:
    #         sys.exit("Error: Error running ./qstar. Make sure that qstar executable is in the same directory.")
    #     try:
    #         subprocess.run("./tree-pop")
    #     except:
    #         sys.exit("Error: Error running ./tree-pop. Make sure that qstar executable is in the same directory.")
    #     qstar_tree = Tree_with_support.get(path = "treefile", schema="newick")
    #     return qstar_tree

