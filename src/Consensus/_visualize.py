
# consensus からimportしないで, newickのstringなどを入力にして関数を書いてもらえると助かります．
## 複数のサポートを可視化する,状況を可視化するにあたってconsensusからのtaxon_nameのimport
# 例
import ete3
import sys
from ete3 import TextFace
import dendropy
import PyQt5
from bitstring import Bits
import numpy as np
from ._consensus import Tree_with_support

def plot_example_func(Tree_with_support):
    """Conversion function for drawing Tree_with_support with ete3

         Parameters
        ----------
        Tree_with_support: consensus.Tree_with_support
    
        Returns
        -------
        t: ete3.Tree()
        rs: ete3.TreeStyle()
    """

    t = ete3.Tree(Tree_with_support.as_string(schema='newick',suppress_rooting=True))
    ts=ete3.TreeStyle()
    color = ["#006BA4", "#FF800E", "#ABABAB", "#595959",
                 "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]
    if Tree_with_support.branch_support is not None:
        #print("compute branch_support")
        _ = get_support(t,Tree_with_support.taxon_namespace,Tree_with_support.branch_support,pos = 0,leaf_support = False)
        ts.legend.add_face(ete3.TextFace("branch_support",fgcolor=color[0]), column=0)
    if Tree_with_support.transfer_support is not None:
        #print("compute transfer_support")
        _ = get_support(t,Tree_with_support.taxon_namespace,Tree_with_support.transfer_support,pos = 1,leaf_support = False)
        ts.legend.add_face(ete3.TextFace("transfer_support",fgcolor=color[1]), column=0)

    return t,ts
    

#ete3.Tree, Tree_with_support.namespace, Tree_with_support.support,  int pos, bool leaf_support
# ete3の木クラス, Tree_with_support.namespace, supportのhash table, Nodeのどこに記述するかのposition, leaf branch にsupportを書くかどうか
def get_support(Node,namespace,support_hashtable,pos = 0,leaf_support = True):
    """A function to add supports(branch support, transfer support) from The Tree of Tree_with_support　class to Nodes in ete3 for visualization.
    
        Parameters
        ---------
        Node: ete3.TreeNode
        namespace: A list of pecies of `dendropy`
        support_hashtable: dict(), key: clade_bit value: support_score
        pos:int, Position information for adding support to a Node using ete3.TreeNode.add_face
        leaf_support: bool
    
        Returns
        -------
        clade_bit: A bit string representing the bipartition of leaf species.
    """
    #https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/tableau-colorblind10.html よりカラーコードを採用
    color = ["#006BA4", "#FF800E", "#ABABAB", "#595959", "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]
    taxonnames_array = np.array([item.label for item in namespace])
    #clade_bool = [False for i in range(consensus.n_taxa)]
    clade_bool = [False for i in range(len(taxonnames_array))]
    if(Node.is_leaf()):
        digits = [np.where( Node.name == taxonnames_array )][0][0][0]
        clade_bool[-(digits+1)] = True
        clade_bit = Bits(clade_bool)
        if int(clade_bit.bin[-1]) == 1:
            clade_bit = (~clade_bit)
    else:
        clade_bit = Bits(clade_bool)
        if int(clade_bit.bin[-1]) == 1:
            clade_bit = (~clade_bit)
        for child in Node.children:
            clade_bit = clade_bit |get_support(child,namespace,support_hashtable,pos,leaf_support)
    if Node.is_leaf() == False or leaf_support == True: 
        if clade_bit.uint in support_hashtable.keys():
            textface=TextFace("{:.3f}".format(support_hashtable[clade_bit.uint]),fgcolor=color[pos])
            if(pos%2 == 0):
                Node.add_face(textface,pos//2,position="branch-top")
            else:
                Node.add_face(textface,pos//2,position="branch-bottom")
    
    return clade_bit

def read_consensus_NeXML(NeXML_path):
    newtree = dendropy.Tree.get(path=NeXML_path,schema="nexml")
    consensus = Tree_with_support(newtree)
    consensus.branch_support = get_support_from_NeXML(newtree,"branch_support")
    consensus.transfer_support = get_support_from_NeXML(newtree,"transfer_support")
    return consensus


def write_consensus_NeXML(Tree_with_support,NexML_path):
    if(Tree_with_support.branch_support != None):
        for edge in Tree_with_support.postorder_edge_iter():
            edge.branch_support = None
            edge.annotations.add_bound_attribute("branch_support")
            edge.branch_support = Tree_with_support.branch_support[int(edge.bipartition)]

    if(Tree_with_support.transfer_support != None):
        for edge in Tree_with_support.postorder_edge_iter():
            edge.transfer_support = None
            edge.annotations.add_bound_attribute("transfer_support")
            edge.transfer_support = Tree_with_support.transfer_support[int(edge.bipartition)]

    Tree_with_support.write(
        path= NexML_path,
        schema='nexml',
        ignore_unrecognized_keyword_arguments=False,
        )

def get_support_from_NeXML(dendropy_Tree_from_NeXML,support_name):
    if(support_name != "branch_support" and support_name != "transfer_support"):
        print("support_name is not correct, please select 'branch_support' or 'transfer_support'")
        sys.exit([1])

    support_added = dict()
    dendropy_Tree_from_NeXML.encode_bipartitions()
    #あんま綺麗じゃないけど...
    if(support_name == "branch_support"): 
        for edge in dendropy_Tree_from_NeXML.postorder_edge_iter():
            if edge.annotations.find(name="branch_support") == None: continue
            support_added[int(edge.bipartition)] = float(str(edge.annotations.find(name="branch_support")).strip("branch_support=").strip("'"))
    if(support_name == "transfer_support"):
        for edge in dendropy_Tree_from_NeXML.postorder_edge_iter():
            if edge.annotations.find(name="transfer_support") == None: continue
            support_added[int(edge.bipartition)] = float(str(edge.annotations.find(name="transfer_support")).strip("transfer_support=").strip("'"))

    return support_added
    
    

