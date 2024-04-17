
# consensus からimportしないで, newickのstringなどを入力にして関数を書いてもらえると助かります．
## 複数のサポートを可視化する,状況を可視化するにあたってconsensusからのtaxon_nameのimport
# 例
import ete3
from ete3 import TextFace
import PyQt5
from bitstring import Bits
import numpy as np
from ._consensus import Tree_with_support

def plot_example_func(newick_str):
    # plot する関数
    pass

#ete3.Tree, Tree_with_support.namespace, Tree_with_support.support,  int pos, bool leaf_support
# ete3の木クラス, Tree_with_support.namespace, supportのhash table, Nodeのどこに記述するかのposition, leaf branch にsupportを書くかどうか
def get_support(Node,namespace,support_hashtable,pos = 0,leaf_support = True):
    #https://viscid-hub.github.io/Viscid-docs/docs/dev/styles/tableau-colorblind10.html よりカラーコードを採用
    color = ["#006BA4", "#FF800E", "#ABABAB", "#595959",
                 "#5F9ED1", "#C85200", "#898989", "#A2C8EC", "#FFBC79", "#CFCFCF"]
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

