
# consensus からimportしないで, newickのstringなどを入力にして関数を書いてもらえると助かります．
## 複数のサポートを可視化する,状況を可視化するにあたってconsensusからのtaxon_nameのimportがひつ
# 例
def plot_example_func(newick_str):
    # plot する関数
    pass

def get_support(Node,consensus):
    clade_bool = [False for i in range(majority.n_taxa)]
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
            clade_bit = clade_bit |get_support(child,majority)
    if clade_bit.uint in consensus.branch_support.keys():
        textface=TextFace(consensus.branch_support[clade_bit.uint])
        Node.add_face(textface,0,position="branch-top")
    print(clade_bit.bin)
    
    return clade_bit

