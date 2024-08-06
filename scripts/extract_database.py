import dendropy
import random
import argparse

def extract_subtrees_from_file(file_path, n, seed):
    # タクソンラベル S1~S100 を生成
    all_taxa_labels = [f"S{i}" for i in range(1, 101)]
    
    # シード値を設定して乱数生成器を初期化
    random.seed(seed)
    
    # ランダムに n 個のタクソンラベルを選択
    if n > len(all_taxa_labels):
        raise ValueError(f"n cannot be greater than the total number of taxa labels ({len(all_taxa_labels)}).")
    
    selected_taxa_labels = random.sample(all_taxa_labels, n)
    
    # Newick形式の系統樹ファイルを行ごとに読み込む
    with open(file_path, "r") as file:
        trees = file.readlines()
    
    # 各行の木に対してサブツリーを抽出し、標準出力に書き込む
    for i, tree_str in enumerate(trees, 1):
        tree = dendropy.Tree.get(data=tree_str.strip(), schema="newick")
        subtree = tree.extract_tree_with_taxa_labels(selected_taxa_labels)
        
        # サブツリーを標準出力に書き込む
        print(subtree.as_string(schema="newick"),end="")

def main():
    parser = argparse.ArgumentParser(description="Extract subtrees from a file containing multiple Newick trees.")
    parser.add_argument("file_path", type=str, help="The path to the Newick tree file.")
    parser.add_argument("n", type=int, help="The number of taxa labels to randomly select (1-100).")
    parser.add_argument("--seed", type=int, default=42, help="The seed for the random number generator (default: 42).")
    
    args = parser.parse_args()
    
    # ファイルからサブツリーを抽出して、標準出力に表示
    extract_subtrees_from_file(args.file_path, args.n, args.seed)

if __name__ == "__main__":
    main()
