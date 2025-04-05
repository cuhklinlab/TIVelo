import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from .start_node import get_start_node
from .prune import graph_prune
from .main_path import get_path, visualize_path


def process_path(adata, group_key, emb_key, njobs=-1, start_mode="stochastic"):
    '''
    :param adata: anndata to use
    :param group_key: cell cluster key to use
    :param emb_key: basis/embedding key to use
    :param filter_subgraph: whether to filter subgraph
    :param njobs: number of parallel jobs
    :return: dictionary of selected main path including start node and following nodes
    '''

    select_info, select_weight = graph_prune(adata, group_key, emb_key)
    start_node = get_start_node(adata, select_info, group_key, njobs=njobs, mode=start_mode)
    start_index = list(select_info.index).index(start_node)

    path_dict = get_path(select_info, select_weight, start_index=start_index)
    ax = visualize_path(select_info, select_weight, path_dict, emb_key, show=False, ax=None)

    adata.uns["path_dict"] = path_dict.copy()
    print("'path_dict' added to adata.uns")

    return path_dict, ax


if __name__ == '__main__':
    # test
    data_path = r"D:\cuhk\project\velocity\dataset\scRNA-seq\BoneMarrow\human_cd34_bone_marrow_processed.h5ad"
    data_name = "bonemarrow"
    group_key = "clusters"
    emb_key = "X_tsne"

    adata = sc.read_h5ad(data_path)
    path_dict, ax = process_path(adata, group_key, emb_key, njobs=-1, start_mode="stochastic")
    print(path_dict)
    print(adata.uns["path_dict"])
    plt.show()










