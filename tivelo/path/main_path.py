import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from .prune import graph_prune


def get_child_index(select_weight, parent_index):
    # given a parent node index i, return child(i)
    child_index = []
    for index in parent_index:
        index_ = list(np.where(select_weight[index, :] > 0)[0])
        child_index += index_
    return list(set(child_index))


def get_child_path(select_weight, parent_index, root_index):
    # given root node index i and its one child node index i' (parent_index), return successor(i)
    child_index_ = get_child_index(select_weight, parent_index)
    child_index_ = list(set(child_index_) - set(root_index))
    child_path = child_index_
    while len(child_index_) > 0:
        child_index = child_index_
        root_index += child_index
        child_index_ = get_child_index(select_weight, child_index)
        child_index_ = list(set(child_index_) - set(root_index))
        child_path += child_index_
    return child_path


def path_selection(select_info, select_weight, start_index, root_list=None):
    '''
    :param select_info: dataframe information for selected clusters and edges
    :param select_weight: weight matrix for selected clusters and edges
    :param start_index: start node index
    :param root_list: list of root cluster excluded from selection
    :return: path of nodes beginning from start node excluding nodes in root list, with most number of cells
    '''
    select_path = [start_index]
    root_index = [start_index] + root_list if root_list is not None else [start_index]
    split_index = list(np.where(select_weight[start_index, :] > 0)[0])

    while len(split_index) > 0:
        ncells_list = []
        split_index = list(np.where(select_weight[start_index, :] > 0)[0])
        split_index = list(set(split_index) - set(root_index))
        if len(split_index) == 0:
            break

        for node_index in split_index:
            child_path = get_child_path(select_weight, [node_index], [start_index, node_index])
            ncells = select_info["ncells"][[node_index] + child_path].sum()
            ncells_list.append(ncells)

        select_index = split_index[np.argmax(ncells_list)]
        select_path += [select_index]
        start_index = select_index
        root_index += [select_index]

    return select_path


def get_path(select_info, select_weight, start_index):
    '''
    :param select_info: dataframe information for selected clusters and edges
    :param select_weight: weight matrix for selected clusters and edges
    :param start_index: start node index
    :return: dictionary of selected main path including start node and following nodes
    '''
    path_dict = {}
    root_list = []

    select_path = path_selection(select_info, select_weight, start_index, root_list=None)
    main_path = list(select_info.index[select_path].values)
    path_dict["start_node"] = select_info.index[start_index]
    path_dict["main_path"] = main_path

    root_list += select_path
    flag = 1
    for index in select_path:
        branch_index = list(set(get_child_index(select_weight, [index])) - set(root_list))
        if len(branch_index) > 0:
            for index_ in branch_index:
                select_path = get_child_path(select_weight, [index_], root_list + [index_])
                select_path = [index] + [index_] + select_path
                path_dict["branch_{}".format(flag)] = list(select_info.index[select_path].values)
                flag += 1

    return path_dict


def visualize_path(select_info, select_weight, path_dict, emb_key="X_umap", show=True, ax=None):
    flag1, flag2 = 0, 0

    # special adjustment for organoid:
    if "Enteroendocrine progenitor" in list(select_info.index):
        stem_idx = list(select_info.index).index("Stem cells")
        select_info["{}1".format(emb_key)][stem_idx] = select_info["{}1".format(emb_key)][stem_idx] + 0.5

        ta_idx = list(select_info.index).index("TA cells")
        select_info["{}1".format(emb_key)][ta_idx] = select_info["{}1".format(emb_key)][ta_idx] + 0.5

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    for i in range(select_weight.shape[0]):
        for j in range(select_weight.shape[1]):
            if select_weight[i, j] > 0:
                x = select_info["{}0".format(emb_key)][i], select_info["{}0".format(emb_key)][j]
                y = select_info["{}1".format(emb_key)][i], select_info["{}1".format(emb_key)][j]
                if select_info.index[i] in path_dict["main_path"] and select_info.index[j] in path_dict["main_path"]:
                    if flag1 == 0:
                        ax.plot(x, y, lw=5 * select_weight[i, j], c="red", zorder=0, label="Main path")
                        flag1 += 1
                    else:
                        ax.plot(x, y, lw=5 * select_weight[i, j], c="red", zorder=0)
                else:
                    if flag2 == 0:
                        ax.plot(x, y, lw=5 * select_weight[i, j], c="green", zorder=0, label="Branch")
                        flag2 += 1
                    else:
                        ax.plot(x, y, lw=5 * select_weight[i, j], c="green", zorder=0)

    start_index = list(select_info.index).index(path_dict["start_node"])
    remain = select_info.drop(path_dict["start_node"], axis=0)
    ax.scatter(remain["{}0".format(emb_key)], remain["{}1".format(emb_key)], s=remain['sizes'], c="tab:blue", zorder=1,
               )
    ax.scatter(select_info["{}0".format(emb_key)][start_index], select_info["{}1".format(emb_key)][start_index],
               s=select_info['sizes'][start_index], c="orange", zorder=1, label="Origin node")

    for count, cluster_name in enumerate(select_info.index):
        ax.text(select_info["{}0".format(emb_key)][count], select_info["{}1".format(emb_key)][count], cluster_name,
                ha="center", zorder=2, fontsize=12)
    if show:
        plt.show()
    else:
        return ax


if __name__ == '__main__':
    # test
    data_path = r"D:\cuhk\project\velocity\dataset\scRNA-seq\BoneMarrow\human_cd34_bone_marrow_processed.h5ad"
    group_key = "clusters"
    emb_key = "X_tsne"

    adata = sc.read_h5ad(data_path)
    select_info, select_weight = graph_prune(adata, group_key, emb_key)

    path_dict = get_path(select_info, select_weight, start_index=5)
    print(path_dict)
    visualize_path(select_info, select_weight, path_dict, emb_key, show=True)




