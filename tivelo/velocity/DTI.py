import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def get_child_dict(adata, group_key, path_dict, threshold=0.1, threshold_trans=1.0, adjust=False):
    '''
    :param adata: anndata to use
    :param group_key: cell cluster key to use
    :param path_dict: path dictionary from recover direction
    :param threshold: threshold for weight of edges
    :param threshold_trans: threshold for product of weights of edges
    :return: child node dictionary for each cluster node, level dictionary for each cluster node
    :return: if conduct adjustment for t_2 in DTI
    '''

    # rerun paga and get full tree
    sc.tl.diffmap(adata)
    sc.tl.paga(adata, groups=group_key)

    weight = adata.uns["paga"]["connectivities"].A
    cluster_name = list(adata.obs[group_key].values.categories)
    start_node = path_dict["start_node"]

    # create level dictionary for data
    thres_list = []
    threshold = threshold * np.max(weight)
    level_dict = {}
    level = 0
    parent_nodes = []
    child_nodes = [start_node]
    child_nodes_ = []
    while len(child_nodes) > 0:
        for node in child_nodes:
            level_dict[node] = level
            i = cluster_name.index(node)
            for i_ in np.where(weight[i, :] > threshold)[0]:
                node_ = cluster_name[i_]
                weight_trans = max([weight[i, j] * weight[j, i_] for j in range(len(cluster_name))])
                thres_list.append(weight[i, i_] / weight_trans)
                if node_ not in parent_nodes + child_nodes and threshold_trans * weight_trans <= weight[i, i_]:
                    child_nodes_.append(node_)

        parent_nodes = parent_nodes + child_nodes
        child_nodes = child_nodes_.copy()
        child_nodes_ = []
        level += 1
        thres_list.sort()

    # get child nodes for every node
    child_dict = {}
    for node in cluster_name:
        child_list = []
        for i in np.where(weight[cluster_name.index(node), :] > threshold)[0]:
            node_ = cluster_name[i]
            try:
                if level_dict[node] == level_dict[node_] - 1:
                    child_list.append(node_)
            except KeyError:
                pass
        child_dict[node] = child_list

    for node in cluster_name:
        if node != start_node and not np.any([node in child_dict[node_] for node_ in cluster_name]):
            child_dict[cluster_name[np.argmax(weight[cluster_name.index(node), :])]] += [node]

    for node in cluster_name:
        parent_nodes = [node_ for node_ in cluster_name if node in child_dict[node_]]
        if len(parent_nodes) > 1:
            i = cluster_name.index(node)
            parent_index = [cluster_name.index(node_) for node_ in parent_nodes]
            weight_node_ = [weight[i, j] for j in parent_index]
            node_ = cluster_name[parent_index[np.argmax(weight_node_)]]
            for node__ in cluster_name:
                if node__ != node_ and node in child_dict[node__]:
                    child_dict[node__].remove(node)

    if not adjust:
        adata.uns["child_dict"] = child_dict.copy()
        adata.uns["level_dict"] = level_dict.copy()
        adata.uns["threshold_list"] = thres_list.copy()
        print("'child_dict' added to adata.uns")
        print("'level_dict' added to adata.uns")
        print("'threshold_list' added to adata.uns")
    return child_dict, level_dict, thres_list


# create directed nearest neighborhood
def get_d_nn(adata, group_key, child_dict, start_node, emb_key=None, root_select="connectivities"):

    start_cell_idx = np.flatnonzero(adata.obs[group_key] == start_node)
    child_cell_idx = []
    for child_node in child_dict[start_node]:
        child_cell_idx += list(np.flatnonzero(adata.obs[group_key] == child_node))

    if root_select == "connectivities":
        root_cell_idx = np.flatnonzero(adata.obs[group_key] == start_node)[0]
    elif root_select == "coordinates" and emb_key is not None:
        coordinates = adata.obsm[emb_key]
        start_xy = coordinates[start_cell_idx, :]
        child_xy = coordinates[child_cell_idx, :]
        dist_xy = np.sum(np.square(start_xy), axis=1, keepdims=True) + np.sum(np.square(child_xy), axis=1, keepdims=True).T \
                  - 2 * np.dot(start_xy, child_xy.T)
        root_cell_idx = start_cell_idx[np.argmin(np.sum(dist_xy, axis=1))]
    elif root_select == "terminal_score":
        root_score = adata.obs["root_cells"]
        start_root_score = root_score[start_cell_idx]
        root_cell_idx = start_cell_idx[np.argmax(start_root_score)]
    elif root_select == "coordinates" and emb_key is None:
        raise KeyError("Embedding key not provided")
    else:
        raise KeyError("No such method for root cell selection")

    adata.uns['iroot'] = root_cell_idx
    sc.tl.dpt(adata)
    dpt_pseudotime = adata.obs["dpt_pseudotime"].values
    cluster_label = np.array(adata.obs[group_key].values)
    cluster_name = list(adata.obs[group_key].values.categories)

    knn = (adata.obsp["distances"].A != 0).astype(float)
    child_c, same_c = np.zeros_like(knn), np.zeros_like(knn)
    for cluster in cluster_name:
        child_node = child_dict[cluster]
        if_child = np.array([j in child_node for j in cluster_label])
        if_same = np.array([j in [cluster] for j in cluster_label])
        child_c[cluster_label == cluster, :] = if_child
        same_c[cluster_label == cluster, :] = if_same

    d_nn = knn * (child_c + same_c) * (dpt_pseudotime[:, None] <= dpt_pseudotime)
    adata.obsp["d_nn"] = d_nn.copy()
    print("'d_nn' added to adata.obsp")

    return d_nn, knn, child_c, same_c


def directed_graph(adata, child_dict, group_key, start_node, emb_key="X_umap", data_name=None, show=False, ax=None):
    # visualize the directed graph
    embedding = adata.obsm[emb_key]
    adata.obs["{}0".format(emb_key)] = embedding[:, 0]
    adata.obs["{}1".format(emb_key)] = embedding[:, 1]
    info = pd.concat([adata.obs.groupby(group_key)["{}0".format(emb_key)].mean(),
                      adata.obs.groupby(group_key)["{}1".format(emb_key)].mean(),
                      pd.Series(adata.obs[group_key].value_counts(), name="ncells")], axis=1)

    groups_sizes = adata.obs[group_key].value_counts()
    median_group_size = np.median(groups_sizes)
    groups_sizes = pd.Series(120 * np.power(groups_sizes / median_group_size, 0.5), name="sizes")
    info = pd.concat([info, groups_sizes], axis=1)

    try:
        colors = pd.Series(adata.uns[group_key + "_colors"], index=info.index, name="colors")
        info = pd.concat([info, colors], axis=1)
    except KeyError:
        pass

    weight = adata.uns["paga"]["connectivities"].A
    weight = weight + weight.T
    cluster_name = list(adata.obs[group_key].values.categories)

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    start_index = list(info.index).index(start_node)
    remain = info.drop(start_node, axis=0)
    ax.scatter(remain["{}0".format(emb_key)], remain["{}1".format(emb_key)], s=remain['sizes'], c="tab:blue",
               zorder=1)
    ax.scatter(info["{}0".format(emb_key)][start_index], info["{}1".format(emb_key)][start_index],
               s=info['sizes'][start_index], c="orange", zorder=0)

    for i in range(weight.shape[0]):
        for j in range(weight.shape[1]):
            if cluster_name[j] in child_dict[cluster_name[i]]:
                x, y = info["{}0".format(emb_key)][i], info["{}1".format(emb_key)][i]
                u, v = info["{}0".format(emb_key)][j], info["{}1".format(emb_key)][j]
                ax.quiver(x, y, u-x, v-y, angles='xy', scale_units='xy', scale=1, zorder=1, color="gray")

    for count, cluster_name in enumerate(info.index):
        ax.text(info["{}0".format(emb_key)][count], info["{}1".format(emb_key)][count], cluster_name,
                ha="left", zorder=2)

    ax.set_xticks([])
    ax.set_yticks([])

    if data_name is not None:
        plt.title(data_name)
    if show:
        plt.show()

    return ax







