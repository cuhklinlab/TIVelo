import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import networkx as nx


def graph_prune(adata, group_key, emb_key, ncells_thre=0.5, weight_thre=0.6):
    '''
    :param adata: anndata to use
    :param group_key: cell cluster key to use
    :param emb_key: basis/embedding key to use
    :param ncells_thre: threshold of percent for number of cells
    :param weight_thre: threshold of weight for cluster edges which are kept
    :param split_weight: threshold of weight for cluster edges which are deleted
    :return: dataframe information for selected clusters and edges; weight matrix for selected clusters and edges
    '''
    sc.tl.diffmap(adata)
    sc.tl.paga(adata, groups=group_key)

    embedding = adata.obsm[emb_key]
    adata.obs["{}0".format(emb_key)] = embedding[:, 0]
    adata.obs["{}1".format(emb_key)] = embedding[:, 1]
    info = pd.concat([adata.obs.groupby(group_key)["{}0".format(emb_key)].mean(),
                      adata.obs.groupby(group_key)["{}1".format(emb_key)].mean(),
                      pd.Series(adata.obs[group_key].value_counts(), name="ncells")], axis=1)

    groups_sizes = adata.obs[group_key].value_counts()
    median_group_size = np.median(groups_sizes)
    groups_sizes = pd.Series(100 * np.power(groups_sizes / median_group_size, 0.5), name="sizes")
    info = pd.concat([info, groups_sizes], axis=1)

    try:
        colors = pd.Series(adata.uns[group_key + "_colors"], index=info.index, name="colors")
        info = pd.concat([info, colors], axis=1)
    except KeyError:
        pass

    # weight matrix from PAGA
    weight = adata.uns["paga"]["connectivities_tree"].A
    weight = weight + weight.T

    select_num = 0
    for i in range(len(info.index)):
        if info.sort_values(by=["ncells"], ascending=False).iloc[:i + 1, 2].sum() > ncells_thre * info["ncells"].sum():
            select_num = i + 1
            break

    if select_num < 5:
        select_num = min(5, len(info.index))
    select_clusters = list(info.sort_values(by=["ncells"], ascending=False).index[:select_num])

    # if a single cluster has over 10% of the total number of cells
    all_clusters = list(info.index)
    for cluster_name in all_clusters:
        pos = all_clusters.index(cluster_name)
        if info.iloc[pos, 2] > 0.1 * info["ncells"].sum():
            if cluster_name not in select_clusters:
                select_clusters.append(cluster_name)

    # if a cluster is connected to clusters just selected with weight over 0.6, also select it
    for cluster_name in select_clusters:
        pos = all_clusters.index(cluster_name)
        indices = np.where(weight[:, pos] > weight_thre * np.max(weight))[0]
        if len(indices) > 0:
            for idx in indices:
                cluster_new = all_clusters[idx]
                if cluster_new not in select_clusters:
                    select_clusters.append(cluster_new)

    select_index = info.index.isin(select_clusters)
    select_info = info.iloc[select_index, :]
    select_weight = weight[select_index, :][:, select_index]

    # connect two or more sub-graphs together
    n_clusters = select_info.shape[0]
    point_list = [i for i in range(n_clusters)]
    link_list = [(i, j) for i in range(n_clusters) for j in range(n_clusters) if select_weight[i, j] > 0 and i < j]
    G = nx.Graph()
    for node in point_list:
        G.add_node(node)
    for link in link_list:
        G.add_edge(link[0], link[1])

    ncells_list = []
    for c in nx.connected_components(G):
        ncells_list.append(select_info.iloc[list(c)]["ncells"].sum())

    sorted_graph = [list(nx.connected_components(G))[i] for i in list(np.argsort(-np.array(ncells_list)))]

    if len(sorted_graph) > 1:
        for p in range(len(sorted_graph)):
            for p_ in range(p + 1, len(sorted_graph)):
                c, c_ = sorted_graph[p], sorted_graph[p_]
                path_list, path_list_ = list(c), list(c_)

                path_list = [list(select_info.index)[i] for i in path_list]
                path_list = [list(info.index).index(i) for i in path_list]
                path_list_ = [list(select_info.index)[i] for i in path_list_]
                path_list_ = [list(info.index).index(i) for i in path_list_]

                path_list.sort()
                path_list_.sort()
                weight_dict = {}

                # select connected edge between c and c_
                for j in path_list:
                    for j_ in path_list_:
                        w = weight[j, j_]
                        weight_dict[(j, j_)] = w

                if max(weight_dict.values()) > 0:
                    for key, value in weight_dict.items():
                        if value == max(weight_dict.values()):
                            j, j_ = key
                            j = list(select_info.index).index(info.index[j])
                            j_ = list(select_info.index).index(info.index[j_])
                            select_weight[j, j_] = value

                # add another node to connect c and c_
                else:
                    not_selected = [i for i in list(info.index) if i not in select_clusters]
                    weight_dict = {}
                    for i in not_selected:
                        idx = list(info.index).index(i)
                        weight_idx = []
                        for j in path_list:
                            for j_ in path_list_:
                                w_c, w_c_ = weight[idx, j], weight[idx, j_]
                                weight_idx.append(w_c * w_c_)
                        weight_idx = np.max(weight_idx)
                        weight_dict[i] = weight_idx

                    if max(weight_dict.values()) > 0:
                        for key, value in weight_dict.items():
                            if value == max(weight_dict.values()):
                                select_clusters.append(key)
                                select_index = info.index.isin(select_clusters)
                                select_info = info.iloc[select_index, :]
                                select_weight = weight[select_index, :][:, select_index]

                    else:
                        raise ValueError("There are more then one unconnected sub-graphs in the dataset."
                                         "Please split the dataset and run TIVelo for them independently.")

    return select_info, select_weight


def visualize(select_info, select_weight, show=False):
    fig, ax = plt.subplots(1, 1)
    for i in range(select_weight.shape[0]):
        for j in range(select_weight.shape[1]):
            if select_weight[i, j] > 0:
                x = select_info["{}0".format(emb_key)][i], select_info["{}0".format(emb_key)][j]
                y = select_info["{}1".format(emb_key)][i], select_info["{}1".format(emb_key)][j]
                ax.plot(x, y, lw=5 * select_weight[i, j], c="black", zorder=0)

    try:
        ax.scatter(select_info["{}0".format(emb_key)], select_info["{}1".format(emb_key)], s=select_info['sizes'],
                   c=select_info["colors"], zorder=1)
    except KeyError:
        ax.scatter(select_info["{}0".format(emb_key)], select_info["{}1".format(emb_key)], s=select_info['sizes'],
                   zorder=1)

    for count, cluster_name in enumerate(select_info.index):
        ax.text(select_info["{}0".format(emb_key)][count], select_info["{}1".format(emb_key)][count], cluster_name,
                ha="center", zorder=2)
    if show:
        plt.show()
    else:
        return ax


if __name__ == '__main__':
    # test
    data_dir = r"D:\cuhk\project\velocity\dataset\atac\HSPCs\3423-MV-2_adata_postpro.h5ad"
    group_key = "leiden"
    emb_key = "X_umap"

    adata_rna = sc.read_h5ad(data_dir)
    select_info, select_weight = graph_prune(adata_rna, group_key, emb_key)
    ax = visualize(select_info, select_weight, show=False)
    plt.show()