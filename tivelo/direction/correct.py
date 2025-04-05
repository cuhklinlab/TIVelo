import numpy as np
from .model import orient_score, visualize
from joblib import Parallel, delayed
import scanpy as sc
from ..utils.velocity_genes import compute_velocity_genes



def recover_direction(adata, order_idx, step=1, vis_idx=None):
    u_path, s_path = adata.layers["Mu"], adata.layers["Ms"]
    u_path, s_path = u_path[order_idx, :], s_path[order_idx, :]

    dict_tree = {}
    tree_dict = {}

    r = Parallel(n_jobs=-1, verbose=10)(delayed(orient_score)(u_path[:, g].squeeze(), s_path[:, g].squeeze(), step=step)
                                        for g in range(adata.n_vars))
    score_list, _, _, u_tree, s_tree, _ = zip(*r)

    for g in range(adata.n_vars):
        tree_dict['u_tree'] = u_tree[g]
        tree_dict['s_tree'] = s_tree[g]
        dict_tree[str(g)] = tree_dict.copy()

    score_list = np.array(score_list)

    mean = np.mean(score_list)
    median = np.median(score_list)
    uq = np.quantile(score_list, 0.75)
    lq = np.quantile(score_list, 0.25)
    minimum = np.min(score_list)
    maximum = np.max(score_list)
    pos_num = np.sum(score_list > 0)

    print("mean: {:.3f}".format(mean), "\nmedian: {:.3f}".format(median), "\nlower quantile: {:.3f}".format(lq),
          "\nupper quantile: {:.3f}".format(uq), "\nminimum: {:.3f}".format(minimum), "\nmaximum: {:.3f}".format(maximum),
          "\nNo. of positive scores: {}".format(pos_num))

    if vis_idx is not None:
        gene_name = adata.var_names[vis_idx]
        score, u_arr, s_arr, u_tree, s_tree, sec_list = orient_score(u_path[:, vis_idx].squeeze(),
                                                                     s_path[:, vis_idx].squeeze())

        fig, (ax1, ax2) = visualize(u_arr, s_arr, u_tree, s_tree, sec_list, fig_size=(12, 9))
        ax1.legend(loc="upper right", ncol=2, prop={'size': 12})
        ax2.legend(loc="upper right", ncol=2, prop={'size': 12})
        ax1.set_title("{}, orientation score: {:.2f}".format(gene_name, score), fontsize=20)
        fig.tight_layout()

        return score_list, dict_tree, (ax1, ax2)

    else:
        return score_list, dict_tree


def correct_path(path_dict, adata, group_key, rev_stat="mean", tree_gene=None, root_select="connectivities",
                 emb_key=None):

    # main path
    start_node = path_dict["start_node"]
    main_path = path_dict["main_path"]

    velocity_genes = compute_velocity_genes(adata, n_top_genes=2000, r2_adjust=True, inplace=False)

    # main_path
    main_path_cell = adata.obs[group_key].isin(main_path)
    adata_new = adata[main_path_cell, :]

    start_cell_idx = np.flatnonzero(adata_new.obs[group_key] == start_node)
    success_node_list = [node for node in main_path if node != start_node]
    success_cell_idx = []
    for success_node in success_node_list:
        success_cell_idx += list(np.flatnonzero(adata_new.obs[group_key] == success_node))

    if root_select == "connectivities":
        root_cell_idx = np.flatnonzero(adata_new.obs[group_key] == start_node)[0]
    elif root_select == "coordinates" and emb_key is not None:
        coordinates = adata_new.obsm[emb_key]
        start_xy = coordinates[start_cell_idx, :]
        success_xy = coordinates[success_cell_idx, :]
        dist_xy = (np.sum(np.square(start_xy), axis=1, keepdims=True) +
                   np.sum(np.square(success_xy), axis=1, keepdims=True).T - 2 * np.dot(start_xy, success_xy.T))
        distances = np.partition(dist_xy, 30, axis=1)[:, :30]
        root_cell_idx = start_cell_idx[np.argmin(np.sum(distances, axis=1))]
    elif root_select == "terminal_score":
        root_score = adata_new.obs["root_cells"]
        start_root_score = root_score[start_cell_idx]
        root_cell_idx = start_cell_idx[np.argmax(start_root_score)]
    elif root_select == "coordinates" and emb_key is None:
        raise KeyError("Embedding key not provided")
    else:
        raise KeyError("No such method for root cell selection")

    adata_new.uns['iroot'] = root_cell_idx
    sc.tl.diffmap(adata_new)
    sc.tl.dpt(adata_new)

    adata_new = adata_new[:, velocity_genes]
    t_paga = adata_new.obs["dpt_pseudotime"].values
    order_idx = np.argsort(t_paga)

    # recover direction
    print("\nmain path:", main_path)
    if tree_gene is not None and tree_gene in adata_new.var_names:
        gene_idx = list(adata_new.var_names).index(tree_gene)
        score_list, dict_all, ax = recover_direction(adata_new, order_idx, vis_idx=gene_idx)
    else:
        score_list, dict_all = recover_direction(adata_new, order_idx)
        ax = None

    adata.uns["score_main"] = score_list.copy()

    # judge if the direction should be reversed
    if rev_stat == "median":
        if np.median(score_list) < 0:
            path_dict["start_node"] = main_path[-1]
            path_dict["main_path"] = main_path[::-1]
    elif rev_stat == "mean":
        if np.mean(score_list) < 0:
            path_dict["start_node"] = main_path[-1]
            path_dict["main_path"] = main_path[::-1]
    else:
        raise KeyError("Unknown statistics for reversing main path: {}".format(rev_stat))

    # branches
    for key, branch in path_dict.items():
        if key.startswith("branch"):
            start_node = branch[0]
            branch_cell = adata.obs[group_key].isin(branch)
            adata_new = adata[branch_cell, :]

            start_cell_idx = np.flatnonzero(adata_new.obs[group_key] == start_node)
            success_node_list = [node for node in branch if node != start_node]
            success_cell_idx = []
            for success_node in success_node_list:
                success_cell_idx += list(np.flatnonzero(adata_new.obs[group_key] == success_node))

            if root_select == "connectivities":
                root_cell_idx = np.flatnonzero(adata_new.obs[group_key] == start_node)[0]
            elif root_select == "coordinates" and emb_key is not None:
                coordinates = adata_new.obsm[emb_key]
                start_xy = coordinates[start_cell_idx, :]
                success_xy = coordinates[success_cell_idx, :]
                dist_xy = (np.sum(np.square(start_xy), axis=1, keepdims=True) +
                           np.sum(np.square(success_xy), axis=1, keepdims=True).T - 2 * np.dot(start_xy, success_xy.T))
                distances = np.partition(dist_xy, 30, axis=1)[:, :30]
                root_cell_idx = start_cell_idx[np.argmin(np.sum(distances, axis=1))]
            elif root_select == "terminal_score":
                root_score = adata_new.obs["root_cells"]
                start_root_score = root_score[start_cell_idx]
                root_cell_idx = start_cell_idx[np.argmax(start_root_score)]
            elif root_select == "coordinates" and emb_key is None:
                raise KeyError("Embedding key not provided")
            else:
                raise KeyError("No such method for root cell selection")

            adata_new.uns['iroot'] = root_cell_idx
            sc.tl.diffmap(adata_new)
            sc.tl.dpt(adata_new)

            adata_new = adata_new[:, velocity_genes]
            t_paga = adata_new.obs["dpt_pseudotime"].values
            order_idx = np.argsort(t_paga)

            print("\n" + key + ":", branch)
            score_list, dict_all = recover_direction(adata_new, order_idx)
            adata.uns["score_" + key] = score_list

            if np.median(score_list) < 0 and np.mean(score_list) < 0:
                path_dict[key] = branch[::-1]

    # further correct the result according to overall result
    main_path = path_dict["main_path"]
    del_key = []

    for key, obj in path_dict.items():
        if key.startswith("branch"):
            if path_dict[key][-1] == main_path[0]:
                path_dict["main_path"] = path_dict[key] + main_path[1:]
                path_dict["start_node"] = path_dict[key][0]
                del_key.append(key)
            elif path_dict[key][0] == main_path[-1]:
                path_dict["main_path"] = main_path + path_dict[key][1:]
                path_dict["start_node"] = main_path[0]
                del_key.append(key)
            elif path_dict[key][0] != main_path[-1] and path_dict[key][-1] != main_path[0]:
                if path_dict[key][-1] in main_path:
                    path_dict[key] = path_dict[key][::-1]

    for branch in del_key:
        path_dict.pop(branch)

    adata.uns["path_dict"] = path_dict.copy()
    print("\n'path_dict' added to adata.uns")
    return path_dict, ax



