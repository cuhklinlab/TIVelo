import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvelo as scv
import pandas as pd
import warnings
warnings.filterwarnings("ignore")


def get_start_node(adata, select_info, group_key, njobs=1, mode="stochastic"):
    '''
    :param adata: anndata to use
    :param select_info: dataframe information for selected clusters and edges
    :param group_key: cell cluster key to use
    :param njobs: number of parallel jobs
    :param mode: start node finding mode
    :return: string for start node cluster
    '''
    # run scvelo and get root/end score
    if "velocity_graph" not in adata.uns.keys():
        if mode == "stochastic":
            scv.tl.velocity(adata, mode='stochastic')
        elif mode == "dynamical":
            scv.tl.recover_dynamics(adata, n_jobs=njobs)
            scv.tl.velocity(adata, mode="dynamical")
        else:
            raise KeyError("mode must be stochastic or dynamical")

    scv.tl.velocity_graph(adata, n_jobs=njobs)
    scv.tl.terminal_states(adata)

    select_clusters = select_info.index.values
    select_index = adata.obs[group_key].isin(select_clusters)
    adata = adata[select_index, :]

    root_score = np.array(adata.obs.groupby(group_key)["root_cells"].mean())
    end_score = np.array(adata.obs.groupby(group_key)["end_points"].mean())

    root_order = abs(np.sort(-root_score))
    end_order = abs(np.sort(-end_score))

    if np.max(root_score) > 0.1:
        flag = 0
        root_index = np.where(root_score == root_order[flag])[0][0]
        while flag < len(root_score) - 1:
            if end_score[root_index] > root_score[root_index] or root_score[root_index] <= 0.1:
                flag += 1
                root_index = np.where(root_score == root_order[flag])[0][0]
            else:
                break

        if end_score[root_index] <= root_score[root_index] and root_score[root_index] > 0.1:
            start_node = select_info.index[root_index]
            return start_node

        elif np.max(end_score) < 0.1:
            root_index = np.where(root_score == root_order[0])[0][0]
            start_node = select_info.index[root_index]
            return start_node

        else:
            flag = 0
            end_index = np.where(end_score == end_order[flag])[0][0]
            while flag < len(end_score) - 1:
                if root_score[end_index] > end_score[end_index] or end_score[end_index] <= 0.1:
                    flag += 1
                    end_index = np.where(end_score == end_order[flag])[0][0]
                else:
                    break
            if root_score[end_index] <= end_score[end_index] and end_score[end_index] > 0.1:
                start_node = select_info.index[end_index]
                return start_node
            else:
                root_index = np.where(root_score == root_order[0])[0][0]
                start_node = select_info.index[root_index]
                return start_node
    else:
        flag = 0
        end_index = np.where(end_score == end_order[flag])[0][0]
        while flag < len(end_score) - 1:
            if root_score[end_index] > end_score[end_index] or end_score[end_index] <= 0.1:
                flag += 1
                end_index = np.where(end_score == end_order[flag])[0][0]
            else:
                break
        if root_score[end_index] <= end_score[end_index] and end_score[end_index] > 0.1:
            start_node = select_info.index[end_index]
            return start_node
        else:
            root_index = np.where(root_score == root_order[0])[0][0]
            start_node = select_info.index[root_index]
            return start_node


if __name__ == '__main__':
    # test
    frame_path = "D:/cuhk/project/velocity/dataset/scRNA-seq/data_frame.csv"
    data_frame = pd.read_csv(frame_path, index_col=0)

    for data_name in data_frame.index:
        data_path = data_frame.loc[data_name]["path"]
        group_key = data_frame.loc[data_name]["group_key"]
        adata = sc.read_h5ad(data_path)

        scv.tl.velocity(adata, mode='stochastic')
        scv.tl.velocity_graph(adata, n_jobs=20)

        scv.tl.terminal_states(adata)
        ax = scv.pl.scatter(adata, color=['root_cells', 'end_points'], show=False)
        plt.show()

        root_score = adata.obs.groupby(group_key)["root_cells"].mean()
        print("\n", data_name, "\nroot score:\n", root_score)

        terminal_score = adata.obs.groupby(group_key)["end_points"].mean()
        print("\nterminal score:\n", terminal_score)






