import numpy as np
import scanpy as sc
import pandas as pd
from tivelo.main import tivelo


if __name__ == '__main__':
    # common setting: batch_size=None, n_epochs=200, alpha_2=0.1, only_s=True

    frame_path = "D:/cuhk/project/velocity/dataset/scRNA-seq/data_frame.csv"
    edge_path = "D:/cuhk/project/velocity/dataset/scRNA-seq/cluster_edges.npy"

    data_frame = pd.read_csv(frame_path, index_col=0)
    edge_dict = np.load(edge_path, allow_pickle=True).item()

    data_name = "retina"
    data_path = data_frame.loc[data_name]["path"]
    adata = sc.read_h5ad(data_path)
    group_key = data_frame.loc[data_name]["group_key"]
    emb_key = data_frame.loc[data_name]["emb_key"]
    cluster_edges = edge_dict[data_name]

    # parameters in step 1&2
    # pancreas: t2=0.9; hindbrain2: t2=2.0; dentategyrus2: 0.39
    t1 = 0.1
    t2 = 0.9
    show_fig = False
    start_mode = "stochastic"
    rev_stat = "mean"
    adjust_DTI = False

    # parameters in step 3
    filter_genes = True
    constrain = True
    loss_fun = "mse"
    only_s = False
    alpha_1 = 1
    alpha_2 = 0.1
    batch_size = 1024
    n_epochs = 100

    velocity_key = "velocity"
    measure_performance = True

    adata = tivelo(adata, group_key, emb_key,
                   data_name=data_name,
                   save_folder="results",
                   njobs=-1,
                   t1=t1,
                   t2=t2,
                   adjust_DTI=adjust_DTI,
                   start_mode=start_mode,
                   rev_stat=rev_stat,
                   show_fig=show_fig,
                   filter_genes=filter_genes,
                   constrain=constrain,
                   loss_fun=loss_fun,
                   only_s=only_s,
                   alpha_1=alpha_1,
                   alpha_2=alpha_2,
                   batch_size=batch_size,
                   n_epochs=n_epochs,
                   velocity_key="velocity",
                   cluster_edges=cluster_edges,
                   measure_performance=measure_performance)

