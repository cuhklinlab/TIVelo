import scanpy as sc
import pandas as pd
import numpy as np
from tivelo.main import tivelo

if __name__ == '__main__':
    frame_path = "D:/cuhk/project/velocity/dataset/atac/data_frame.csv"
    data_frame = pd.read_csv(frame_path, index_col=0)

    data_name = "HSPCs"
    rna_path = data_frame.loc[data_name]["rna"]
    atac_path = data_frame.loc[data_name]["atac"]
    group_key = data_frame.loc[data_name]["group_key"]
    emb_key = data_frame.loc[data_name]["emb_key"]

    adata_rna = sc.read_h5ad(rna_path)
    adata_atac = sc.read_h5ad(atac_path)

    edge_path = "D:/cuhk/project/velocity/dataset/atac/cluster_edges.npy"
    edge_dict = np.load(edge_path, allow_pickle=True).item()
    cluster_edges = edge_dict[data_name]

    # step 2 parameters
    tree_folder = "results"
    save_folder = "results"
    rev_stat = "mean"

    # step 3 parameters
    t1 = 0.1
    t2 = 1.0
    show_fig = True
    start_mode = "stochastic"
    # start_mode = "dynamical" # for MouseSkin

    # network configurations
    filter_genes = True
    save_coeff = True
    constrain = True
    loss_fun = "mse"
    only_s = False
    alpha_1 = 1
    alpha_2 = 0.1
    batch_size = 1024
    n_epochs = 100

    velocity_key = "velocity"
    adjust_DTI = False
    measure_performance = True

    adata = tivelo(adata_rna, group_key, emb_key,
                   data_name=data_name,
                   save_folder=save_folder,
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


