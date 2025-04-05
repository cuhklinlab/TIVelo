import scanpy as sc
import pandas as pd
import numpy as np

import sys
sys.path.append("..")
from tivelo.main import tivelo


if __name__ == '__main__':
    data_name = "HSPCs"
    rna_path = "/lustre/project/Stat/s1155184322/datasets/velocity/multi-omics/3423-MV-2_adata_postpro.h5ad"
    atac_path = "/lustre/project/Stat/s1155184322/datasets/velocity/multi-omics/3423-MV-2_adata_atac_postpro.h5ad"
    adata_rna = sc.read(rna_path)
    adata_atac = sc.read(atac_path)

    group_key = "leiden"
    emb_key = "X_umap"
    cluster_edges = [("HSC", "MPP"), ("MPP", "LMPP"),  ("LMPP", "GMP"), ("HSC", "MEP"), ("MEP", "Erythrocyte"),
                    ("MEP", "Prog MK"), ("Prog MK", "Platele")]

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
                   njobs=32,
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


