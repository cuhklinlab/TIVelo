import numpy as np
import scanpy as sc
import pandas as pd

import sys
sys.path.append("..")
from tivelo.main import tivelo


if __name__ == '__main__':
    # data_name = "rpe1"
    data_name = "u2os"
    data_path = "/lustre/project/Stat/s1155184322/datasets/velocity/FUCCI/adata_fucci_u2os_processed.h5ad"
    adata = sc.read_h5ad(data_path)

    # velocity
    group_key = "leiden"
    emb_key = "X_umap"
    res = 0.6

    # parameters in step 1&2
    # pancreas: t2=0.9; hindbrain2: t2=2.0; rpe1: t1=0.4; u2os: t2=3.2
    t1 = 0.1
    t2 = 1.0
    show_fig = True
    start_mode = "stochastic"
    rev_stat = "mean"
    adjust_DTI = False

    # parameters in step 3
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
                   )



