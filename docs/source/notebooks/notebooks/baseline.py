import torch
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv
import anndata as ad
from sklearn.neighbors import NearestNeighbors
import unitvelo as utv
import multivelo as mv
from velovi import preprocess_data, VELOVI
from tivelo.utils.metrics import inner_cluster_coh, cross_boundary_correctness, cross_boundary_scvelo_probs, \
    cross_boundary_correctness2, inner_cluster_coh2, velo_coh


# copy from celldancer package:
# https://github.com/GuangyuWangLab2021/cellDancer/blob/main/src/celldancer/utilities.py#L233
def to_dynamo(cellDancer_df):
    cellDancer_df = cellDancer_df.sort_values('cellID')

    spliced = cellDancer_df.pivot(index='cellID', columns='gene_name', values='splice')
    unspliced = cellDancer_df.pivot(index='cellID', columns='gene_name', values='unsplice')

    spliced_predict = cellDancer_df.pivot(index='cellID', columns='gene_name', values='splice_predict')
    unspliced_predict = cellDancer_df.pivot(index='cellID', columns='gene_name', values='unsplice_predict')

    alpha = cellDancer_df.pivot(index='cellID', columns='gene_name', values='alpha')
    beta = cellDancer_df.pivot(index='cellID', columns='gene_name', values='beta')
    gamma = cellDancer_df.pivot(index='cellID', columns='gene_name', values='gamma')

    one_gene = cellDancer_df['gene_name'].iloc[0]
    one_cell = cellDancer_df['cellID'].iloc[0]

    adata1 = ad.AnnData(spliced)

    # var
    adata1.var['highly_variable_genes'] = True
    # adata1.var['loss'] = (cellDancer_df[cellDancer_df['cellID'] == one_cell]['loss']).tolist()
    loss = cellDancer_df.pivot(index='gene_name', columns='cellID', values='loss').iloc[:, 0]
    loss.index = loss.index.astype(str)
    adata1.var['loss'] = loss
    # celldancer uses all genes (high variable) for dynamics and transition.
    adata1.var['use_for_dynamics'] = True
    adata1.var['use_for_transition'] = True

    # obs
    if 'clusters' in cellDancer_df:
        clusters = cellDancer_df.pivot(index='cellID', columns='gene_name', values='clusters').iloc[:, 0]
        clusters.index = clusters.index.astype(str)
        adata1.obs['clusters'] = clusters
    #  layers
    adata1.layers['X_spliced'] = spliced
    adata1.layers['X_unspliced'] = unspliced

    adata1.layers['M_s'] = spliced
    adata1.layers['M_u'] = unspliced
    adata1.layers['velocity_S'] = spliced_predict - spliced

    adata1.layers['velocity_U'] = unspliced_predict - unspliced
    adata1.layers['alpha'] = alpha
    adata1.layers['beta'] = beta
    adata1.layers['gamma'] = gamma

    # obsm
    adata1.obsm['X_cdr'] = cellDancer_df[cellDancer_df['gene_name'] == one_gene][['embedding1', 'embedding2']].values
    # assuming no downsampling is used for the cell velocities in the cellDancer_df
    if 'velocity1' in cellDancer_df:
        adata1.obsm['velocity_cdr'] = cellDancer_df[cellDancer_df['gene_name'] == one_gene][
            ['velocity1', 'velocity2']].values

    # obsp
    n_neighbors = 20
    nn = NearestNeighbors(n_neighbors=n_neighbors)
    nn.fit(adata1.obsm['X_cdr'])
    connect_knn = nn.kneighbors_graph(mode='connectivity')
    distance_knn = nn.kneighbors_graph(mode='distance')
    adata1.obsp['connectivities'] = connect_knn
    adata1.obsp['distances'] = distance_knn

    # uns
    dynamics_info = {'filter_gene_mode': 'final',
                     't': None,
                     'group': None,
                     'X_data': None,
                     'X_fit_data': None,
                     'asspt_mRNA': 'ss',
                     'experiment_type': 'conventional',
                     'normalized': True,
                     'model': 'static',
                     'est_method': 'ols',
                     'has_splicing': True,
                     'has_labeling': False,
                     'splicing_labeling': False,
                     'has_protein': False,
                     'use_smoothed': True,
                     'NTR_vel': False,
                     'log_unnormalized': False,
                     'fraction_for_deg': False}

    adata1.uns['dynamics'] = dynamics_info

    return adata1


def df_to_adata(adata, cellDancer_df, group_key, emb_key):
    # transform celldancer dataframe to anndata
    velocity_key = "velocity_S"
    adata_from_dancer = to_dynamo(cellDancer_df)
    for key in adata.obs.keys():
        adata_from_dancer.obs[key] = adata.obs[key]
    adata_from_dancer.uns["neighbors"] = adata.uns["neighbors"]
    adata_from_dancer.layers["Mu"] = adata_from_dancer.layers["M_u"]
    adata_from_dancer.layers["Ms"] = adata_from_dancer.layers["M_s"]
    adata_from_dancer.obsm["{}_{}".format(velocity_key, emb_key.split("_")[-1])] = adata_from_dancer.obsm["velocity_cdr"]
    adata_from_dancer.obsm[emb_key] = adata_from_dancer.obsm["X_cdr"]
    try:
        adata_from_dancer.uns["{}_colors".format(group_key)] = adata.uns["{}_colors".format(group_key)]
    except KeyError:
        pass
    del adata_from_dancer.layers["M_u"], adata_from_dancer.layers["M_s"], adata_from_dancer.obsm["velocity_cdr"], \
        adata_from_dancer.obsm["X_cdr"]

    return adata_from_dancer


def add_velovi_outputs_to_adata(adata, vae):
    latent_time = vae.get_latent_time(n_samples=25)
    velocities = vae.get_velocity(n_samples=25, velo_statistic="mean")

    t = latent_time
    scaling = 20 / t.max(0)

    adata.layers["velocity"] = velocities / scaling
    adata.layers["latent_time_velovi"] = latent_time

    adata.var["fit_alpha"] = vae.get_rates()["alpha"] / scaling
    adata.var["fit_beta"] = vae.get_rates()["beta"] / scaling
    adata.var["fit_gamma"] = vae.get_rates()["gamma"] / scaling
    adata.var["fit_t_"] = (torch.nn.functional.softplus(vae.module.switch_time_unconstr).detach().cpu().numpy()) * scaling
    adata.layers["fit_t"] = latent_time.values * np.array(scaling)[None, :]
    adata.var['fit_scaling'] = 1.0


def run_velovi(adata):
    adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train(max_epochs=500)

    add_velovi_outputs_to_adata(adata, vae)
    scv.tl.velocity_graph(adata, n_jobs=32)
    scv.tl.velocity_embedding(adata)
    return adata


def run_baseline(adata, model, data_name, group_key, emb_key, cluster_edges=None, adata_atac=None, save_folder="results",
                 show_fig=False, measure_performance=True, unitvelo_mode="1"):
    if model == "celldancer":
        velocity_key = "velocity_S"
    elif model == "multivelo":
        velocity_key = "velo_s_norm"
    else:
        velocity_key = "velocity"

    result_path = save_folder + "/{}/".format(data_name)
    figs_path = save_folder + "/{}/figs/".format(data_name)

    if not os.path.exists(result_path):
        os.makedirs(result_path)
    if not os.path.exists(figs_path):
        os.makedirs(figs_path)

    # preprocess
    if group_key is None or group_key not in adata.obs.keys():
        sc.tl.leiden(adata, resolution=0.6)
        group_key = "leiden"
    if emb_key is None or emb_key not in adata.obsm.keys():
        scv.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
        scv.tl.umap(adata)
        emb_key = "umap"
    if "neighbors" not in adata.uns.keys() or adata.uns["neighbors"]["params"]["n_neighbors"] != 30:
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
        sorted_indices = np.argsort(adata.obsp["distances"].A + np.identity(adata.n_obs), axis=1)
        sorted_indices = np.fliplr(sorted_indices)
        indices = sorted_indices[:, 0: 30]
        adata.uns["neighbors"]["indices"] = indices
    if "indices" not in adata.uns["neighbors"].keys():
        sorted_indices = np.argsort(adata.obsp["distances"].A + np.identity(adata.n_obs), axis=1)
        sorted_indices = np.fliplr(sorted_indices)
        indices = sorted_indices[:, 0: 30]
        adata.uns["neighbors"]["indices"] = indices

    try:
        adata_ = sc.read_h5ad(result_path + "{}.h5ad".format(model))
    except FileNotFoundError:
        if model == "scvelo":
            adata_ = scv.tl.velocity(adata, mode='stochastic', copy=True, vkey=velocity_key)
            scv.tl.velocity_graph(adata_, vkey=velocity_key, n_jobs=32)
            scv.tl.velocity_embedding(adata_, vkey=velocity_key)
            adata_.write_h5ad(result_path + "{}.h5ad".format(model))
        elif model == "scvelo2":
            scv.tl.recover_dynamics(adata, n_jobs=32)
            adata_ = scv.tl.velocity(adata, mode='dynamical', copy=True, vkey=velocity_key)
            scv.tl.velocity_graph(adata_, vkey=velocity_key, n_jobs=32)
            scv.tl.velocity_embedding(adata_, vkey=velocity_key)
            adata_.write_h5ad(result_path + "{}.h5ad".format(model))
        elif model == "unitvelo":
            velo_config = utv.config.Configuration()
            velo_config.MAX_ITER = 10000
            if unitvelo_mode == "2":
                velo_config.FIT_OPTION = '2'
            adata_ = utv.run_model(adata, group_key, config_file=velo_config)
            adata_.write_h5ad(result_path + "{}.h5ad".format(model))
        elif model == "velovi":
            adata_ = run_velovi(adata)
            adata_.write_h5ad(result_path + "{}.h5ad".format(model))
        elif model == "celldancer":
            try:
                # save the output of celldancer here: save_folder/celldancer/data_name/data_name.csv
                cellDancer_df = pd.read_csv(save_folder + "/{}/celldancer/{}.csv".format(data_name, data_name))
            except FileNotFoundError:
                raise FileNotFoundError("No celldancer output file found at {}".format(save_folder))
            adata_ = df_to_adata(adata, cellDancer_df, group_key, emb_key)
            scv.tl.velocity_graph(adata_, vkey=velocity_key, xkey="Ms", n_jobs=32)
            adata_.write_h5ad(result_path + "{}.h5ad".format(model))
        elif model == "multivelo" and adata_atac is not None:
            adata_ = mv.recover_dynamics_chrom(adata, adata_atac, max_iter=5, init_mode="invert", parallel=False,
                                               n_jobs=10, save_plot=False, rna_only=False, fit=True, n_anchors=500,
                                               extra_color_key=group_key)
            adata_.write_h5ad(result_path + "{}.h5ad".format(model))
        else:
            raise KeyError("No such model for comparison.")

    # velocity stream plot
    if model == "multivelo":
        ax = mv.velocity_embedding_stream(adata_, basis=emb_key, color=group_key, show=False, title="")
    else:
        ax = scv.pl.velocity_embedding_stream(adata_, vkey=velocity_key, color=group_key, title="", show=False)

    plt.tight_layout()
    plt.savefig(figs_path + "{}_velo.png".format(model))
    if show_fig:
        plt.show()
    plt.close()

    # if model == "scvelo2" or model == "unitvelo":
    #     adata_ = adata_[:, adata_.var.loc[adata_.var['{}_genes'.format(velocity_key)] == True].index]

    # measure performance
    if measure_performance:
        if model == "scvelo2" or model == "unitvelo":
            adata_ = adata_[:, adata_.var.loc[adata_.var['{}_genes'.format(velocity_key)] == True].index]
        
        if cluster_edges is not None:
            _, cbdir = cross_boundary_correctness(adata_, cluster_key=group_key, velocity_key=velocity_key,
                                                  cluster_edges=cluster_edges, x_emb=emb_key)
            _, cbdir2 = cross_boundary_correctness2(adata_, cluster_key=group_key, velocity_key=velocity_key,
                                                    cluster_edges=cluster_edges)
            _, trans_probs = cross_boundary_scvelo_probs(adata_, cluster_key=group_key, cluster_edges=cluster_edges,
                                                         trans_g_key="{}_graph".format(velocity_key))
            _, icvcoh = inner_cluster_coh(adata_, cluster_key=group_key, velocity_key=velocity_key)
            _, icvcoh2 = inner_cluster_coh2(adata_, cluster_key=group_key, velocity_key=velocity_key, x_emb=emb_key)
            velocoh = velo_coh(adata_, velocity_key=velocity_key, trans_g_key="{}_graph".format(velocity_key))

            print("{}:\n".format(model), "CBDir:", "%.4f" % cbdir, "ICVCoh:", "%.4f" % icvcoh, "\n",
                  "CBDir2:", "%.4f" % cbdir2, "ICVCoh2:", "%.4f" % icvcoh2, "\n",
                  "TransProbs:", "%.4f" % trans_probs, "VeloCoh:", "%.4f" % velocoh)

        else:
            _, icvcoh = inner_cluster_coh(adata_, cluster_key=group_key, velocity_key=velocity_key)
            _, icvcoh2 = inner_cluster_coh2(adata_, cluster_key=group_key, velocity_key=velocity_key, x_emb=emb_key)
            velocoh = velo_coh(adata_, velocity_key=velocity_key, trans_g_key="{}_graph".format(velocity_key))

            print("{}:\n".format(model), "ICVCoh:", "%.4f" % icvcoh, "\n", "ICVCoh2:", "%.4f" % icvcoh2, "\n",
                  "VeloCoh:", "%.4f" % velocoh)

    return adata_





