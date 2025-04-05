import os
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvelo as scv
from .path.process import process_path
from .direction.correct import correct_path
from .velocity.DTI import get_child_dict, get_d_nn, directed_graph
from .velocity.model import get_velocity
from .velocity.model_rate import get_velocity_rate
from .utils.metrics import inner_cluster_coh, cross_boundary_correctness, cross_boundary_scvelo_probs,\
    cross_boundary_correctness2, inner_cluster_coh2, velo_coh


def tivelo(adata, group_key, emb_key, res=0.6, data_name="data", save_folder="results", njobs=-1,
           start_mode="stochastic", rev_stat="mean", tree_gene=None, t1=0.1, t2=1, show_fig=True, filter_genes=True,
           constrain=True, loss_fun="mse", only_s=True, alpha_1=1, alpha_2=0.1, batch_size=1024, n_epochs=100,
           adjust_DTI=False, show_DTI=False, velocity_key="velocity", cluster_edges=None, measure_performance=True,
           rate_mode=False):

    # create path
    result_path = save_folder + "/{}/".format(data_name)
    figs_path = save_folder + "/{}/figs/".format(data_name)
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    if not os.path.exists(figs_path):
        os.makedirs(figs_path)

    # load existing adata
    try:
        adata_ = sc.read_h5ad(result_path + "tivelo.h5ad")

        # load dictionary
        child_dict_ = adata_.uns["child_dict"]
        path_dict = adata_.uns["path_dict"]
        start_node = adata_.uns["path_dict"]["start_node"]
        thres_list = adata_.uns["threshold_list"]

        # add PAGA result to original adata
        if "paga" not in adata_.uns.keys():
            sc.tl.diffmap(adata_)
            sc.tl.paga(adata_, groups=group_key)

    except FileNotFoundError:
        # set group_key, emb_key and neighbourhood
        if group_key is None or group_key not in adata.obs.keys():
            sc.tl.leiden(adata, resolution=res)
            group_key = "leiden"
        if "neighbors" not in adata.uns.keys():
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
        if emb_key is None or emb_key not in adata.obsm.keys():
            scv.tl.umap(adata)
            emb_key = "X_umap"

        # step 1
        path_dict, ax = process_path(adata, group_key, emb_key, njobs=njobs, start_mode=start_mode)

        ax.axis('off')
        ax.legend(loc="best", ncol=1, prop={'size': 10}, frameon=False)
        plt.savefig(figs_path + "path.pdf")
        if show_fig:
            plt.show()
        plt.close()

        # step 2
        path_dict, ax = correct_path(path_dict, adata, group_key, rev_stat=rev_stat, tree_gene=tree_gene,
                                     root_select="connectivities")

        if ax is not None:
            plt.savefig(figs_path + "tree.pdf")
            if show_fig:
                plt.show()
            plt.close()

        # step 3
        # get DNN
        child_dict_, level_dict_, thres_list = get_child_dict(adata, group_key, path_dict, threshold=t1,
                                                              threshold_trans=t2, adjust=False)
        start_node = list(level_dict_.keys())[0]
        d_nn_, knn_, child_c_, same_c_ = get_d_nn(adata, group_key, child_dict_, start_node, emb_key=emb_key,
                                                  root_select="connectivities")

        # velocity inference
        if not rate_mode:
            adata_, v_u, v_s = get_velocity(adata, d_nn_, knn_, same_c=same_c_, n_epochs=n_epochs, loss_fun=loss_fun,
                                            only_s=only_s, constrain=constrain, alpha_1=alpha_1, alpha_2=alpha_2,
                                            batch_size=batch_size, filter_genes=filter_genes)
        else:
            adata_, v_u, v_s = get_velocity_rate(adata, d_nn_, knn_, same_c=same_c_, n_epochs=n_epochs, loss_fun=loss_fun,
                                                 only_s=only_s, constrain=constrain, alpha_1=alpha_1, alpha_2=alpha_2,
                                                 batch_size=batch_size, filter_genes=filter_genes)

        scv.tl.velocity_graph(adata_, vkey=velocity_key, n_jobs=njobs)
        scv.tl.velocity_embedding(adata_)
        adata_.write_h5ad(result_path + "tivelo.h5ad")

    # velocity stream plot
    ax = scv.pl.velocity_embedding_stream(adata_, vkey=velocity_key, color=group_key, title='', cutoff_perc=0,
                                          show=False)
    plt.savefig(figs_path + "tivelo.png")
    if show_fig:
        plt.show()
    plt.close()

    # DTI
    if show_DTI:
        ax = directed_graph(adata_, child_dict_, group_key, start_node, emb_key=emb_key, data_name=data_name)
        ax.axis('off')
        plt.savefig(figs_path + "DTI.pdf")
        if show_fig:
            plt.show()
        plt.close()

    # adjust DTI
    if adjust_DTI:
        for t_ in thres_list:
            child_dict_new, _, _ = get_child_dict(adata_, group_key, path_dict, threshold=t1, threshold_trans=t_,
                                                  adjust=False)
            if child_dict_new != child_dict_:
                child_dict_ = child_dict_new
                ax = directed_graph(adata_, child_dict_, group_key, start_node, emb_key=emb_key, data_name=data_name)
                ax.set_title(data_name + ", " + r"$t_2$={}".format(np.round(t_, 2)))
                plt.savefig(figs_path + "/DTI_t2={}.png".format(np.round(t_, 2)))
                if show_fig:
                    plt.show()

    # measure performance
    if measure_performance:
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

            print("TIVelo:\n", "CBDir:", "%.4f" % cbdir, "ICVCoh:", "%.4f" % icvcoh, "\n",
                  "CBDir2:", "%.4f" % cbdir2, "ICVCoh2:", "%.4f" % icvcoh2, "\n",
                  "TransProbs:", "%.4f" % trans_probs, "VeloCoh:", "%.4f" % velocoh)

        else:
            _, icvcoh = inner_cluster_coh(adata_, cluster_key=group_key, velocity_key=velocity_key)
            _, icvcoh2 = inner_cluster_coh2(adata_, cluster_key=group_key, velocity_key=velocity_key, x_emb=emb_key)
            velocoh = velo_coh(adata_, velocity_key=velocity_key, trans_g_key="{}_graph".format(velocity_key))

            print("TIVelo:\n", "ICVCoh:", "%.4f" % icvcoh, "\n", "ICVCoh2:", "%.4f" % icvcoh2, "\n",
                  "VeloCoh:", "%.4f" % velocoh)

    return adata_



