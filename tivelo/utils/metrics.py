import numpy as np
from sklearn.metrics.pairwise import cosine_similarity


def acc(y_true, y_pred):
    """
    Calculate clustering accuracy. Require scikit-learn installed
    # Arguments
        y: true labels, numpy.array with shape `(n_samples,)`
        y_pred: predicted labels, numpy.array with shape `(n_samples,)`
    # Return
        accuracy, in [0,1]
    """
    y_true = y_true.astype(np.int64)
    assert y_pred.size == y_true.size
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    from scipy.optimize import linear_sum_assignment
    ind = linear_sum_assignment(w.max() - w)
    ind = np.asarray(ind)
    ind = np.transpose(ind)
    return sum([w[i, j] for i, j in ind]) * 1.0 / y_pred.size


def keep_type(adata, nodes, target, cluster_key):

    return nodes[adata.obs[cluster_key][nodes].values == target]


def cross_boundary_correctness(adata, cluster_key, velocity_key, cluster_edges, return_raw=False, x_emb="X_umap"):

    scores = {}
    all_scores = {}

    if x_emb == "X_umap":
        v_emb = adata.obsm['{}_umap'.format(velocity_key)]
    else:
        v_emb = adata.obsm[[key for key in adata.obsm if key.startswith(velocity_key)][0]]

    x_emb = adata.obsm[x_emb]

    for u, v in cluster_edges:
        sel = adata.obs[cluster_key] == u
        nbs = adata.uns['neighbors']['indices'][sel]

        boundary_nodes = map(lambda nodes: keep_type(adata, nodes, v, cluster_key), nbs)
        x_points = x_emb[sel]
        x_velocities = v_emb[sel]

        type_score = []
        for x_pos, x_vel, nodes in zip(x_points, x_velocities, boundary_nodes):
            if len(nodes) == 0:
                continue

            position_dif = x_emb[nodes] - x_pos
            dir_scores = cosine_similarity(position_dif, x_vel.reshape(1, -1)).flatten()
            type_score.append(np.mean(dir_scores))

        scores[(u, v)] = np.mean(type_score)
        all_scores[(u, v)] = type_score

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values() if not np.isnan(sc)])


def cross_boundary_correctness2(adata, cluster_key, velocity_key, cluster_edges, return_raw=False):

    scores = {}
    all_scores = {}

    velocities = adata.layers[velocity_key]
    x = adata.layers["Ms"]

    for u, v in cluster_edges:
        sel = adata.obs[cluster_key] == u
        nbs = adata.uns['neighbors']['indices'][sel]

        boundary_nodes = map(lambda nodes: keep_type(adata, nodes, v, cluster_key), nbs)
        x_points = x[sel]
        x_velocities = velocities[sel]

        type_score = []
        for x_pos, x_vel, nodes in zip(x_points, x_velocities, boundary_nodes):
            if len(nodes) == 0:
                continue

            position_dif = x[nodes] - x_pos
            dir_scores = cosine_similarity(position_dif, x_vel.reshape(1, -1)).flatten()
            type_score.append(np.mean(dir_scores))

        scores[(u, v)] = np.mean(type_score)
        all_scores[(u, v)] = type_score

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values() if not np.isnan(sc)])


def inner_cluster_coh(adata, cluster_key, velocity_key, return_raw=False):

    clusters = np.unique(adata.obs[cluster_key])
    scores = {}
    all_scores = {}
    for cat in clusters:
        sel = adata.obs[cluster_key] == cat
        nbs = adata.uns['neighbors']['indices'][sel]
        same_cat_nodes = map(lambda nodes: keep_type(adata, nodes, cat, cluster_key), nbs)
        velocities = adata.layers[velocity_key]
        cat_vels = velocities[sel]
        cat_score = [cosine_similarity(cat_vels[[ith]], velocities[nodes]).mean()
                     for ith, nodes in enumerate(same_cat_nodes)
                     if len(nodes) > 0]
        all_scores[cat] = cat_score
        scores[cat] = np.mean(cat_score)

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values() if not np.isnan(sc)])


def inner_cluster_coh2(adata, cluster_key, velocity_key, return_raw=False, x_emb="X_umap"):
    scores = {}
    all_scores = {}

    if x_emb == "X_umap":
        v_emb = adata.obsm['{}_umap'.format(velocity_key)]
    else:
        v_emb = adata.obsm[[key for key in adata.obsm if key.startswith(velocity_key)][0]]

    clusters = np.unique(adata.obs[cluster_key])
    for cat in clusters:
        sel = adata.obs[cluster_key] == cat
        nbs = adata.uns['neighbors']['indices'][sel]
        same_cat_nodes = map(lambda nodes: keep_type(adata, nodes, cat, cluster_key), nbs)
        cat_vels = v_emb[sel]
        cat_score = [cosine_similarity(cat_vels[[ith]], v_emb[nodes]).mean()
                     for ith, nodes in enumerate(same_cat_nodes)
                     if len(nodes) > 0]
        all_scores[cat] = cat_score
        scores[cat] = np.mean(cat_score)

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values() if not np.isnan(sc)])


def cross_boundary_scvelo_probs(adata, cluster_key, cluster_edges, trans_g_key, return_raw=False):
    scores = {}
    all_scores = {}

    for u, v in cluster_edges:
        sel = adata.obs[cluster_key] == u
        nbs = adata.uns['neighbors']['indices'][sel]
        boundary_nodes = map(lambda nodes: keep_type(adata, nodes, v, cluster_key), nbs)
        type_score = [trans_probs.toarray()[:, nodes].mean()
                      for trans_probs, nodes in zip(adata.uns[trans_g_key][sel], boundary_nodes)
                      if len(nodes) > 0]
        scores[(u, v)] = np.mean(type_score)
        all_scores[(u, v)] = type_score
    if return_raw:
        return all_scores
    return scores, np.mean([sc for sc in scores.values() if not np.isnan(sc)])


def cross_boundary_coh(adata, cluster_key, velocity_key, cluster_edges, return_raw=False):
    scores = {}
    all_scores = {}
    for u, v in cluster_edges:
        sel = adata.obs[cluster_key] == u
        nbs = adata.uns['neighbors']['indices'][sel]
        boundary_nodes = map(lambda nodes: keep_type(adata, nodes, v, cluster_key), nbs)

        velocities = adata.layers[velocity_key]
        v_us = velocities[sel]
        type_score = [cosine_similarity(v_us[[ith]], velocities[nodes]).mean()
                      for ith, nodes in enumerate(boundary_nodes)
                      if len(nodes) > 0]
        scores[(u, v)] = np.mean(type_score)
        all_scores[(u, v)] = type_score

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values() if not np.isnan(sc)])


def velo_coh(adata, velocity_key, trans_g_key):
    s = adata.layers["Ms"]
    T = adata.uns[trans_g_key].toarray()
    T = np.exp(T) / np.sum(np.exp(T), axis=1, keepdims=True)
    Ts = np.dot(T, s)
    delta = Ts - s
    velocities = adata.layers[velocity_key]

    velo_coh = np.sum(delta * velocities, axis=1) / (np.linalg.norm(delta, axis=1) * np.linalg.norm(velocities, axis=1) + 1e-8)
    return np.mean(velo_coh)


