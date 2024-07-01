import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors


def gaussian_kernel(X, mu=0, sigma=1):
    return np.exp(-(X - mu) ** 2 / (2 * sigma ** 2)) / np.sqrt(2 * np.pi * sigma ** 2)


def sampling_neighbors(gene_unsplice_splice, step=(30, 30), percentile=25):
    grs = []
    for dim_i in range(gene_unsplice_splice.shape[1]):
        m, M = np.min(gene_unsplice_splice[:, dim_i]), np.max(gene_unsplice_splice[:, dim_i])
        m = m - 0.025 * np.abs(M - m)
        M = M + 0.025 * np.abs(M - m)
        gr = np.linspace(m, M, step[dim_i])
        grs.append(gr)
    meshes_tuple = np.meshgrid(*grs)
    gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T
    gridpoints_coordinates = gridpoints_coordinates + norm.rvs(loc=0, scale=0.15, size=gridpoints_coordinates.shape)

    np.random.seed(42)  # set random seed

    nn = NearestNeighbors()

    neighbors_1 = min((gene_unsplice_splice[:, 0:2].shape[0] - 1), 20)
    nn.fit(gene_unsplice_splice[:, 0:2])
    dist, ixs = nn.kneighbors(gridpoints_coordinates, neighbors_1)

    ix_choice = ixs[:, 0].flat[:]
    ix_choice = np.unique(ix_choice)

    nn = NearestNeighbors()

    neighbors_2 = min((gene_unsplice_splice[:, 0:2].shape[0] - 1), 20)
    nn.fit(gene_unsplice_splice[:, 0:2])
    dist, ixs = nn.kneighbors(gene_unsplice_splice[ix_choice, 0:2], neighbors_2)

    density_extimate = gaussian_kernel(dist, mu=0, sigma=0.5).sum(1)
    bool_density = density_extimate > np.percentile(density_extimate, percentile)
    ix_choice = ix_choice[bool_density]

    return (ix_choice)


def scatter_gene(adata, velocity_key="velocity", gene_name=None, ax=None, color=None, custom_xlim=None,
                 custom_ylim=None, vmin=None, vmax=None, alpha=1, s=10, velocity=True, arrow_grid=(15, 15),
                 key_cluster=None, length_scale=5):

    if color in adata.obs.keys() and color + "_colors" in adata.uns.keys():
        group_list = list(adata.obs[color].cat.categories)
        color_list = adata.uns[color + "_colors"]
        color_dict = dict(zip(group_list, color_list))
        if key_cluster not in list(adata.obs[color].cat.categories):
            c = [color_dict[category] for category in adata.obs[color]]
        else:
            c = ["#FF0000" if category == key_cluster else "#95D9EF" for category in adata.obs[color]]

    elif color is None:
        c = "#95D9EF"

    assert gene_name, '\nError! gene is required!\n'

    gene_idx = list(adata.var_names).index(gene_name)
    Mu, Ms = adata.layers['Mu'][:, gene_idx], adata.layers['Ms'][:, gene_idx]
    ax.scatter(Ms, Mu, c=c, s=s, alpha=alpha, vmin=vmin, vmax=vmax, edgecolor="none")

    if custom_xlim is not None:
        ax.set_xlim(custom_xlim[0], custom_xlim[1])
    if custom_ylim is not None:
        ax.set_ylim(custom_ylim[0], custom_ylim[1])

    if velocity:
        u_s = np.concatenate((Mu[:, None], Ms[:, None]), axis=1)
        u_max, s_max = np.max(Mu), np.max(Ms)

        velocity = adata.layers[velocity_key][:, gene_idx]
        velocity_u = adata.layers[velocity_key + "_u"][:, gene_idx]
        # velocity_u, velocity = np.clip(velocity_u, -u_max/20, u_max/20), np.clip(velocity, -s_max/10, s_max/10)

        # preprocess velocity vector using sigmoid function
        velocity_u = u_max / length_scale * (2 / (1 + np.exp(-velocity_u)) - 1)
        velocity = s_max / length_scale * (2 / (1 + np.exp(-velocity)) - 1)

        sampling_idx = sampling_neighbors(u_s, step=arrow_grid, percentile=5)  # Sampling
        u_s_downsample = u_s[sampling_idx, :]
        velocity_downsample = velocity[sampling_idx]
        velocity_u_downsample = velocity_u[sampling_idx]

        plt.scatter(u_s_downsample[:, 1], u_s_downsample[:, 0], color="none", s=s, edgecolor="k")
        plt.quiver(u_s_downsample[:, 1], u_s_downsample[:, 0],
                   velocity_downsample, velocity_u_downsample,
                   angles='xy', scale_units='xy', scale=1)

    return ax