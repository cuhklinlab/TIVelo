import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import BatchSampler, RandomSampler
from tqdm import tqdm
import numpy as np
from ..utils.velocity_genes import compute_velocity_genes


class model_velo(object):
    def __init__(self, x_u, x_s, d_nn, knn):
        self.x_u = torch.from_numpy(x_u)
        self.x_s = torch.from_numpy(x_s)
        self.n_obs = self.x_u.shape[0]
        self.n_vars = self.x_u.shape[1]
        self.knn = torch.from_numpy(knn).to(torch.float32)
        self.d_nn = torch.from_numpy(d_nn).to(torch.float32)

        self.n_neighs_d = torch.sum(self.d_nn, dim=1)
        self.n_neighs_d[torch.where(self.n_neighs_d == 0)] = 1
        self.d_u = torch.mm(self.d_nn, self.x_u) / self.n_neighs_d.unsqueeze(1)
        self.d_s = torch.mm(self.d_nn, self.x_s) / self.n_neighs_d.unsqueeze(1)

        v_u, v_s = self.d_u - self.x_u, self.d_s - self.x_s
        self.v_u = nn.Parameter(v_u, requires_grad=True)
        self.v_s = nn.Parameter(v_s, requires_grad=True)

    def simple_fit(self):
        v_u, v_s = self.d_u - self.x_u, self.d_s - self.x_s
        v_u, v_s = v_u.clone().detach().numpy(), v_s.clone().detach().numpy()
        return v_u, v_s

    def fit(self, n_epochs=100, lr=0.001, alpha_1=1, alpha_2=0.1, batch_size=None, loss_fun="mse", only_s=False,
            same_c=None):
        same_c = torch.from_numpy(same_c).to(torch.float32) if same_c is not None else None
        mse = nn.MSELoss()
        cos_sim = nn.CosineSimilarity(dim=1)
        optimizer = optim.Adam([{'params': self.v_u}, {'params': self.v_s}], amsgrad=True, lr=lr)

        # compute loss
        pbar = tqdm(range(n_epochs))
        for epoch in pbar:
            if batch_size is None:
                optimizer.zero_grad()

                mse_u = mse(self.x_u + self.v_u, self.d_u)
                mse_s = mse(self.x_s + self.v_s, self.d_s)

                cos_xu = torch.mean(cos_sim(self.v_u, self.d_u - self.x_u))
                cos_xs = torch.mean(cos_sim(self.v_s, self.d_s - self.x_s))

                # cosine similarity with velocity
                norm_v_u = torch.norm(self.v_u, p=2, dim=1, keepdim=True)
                norm_v_s = torch.norm(self.v_s, p=2, dim=1, keepdim=True)

                cos_u = torch.mm(self.v_u, self.v_u.t()) / (torch.mm(norm_v_u, norm_v_u.t()) + 1e-6)
                cos_s = torch.mm(self.v_s, self.v_s.t()) / (torch.mm(norm_v_s, norm_v_s.t()) + 1e-6)

                if same_c is None:
                    n_neighs = torch.sum(self.knn, dim=1)
                    cos_u = torch.mean(torch.sum(cos_u * self.knn, dim=1) / (n_neighs + 1e-6))
                    cos_s = torch.mean(torch.sum(cos_s * self.knn, dim=1) / (n_neighs + 1e-6))
                else:
                    n_neighs = torch.sum(same_c * self.knn, dim=1)
                    cos_u = torch.mean(torch.sum(cos_u * same_c * self.knn, dim=1) / (n_neighs + 1e-6))
                    cos_s = torch.mean(torch.sum(cos_s * same_c * self.knn, dim=1) / (n_neighs + 1e-6))

                if loss_fun == "mse":
                    if not only_s:
                        loss = mse_u + alpha_1 * mse_s - alpha_2 * cos_u - alpha_2 * cos_s
                        loss.backward()
                        optimizer.step()

                        pbar.set_description("Model training")
                        pbar.set_postfix(mse_u=np.around(mse_u.clone().detach().numpy(), 3),
                                         mse_s=np.around(mse_s.clone().detach().numpy(), 3),
                                         cos_u=np.around(cos_u.clone().detach().numpy(), 3),
                                         cos_s=np.around(cos_s.clone().detach().numpy(), 3))
                    else:
                        loss = mse_s - alpha_2 * cos_s
                        loss.backward()
                        optimizer.step()

                        pbar.set_description("Model training")
                        pbar.set_postfix(mse_s=np.around(mse_s.clone().detach().numpy(), 3),
                                         cos_s=np.around(cos_s.clone().detach().numpy(), 3))
                elif loss_fun == "cos":
                    if not only_s:
                        loss = -cos_xu - alpha_1 * cos_xs - alpha_2 * cos_u - alpha_2 * cos_s
                        loss.backward()
                        optimizer.step()

                        pbar.set_description("Model training")
                        pbar.set_postfix(cos_xu=np.around(cos_xu.clone().detach().numpy(), 3),
                                         cos_xs=np.around(cos_xs.clone().detach().numpy(), 3),
                                         cos_u=np.around(cos_u.clone().detach().numpy(), 3),
                                         cos_s=np.around(cos_s.clone().detach().numpy(), 3))
                    else:
                        loss = -cos_xs - alpha_2 * cos_s
                        loss.backward()
                        optimizer.step()

                        pbar.set_description("Model training")
                        pbar.set_postfix(cos_xs=np.around(cos_xs.clone().detach().numpy(), 3),
                                         cos_s=np.around(cos_s.clone().detach().numpy(), 3))

            else:
                train_loader = list(BatchSampler(RandomSampler(range(self.n_obs)), batch_size=batch_size, drop_last=False))
                mse_u_, mse_s_, cos_xu_, cos_xs_, cos_u_, cos_s_ = [], [], [], [], [], []

                for step, indices in enumerate(train_loader):
                    optimizer.zero_grad()

                    x_u, x_s = self.x_u[indices, :], self.x_s[indices, :]
                    d_u, d_s = self.d_u[indices, :], self.d_s[indices, :]

                    mse_u = mse(x_u + self.v_u[indices, :], d_u)
                    mse_s = mse(x_s + self.v_s[indices, :], d_s)

                    cos_xu = torch.mean(cos_sim(self.v_u[indices, :], d_u - x_u))
                    cos_xs = torch.mean(cos_sim(self.v_s[indices, :], d_s - x_s))

                    norm_v_u = torch.norm(self.v_u[indices, :], p=2, dim=1, keepdim=True)
                    norm_v_s = torch.norm(self.v_s[indices, :], p=2, dim=1, keepdim=True)

                    cos_u = torch.mm(self.v_u[indices, :], self.v_u[indices, :].t()) / (torch.mm(norm_v_u, norm_v_u.t()) + 1e-6)
                    cos_s = torch.mm(self.v_s[indices, :], self.v_s[indices, :].t()) / (torch.mm(norm_v_s, norm_v_s.t()) + 1e-6)

                    n_neighs = torch.mean(torch.sum(self.knn[indices, :][:, indices], dim=1))
                    cos_u = torch.sum(torch.mean(cos_u * self.knn[indices, :][:, indices], dim=1)) / (n_neighs + 1e-6)
                    cos_s = torch.sum(torch.mean(cos_s * self.knn[indices, :][:, indices], dim=1)) / (n_neighs + 1e-6)

                    if loss_fun == "mse":
                        if not only_s:
                            loss = mse_u + alpha_1 * mse_s - alpha_2 * cos_u - alpha_2 * cos_s
                            loss.backward()
                            optimizer.step()
                        else:
                            loss = mse_s - alpha_2 * cos_s
                            loss.backward()
                            optimizer.step()

                    elif loss_fun == "cos":
                        if not only_s:
                            loss = -cos_xu - alpha_1 * cos_xs - alpha_2 * cos_u - alpha_2 * cos_s
                            loss.backward()
                            optimizer.step()
                        else:
                            loss = -cos_xs - alpha_2 * cos_s
                            loss.backward()
                            optimizer.step()
                    else:
                        raise KeyError("No such loss function for simple method")

                    mse_u_.append(mse_u.unsqueeze(0).clone().detach())
                    mse_s_.append(mse_s.unsqueeze(0).clone().detach())
                    cos_xu_.append(cos_xu.unsqueeze(0).clone().detach())
                    cos_xs_.append(cos_xs.unsqueeze(0).clone().detach())
                    cos_u_.append(cos_u.unsqueeze(0).clone().detach())
                    cos_s_.append(cos_s.unsqueeze(0).clone().detach())

                mse_u = np.mean(torch.cat(mse_u_).numpy())
                mse_s = np.mean(torch.cat(mse_s_).numpy())
                cos_xu = np.mean(torch.cat(cos_xu_).numpy())
                cos_xs = np.mean(torch.cat(cos_xs_).numpy())
                cos_u = np.mean(torch.cat(cos_u_).numpy())
                cos_s = np.mean(torch.cat(cos_s_).numpy())

                pbar.set_description("Model training")

                if loss_fun == "mse":
                    if not only_s:
                        pbar.set_postfix(mse_u=np.around(mse_u, 3), mse_s=np.around(mse_s, 3),
                                         cos_u=np.around(cos_u, 3), cos_s=np.around(cos_s, 3))
                    else:
                        pbar.set_postfix(mse_s=np.around(mse_s, 3), cos_s=np.around(cos_s, 3))
                elif loss_fun == "cos":
                    if not only_s:
                        pbar.set_postfix(cos_xu=np.around(cos_xu, 3), cos_xs=np.around(mse_s, 3),
                                         cos_u=np.around(cos_u, 3), cos_s=np.around(cos_s, 3))
                    else:
                        pbar.set_postfix(cos_xs=np.around(cos_xs, 3), cos_s=np.around(cos_s, 3))

        v_u, v_s = self.v_u.clone().detach().numpy(), self.v_s.clone().detach().numpy()
        return v_u, v_s


def get_velocity(adata, d_nn, knn, same_c=None, n_epochs=100, loss_fun="mse", only_s=False, constrain=True, alpha_1=1,
                 alpha_2=0.1, batch_size=512, filter_genes=False):

    if filter_genes:
        velocity_genes = compute_velocity_genes(adata, n_top_genes=2000, inplace=False)
        adata = adata[:, velocity_genes]

    x_u = adata.layers["Mu"]
    x_s = adata.layers["Ms"]

    model = model_velo(x_u, x_s, d_nn, knn)
    if constrain:
        v_u, v_s = model.fit(batch_size=batch_size, n_epochs=n_epochs, loss_fun=loss_fun, only_s=only_s,
                             alpha_1=alpha_1, alpha_2=alpha_2, same_c=same_c)
    else:
        v_u, v_s = model.simple_fit()

    adata.layers["velocity"] = v_s
    adata.layers["velocity_u"] = v_u

    return adata, v_u, v_s




