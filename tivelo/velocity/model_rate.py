import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils import data
import os
from tqdm import tqdm
import numpy as np
from ..utils.velocity_genes import compute_velocity_genes


class DNNLayer(nn.Module):
    def __init__(self, n_genes, n_dims=[256, 64]):
        super(DNNLayer, self).__init__()
        self.n_genes = n_genes

        self.alpha, self.beta, self.gamma = None, None, None
        
        self.layers = nn.ModuleList()
        # input layer
        self.layers.append(nn.Sequential(nn.Linear(self.n_genes * 2, n_dims[0]), nn.ReLU(True)))
        # hidden layers
        for i in range(len(n_dims) - 1):
            self.layers.append(nn.Sequential(nn.Linear(n_dims[i], n_dims[i + 1]), nn.ReLU(True)))

        out_layer_dim = n_genes * 3
        self.layers.append(nn.Sequential(nn.Linear(n_dims[-1], out_layer_dim), nn.ReLU(True)))

    def forward(self, x_u, x_s):
        h = torch.cat([x_u, x_s], dim=1)
        for i, layer in enumerate(self.layers):
            h = layer(h)
        x = h

        self.alpha = x[:, 0: self.n_genes]
        self.beta = x[:, self.n_genes: 2 * self.n_genes]
        self.gamma = x[:, 2 * self.n_genes: 3 * self.n_genes]

        v_u = self.alpha - self.beta * x_u
        v_s = self.beta * x_u - self.gamma * x_s

        return v_u, v_s

    def get_current_batch_kinetic_rates(self):
        return self.current_kinetic_rates


class DNN(object):
    def __init__(self, x_u, x_s, d_nn, knn, n_dims, normalize=False):

        assert x_u.shape[0] == x_s.shape[0] and x_u.shape[1] == x_s.shape[1],\
            "Data dimensions are not consistent"
        # configurations about inputs u and s
        self.x_u = torch.from_numpy(x_u)
        self.x_s = torch.from_numpy(x_s)

        self.x_u_norm = self.x_u / torch.nan_to_num(torch.std(self.x_u, dim=0, keepdim=True), 1)
        self.x_s_norm = self.x_s / torch.nan_to_num(torch.std(self.x_s, dim=0, keepdim=True), 1)

        self.d_nn = torch.from_numpy(d_nn).to(torch.float32)
        self.knn = torch.from_numpy(knn)

        self.n_cells, self.n_genes = x_u.shape[0], x_s.shape[1]
        self.n_dims = n_dims

        # initialize DNN
        self.dnn = DNNLayer(self.n_genes, n_dims)
        self.normalize = normalize

        # d_u, d_s
        n_neighs_d = torch.sum(self.d_nn, dim=1)
        n_neighs_d[torch.where(n_neighs_d == 0)] = 1
        self.d_u = torch.mm(self.d_nn, self.x_u) / n_neighs_d.unsqueeze(1)
        self.d_s = torch.mm(self.d_nn, self.x_s) / n_neighs_d.unsqueeze(1)
        self.v_u, self.v_s = self.d_u - self.x_u, self.d_s - self.x_s

    def pretrain(self, n_epochs=100, lr=1e-3, batch_size=128, loss_fun="mse", only_s=False):
        mse = nn.MSELoss()
        cos_sim = nn.CosineSimilarity(dim=1)
        parameter = self.dnn.parameters()
        optimizer = optim.Adam(parameter, amsgrad=True, lr=lr)
        batch_size = self.n_cells if batch_size is None else batch_size
        self.dnn.train()

        if self.normalize:
            train_data = data.TensorDataset(self.x_u_norm, self.x_s_norm, torch.arange(self.n_cells))
            train_loader = data.DataLoader(train_data, batch_size=batch_size, shuffle=False)
        else:
            train_data = data.TensorDataset(self.x_u, self.x_s, torch.arange(self.n_cells))
            train_loader = data.DataLoader(train_data, batch_size=batch_size, shuffle=False)

        pbar = tqdm(range(n_epochs))
        for epoch in pbar:
            loss_var = []
            for u, s, idx in train_loader:
                optimizer.zero_grad()
                
                v_u, v_s = self.dnn(u, s)

                if loss_fun == "mse":
                    if only_s:
                        loss = mse(v_s, self.v_s[idx, :])
                    else:
                        loss = mse(v_u, self.v_u[idx, :]) + mse(v_s, self.v_s[idx, :])
                elif loss_fun == "cos":
                    if only_s:
                        loss = -torch.mean(cos_sim(v_s, self.v_s[idx, :]))
                    else:
                        loss = -torch.mean(cos_sim(v_u, self.v_u[idx, :])) - torch.mean(cos_sim(v_s, self.v_s[idx, :]))
                else:
                    raise KeyError("No such loss function for pretraining")

                loss.backward()
                optimizer.step()
                loss_var.append(loss.clone().detach().numpy())

            pbar.set_description("Pretraining")
            if loss_fun == "mse":
                pbar.set_postfix(loss=np.around(np.mean(loss_var), 3))
            elif loss_fun == "cos":
                pbar.set_postfix(loss=np.around(-np.mean(loss_var), 3))

    def fit(self, n_epochs=100, lr=1e-3, batch_size=64, constrain=True, loss_fun="mse", alpha_1=1, alpha_2=1,
            only_s=False, same_c=None):
        mse = nn.MSELoss()
        cos_sim = nn.CosineSimilarity(dim=1)
        same_c = torch.from_numpy(same_c).to(torch.float32) if same_c is not None else None
        parameter = self.dnn.parameters()
        batch_size = self.n_cells if batch_size is None else batch_size
        optimizer = optim.Adam(parameter, amsgrad=True, lr=lr)
        self.dnn.train()

        if self.normalize:
            train_data = data.TensorDataset(self.x_u_norm, self.x_s_norm, self.d_u, self.d_s, torch.arange(self.n_cells))
            train_loader = data.DataLoader(train_data, batch_size=batch_size, shuffle=False)
        else:
            train_data = data.TensorDataset(self.x_u, self.x_s, self.d_u, self.d_s, torch.arange(self.n_cells))
            train_loader = data.DataLoader(train_data, batch_size=batch_size, shuffle=False)

        pbar = tqdm(range(n_epochs))
        for epoch in pbar:
            mse_u_, mse_s_, cos_u_, cos_s_, cos_xu_, cos_xs_ = [], [], [], [], [], []

            for u, s, d_u, d_s, idx in train_loader:
                optimizer.zero_grad()
                v_u, v_s = self.dnn(u, s)

                # cosine similarity for velocity
                norm_v_u = torch.norm(v_u, p=2, dim=1, keepdim=True)
                norm_v_s = torch.norm(v_s, p=2, dim=1, keepdim=True)
                cos_u = torch.mm(v_u, v_u.t()) / (torch.mm(norm_v_u, norm_v_u.t()) + 1e-6)
                cos_s = torch.mm(v_s, v_s.t()) / (torch.mm(norm_v_s, norm_v_s.t()) + 1e-6)

                if same_c is None:
                    n_neighs = torch.mean(torch.sum(self.knn[idx, :][:, idx], dim=1))
                    cos_u = torch.sum(torch.mean(cos_u * self.knn[idx, :][:, idx], dim=1)) / (n_neighs + 1e-6)
                    cos_s = torch.sum(torch.mean(cos_s * self.knn[idx, :][:, idx], dim=1)) / (n_neighs + 1e-6)
                else:
                    n_neighs = torch.mean(torch.sum(same_c * self.knn[idx, :][:, idx], dim=1))
                    cos_u = torch.sum(torch.mean(cos_u * same_c * self.knn[idx, :][:, idx], dim=1)) / (n_neighs + 1e-6)
                    cos_s = torch.sum(torch.mean(cos_s * same_c * self.knn[idx, :][:, idx], dim=1)) / (n_neighs + 1e-6)

                # loss
                mse_u, mse_s = mse(u + v_u, d_u), mse(s + v_s, d_s)
                cos_xu, cos_xs = torch.mean(cos_sim(v_u, d_u - u)), torch.mean(cos_sim(v_s, d_s - s))

                mse_u_.append(mse_u.clone().detach().numpy())
                mse_s_.append(mse_s.clone().detach().numpy())
                cos_u_.append(cos_u.clone().detach().numpy())
                cos_s_.append(cos_s.clone().detach().numpy())
                cos_xu_.append(cos_xu.clone().detach().numpy())
                cos_xs_.append(cos_xs.clone().detach().numpy())

                if loss_fun == "mse":
                    if not only_s:
                        if constrain:
                            loss = mse_u + alpha_1 * mse_s - alpha_2 * cos_u - alpha_2 * cos_s
                        else:
                            loss = mse_u + alpha_1 * mse_s
                        loss.backward()
                        optimizer.step()

                    else:
                        if constrain:
                            loss = mse_s - alpha_2 * cos_s
                        else:
                            loss = mse_s
                        loss.backward()
                        optimizer.step()

                elif loss_fun == "cos":
                    if not only_s:
                        if constrain:
                            loss = -cos_xu - alpha_1 * cos_xs - alpha_2 * cos_u - alpha_2 * cos_s
                        else:
                            loss = -cos_xu - alpha_1 * cos_xs
                        loss.backward()
                        optimizer.step()

                    else:
                        if constrain:
                            loss = -cos_xs - alpha_2 * cos_s
                        else:
                            loss = -cos_xs
                        loss.backward()
                        optimizer.step()

            pbar.set_description("DNN training")
            pbar.set_postfix(mse_u=np.around(np.mean(mse_u_), 3), mse_s=np.around(np.mean(mse_s_), 3),
                             cos_u=np.around(np.mean(cos_u_), 3), cos_s=np.around(np.mean(cos_s_), 3),
                             cos_xu=np.around(np.mean(cos_xu_), 3), cos_xs=np.around(np.mean(cos_xs_), 3))

        # coefficient and velocity
        self.dnn.eval()

        if self.normalize:
            v_u, v_s = self.dnn.forward(self.x_u_norm, self.x_s_norm)
        else:
            v_u, v_s = self.dnn.forward(self.x_u, self.x_s)
                
        alpha, beta, gamma = self.dnn.alpha, self.dnn.beta, self.dnn.gamma
        v_u, v_s = alpha - beta * self.x_u, beta * self.x_u - gamma * self.x_s

        alpha, beta, gamma = alpha.clone().detach().numpy(), beta.clone().detach().numpy(), gamma.clone().detach().numpy()
        v_u, v_s = v_u.clone().detach().numpy(), v_s.clone().detach().numpy()
            
        return v_u, v_s, alpha, beta, gamma
    
    
def get_velocity_rate(adata, d_nn, knn, same_c=None, n_dims=[128, 64], n_epochs=100, method="simple", loss_fun="mse",
                      only_s=False, normalize=True, constrain=True, pretrain=False, lr=1e-3, alpha_1=1, alpha_2=0.1,
                      batch_size=512, save_coeff=True, filter_genes=False):

    if filter_genes:
        velocity_genes = compute_velocity_genes(adata, n_top_genes=2000, inplace=False)
        adata = adata[:, velocity_genes]

    x_u = adata.layers["Mu"]
    x_s = adata.layers["Ms"]


    # DNN implementation
    dnn = DNN(x_u, x_s, d_nn, knn, n_dims=n_dims, normalize=normalize)
    if pretrain:
        dnn.pretrain(n_epochs=n_epochs, batch_size=batch_size, loss_fun="cos", only_s=only_s)
    v_u, v_s, alpha, beta, gamma = dnn.fit(n_epochs=n_epochs, batch_size=batch_size, alpha_1=alpha_1,
                                            alpha_2=alpha_2, lr=lr, constrain=constrain, loss_fun="mse",
                                            only_s=only_s, same_c=same_c)
    if save_coeff:
        adata.layers["alpha"] = alpha
        adata.layers["beta"] = beta
        adata.layers["gamma"] = gamma

    adata.layers["velocity"] = v_s
    adata.layers["velocity_u"] = v_u

    return adata, v_u, v_s

