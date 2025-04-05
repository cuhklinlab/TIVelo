import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from .linear_tree import create_tree, create_forecast



def get_slope(tree, data):
    if not isinstance(tree, dict):
        return tree[1]
    if data[tree['spInd']] > tree['spVal']:
        if isinstance(tree, dict):
            return get_slope(tree['left'], data)
        else:
            return tree[1]
    else:
        if isinstance(tree, dict):
            return get_slope(tree['right'], data)
        else:
            return tree[1]


def rle(inarray):
    # for a sequence, return the corresponding run-length, position and element
    ia = np.asarray(inarray)
    n = len(ia)
    if n == 0:
        return None, None, None
    else:
        y = np.array(ia[1:] != ia[:-1])
        i = np.append(np.where(y), n - 1)
        z = np.diff(np.append(-1, i))
        p = np.cumsum(np.append(0, z))[:-1]
        return z, p, ia[i]


def orient_score(u, s, rate=1e-10, dur=100, kernel_size=100, step=1):
    n_cells = len(u)

    # preprocess
    u = np.convolve(u, np.ones(kernel_size) / kernel_size, mode="same")
    s = np.convolve(s, np.ones(kernel_size) / kernel_size, mode="same")
    if np.sum(u) > 0:
        u = u / np.sum(u)
    if np.sum(s) > 0:
        s = s / np.sum(s)

    x = np.arange(n_cells)
    u_arr = np.concatenate((x[:, None], u[:, None]), axis=1)
    s_arr = np.concatenate((x[:, None], s[:, None]), axis=1)

    # create tree
    u_tree = create_tree(u_arr, rate, dur, step)
    s_tree = create_tree(s_arr, rate, dur, step)

    # calculate slope
    u_slope, s_slope = np.zeros(n_cells), np.zeros(n_cells)
    for i in range(n_cells):
        u_slope[i] = get_slope(u_tree, u_arr[i, :])
        s_slope[i] = get_slope(s_tree, s_arr[i, :])

    # select sections
    u_rise, s_rise = u_slope > 0, s_slope > 0
    u_fall, s_fall = u_slope < 0, s_slope < 0
    rise_sec = u_rise & s_rise
    fall_sec = u_fall & s_fall

    # calculate run length, normalized reads, weight of location for rising section
    zr, pr, er = rle(rise_sec)
    pr_ = np.append(pr, n_cells)
    rise_l, rise_w, u_norm, s_norm = np.zeros(n_cells), np.zeros(n_cells), np.zeros(n_cells), np.zeros(n_cells)
    for i in range(len(pr)):
        rise_l[pr_[i]: pr_[i+1]] = zr[i]
        rise_w[pr_[i]: pr_[i+1]] = pr[i]
        u_norm[pr_[i]: pr_[i+1]] = u[pr_[i]: pr_[i+1]] - u[pr_[i]]
        s_norm[pr_[i]: pr_[i + 1]] = s[pr_[i]: pr_[i + 1]] - s[pr_[i]]

    rise_diff = u - s

    # calculate run length, normalized reads, weight of location for falling section
    zf, pf, ef = rle(fall_sec)
    pf_ = np.append(pf, n_cells)
    fall_l, fall_w, u_norm, s_norm = np.zeros(n_cells), np.zeros(n_cells), np.zeros(n_cells), np.zeros(n_cells)
    for i in range(len(pf)):
        fall_l[pf_[i]: pf_[i+1]] = zf[i]
        fall_w[pf_[i]: pf_[i+1]] = pf[i]
        u_norm[pf_[i]: pf_[i+1]] = u[pf_[i]: pf_[i+1]] - u[pf_[i]]
        s_norm[pf_[i]: pf_[i + 1]] = s[pf_[i]: pf_[i + 1]] - s[pf_[i]]

    fall_diff = u - s

    # calculate orientation scores
    pos_score = np.sum(rise_sec * rise_l * rise_diff)
    neg_score = np.sum(fall_sec * fall_l * fall_diff)
    score = pos_score - neg_score

    # calculate section score
    rise_p_ = [pr[i] for i in np.arange(len(pr)) if er[i] == 1]
    rise_l_ = [zr[i] for i in np.arange(len(zr)) if er[i] == 1]
    pos_score_ = [np.sum(((u_rise & s_rise) * rise_l * rise_diff)[pr_[i]: pr_[i+1]]) for i in
                  np.arange(len(pr)) if er[i] == 1]
    fall_p_ = [pf[i] for i in np.arange(len(pf)) if ef[i] == 1]
    fall_l_ = [zf[i] for i in np.arange(len(zf)) if ef[i] == 1]
    neg_score_ = [np.sum(((u_fall & s_fall) * fall_l * fall_diff)[pf_[i]: pf_[i + 1]]) for i in
                   np.arange(len(pf)) if ef[i] == 1]

    sec_list = [rise_sec, fall_sec, rise_p_, fall_p_, rise_l_, fall_l_, pos_score_, neg_score_]

    return score, u_arr, s_arr, u_tree, s_tree, sec_list


def visualize(u_arr, s_arr, u_tree, s_tree, sec_list, fig_size=(12, 9)):
    u_hat = create_forecast(u_tree, u_arr[:, 0])
    s_hat = create_forecast(s_tree, s_arr[:, 0])
    rise_sec, fall_sec, rise_p_, fall_p_, rise_l_, fall_l_, pos_score_, neg_score_ = sec_list

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=fig_size, sharex=True)

    ax1.plot(u_arr[:, 0], u_arr[:, 1], c="blue", label="Mu")
    ax1.plot(s_arr[:, 0], s_arr[:, 1], c="red", label="Ms")

    # search for the beginning points
    n_cells = len(u_arr)
    u_slope, s_slope = np.zeros(n_cells), np.zeros(n_cells)
    for i in range(n_cells):
        u_slope[i] = get_slope(u_tree, u_arr[i, :])
        s_slope[i] = get_slope(s_tree, s_arr[i, :])
    zu, pu, eu = rle(u_slope)
    zs, ps, es = rle(s_slope)
    pu_ = np.append(pu, n_cells)
    ps_ = np.append(ps, n_cells)
    # plot
    for i in np.arange(len(pu)):
        if i == 0:
            ax2.plot(u_arr[pu_[i]:pu_[i + 1], 0], np.maximum(u_hat[pu_[i]:pu_[i + 1]], 0), c="blue", label="Fitted Mu")
        else:
            ax2.plot(u_arr[pu_[i]:pu_[i + 1], 0], np.maximum(u_hat[pu_[i]:pu_[i + 1]], 0), c="blue")
        ax2.plot(u_arr[pu[i]-1:pu[i]+1, 0], np.maximum(u_hat[pu[i]-1:pu[i]+1], 0), c="blue", ls="dotted")
    for i in np.arange(len(ps)):
        if i == 0:
            ax2.plot(s_arr[ps_[i]:ps_[i + 1], 0], np.maximum(s_hat[ps_[i]:ps_[i + 1]], 0), c="red", label="Fitted Ms")
        else:
            ax2.plot(s_arr[ps_[i]:ps_[i + 1], 0], np.maximum(s_hat[ps_[i]:ps_[i + 1]], 0), c="red")
        ax2.plot(s_arr[ps[i]-1:ps[i]+1, 0], np.maximum(s_hat[ps[i] - 1:ps[i]+1], 0), c="red", ls="dotted")

    max_val = max(np.max(u_arr[:, 1]), np.max(s_arr[:, 1]))
    ax1.fill_between(np.arange(len(rise_sec)), 0, max_val, where=rise_sec, color="orange", alpha=0.1, label="Rising")
    ax1.fill_between(np.arange(len(fall_sec)), 0, max_val, where=fall_sec, color="green", alpha=0.1, label="Falling")
    # ax1.text(0, max_val, "Total Score: {:.2f}".format(score))
    for i in range(len(rise_p_)):
        ax2.text(rise_p_[i] + rise_l_[i] / 2, 0, "{:.2f}".format(pos_score_[i]), ha="center", va="bottom", fontsize=12)
    for i in range(len(fall_p_)):
        ax2.text(fall_p_[i] + fall_l_[i] / 2, 0, "{:.2f}".format(-neg_score_[i]), ha="center", va="bottom", fontsize=12)

    max_val = max(np.max(u_hat), np.max(s_hat))
    ax2.fill_between(np.arange(len(rise_sec)), 0, max_val, where=rise_sec, color="orange", alpha=0.1)
    ax2.fill_between(np.arange(len(fall_sec)), 0, max_val, where=fall_sec, color="green", alpha=0.1)

    return fig, (ax1, ax2)











