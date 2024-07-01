import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)


@jit(['Tuple((f8[:,:], f8[:,:]))(f8[:,:], i4, f8)'], nopython=True, fastmath=True, cache=True)
def bin_split_dataset(dataSet, feature, value):
    mat0 = dataSet[np.where(dataSet[:, feature] > value)[0], :]
    mat1 = dataSet[np.where(dataSet[:, feature] <= value)[0], :]
    return mat0, mat1


@jit(nopython=True, fastmath=True, cache=True)
def linear_solve(dataSet):
    m, n = np.shape(dataSet)
    X = np.ones((m, n))
    X[:, 1:n] = dataSet[:, 0:n-1]
    Y = dataSet[:, -1]
    xTx = np.dot(X.T, X)
    if np.linalg.det(xTx) == 0.0:
        raise NameError('This matrix is singular, cannot do inverse, try increasing dur')
    ws = np.dot(np.linalg.inv(xTx), np.dot(X.T, np.copy(Y)))
    return ws, X, Y


@jit(['f8[:](f8[:,:])'], nopython=True, fastmath=True, cache=True)
def model_leaf(dataSet):
    ws, X, Y = linear_solve(dataSet)
    return ws


@jit(['f8(f8[:,:])'], nopython=True, fastmath=True, cache=True)
def model_error(dataSet):
    ws, X, Y = linear_solve(dataSet)
    yHat = np.dot(np.copy(X), np.copy(ws))
    return np.sum(np.power(Y - yHat, 2))


@jit(['Tuple((i4, f8, f8[:], f8, f8))(f8[:,:], i4, i4)'], nopython=True, fastmath=True, cache=True)
def choose_best_split(dataSet, dur, step):
    m, n = np.shape(dataSet)
    S = model_error(dataSet)
    bestS = np.inf
    bestIndex = 0
    bestValue = 0
    for featIndex in range(n-1):
        l, u = int(np.min(dataSet[:, featIndex])), int(np.max(dataSet[:, featIndex]))
        for splitVal in range(l, u+1, step):
            mat0, mat1 = bin_split_dataset(dataSet, featIndex, splitVal)
            if (np.shape(mat0)[0] < dur) or (np.shape(mat1)[0] < dur):
                continue
            newS = model_error(mat0) + model_error(mat1)
            if newS < bestS:
                bestIndex = featIndex
                bestValue = splitVal
                bestS = newS

    return bestIndex, bestValue, model_leaf(dataSet), S, bestS


def is_tree(obj):
    return isinstance(obj, dict)


def model_tree_eval(model, inDat):
    n = np.shape(inDat)[1]
    X = np.ones((1, n+1))
    X[:, 1:n+1] = inDat
    return float(np.dot(X, model))


def tree_forecast(tree, inData):
    if not is_tree(tree):
        return model_tree_eval(tree, inData)
    if inData[tree['spInd']] > tree['spVal']:
        if is_tree(tree['left']):
            return tree_forecast(tree['left'], inData)
        else:
            return model_tree_eval(tree['left'], inData)
    else:
        if is_tree(tree['right']):
            return tree_forecast(tree['right'], inData)
        else:
            return model_tree_eval(tree['right'], inData)


def create_forecast(tree, testData):
    m = len(testData)
    yHat = np.zeros(m)
    for i in range(m):
        yHat[i] = tree_forecast(tree, np.mat(testData[i]))
    return yHat


def draw(dataSet, tree):
    plt.scatter(dataSet[:, 0], dataSet[:, 1], s=5)
    yHat = create_forecast(tree, dataSet[:, 0])
    plt.plot(dataSet[:, 0], yHat, linewidth=2.0, color='red')
    plt.show()


def create_tree(dataSet, rate, dur, step):
    if len(set(dataSet[:, -1])) == 1:
        return model_leaf(dataSet)
    feat, val, ws, S, bestS = choose_best_split(dataSet, dur, step)
    if (S - bestS) < rate:
        return ws
    retTree = {}
    retTree['spInd'] = feat
    retTree['spVal'] = val
    lSet, rSet = bin_split_dataset(dataSet, feat, val)
    retTree['left'] = create_tree(lSet, rate, dur, step)
    retTree['right'] = create_tree(rSet, rate, dur, step)
    return retTree


if __name__ == '__main__':
    import tushare as ts

    df = ts.get_k_data(code='002230', start='2000-01-01')
    e = pd.DataFrame()
    e['idx'] = df.index
    e['close'] = df['close']
    arr = np.array(e)

    tree = create_tree(arr, 100, 10, 1)
    draw(arr, tree)


