import numpy as np
import torch
from torch.utils.data.dataset import TensorDataset
from torch.utils.data import DataLoader
import torch.nn.functional as F
import torch.nn as nn
from torch.autograd import Variable
import numpy as np
import scipy.linalg as slin
import scipy.sparse as sp
import networkx as nx
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import os
import glob
import re
import pickle
import math
from torch.optim.adam import Adam
from sklearn.linear_model import LinearRegression
import sys

# data generating functions

def simulate_random_dag(d: int,
                        degree: float,
                        graph_type: str,
                        w_range: tuple = (0.5, 2.0)) -> nx.DiGraph:
    """Simulate random DAG with some expected degree.

    Args:
        d: number of nodes
        degree: expected node degree, in + out
        graph_type: {erdos-renyi, barabasi-albert, full}
        w_range: weight range +/- (low, high)

    Returns:
        G: weighted DAG
    """
    if graph_type == 'erdos-renyi':
        prob = float(degree) / (d - 1)
        B = np.tril((np.random.rand(d, d) < prob).astype(float), k=-1)
    elif graph_type == 'barabasi-albert':
        m = int(round(degree / 2))
        B = np.zeros([d, d])
        bag = [0]
        for ii in range(1, d):
            dest = np.random.choice(bag, size=m)
            for jj in dest:
                B[ii, jj] = 1
            bag.append(ii)
            bag.extend(dest)
    elif graph_type == 'full':  # ignore degree, only for experimental use
        B = np.tril(np.ones([d, d]), k=-1)
    else:
        raise ValueError('unknown graph type')
    # random permutation
    P = np.random.permutation(np.eye(d, d))  # permutes first axis only
    B_perm = P.T.dot(B).dot(P)
    U = np.random.uniform(low=w_range[0], high=w_range[1], size=[d, d])
    U[np.random.rand(d, d) < 0.5] *= -1
    W = (B_perm != 0).astype(float) * U
    G = nx.DiGraph(W)
    return G


def simulate_sem(G: nx.DiGraph,
                 n: int, x_dims: int,
                 sem_type: str,
                 noise_scale: float = 1.0) -> np.ndarray:
    """Simulate samples from SEM with specified type of noise.

    Args:
        G: weigthed DAG
        n: number of samples
        sem_type: {linear-gauss,linear-exp,linear-gumbel}
        noise_scale: scale parameter of noise distribution in linear SEM

    Returns:
        X: [n,d] sample matrix
    """
    W = nx.to_numpy_array(G)
    d = W.shape[0]
    X = np.zeros([n, d, x_dims])
    ordered_vertices = list(nx.topological_sort(G))
    assert len(ordered_vertices) == d
    for j in ordered_vertices:
        parents = list(G.predecessors(j))
        # eta = (np.sin(X[:, parents])+1.).dot(W[parents, j])  # [n,]
        eta = X[:, parents, 0].dot(W[parents, j])
        if sem_type == 'linear-gauss':
            X[:, j, 0] = eta + np.random.normal(scale=noise_scale, size=n)
        elif sem_type == 'linear-exp':
            X[:, j, 0] = eta + np.random.exponential(scale=noise_scale, size=n)
        elif sem_type == 'linear-gumbel':
            X[:, j, 0] = eta + np.random.gumbel(scale=noise_scale, size=n)
        else:
            raise ValueError('unknown sem type')
    for i in range(x_dims-1):
        X[:, :, i+1] = np.random.normal(scale=noise_scale, size=1)*X[:, :, 0] + np.random.normal(scale=noise_scale, size=1) + np.random.normal(scale=noise_scale, size=(n, d))
    X[:, :, 0] = np.random.normal(scale=noise_scale, size=1) * X[:, :, 0] + np.random.normal(scale=noise_scale, size=1) + np.random.normal(scale=noise_scale, size=(n, d))
 #       for j in ordered_vertices:
 #           parents = list(G.predecessors(j))
#            eta = X[:, parents, i].dot(W[parents, j])
#            if sem_type == 'linear-gauss':
#                X[:, j, i] = eta + np.random.normal(scale=noise_scale, size=n)
#            elif sem_type == 'linear-exp':
#                X[:, j, i] = eta + np.random.exponential(scale=noise_scale, size=n)
#            elif sem_type == 'linear-gumbel':
#                X[:, j, i] = eta + np.random.gumbel(scale=noise_scale, size=n)
#            else:
#                raise ValueError('unknown sem type')
    return X

def import_sem(G: nx.DiGraph,
                 n: int,
                 sem_type: str,
                 noise_scale: float = 1.0) -> np.ndarray:
    """Simulate samples from SEM with specified type of noise.

    Args:
        G: weigthed DAG
        n: number of samples
        sem_type: {linear-gauss,linear-exp,linear-gumbel}
        noise_scale: scale parameter of noise distribution in linear SEM

    Returns:
        X: [n,d] sample matrix
    """
    W = nx.to_numpy_array(G)
    d = W.shape[0]
    X = np.zeros([n, d])


    df1 = pd.read_excel('../Sachs_data/1.xls', header=None)
    d_train1 = df1.as_matrix()
    df2 = pd.read_excel('../Sachs_data/2.xls', header=None)
    d_train2 = df2.as_matrix()
    df3 = pd.read_excel('../Sachs_data/3.xls', header=None)
    d_train3 = df3.as_matrix()
    df4 = pd.read_excel('../Sachs_data/4.xls', header=None)
    d_train4 = df4.as_matrix()
    df5 = pd.read_excel('../Sachs_data/5.xls', header=None)
    d_train5 = df5.as_matrix()
    df6 = pd.read_excel('../Sachs_data/6.xls', header=None)
    d_train6 = df6.as_matrix()
    df7 = pd.read_excel('../Sachs_data/7.xls', header=None)
    d_train7 = df7.as_matrix()
    df8 = pd.read_excel('../Sachs_data/8.xls', header=None)
    d_train8 = df8.as_matrix()
    df9 = pd.read_excel('../Sachs_data/9.xls', header=None)
    d_train9 = df9.as_matrix()
    df10 = pd.read_excel('../Sachs_data/10.xls', header=None)
    d_train10 = df10.as_matrix()
    df11 = pd.read_excel('../Sachs_data/11.xls', header=None)
    d_train11 = df11.as_matrix()
    df12 = pd.read_excel('../Sachs_data/12.xls', header=None)
    d_train12 = df12.as_matrix()
    df13 = pd.read_excel('../Sachs_data/13.xls', header=None)
    d_train13 = df13.as_matrix()
    df14 = pd.read_excel('../Sachs_data/14.xls', header=None)
    d_train14 = df14.as_matrix()

    d_train=np.vstack((d_train1[1:,:], d_train2[1:,:], d_train3[1:,:], d_train4[1:,:], d_train5[1:,:], d_train6[1:,:], d_train7[1:,:], d_train8[1:,:], d_train9[1:,:]))#, d_train10[1:,:], d_train11[1:,:], d_train12[1:,:], d_train13[1:,:], d_train14[1:,:]))

    #return (d_train[0:11650,:]).astype(np.float64)
    return (d_train[0:7464, :]).astype(np.float64)


def simulate_population_sample(W: np.ndarray,
                               Omega: np.ndarray) -> np.ndarray:
    """Simulate data matrix X that matches population least squares.

    Args:
        W: [d,d] adjacency matrix
        Omega: [d,d] noise covariance matrix

    Returns:
        X: [d,d] sample matrix
    """
    d = W.shape[0]
    X = np.sqrt(d) * slin.sqrtm(Omega).dot(np.linalg.pinv(np.eye(d) - W))
    return X


def count_accuracy(G_true: nx.DiGraph,
                   G: nx.DiGraph,
                   G_und: nx.DiGraph = None) -> tuple:
    """Compute FDR, TPR, and FPR for B, or optionally for CPDAG B + B_und.

    Args:
        G_true: ground truth graph
        G: predicted graph
        G_und: predicted undirected edges in CPDAG, asymmetric

    Returns:
        fdr: (reverse + false positive) / prediction positive
        tpr: (true positive) / condition positive
        fpr: (reverse + false positive) / condition negative
        shd: undirected extra + undirected missing + reverse
        nnz: prediction positive
    """
    B_true = nx.to_numpy_array(G_true) != 0
    B = nx.to_numpy_array(G) != 0
    B_und = None if G_und is None else nx.to_numpy_array(G_und)
    d = B.shape[0]
    # linear index of nonzeros
    if B_und is not None:
        pred_und = np.flatnonzero(B_und)
    pred = np.flatnonzero(B)
    cond = np.flatnonzero(B_true)
    cond_reversed = np.flatnonzero(B_true.T)
    cond_skeleton = np.concatenate([cond, cond_reversed])
    # true pos
    true_pos = np.intersect1d(pred, cond, assume_unique=True)
    if B_und is not None:
        # treat undirected edge favorably
        true_pos_und = np.intersect1d(pred_und, cond_skeleton, assume_unique=True)
        true_pos = np.concatenate([true_pos, true_pos_und])
    # false pos
    false_pos = np.setdiff1d(pred, cond_skeleton, assume_unique=True)
    if B_und is not None:
        false_pos_und = np.setdiff1d(pred_und, cond_skeleton, assume_unique=True)
        false_pos = np.concatenate([false_pos, false_pos_und])
    # reverse
    extra = np.setdiff1d(pred, cond, assume_unique=True)
    reverse = np.intersect1d(extra, cond_reversed, assume_unique=True)
    # compute ratio
    pred_size = len(pred)
    if B_und is not None:
        pred_size += len(pred_und)
    cond_neg_size = 0.5 * d * (d - 1) - len(cond)
    fdr = float(len(reverse) + len(false_pos)) / max(pred_size, 1)
    tpr = float(len(true_pos)) / max(len(cond), 1)
    fpr = float(len(reverse) + len(false_pos)) / max(cond_neg_size, 1)
    # structural hamming distance
    B_lower = np.tril(B + B.T)
    if B_und is not None:
        B_lower += np.tril(B_und + B_und.T)
    pred_lower = np.flatnonzero(B_lower)
    cond_lower = np.flatnonzero(np.tril(B_true + B_true.T))
    extra_lower = np.setdiff1d(pred_lower, cond_lower, assume_unique=True)
    missing_lower = np.setdiff1d(cond_lower, pred_lower, assume_unique=True)
    shd = len(extra_lower) + len(missing_lower) + len(reverse)
    return fdr, tpr, fpr, shd, pred_size





'''
COMPUTE SCORES FOR BN
'''
def compute_BiCScore(G, D, is_discrete=True):
    '''compute the bic score'''
    # score = gm.estimators.BicScore(self.data).score(self.model)
    origin_score = []
    num_var = G.shape[0]
    for i in range(num_var):
        parents = np.where(G[:,i] !=0)
        if is_discrete:
            score_one = compute_local_BiCScore(D, i, parents)
        else:
            score_one = compute_local_BiCScore_gauss(D, i, parents)
        origin_score.append(score_one)

    score = sum(origin_score)

    return score


def compute_local_BiCScore_gauss(np_data, target, parents):
    # use dictionary
    sample_size = np_data.shape[0]
    var_size = np_data.shape[1]

    # build dictionary and populate
    count_d = dict()
    if len(parents) < 1:
        a = 1

    if parents == []:
        sigma_sq = np.var(np_data[:, target]) ** 2
        mu = np.mean(np_data[:, target])
        loglik = - sample_size / 2 * np.log(2 * np.pi) - sample_size * np.log(np.sqrt(sigma_sq)) \
                 - 1 / (2 * sigma_sq) * np.dot(
                (np_data[:, target] - mu), \
                (np_data[:, target] - mu))

    else:
        reg = LinearRegression().fit(np_data[:,parents], np_data[:,target])
        w, w_0 = reg.coef_, reg.intercept_
        sigma_sq = 1/ sample_size * (
                np.dot((np_data[:, target] - (w_0 + np.dot(np_data[:, parents], w))),(np_data[:, target] - (w_0 + np.dot(np_data[:, parents], w))))
                )
        loglik = - sample_size / 2 * np.log(2 * np.pi) - sample_size * np.log(np.sqrt(sigma_sq)) \
                 - 1 / (2 * sigma_sq) * np.dot(
            (np_data[:, target] - (w_0 + np.dot(np_data[:, parents], w))), \
            (np_data[:, target] - (w_0 + np.dot(np_data[:, parents], w))))

    if sigma_sq == 0:
        print('error in sigma')

    # compute likelihood


    # penality
    num_param = np.size(parents) + 2  # count_faster(count_d) - len(count_d) - 1 # minus top level and minus one
    bic = loglik - 0.5 * math.log(sample_size) * num_param

    return bic


def compute_local_BiCScore(np_data, target, parents):
    # use dictionary
    sample_size = np_data.shape[0]
    var_size = np_data.shape[1]

    # build dictionary and populate
    count_d = dict()
    if len(parents) < 1:
        a = 1

    # unique_rows = np.unique(self.np_data, axis=0)
    # for data_ind in range(unique_rows.shape[0]):
    #     parent_combination = tuple(unique_rows[data_ind,:].reshape(1,-1)[0])
    #     count_d[parent_combination] = dict()
    #
    #     # build children
    #     self_value = tuple(self.np_data[data_ind, target].reshape(1,-1)[0])
    #     if parent_combination in count_d:
    #         if self_value in count_d[parent_combination]:
    #             count_d[parent_combination][self_value] += 1.0
    #         else:
    #             count_d[parent_combination][self_value] = 1.0
    #     else:
    #         count_d[parent_combination] = dict()
    #         count_d

    # slower implementation
    for data_ind in range(sample_size):
        parent_combination = tuple(np_data[data_ind, parents].reshape(1, -1)[0])
        self_value = tuple(np_data[data_ind, target].reshape(1, -1)[0])
        if parent_combination in count_d:
            if self_value in count_d[parent_combination]:
                count_d[parent_combination][self_value] += 1.0
            else:
                count_d[parent_combination][self_value] = 1.0
        else:
            count_d[parent_combination] = dict()
            count_d[parent_combination][self_value] = 1.0

    # compute likelihood
    loglik = 0.0
    # for data_ind in range(sample_size):
    # if len(parents) > 0:
    num_parent_state = np.prod(np.amax(np_data[:, parents], axis=0) + 1)
    # else:
    #    num_parent_state = 0
    num_self_state = np.amax(np_data[:, target], axis=0) + 1

    for parents_state in count_d:
        local_count = sum(count_d[parents_state].values())
        for self_state in count_d[parents_state]:
            loglik += count_d[parents_state][self_state] * (
                        math.log(count_d[parents_state][self_state] + 0.1) - math.log(local_count))

    # penality
    num_param = num_parent_state * (
                num_self_state - 1)  # count_faster(count_d) - len(count_d) - 1 # minus top level and minus one
    bic = loglik - 0.5 * math.log(sample_size) * num_param

    return bic

#############
#latent variable graph
def compute_BiC_c_comps_sets(D, c_components, parent_sets, bidirect_tuple):
    '''compute the bic score'''
    # TODO fix this
    nVars = D.shape[1]
    nSamples = D.shape[0]

    # find all c components in graph
    # c_components, sizes, members = find_c_components(G)
    isParent = np.zeros(nVars)

    scores = np.zeros( (1, len(c_components)))

    covMat = np.cov(D.T)
    tol = 10e-2

    # for iComp in c_components:
        # vars = c_components[iComp]

        # convert to mag
    [compMag, district, parents] = find_componenet_mag_sets(c_components, nVars, parent_sets, bidirect_tuple)

    district = np.where(district)[0] # list(set(c_components) - set(parent_sets)) + list(set(parent_sets))
    # scores(iComp) = logdet(2 * pi * covMat(district, district)) + (nSamples - 1) / nSamples;
    scores = compute_local_BiC_comp_sets(compMag, c_components, district, len(c_components), covMat, nSamples, tol)

    # score_one = compute_local_BiC_c_component(D, iComp, vars)
    # origin_score.append(score_one)

    score = np.sum(scores)

    return score

def compute_local_BiC_comp_sets(compMag, component, district, compSize, CovMat, nSamples,  tol):

    if district.size == 1:
        comCovMat = CovMat[district, district]
        sign, logdet = 1, np.log(2 * math.pi * comCovMat) #  np.linalg.slogdet(2 * math.pi * comCovMat)
        sc = sign * logdet + (nSamples - 1) / nSamples
    else:
        nVar = CovMat.shape[0]
        curCovMat = CovMat[np.ix_(district, district)]
        compMag = compMag[np.ix_(district, district)]
        _, _, curHatCovMat, _ = RICF_fit(compMag, curCovMat, tol)

        remParents = np.zeros(nVar)
        remParents[district] = 1
        remParents[list(component)] = 0
        parInds = np.where(remParents[district].astype(int))[0]

        sign, logdet = np.linalg.slogdet(curHatCovMat)
        logdet_signed = sign* logdet
        if np.sum(remParents) > 0:
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed - math.log(np.prod(np.diag(curHatCovMat[np.ix_(parInds, parInds)])))
            l3 = (nSamples - 1) / nSamples * (np.trace( matrix_division_left(curHatCovMat,curCovMat)) - np.sum(remParents))
            sc = l1 + l2 +l3
        else:
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed
            l3 = (nSamples - 1) / nSamples * np.trace(matrix_division_left(curHatCovMat,curCovMat))
            sc = l1 + l2 + l3

    # complexity
    # nsf = -nSamples / 2
    # curScore = -2 * tmpSll + log(nSamples) * (nVars + nEdges)

    nVars  = np.size(component)
    # nEdges = 0
    nEdges = np.where(compMag)[0].size/2
    bp = math.log(nSamples)/2 * (2 * nVars + nEdges)
    sc = - sc * nSamples/2 - bp

    return sc


def find_componenet_mag_sets(component, nVars, parent_sets, bidirect_tuple):
    compMag = np.zeros((nVars,nVars))
    if len(component) > 1:  # fix this for more than size 2, per edges
        # compMag[component[0], component[1]] =  2 # assume 1 here - mag[component, component]
        # compMag[component[1], component[0]] = 2
        for edge_ind, (parent, child) in enumerate(bidirect_tuple):
            compMag[parent, child] = 2
            compMag[child, parent] = 2

    district, parents = np.zeros(( nVars,)), np.zeros((nVars,))
    district[list(component)] = True
    for index, iVar in enumerate(component):
        iParents = list(parent_sets[index]) #isParent[:, iVar]
        compMag[iParents, iVar] = 2
        compMag[iVar, iParents] = 3
        district[iParents] = True
        parents[iParents] = True

    return compMag, district, parents

def compute_BiC_c_component(G, D):
    '''compute the bic score'''
    # score = gm.estimators.BicScore(self.data).score(self.model)
    origin_score = []
    nVars = G.shape[0]
    nSamples = D.size()[0]

    # find all c components in graph
    c_components, sizes, members = find_c_components(G)
    isParent = np.zeros(nVars)

    scores = np.zeros( (1, len(c_components)))

    covMat = np.cov(D.T)
    tol = 10e-3

    for iComp in c_components:
        vars = c_components[iComp]

        # convert to mag
        [compMag, district] = find_componenet_mag(vars, nVars, G, isParent)
        # scores(iComp) = logdet(2 * pi * covMat(district, district)) + (nSamples - 1) / nSamples;
        scores[iComp] = compute_local_BiC_c_component(compMag, vars, district, sizes(iComp), covMat, nSamples, tol)

        # score_one = compute_local_BiC_c_component(D, iComp, vars)
        # origin_score.append(score_one)

    score = sum(scores)

    return score

def find_componenet_mag(component, nVars, mag, isParent):
    compMag = np.zeros((nVars,nVars))
    compMag[component, component] = mag[component, component]
    district, parents = np.zeros(nVars), np.zeros( nVars)
    district[component] = True
    for iVar in component:
        iParents = isParent[:, iVar]
        compMag[iParents, iVar] = 2
        compMag[iVar, iParents] = 3
        district[iParents] = True
        parents[iParents] = True

    return compMag, district, parents

def find_c_components(A):
    #
    # Number of nodes
    N = A.size()[0] # size(A, 1)
    # Remove diagonals
    for i in range(N):
        A[i, i] = 0 # A(1: N + 1:end) = 0
    # make symmetric, just in case i isn't
    A = A + np.transpose(A)
    # Have we visited a particular node yet?
    isDiscovered = np.zeros(N)
    inComponent = np.zeros( N)
    # Empty members cell
    members = dict.fromkeys(range(N))
    # check every node
    for n in range(N):
        if ~isDiscovered[n]:
            # started a new group so add  it to  members
            members[n]= n
            inComponent[n] = n # members[n]  # size(members, 2)
            # account for discovering n
            isDiscovered[n] = 1
            # set the ptr to 1
            ptr = 0
            while (ptr < len(members[n])):
                # find neighbors
                nbrs = A[:, members[n][ptr]] # np.where(A[:, members[n][ptr]])
                # here are the neighbors that are undiscovered
                newNbrs = nbrs[isDiscovered[nbrs] == 0]
                # we can now mark them as discovered
                isDiscovered[newNbrs] = 1
                # add them to member list
                members[n][n + 1: n + len(newNbrs)] = newNbrs
                inComponent[nbrs] = n # size(members, 2)
                # increment ptr so  we check the next member of this component
                ptr = ptr + 1

    # number of  components
    nComponents = len(members)
    sizes = np.zeros(nComponents)
    for n in range(nComponents):
        # compute sizes of components
        sizes[n] = len(members[n])

    #[sizes, idx] = sort(sizes, 'descend');
    # members = members(idx);


    return nComponents, sizes, members

def compute_local_BiC_c_component(compMag, component, district, compSize, covMat, nSamples,  tol):

    if sum(district) == 1:
        sc = np.linalg.slogdet(2 * math.pi * covMat[district, district]) + (nSamples - 1) / nSamples
    else:
        curCovMat = covMat[district, district]
        compMag = compMag[district, district]
        _, _, curHatCovMat, _ = RICF_fit(compMag, curCovMat, tol)
        remParents = district
        remParents[component] = False
        parInds = remParents[district]

        sign, logdet = np.linalg.slogdet(curHatCovMat)
        logdet_signed = sign* np.exp(logdet)
        if any(remParents):
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed - math.log(np.prod(np.diag(curHatCovMat[parInds, parInds])))
            l3 = (nSamples - 1) / nSamples * (np.trace( matrix_division_left(curHatCovMat,curCovMat) - np.sum(remParents)))
            sc = l1 + l2 +l3
        else:
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed
            l3 = (nSamples - 1) / nSamples * np.trace(matrix_division_left(curHatCovMat,curCovMat))
            sc = l1 + l2 + l3

    return sc

def compute_local_BiC_c_component_node(compMag, component, district, compSize, covMat, nSamples,  tol):

    if sum(district) == 1:
        sc = np.linalg.slogdet(2 * math.pi * covMat[district, district]) + (nSamples - 1) / nSamples
    else:
        curCovMat = covMat[district, district]
        compMag = compMag[district, district]
        _, _, curHatCovMat, _ = RICF_fit(compMag, curCovMat, tol)
        remParents = district
        remParents[component] = False
        parInds = remParents[district]

        sign, logdet = np.linalg.slogdet(curHatCovMat)
        logdet_signed = sign* np.exp(logdet)
        if any(remParents):
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed - math.log(np.prod(np.diag(curHatCovMat[parInds, parInds])))
            l3 = (nSamples - 1) / nSamples * (np.trace( matrix_division_left(curHatCovMat,curCovMat) - np.sum(remParents)))
            sc = l1 + l2 +l3
        else:
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed
            l3 = (nSamples - 1) / nSamples * np.trace(matrix_division_left(curHatCovMat,curCovMat))
            sc = l1 + l2 + l3

    return sc

def RICF_fit(smm, covMat, tol):

    # % if any(eig(covMat)<= 0);
    # %     errprintf('Covariance matrix is not positive definite\n')
    # % end
    if (smm==4).any():
        sys.exit('Graph includes bows\n')

    # a=0; b=0; c=0; d=0; e=0;

    nVars = covMat.shape[1]
    dg = np.multiply(smm == 2, smm.T == 3)
    bg = np.multiply(smm == 2, smm.T == 2)
    # %
    # % dg = (smm ==2 & smm'==3) | (smm==2 & smm'==4);
    # % bg = (smm ==2 & smm'==4) | (smm==2 & smm'==2) | (smm==4 & smm'==2);

    # starting values
    omega = np.diag(np.diag(covMat),0)
    beta = np.eye(nVars)

    # list of parents, spouses etc
    par, sp = np.zeros((nVars, nVars)), np.zeros((nVars, nVars)) #deal(false(nVars));
    for iVar in range(nVars):
        par[iVar, dg[:, iVar]] = True
        sp[iVar, bg[:, iVar]] = True
    iter=0

    ricf = {}
    while True:
        iter = iter+1
        ricf[iter] = {}
        ricf[iter]['omega'] = omega.copy()
        ricf[iter]['beta'] = beta.copy()

        omega = omega.copy()
        beta = beta.copy()
        # for each variable
        for iVar in range(nVars):

           vcomp = list(range(0,iVar)) + list(range(iVar+1, nVars))
           iPar= np.where(par[iVar,:])[0]
           iSp = np.where(sp[iVar, :])[0]

           if iSp.size == 0:
               if iPar.size > 0:
                   if iter==1:
                       if len(vcomp) == 1:
                           inv_mat = 1 /  (covMat[iPar, iPar])
                       else:
                           inv_mat =   np.linalg.inv(covMat[np.ix_(iPar, iPar)])
                       beta[iVar, iPar] = np.matmul(-covMat[iVar, iPar],inv_mat) # TODO check this
                       omega[iVar, iVar] = covMat[iVar, iVar] + np.matmul(beta[iVar, iPar],covMat[iPar, iVar])

           elif iPar.size > 0: # spouses and parents
               oInv = np.zeros((nVars, nVars))
               if len(vcomp) == 1:
                   oInv[vcomp, vcomp] = 1 / oInv[vcomp, vcomp]
               else:
                   oInv[np.ix_(vcomp, vcomp)] = np.linalg.inv(omega[np.ix_(vcomp, vcomp)]) #; %\Omega_{-i, -i}^-1
               Z = np.matmul(oInv[np.ix_(iSp, vcomp)], beta[vcomp, :])
               if Z.ndim == 1:
                   Z = Z.reshape(1, -1)
               nPar = iPar.size
               nSp = iSp.size
               range1 = list(range(nPar))
               range2= list(range(nPar, nPar+nSp))
               # % XX
               XX = np.zeros( (nPar+nSp,(nPar+nSp)) )
               XX[:] = np.nan
               #% Upper left quadrant
               XX[np.ix_(range1, range1)] = covMat[np.ix_(iPar, iPar)]
               #% Upper right quadrant
               XX[np.ix_(range1, range2)] = np.matmul(covMat[iPar, :], Z.T)
               #% Lower left quadrant
               XX[np.ix_(range2, range1)] = XX[np.ix_(range(nPar), range(nPar,nPar+nSp))].T
               #% Lower right quadrant
               XX[np.ix_(range2, range2)] = np.matmul(np.matmul(Z,covMat),Z.T)
               #% YX <- c(S[v,parv], S[v,]%*%t(Z))
               #% temp <- YX %*% solve(XX)
               YX = np.hstack((covMat[iVar, iPar].reshape(1,-1), np.array([np.matmul(covMat[iVar, :],Z.T)]).reshape(1,-1)
                                    )).reshape((-1,1))
               # YX = np.array( [covMat[iVar, iPar], np.matmul(covMat[iVar, :], Z.T)]) #.T

               if YX.size == 1:
                    temp = YX / XX
                    beta[iVar, iPar] = -temp[range1]
                    omega[iVar, iSp] = temp[range2]
                    omega[iSp, iVar] = omega[iVar, iSp]
                    tempVar = covMat[iVar, iVar] - temp * YX
                    omega[iVar, iVar] = tempVar + np.matmul(omega[iVar, iSp]/ omega[np.ix_(iSp, iSp)] \
                                                            , omega[iSp, iVar])
               else:
                   temp = matrix_division_right(YX.T, XX) #np.lin.solve(YX.dot(YX.T), YX).dot(XX)  # YX.T/XX # TODO check
                   # % update beta, omega
                   beta[iVar, iPar] = -temp[0,range1]
                   omega[iVar, iSp] = temp[0, range2]
                   omega[iSp, iVar] = omega[iVar, iSp]

                   tempVar = covMat[iVar, iVar] - np.matmul(temp, YX)

                   if iSp.size == 1:
                       omega[iVar, iVar] = tempVar + omega[iVar, iSp]/ omega[np.ix_(iSp, iSp)] * omega[iSp, iVar]
                   else:
                        omega[iVar, iVar] = tempVar + np.matmul(matrix_division_right(omega[iVar, iSp], omega[np.ix_(iSp, iSp)])\
                                       ,omega[iSp, iVar])

           else:
               oInv = np.zeros((nVars, nVars))
               if len(vcomp) == 1:
                    oInv[vcomp, vcomp] = 1/omega[vcomp, vcomp]
               else:
                    oInv[np.ix_(vcomp, vcomp)] = np.linalg.inv(omega[ np.ix_(vcomp, vcomp)]) # %\Omega_{-i, -i}^-1
               Z = np.matmul(oInv[np.ix_(iSp, vcomp)],beta[vcomp, :])
               XX = np.matmul( np.matmul(Z,covMat), Z.T)
               YX = np.matmul(covMat[iVar, :], Z.T).T

               if YX.size== 1:
                   omega[iVar, iSp] =  YX / XX
                   omega[iSp, iVar] = omega[iVar, iSp]
                   tempVar = covMat[iVar, iVar] -  omega[[iVar],[iSp]] * YX
                             # omega[iVar, iSp] * YX
                   omega[iVar, iVar] = tempVar + \
                                       np.matmul(np.matmul(omega[[iVar], [iSp]], oInv[np.ix_(iSp,iSp)]), omega[iSp, [iVar]])
                                       # omega[iVar, iSp] * oInv[iSp, iSp] * omega[iSp, iVar]

               else:
                   omega[iVar, iSp] =  matrix_division_right(YX.T, XX) # YX'/XX;
                   omega[iSp, iVar] = omega[iVar, iSp]
                   tempVar = covMat[iVar, iVar] - np.matmul(omega[iVar,iSp],YX)
                   omega[iVar, iVar] = tempVar+ \
                                       np.matmul(np.matmul(omega[iVar, iSp],oInv[np.ix_(iSp,iSp)]),omega[iSp, iVar])

        if (sum(sum(abs(ricf[iter]['omega']-omega))) + sum(sum(abs(ricf[iter]['beta']-beta))) < tol):
           break

        elif iter>50:
           print('RICF fit did not converge')
           break

    # %    if iter>20
    # %        display('20')
    # %        iter;
    # %    end
    hatCovMat = np.matmul( np.matmul(np.linalg.inv(beta), omega), np.linalg.inv(beta.T))
    # % load('maria.mat')
    # % kk = vertcat(kk,[a b c d e]);
    # % [a b c d]
    # % save maria.mat kk

    return beta, omega, hatCovMat, ricf

def matrix_division_right(YX, XX):
    # solve for YX/XX or x XX = YX

    x = np.matmul(np.matmul(YX, XX.T), np.linalg.inv(np.matmul(XX, XX.T)))

    return x


def matrix_division_left(B, b):
    # solves for B\b or Bx = b
    x, resid, rank, s = np.linalg.lstsq(B,b)

    # x = np.matmul(np.matmul(YX, XX.T), np.linalg.inv(np.matmul(XX, XX.T)))

    return x


def score_write_to_file(scores, file_name):
    # write scores to find
    with open(file_name, 'wb') as handle:
        pickle.dump(scores, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return

def check_connected_component(set, node_ids):
    # set is a list of tuples, which consists a tuple of edge tupule
    new_set = []

    # single set
    if len(set)  < 2:
        return set

    for edge_lists in set:
        # check connected
        flag = True
        for node in node_ids:
            if not (any(node in i for i in edge_lists)):
                flag = False
        if flag:
            new_set.append(edge_lists)

    return new_set

def find_bi_connected_node(iVar, comp_edges):
    # find nodes that bi-direct connect to iVar in comp edges

    nodes = []

    if len(comp_edges)  < 2:
        return nodes

    for edge_lists in comp_edges:
        if iVar in  edge_lists:
            bi_node = set( list(edge_lists)) - set([iVar])
            nodes.append( list(bi_node)[0])

    return nodes

def has_cycle_in_c_comp(parent_set_config, vars_in_comp):
    # check if parent set has a directed cycle of length 2 only inside a c-comp
    num_var_parent = len(parent_set_config)
    vars_in_comp_list = list(vars_in_comp)
    flag = False

    for index_var, iVar in enumerate(vars_in_comp_list): # range(num_var):
        iParents = list(parent_set_config[index_var])
        for each_parent in iParents:
            if each_parent in vars_in_comp_list:
                each_parent_index = vars_in_comp_list.index(each_parent)
                # if num_var_parent > each_parent:  # has data for that parent
                par_parents = list(parent_set_config[each_parent_index])
                if iVar in par_parents:
                    flag = True
                    break

        if flag == True:
            break

    return flag