import dynamo as dyn
import scvelo as scv
import numpy as np
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import euclidean_distances
from scipy.spatial import distance_matrix
from collections import Counter
import re
from sklearn.decomposition import PCA
from sklearn import metrics
import pandas as pd
import scvelo as scv
from kneed import KneeLocator
import matplotlib.pyplot as plt
from scipy import interpolate
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
import sklearn
import sklearn.covariance
import scipy
import scipy.sparse.linalg
import numpy.linalg
import sklearn.metrics
import sklearn.model_selection
import config
import os

def get_optimal_K(X, kmin=1, kmax=21):
    distortions = []
    inertias = []
    mapping1 = {}
    k2inertia = {}
    K = range(kmin, kmax)
    for k in K:
        # Building and fitting the model
        kmeanModel = KMeans(n_clusters=k, random_state=0).fit(X)
        kmeanModel.fit(X)

        distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_,
                                            'euclidean'), axis=1)) / X.shape[0])
        inertias.append(kmeanModel.inertia_)

        mapping1[k] = sum(np.min(cdist(X, kmeanModel.cluster_centers_,
                                       'euclidean'), axis=1)) / X.shape[0]
        k2inertia[k] = kmeanModel.inertia_
    kn = KneeLocator(K, distortions, curve='convex', direction='decreasing')
    return kn.knee


def read_list(path):
    with open(path, 'r') as f:
        data = []
        for line in f:
            data.append(line.replace('\n', ''))
        return data

    
def make_dir(path, abort=True):
    if os.path.exists(path):
        print(path + ' : exists')
        if abort:
            exit(0)
        elif os.path.isdir(path):
            print(path + ' : is a directory, continue using the old one')
            return False
        else:
            print(path + ' : is not a directory, creating one')
            os.makedirs(path)
            return True
    else:
        os.makedirs(path)
        return True

    
def gen_two_gene_graph(adata, g1, g2, save_dir):
    print([g1, g2])
    sub_data = adata[:, [g1, g2]]
    # scv.pp.filter_and_normalize(sub_data)
    # scv.pp.moments(sub_data)
    # scv.tl.velocity(sub_data)
    # scv.tl.velocity_graph(sub_data)
    filename = g1 + '_' + g2 + '_vel_stream.png'
    path = os.path.join(save_dir,
                        filename)
    scv.pl.velocity_embedding_stream(sub_data, show=False, save=filename)
    
    

def gen_all_gene_pair_vector_field(adata, gene_list, save_dir=config.two_gene_graph_dir):
    make_dir(save_dir, abort=False)
    for i in range(len(gene_list)):
        for j in range(i+1, len(gene_list)):
            g1, g2 = gene_list[i], gene_list[j]
            if g1 in adata.var_names and g2 in adata.var_names:
                gen_two_gene_graph(adata, g1, g2, save_dir)
            else:
                print('both or one of (%s, %s) not in adata gene list' % (g1, g2)) 


def analyze_specific_cluster(adata, indices,
                             predicted_velocities, errors,
                             model, cluster_label,
                             pca_model=None, emt_genes=None,
                             selected_genes=config.selected_genes_jacobian):
    '''
    adata: annData obj
    indices: true or false indicating a sample belonging to a cluster or not
    model: lr model
    cluster label: cluster id/label
    errors: (true-pred)^2, a vector for each observation
    '''
    cluster_data = adata[indices, :]

    # use velocities to define centroids
    vel_norms = []
    for i in range(len(cluster_data)):
        # calculate norm of velocities by sparse matrix API or numpy API
        # Note: somehow velocity is not sparse by default
        # vel_norm = scipy.sparse.linalg.norm(cluster_data.layers['velocity'][i])
        vel_norm = numpy.linalg.norm(cluster_data.layers['velocity'][i])
        vel_norms.append(vel_norm)
    cluster_data.obs['vel_norms'] = vel_norms

    # use MSE and distance to define centroids
    N = len(cluster_data)
    dist = np.zeros((N, N))
    for i in range(len(cluster_data)):
        for j in range(i + 1, len(cluster_data)):
            diff = cluster_data.layers['velocity'][i] \
                   - cluster_data.layers['velocity'][j]
            diff = numpy.linalg.norm(diff)
            # diff = scipy.sparse.linalg.norm(diff)
            dist[i, j] = diff
            dist[j, i] = diff

    '''
    calculate neighbor MSE sum
    '''
    num_neighbors = 10
    neighbor_errors = []
    for i in range(N):
        argsorted = np.argsort(dist[i, :])
        total_error = np.sum(errors[argsorted[:num_neighbors]])
        neighbor_errors.append(total_error)

    cluster_data.obs['neighbor_errors'] = neighbor_errors
    suptitle = 'cluster:' + str(cluster_label)
    scv.pl.scatter(cluster_data,
                   color=['neighbor_errors', 'vel_norms'],
                   save='cluster%s_centroids.png' % str(cluster_label))

    # sometimes fail because of lack of samples
    try:
        scv.pl.velocity_embedding_stream(
            cluster_data,
            color=[
                'neighbor_errors',
                'vel_norms'],
            save='cluster%s_vel_centroids.png' %
            str(cluster_label),
            figsize=(
                14,
                10))
    except Exception as e:
        print('failed to generate velocity embedding stream')
        print(e)

    jacob = model.coef_
    # selected_genes = ['FN1', 'SNAI2', 'ZEB2', 'TWIST1']
    # if PCA space, we need to transform coefficient to true Jacobian.
    if pca_model:
        Q = pca_model.components_  # n components x n features
        jacob = Q.T @ jacob @ Q
    analyze_jacobian(cluster_data, jacob, selected_genes, emt_genes)


def analyze_jacobian(adata, jacob, selected_genes, emt_genes=None, topk=5):
    genes = adata.var_names
    # adata.var_names.get_loc('FN1')
    emt_genes = set(emt_genes)
    for gene in selected_genes:
        if not (gene in genes):
            print(gene, 'is not in top gene list')
            continue
        idx = adata.var_names.get_loc(gene)
        row_coef = jacob[idx, :]
        args = np.argsort(-np.abs(row_coef))
        top_genes = genes[args[:topk]]
        print('for gene:', gene)
        print('top inhib/exhibit genes:', top_genes)
        print('top inhib/exhibit genes coefs:', row_coef[args[:topk]])

        is_in_emt = [1 if gene in emt_genes else 0 for gene in top_genes ]
        print('Whether top genes in gene list (1-yes, 0-no)', is_in_emt)
        print('number of top5 genes in known emt list:', sum(is_in_emt))
        

def filter_a549_MET_samples(adata, meta, day0_only=config.day0_only):
    cell_ids = np.array(meta["Unnamed: 0"]) + 'x'
    for ID in range(len(cell_ids)):
        # This is needed to make the cell ids have the same syntax as the loom
        # files
        cell_ids[ID] = re.sub('x', "x-", cell_ids[ID], count=1)
        cell_ids[ID] = re.sub('_', ":", cell_ids[ID])

    meta['Unnamed: 0'] = cell_ids
    cells = meta['Unnamed: 0'].to_numpy()
    # time_raw = [meta['Time'][cells==cell][cell] for cell in adata.obs_names]
    time_raw = np.array([[meta['Time'][np.squeeze(
        np.argwhere(cells == cell))]][0] for cell in adata.obs_names])
    adata.obs['Time'] = time_raw

    time_raw = adata.obs['Time']
    # no _rm in time means EMT
    if day0_only:
        is_day0 = np.array([time == '0d' for time in time_raw])
        adata = adata[is_day0]
    else:
        is_EMT = np.array([time.find('_rm') == -1 for time in time_raw])
        adata = adata[is_EMT]
    return adata
