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
import scipy.spatial.distance

def calc_distance_matrix(data):
    dist_mat = scipy.spatial.distance.pdist(data, metric='euclidean')
    dist_mat = scipy.spatial.distance.squareform(dist_mat)
    return dist_mat

        
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
    print(sub_data.var_names)
    scv.pp.filter_and_normalize(sub_data)
    scv.pp.moments(sub_data)
    scv.tl.velocity(sub_data, var_names=[g1, g2])
    scv.tl.velocity_graph(sub_data)
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
    # dist = np.zeros((N, N))
    # for i in range(len(cluster_data)):
    #     for j in range(i + 1, len(cluster_data)):
    #         diff = cluster_data.layers['velocity'][i] \
    #                - cluster_data.layers['velocity'][j]
    #         diff = numpy.linalg.norm(diff)
    #         # diff = scipy.sparse.linalg.norm(diff)
    #         dist[i, j] = diff
    #         dist[j, i] = diff
    dist = calc_distance_matrix(cluster_data.layers['velocity'])
    
    '''
    calculate neighbor MSE sum
    '''
    num_neighbors = config.MAR_neighbor_num
    neighbor_errors = []
    for i in range(N):
        argsorted = np.argsort(dist[i, :])
        total_error = np.sum(errors[argsorted[:num_neighbors]])
        neighbor_errors.append(total_error)

    cluster_data.obs['whole_cluster_neighbor_errors'] = neighbor_errors
    suptitle = 'cluster:' + str(cluster_label)
    scv.pl.scatter(cluster_data,
                   color=['whole_cluster_neighbor_errors', 'vel_norms'],
                   save='cluster%s_centroids.png' % str(cluster_label),
                   show=False)

    # sometimes fail because of lack of samples
    try:
        scv.pl.velocity_embedding_stream(
            cluster_data,
            color=[
                'whole_cluster_neighbor_errors',
                'vel_norms'],
            save='cluster%s_vel_centroids.png' %
            str(cluster_label),
            figsize=(
                14,
                10),
            show=False)
        pass
    except Exception as e:
        print('failed to generate velocity embedding stream')
        print(e)

    jacob = model.coef_
    # selected_genes = ['FN1', 'SNAI2', 'ZEB2', 'TWIST1']
    # if PCA space, we need to transform coefficient to true Jacobian.
    if pca_model:
        Q = pca_model.components_  # n components x n features
        jacob = Q.T @ jacob @ Q

    analyze_adata_jacobian_genes(cluster_data, jacob, selected_genes, emt_genes)


def calc_eigen(jacob):
    eig_vals, eig_vectors = np.linalg.eig(jacob)
    return eig_vals, eig_vectors


def analyze_jacob_eigen(jacob):
    eig_vals, eig_vectors = calc_eigen(jacob)
    print('jacobian shape:', jacob.shape)
    
    reals, imgs = [num.real for num in eig_vals], [num.imag for num in eig_vals]
    plt.scatter(reals, imgs)
    plt.xlabel('real')
    plt.ylabel('image')
    # plt.show()
    plt.savefig('./figures/jacob_eigen_complex_plane.png')
    sorted_reals = sorted(reals, reverse=True)
    print(sorted_reals)
    sorted_imgs = sorted(imgs, reverse=True)
    print(sorted_imgs)
    
    pass


def analyze_adata_jacobian_genes(adata, jacob, selected_genes, emt_genes=None, topk=5):
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

        is_in_emt = [1 if gene in emt_genes else 0 for gene in top_genes ]

        print('########################################')
        print('for gene:', gene)
        print('top inhib/exhibit genes:', top_genes)
        print('top inhib/exhibit genes coefs:', row_coef[args[:topk]])        
        print('Whether top genes in emt gene list (1-yes, 0-no)', is_in_emt)
        print('number of top5 genes in known emt list:', sum(is_in_emt))
        print('########################################')

        
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

    

def neighbor_MAR(data_mat, labels, neighbor_num=100, dist_mat=None, pca_model=None):
    if dist_mat is None:
        dist_mat = calc_distance_matrix(data_mat)

    num_feature = neighbor_num // 10
    if pca_model:
        ratio = sum(pca_model.explained_variance_ratio_[:num_feature])
        print('neighborhood pca explained variance:', ratio, '#feature:', num_feature)
    data_mat = data_mat[:, :num_feature]
    sub_labels = labels[:, :num_feature]
    # sub_labels = labels
    mses, r2s, max_eigenvals, feature_norms = [], [], [], []
    print('dist shape:',  dist_mat.shape)
    for i in range(len(data_mat)):
        neighbors = np.argsort(dist_mat[i, :])[:neighbor_num]        
        specific_mat = data_mat[neighbors, :]
        neighbor_labels = sub_labels[neighbors, :]
        model = LinearRegression().fit(specific_mat, neighbor_labels)
        predicted_vals = model.predict(specific_mat)

        r2_score = model.score(specific_mat, neighbor_labels)
        mse = sklearn.metrics\
                     .mean_squared_error(neighbor_labels,
                                         predicted_vals)
        # print('center:%d, r^2 score: %.5f, mse:%.5f' % (i, r2_score, mse))
        # print('model.coef_ shape:', model.coef_.shape)
        eigen_vals, eigen_vectors = calc_eigen(model.coef_)
        eigen_reals = [num.real for num in eigen_vals]
        feature_norm = np.mean(np.linalg.norm(neighbor_labels, axis=1))
        max_eigenvals.append(max(eigen_reals))
        mses.append(mse)
        r2s.append(r2_score)
        feature_norms.append(feature_norm)
    return mses, r2s, max_eigenvals, feature_norms


def centroid_neighbor_MAR(data_mat, labels, neighbor_num, dist_mat=None, pca_model=None):
    if dist_mat is None:
        dist_mat = calc_distance_matrix(data_mat)

    num_feature = neighbor_num // 10
    center = np.mean(data_mat, axis=0)
    dist_vec = np.linalg.norm(data_mat - center, axis=1)
    # print('dist_vec shape:', dist_vec.shape)
    neighbors = np.argsort(dist_vec)[:neighbor_num]
    is_center_neighbors = np.full(len(data_mat), False)
    is_center_neighbors[neighbors] = True
    # print('neighbor len:', neighbors)
    
    specific_mat = data_mat[neighbors, :num_feature]
    neighbor_labels = labels[neighbors, :num_feature]
    # neighbor_labels = labels[neighbors, :] # predict all velocities
    model = LinearRegression().fit(specific_mat, neighbor_labels)
    predicted_vals = model.predict(specific_mat)
    sample_mses = np.sum((neighbor_labels - predicted_vals)**2, axis=1)

    min_id = np.argsort(dist_vec)[0]

    print("analyze mean centroid jacobian")
    jacob = model.coef_
    if not (pca_model is None):
        pca_jacob = jacob
        Q = pca_model.components_[:num_feature, :]  # n principle components x n original features
        gene_jacob = Q.T @ jacob @ Q

        print('PCA space jacob:')
        analyze_jacob_eigen(pca_jacob)
        # print('gene space jacob:')
        # analyze_jacob_eigen(gene_jacob)
    return sample_mses, is_center_neighbors, min_id
