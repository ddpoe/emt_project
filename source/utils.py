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
from sklearn.manifold import TSNE
from sklearn import metrics
import pandas as pd
import scvelo as scv
from kneed import KneeLocator
import matplotlib.pyplot as plt
from scipy import interpolate
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression, Lasso
import sklearn
import sklearn.covariance
import scipy
import scipy.sparse.linalg
import numpy.linalg
import sklearn.metrics
import sklearn.model_selection
import config
import os
import sys
import scipy.spatial.distance
import matplotlib
matplotlib.use('Agg')

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

    
def filter_a549_MET_samples(adata, meta, include_a549_days=config.include_a549_days):
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
    if include_a549_days != 'all':
        is_day0 = np.array([time in include_a549_days for time in time_raw])
        adata = adata[is_day0]
        if len(include_a549_days) == 1:
            adata.obs['Time'] = 'r' # prevent SCVELO plot bug
    else:
        is_EMT = np.array([time.find('_rm') == -1 for time in time_raw])
        adata = adata[is_EMT]
    return adata


def read_CBC_GBC_mapping(path):
    def get_mapping(frame):
        res = {}
        for i in range(len(frame)):
            cbc = frame['CBC'][i]
            gbc = frame['GBC'][i]
            res[cbc] = gbc
        return res
    frame = pd.read_csv(path, delimiter='\t')
    return get_mapping(frame)


def read_gbc2dosage_mapping(path):
    def get_mapping(frame):
        res = {}
        for i in range(len(frame)):
            gbc = frame['Name'][i]
            dosage = float(frame['TGFb'][i])
            res[gbc] = dosage
        return res
    frame = pd.read_csv(path, delimiter='\t')
    return get_mapping(frame)


def dosage_count(adata):
    all_dosage_count = {}
    for i in range(len(adata)):
        dosage = adata.obs['dosage'][i]
        if not (dosage in all_dosage_count):
            all_dosage_count[dosage] = 0
        all_dosage_count[dosage] += 1
    for dosage in all_dosage_count:
        print('dosage:', dosage, '#samples:', all_dosage_count[dosage])

        
def process_kazu_loom_data(adata, cbc2gbc_path, gbc_info_path, dosage_range):
    '''
    format kazu loom data: 
    1) modify observations name 
    2) add dosage information 
    3) discard any samples without meta dosage info
    '''
    adata_cbc_codes = [x[x.find(':')+1:] for x in list(adata.obs_names)]
    adata.obs_names = adata_cbc_codes
    adata.obs['dosage'] = np.zeros(len(adata), dtype=np.float)
    cbc2gbc = read_CBC_GBC_mapping(cbc2gbc_path)
    gbc2dosage = read_gbc2dosage_mapping(gbc_info_path)
    is_cbc_in_mapping = []
    no_meta_info_count = 0
    for cbc in adata_cbc_codes:
        if (not (cbc in cbc2gbc)) or (not (cbc2gbc[cbc] in gbc2dosage)):
            is_cbc_in_mapping.append(False)
            no_meta_info_count += 1
            continue
        else:
            is_cbc_in_mapping.append(True)
            gbc = cbc2gbc[cbc]
            adata.obs['dosage'][cbc] = gbc2dosage[gbc]
    adata = adata[is_cbc_in_mapping]
    print('total observation discarded due to no gbc/cbc info:', no_meta_info_count)
    print('new #observations after processing kazu loom data:', len(adata))
    # dosage statistics
    dosage_count(adata)
    # filter by dosage
    adata = adata[ np.logical_and((dosage_range[0] <= adata.obs['dosage']), (adata.obs['dosage'] <= dosage_range[1]))]
    # dosage stats after filtering
    print('all dosage counting after filtering by dosage range passed in config parameter:')
    dosage_count(adata)
    return adata


def gen_analysis_figures(adata, prefix):
    whole_neighbor_obs = ['whole_neighbor_MAR_r2',
                          'whole_neighbor_MAR_mse',
                          'whole_neighbor_bias_norms',
                          'whole_neighbor_avg_vel_norms',
                          'whole_neighbor_max_eigenVal_real',
                          'whole_neighbor_MAR_is_eigenstable',
                          'whole_data_MAR_squared_error']
    if config.use_dataset == 'kazu_mcf10a':
        whole_neighbor_obs += ['dosage']
    # debug
    # print(whole_neighbor_obs)
    # print(adata.obs_keys())
    if not config.only_whole_data:
        scv.pl.scatter(
            adata,
            color=[
                'root_cells',
                'end_points',
                'cluster_squared_error',
                'whole_data_MAR_squared_error',
                'Clusters'],
            save=prefix + 'error_root_end_points.png',
            show=False)
        scv.pl.velocity_embedding_stream(
            adata,
            color=whole_neighbor_obs +  ['clusters_neighbor_MAR_r2',
                                         'clusters_neighbor_MAR_mse',
                                         'clusters_centroid_neighborhood_sample_mses',
                                         'is_centroid_neighbor_indicator',
                                         'Clusters',
                                         'vel_norms'],
            colorbar=True,
            save=prefix + '_neighbor_MAR_stats.png',
            show=False)

    else:
        scv.pl.scatter(
            adata,
            color=[
                'root_cells',
                'end_points',
                'whole_data_MAR_squared_error',
                'Clusters'],
            save=prefix + '_error_root_end_points.png',
            show=False)
        scv.pl.velocity_embedding_stream(adata,
                                         color=whole_neighbor_obs + ['Time',
                                                                     'Clusters',
                                                                     'vel_norms'],
                                         save=prefix + '_neighbor_MAR_stats.png',
                                         show=False)    
