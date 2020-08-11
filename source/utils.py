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

    
