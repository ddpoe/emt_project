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

def get_optimal_K(X, kmin=1, kmax=21):
    distortions = [] 
    inertias = [] 
    mapping1 = {} 
    k2inertia = {} 
    K = range(kmin, kmax)   
    for k in K: 
    #Building and fitting the model 
        kmeanModel = KMeans(n_clusters=k,random_state=0).fit(X) 
        kmeanModel.fit(X)     
        
        distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_, 
                                            'euclidean'),axis=1)) / X.shape[0]) 
        inertias.append(kmeanModel.inertia_) 
    
        mapping1[k] = sum(np.min(cdist(X, kmeanModel.cluster_centers_, 
                                       'euclidean'),axis=1)) / X.shape[0] 
        k2inertia[k] = kmeanModel.inertia_ 
    kn = KneeLocator(K, distortions, curve='convex', direction='decreasing')
    return kn

