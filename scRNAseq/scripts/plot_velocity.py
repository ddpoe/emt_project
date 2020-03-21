import sys
import dynamo as dyn
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns

def plot_X(X, dim1=0, dim2=1, create_figure=False, figsize=(6, 6), **kwargs):
    if create_figure:
        plt.figure(figsize=figsize)
    plt.scatter(X[:, dim1], X[:, dim2], **kwargs)

def plot_V(X, V, dim1=0, dim2=1, create_figure=False, figsize=(6, 6), **kwargs):
    if create_figure:
        plt.figure(figsize=figsize)
    plt.quiver(X[:, dim1], X[:, dim2], V[:, dim1], V[:, dim2])

def color_set(data,gene):
    gene_idx=np.where(data.var_names==gene)[0][0]
    c=np.squeeze(np.asarray(data.X[:,gene_idx].todense()))
    return c

dataset=sys.argv[1]

embs=['X_pca','X_umap']
#3 is the stationary distribution, it will used later
plots=['FN1','VIM',3]

P=np.loadtxt(dataset+'_P.npy')
Idx=np.load(dataset+'_idx.npy',allow_pickle=True)
adata=dyn.read_h5ad(dataset+'_processed.h5ad')

kmc=dyn.tl.KernelMarkovChain(P=P,Idx=Idx)

X = adata.layers['M_s'][:, adata.var['use_for_velocity']]
V = adata.layers['velocity_S'][:, adata.var['use_for_velocity']]
sd = kmc.compute_stationary_distribution()

for plot in plots:
    for emb in embs:
        X_emb=adata.obsm[emb]
        Uc = kmc.compute_density_corrected_drift(X_emb, normalize_vector=True)
        U_grid, X_grid = dyn.tl.smoothen_drift_on_grid(X_emb[:, :2], Uc[:, :2], 30, k=50, smoothness=0.5)
        if type(plot)==str:
            c=color_set(adata,plot)
            title=dataset+': '+plot+" | "+emb.split("_")[1]
        else: 
            c=sd
            title=dataset+': stationary distribution | '+emb.split("_")[1]

        plot_X(X_emb, create_figure=True, figsize=(12, 6),c=c)
        plt.colorbar()
        plot_V(X_grid, U_grid, facecolor='k')
        plt.title(title)
        plt.savefig(title+'.png') 
