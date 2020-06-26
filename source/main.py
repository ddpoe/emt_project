from utils import *

def run_graphlasso(X, lm=0.00001, prefix=''):
    model = sklearn.covariance.GraphicalLasso(alpha=lm)
    model.fit(X)
    print('mean:', model.location_)
    print('covariance:\n', model.covariance_)
    print('precision:\n', model.precision_)
    print('iters:', model.n_iter_)
    np.save(prefix + "_precision_matrix.csv", model.precision_)
    pass


def main():
    loom_data_path = '../data/a549_tgfb1.loom'
    meta_path = '../data/a549_tgfb1_meta.csv'
    adata=dyn.read_loom(loom_data_path)
    meta=pd.read_csv(meta_path)
    n_top_genes = 2000
    print('data read complete', flush=True)
    cells=meta['Unnamed: 0'].to_numpy()
    
    # print('flag1')
    # treatment=np.array([[meta['Time'][np.squeeze(np.argwhere(cells==cell))]][0] for cell in adata.obs_names])
    # print('flag2')
    # adata.obs['treatment']=treatment
    # print('flag3')
    count_matrix=adata.X.todense()[:, ...]
    
    scv.pp.filter_and_normalize(adata, n_top_genes=n_top_genes)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    velocities=adata.layers['velocity']

    scv.tl.terminal_states(adata)
    root_cells = adata.obs['root_cells']
    end_cells = adata.obs['end_points']
    print('root cell shape:', root_cells.shape)
    print(root_cells[:5])

    # gen figures
    scv.pl.velocity_embedding_stream(adata, save='vel_stream.png')
    scv.pl.scatter(adata, color=[ 'root_cells', 'end_points'], save='root_end_points.png')
    
    '''
    MAR
    '''
    pca=PCA(n_components=n_top_genes,random_state=0).fit(count_matrix)
    pca_count_matrix=pca.transform(count_matrix)
    # knee = get_optimal_K(pca_count_matrix, kmin=1, kmax=21)
    knee = 7 # calculated
    print('knee of kmeans graph:', knee)
    kmeans = KMeans(n_clusters=knee ,random_state=0).fit(pca_count_matrix)
    cluster_centers=kmeans.cluster_centers_
    labels=kmeans.labels_
    print('computing  MAR')
    label_set = set(labels)
    errors = np.zeros(len(labels))
    adata.obs['kmeans_labels'] = labels
    for label in label_set:
        print('label:', label)
        indices = labels == label
        label_count_matrix = pca_count_matrix[indices, ...]
        label_velocities = velocities[indices, ...]
        model=LinearRegression().fit(label_count_matrix, label_velocities)
        predicted_velocities = model.predict(label_count_matrix)
        
        diff = predicted_velocities - label_velocities
        diff = np.sum(np.abs(diff)**2,axis=-1)**.5
        errors[indices] = diff
        
    adata.obs['mar_mse'] = errors
    # scv.pl.scatter(adata, color=[ 'root_cells', 'end_points', 'errors', 'kmeans_labels'], save='error_root_end_points.png')
    scv.pl.scatter(adata, color=[ 'root_cells', 'end_points', 'mar_mse', 'kmeans_labels'], save='error_root_end_points.png')
    

    
    '''
    graph lasso
    '''
    print('count matrix dimension:', count_matrix.shape)
    print('velocities dimension:', velocities.shape)
    count_matrix=adata.X.A
    print('computing graph lassos')
    run_graphlasso(count_matrix, prefix='count_matrix')
    run_graphlasso(velocities, prefix='velocity')

if __name__ == '__main__':
    main()
