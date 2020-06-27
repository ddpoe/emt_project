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
    # adata = scv.datasets.pancreas()
    
    n_top_genes = 2000
    print('data read complete', flush=True)
    # cells = meta['Unnamed: 0'].to_numpy()
    # print('flag1')
    # treatment=np.array([[meta['Time'][np.squeeze(np.argwhere(cells==cell))]][0] for cell in adata.obs_names])
    # print('flag2')
    # adata.obs['treatment']=treatment
    # print('flag3')
    count_matrix = adata.X.todense()[:, ...]
    print('count matrix nonzerorate:', np.count_nonzero(count_matrix) / np.prod(count_matrix.shape))

    scv.pp.filter_and_normalize(adata, n_top_genes=n_top_genes)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    velocities=adata.layers['velocity']
    print('velocities matrix nonzero rate:', np.count_nonzero(velocities) / np.prod(velocities.shape))
    scv.tl.terminal_states(adata)
    root_cells = adata.obs['root_cells']
    end_cells = adata.obs['end_points']
    print('root cell shape:', root_cells.shape)
    print(root_cells[:5])

    # gen figures
    scv.pl.velocity_embedding_stream(adata, save='vel_stream.png')


    def main_MAR():
        '''
        MAR
        '''
        pca=PCA(n_components=n_top_genes,random_state=7).fit(count_matrix)
        pca_count_matrix=pca.transform(count_matrix)
        # knee = get_optimal_K(pca_count_matrix, kmin=1, kmax=21)
        knee = 7 # calculated
        print('knee of kmeans graph:', knee)
        kmeans = KMeans(n_clusters=knee ,random_state=0).fit(pca_count_matrix)
        cluster_centers=kmeans.cluster_centers_
        kmeans_labels=kmeans.labels_
        adata.obs['kmeans_labels'] = kmeans_labels
        
        print('computing  MAR')

        scv_labels = adata.obs['Clusters']

        labels = scv_labels
        label_set = set(labels)
        errors = np.zeros(len(labels))
        # model=LinearRegression().fit(pca_count_matrix, velocities)
        for label in label_set:
            print('label:', label)
            indices = labels == label
            label_count_matrix = pca_count_matrix[indices, ...]
            label_velocities = velocities[indices, ...]
            model=LinearRegression().fit(label_count_matrix, label_velocities)
            predicted_velocities = model.predict(label_count_matrix)

            diff = predicted_velocities - label_velocities
            diff = np.sum(diff**2, axis=-1)
            errors[indices] = diff
            
            adata.obs['sample_squared_error'] = errors

            score = model.score(label_count_matrix, label_velocities)
            print('label:%d, r^2 score: %f.5' % (label, score))
            # scv.pl.scatter(adata, color=[ 'root_cells', 'end_points', 'errors', 'kmeans_labels'], save='error_root_end_points.png')
        scv.pl.scatter(adata, color=[ 'root_cells', 'end_points', 'sample_squared_error', 'kmeans_labels', 'Clusters'], save='error_root_end_points.png')
    


    def main_graphlasso():
        '''
        graph lasso
        '''
        # count_matrix=adata.X.A
        print('count matrix dimension:', count_matrix.shape)
        print('velocities dimension:', velocities.shape)
        
        print('computing graph lassos')
        # run_graphlasso(count_matrix[:, :1000], prefix='count_matrix')
        run_graphlasso(velocities[:, :100], prefix='velocity')

    main_MAR()
    # main_graphlasso()
    
if __name__ == '__main__':
    main()
