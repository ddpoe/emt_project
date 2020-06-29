from utils import *

def run_graphlasso(X, lm=0.001, prefix=''):
    model = sklearn.covariance.GraphicalLasso(alpha=lm)
    model.fit(X)
    print('mean:', model.location_)
    print('covariance:\n', model.covariance_)
    print('precision:\n', model.precision_)
    print('iters:', model.n_iter_)
    np.save(prefix + "_precision_matrix.csv", model.precision_)
    pass



def analyze_specific_cluster(adata, indices,
                             predicted_velocities, errors,
                             model, cluster_label):
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
        for j in range(i+1, len(cluster_data)):
            diff = cluster_data.layers['velocity'][i] \
                   - cluster_data.layers['velocity'][j]
            diff = numpy.linalg.norm(diff)
            # diff = scipy.sparse.linalg.norm(diff)
            dist[i, j] = diff
            dist[j, i] = diff


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
        scv.pl.velocity_embedding_stream(cluster_data,
                                         color=['neighbor_errors', 'vel_norms'],
                                         save='cluster%s_vel_centroids.png' % str(cluster_label))
    except Exception as e:
        print('failed to generate velocity embedding stream')
        print(e)
        
    jacob = model.coef_
    selected_genes = ['FN1']

    # todo : if PCA space, we need to transform coefficient to true Jacobian.
    analyze_jacobian(cluster_data, jacob, selected_genes)

    
def analyze_jacobian(adata, jacob, selected_genes, topk=5):
    genes = adata.var_names
    # adata.var_names.get_loc('FN1')

    for gene in selected_genes:
        idx = adata.var_names.get_loc(gene)
        row_coef = jacob[idx, :]
        args = np.argsort(-np.abs(row_coef))
        print('gene:', gene)
        print('top inhib/exhibit genes:', genes[args[:topk]])
        print('top inhib/exhibit genes coefs:', row_coef[args[:topk]])
        
    
def main():
    loom_data_path = '../data/a549_tgfb1.loom'
    meta_path = '../data/a549_tgfb1_meta.csv'
    # adata=dyn.read_loom(loom_data_path)
    adata=scv.read_loom(loom_data_path)
    meta=pd.read_csv(meta_path)
    # adata = scv.datasets.pancreas()
    emt_gene_path = '../data/gene_lists/emt_genes_weikang.txt'
    emt_genes = read_list(emt_gene_path)

    
    # filter by emt genes
    print('filtering genes by only using known emt genes')
    intersection_genes = set(adata.var_names).intersection(emt_genes)
    adata = adata[:, list(intersection_genes)]
    
    n_top_genes = 50 # not 2000 because of # observations
    print('data read complete', flush=True)

    CellIDs=np.array(meta["Unnamed: 0"])+'x'
    for ID in range(len(CellIDs)):
        #This is needed to make the cell ids have the same syntax as the loom files 
        CellIDs[ID]=re.sub('x',"x-",CellIDs[ID],count=1)
        CellIDs[ID]=re.sub('_',":",CellIDs[ID])
    
    meta['Unnamed: 0']=CellIDs
    cells = meta['Unnamed: 0'].to_numpy()
    time_raw = np.array([[meta['Time'][np.squeeze(np.argwhere(cells==cell))]][0] for cell in adata.obs_names])
    adata.obs['Time'] = time_raw

    time_raw = adata.obs['Time']
    # no _rm in time means EMT
    is_EMT = np.array([time.find('_rm')==-1 for time in time_raw])
    adata = adata[is_EMT]
    treatment = meta['Treatment']
    adata.obs['treatment']=treatment
    # print('flag3')
    
    count_matrix = adata.X.todense()[:, ...]
    print('count matrix nonzerorate:', np.count_nonzero(count_matrix) / np.prod(count_matrix.shape))

    scv.pp.filter_and_normalize(adata, n_top_genes=n_top_genes)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    velocities = adata.layers['velocity']
    print('velocities matrix nonzero rate:', np.count_nonzero(velocities) / np.prod(velocities.shape))
    scv.tl.terminal_states(adata)
    root_cells = adata.obs['root_cells']
    end_cells = adata.obs['end_points']
    print('root cell shape:', root_cells.shape)

    # gen figures
    scv.pl.velocity_embedding_stream(adata, save='vel_stream.png')


    # gen ranked genes
    scv.tl.rank_velocity_genes(adata, groupby='Clusters')
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.head().to_csv('rank_genes_vf.csv')

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

        # adata.obs['Clusters'] = kmeans_labels
        print('computing  MAR')
        
        scv_labels = adata.obs['Clusters']

        labels = scv_labels
        label_set = set(labels)
        errors = np.zeros(len(labels))
        # model=LinearRegression().fit(pca_count_matrix, velocities)
        whole_data_label = -1
        label_set.add(whole_data_label) # -1 denote for whole dataset
        
        for label in label_set:
            print('label:', label)
            if label == whole_data_label:
                indices = np.full(len(labels), True)
            else:
                indices = labels == label

            print('#samples in this cluster:', np.sum(indices))

            # choose: PCA reduced by sklearn or reduced by packages?
            # label_count_matrix = pca_count_matrix[indices, ...]
            label_count_matrix = adata.X.todense()[indices, ...]
    
            label_velocities = velocities[indices, ...]
            
            model=LinearRegression().fit(label_count_matrix, label_velocities)
            predicted_velocities = model.predict(label_count_matrix)
            diff = predicted_velocities - label_velocities
            diff = np.sum(diff**2, axis=-1)
            errors[indices] = diff            
            analyze_specific_cluster(adata, indices, predicted_velocities, diff, model, label)

            r2_score = model.score(label_count_matrix, label_velocities)
            # explained_variance = sklearn.metrics\
            #                             .explained_variance_score(label_velocities,
            #                                                       predicted_velocities)
            mse = sklearn.metrics\
                         .mean_squared_error(label_velocities,
                                             predicted_velocities)
            print('label:%d, r^2 score: %.5f, mse:%.5f'\
                  % (label, r2_score, mse))
            # scv.pl.scatter(adata, color=[ 'root_cells', 'end_points', 'errors', 'kmeans_labels'], save='error_root_end_points.png')
            if label == whole_data_label:
                adata.obs['whole_data_squared_error'] = errors
                pass
            else:
                adata.obs['cluster_squared_error'] = errors
                
        scv.pl.scatter(adata, color=[ 'root_cells', 'end_points', 'cluster_squared_error', 'whole_data_squared_error', 'kmeans_labels', 'Clusters'], save='error_root_end_points.png')
    


    def main_graphlasso():
        '''
        graph lasso
        '''
        # count_matrix=adata.X.A
        print('count matrix dimension:', count_matrix.shape)
        print('velocities dimension:', velocities.shape)
        
        print('computing graph lassos')
        # run_graphlasso(count_matrix[:, :1000], prefix='count_matrix')
        print('velocities shape:', velocities.shape)
        run_graphlasso(velocities * 100, prefix='velocity')

    main_MAR()
    # main_graphlasso()
    
if __name__ == '__main__':
    main()
