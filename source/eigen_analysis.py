from utils import *

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
        neighbor_errors.append(total_error/num_neighbors)

    cluster_data.obs['whole_cluster_avg_neighbor_errors'] = neighbor_errors
    suptitle = 'cluster:' + str(cluster_label)
    scv.pl.scatter(cluster_data,
                   color=['whole_cluster_avg_neighbor_errors', 'vel_norms'],
                   save='cluster%s_centroids.png' % str(cluster_label),
                   show=False)

    # sometimes fail because of lack of samples
    try:
        scv.pl.velocity_embedding_stream(
            cluster_data,
            color=[
                'whole_cluster_avg_neighbor_errors',
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
    return np.array(eig_vals), np.array(eig_vectors)


def analyze_jacob_eigen_complex_plane(jacob):
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
    else:
        is_EMT = np.array([time.find('_rm') == -1 for time in time_raw])
        adata = adata[is_EMT]
    return adata

    

def neighbor_MAR(data_mat, labels, neighbor_num=100, dist_mat=None, pca_model=None):
    '''
    labels: RNA velocities
    '''
    if dist_mat is None:
        dist_mat = calc_distance_matrix(data_mat)

    num_feature = neighbor_num // 10
    if pca_model:
        ratio = sum(pca_model.explained_variance_ratio_[:num_feature])
        print('neighborhood pca explained variance:', ratio, '#feature:', num_feature)
    data_mat = data_mat[:, :num_feature]
    sub_labels = labels[:, :num_feature]
    # sub_labels = labels
    mses, r2s, max_eigenvals, feature_norms, bias_norms, jacobs, models = [], [], [], [], [], [], []
    print('dist shape:',  dist_mat.shape)
    for i in range(len(data_mat)):
        neighbors = np.argsort(dist_mat[i, :])[:neighbor_num]        
        specific_mat = data_mat[neighbors, :]
        neighbor_labels = sub_labels[neighbors, :]
        model = LinearRegression().fit(specific_mat, neighbor_labels)
        # model = Lasso(alpha=config.lasso_alpha).fit(specific_mat, neighbor_labels)
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
        bias_norm = np.linalg.norm(model.intercept_)
        bias_norms.append(bias_norm)
        jacobs.append(model.coef_)
        models.append(model)
    return mses, r2s, max_eigenvals, feature_norms, bias_norms, jacobs, models


def centroid_neighbor_MAR(data_mat, labels, neighbor_num, dist_mat=None, pca_model=None, center=None):
    if dist_mat is None:
        dist_mat = calc_distance_matrix(data_mat)

    num_feature = neighbor_num // 10
    if center is None:
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
        analyze_jacob_eigen_complex_plane(pca_jacob)
        # print('gene space jacob:')
        # analyze_jacob_eigen(gene_jacob)
    return sample_mses, is_center_neighbors, min_id


def analyze_group_jacobian(adata, jacobs, topk_eigen=4, topk_gene=100, pca_model=None, group_name='someGroup', top_genes_for_plot=50, gene_count_df=None):
    '''
    analyze one group's jacob and gene
    '''
    if pca_model:
        Q = pca_model.components_ # n components x n features

    topk_gene_lists = []
    adata.uns['top_eig_genes'] = [[] for _ in range(len(adata))]
    for i, jacob in enumerate(jacobs):
        eig_vals, eig_vectors = calc_eigen(jacob)
        reals, imgs = [num.real for num in eig_vals], [num.imag for num in eig_vals]
        args = np.argsort(-np.abs(reals))
        
        topk_eigenvals = np.array(eig_vals[args[:topk_eigen]])
        # topk_eigenvectors = eig_vectors[args[:topk_eigen], ...]
        if pca_model:
            topk_eigenvectors = (Q.T[:, :config.MAR_neighbor_num//10]) @ np.array(eig_vectors) # columns are eigenvectors by math convention
            topk_eigenvectors = topk_eigenvectors.T
            # topk_eigenvectors = topk_eigenvectors[..., :topk_eigen]
        # print('gene eigen vector shape:', topk_eigenvectors.shape)
        significant_genes = set()
        for j in range(len(topk_eigenvectors)):
            topk_genes_args = np.argsort(-np.abs(topk_eigenvectors[j]))[:topk_gene] # based on |c| if complex
            significant_genes.update(topk_genes_args)
        # print('significant genes:', significant_genes)
        adata.uns['top_eig_genes'][i] = list(significant_genes)

    top_eig_genes = adata.uns['top_eig_genes']
    gene_names = adata.var_names
    count_map = {name:0 for name in gene_names}
    for i in range(len(top_eig_genes)):
        for j in range(len(top_eig_genes[i])):
            name = gene_names[top_eig_genes[i][j]]
            count_map[name] += 1
        items = list(count_map.items())
    # choose topk genes to draw
    items = np.array(sorted(items, key=lambda x: x[1], reverse=True))
    # print('top eigen genes:', items[:top_genes_for_plot])

    if gene_count_df is None:
        df = pd.DataFrame(data=items[:, 1], index=items[:, 0])
        df.to_csv(os.path.join('./figures', group_name + '_geneCountInMaxEigenAnalysis.csv'))
    else:
        gene_count_df[group_name] = pd.Series(data=items[:, 1], index=items[:, 0])
        gene_count_df.to_csv(os.path.join('./', group_name + 'geneCountInMaxEigenAnalysis.csv'))


def analyze_MAR_biases(adata, models):
    X = adata.X
    biases_v = []
    biases_e = []
    for sample_i in range(len(adata)):
        x = X[sample_i]
        bias = models[sample_i].intercept_
        A = models[sample_i].coef_
        biases_v.append(bias)
        bias_e, residuals, rank, s = numpy.linalg.lstsq(A, bias)
        biases_e.append(bias_e)

    # note only 2 dim of bias_v will be drawn: need to transform
    biases_v_emb = PCA(n_components=2,
                       random_state=config.random_state).fit_transform(biases_v)
    biases_v_emb_tsne = TSNE(n_components=2,
                       random_state=config.random_state).fit_transform(biases_v)    
    scv.pl.velocity_embedding_stream(adata,
                                     V = biases_v_emb,
                                     save='bias_v_pca_stream.png', show=False)
    scv.pl.velocity_embedding_stream(adata,
                                     V = biases_v_emb_tsne,
                                     save='bias_v_tsne_stream.png', show=False)
    biases_e_emb = PCA(n_components=2,
                       random_state=config.random_state).fit_transform(biases_e)
    scv.pl.velocity_embedding_stream(adata,
                                     X = biases_e_emb,
                                     V = biases_v_emb,
                                     save='bias_e_stream.png', show=False)
