from utils import *
from eigen_analysis import *
import config

def run_graphlasso(X, lm=0.001, prefix=''):
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
    
    adata = None
    if not config.use_dataset == 'pancreas':
        adata = scv.read_loom(loom_data_path) # adata=dyn.read_loom(loom_data_path)
        meta = pd.read_csv(meta_path)
        # gen_all_gene_pair_vector_field(adata, config.selected_genes_jacobian)
    else:
        adata = scv.datasets.pancreas()

    emt_gene_path = '../data/gene_lists/emt_genes_weikang.txt'
    emt_genes = read_list(emt_gene_path)

    '''
    filter by emt genes?
    '''
    if config.use_emt_gene_filter:
        print('filtering genes by only using known emt genes')
        intersection_genes = set(adata.var_names).intersection(emt_genes)
        adata = adata[:, list(intersection_genes)]

        print('intersection genes:', intersection_genes)

    # n_top_genes = 50  # not 2000 because of # observations
    n_top_genes = 2000
    print('data read complete', flush=True)
    print('using %d top dispersed genes' % n_top_genes)

    # print('flag3')
    count_matrix = adata.X.todense()[:, ...]
    print('count matrix nonzerorate:', np.count_nonzero(
        count_matrix) / np.prod(count_matrix.shape))

    scv.pp.filter_and_normalize(adata, n_top_genes=n_top_genes)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)

    scv.tl.terminal_states(adata)
    scv.pl.velocity_embedding_stream(adata, save='vel_stream.png', show=False)

    if config.use_dataset == 'a549':
        adata = filter_a549_MET_samples(adata, meta, include_a549_days=config.include_a549_days)
    else:
        adata.obs['Time'] = np.zeros(len(adata), dtype=np.int)
        
    root_cells = adata.obs['root_cells']
    end_cells = adata.obs['end_points']
    print('root cell shape:', root_cells.shape)

    # gen figures

    
    count_matrix = adata.X.todense()[:, ...]
    velocities = adata.layers['velocity']
    print(
        'velocities matrix nonzero rate:',
        np.count_nonzero(velocities) /
        np.prod(
            velocities.shape))
    
    def main_MAR(only_whole_data=False, use_pca=True):
        '''
        MAR main interface
        '''

        '''
        do whole pca once
        '''
        # if use_pca:
        #     print('applying PCA model to gene space')
        #     num_pc = 100
        #     pca_model = PCA(
        #         n_components=num_pc,
        #         random_state=7).fit(count_matrix)
        #     # pca=PCA(n_components=n_top_genes, random_state=7).fit(count_matrix)
        #     pca_count_matrix = pca_model.transform(count_matrix)

        # else:
        #     num_pc = n_top_genes
        #     pca_model = None


        '''
        determine clusters
        '''
        # knee = get_optimal_K(pca_count_matrix, kmin=1, kmax=21)
        # knee = get_optimal_K(count_matrix, kmin=1, kmax=21)
        knee = 4 # calculated
        if config.use_dataset == 'a549':
            if not config.include_a549_days:
                knee = 8
            else:
                knee = 4
        elif config.use_dataset == 'pancreas':
            knee = 7
        print('knee of kmeans graph:', knee)
        kmeans = KMeans(n_clusters=knee, random_state=0).fit(count_matrix)
        cluster_centers = kmeans.cluster_centers_
        kmeans_labels = kmeans.labels_
        adata.obs['kmeans_labels'] = kmeans_labels
        adata.obs['Clusters'] = kmeans_labels
        print('computing  MAR')
        # scv_labels is kmeans label set earlier
        scv_labels = adata.obs['Clusters']
        labels = scv_labels
        label_set = set(labels)
        errors = np.zeros(len(labels))        

        '''
        gen ranked genes
        '''
        if config.use_dataset == 'pancreas':
            # due to package issue, somehow pancreas dataset has clusters not capitalized
            scv.tl.rank_velocity_genes(adata, groupby='clusters')
        else:
            scv.tl.rank_velocity_genes(adata, groupby='Clusters')
        df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
        df.head().to_csv('./figures/rank_genes_vf.csv')


        '''
        set AnnData Obj layer/obs
        '''
        adata.obs['clusters_neighbor_MAR_mse'] = np.zeros(len(adata))
        adata.obs['clusters_neighbor_MAR_r2'] = np.zeros(len(adata))
        adata.obs['whole_neighbor_MAR_mse'] = np.zeros(len(adata))
        adata.obs['whole_neighbor_MAR_r2'] = np.zeros(len(adata))
        adata.obs['whole_neighbor_bias_norms'] = np.zeros(len(adata))
        adata.obs['whole_neighbor_max_eigenVal_real'] = np.zeros(len(adata))
        adata.obs['whole_neighbor_avg_vel_norms'] = np.zeros(len(adata))
        
        adata.obs['clusters_centroid_neighborhood_sample_mses'] = np.full(len(adata), -1)
        adata.obs['is_centroid_neighbor_indicator'] = np.zeros(len(adata), dtype=np.int)
        adata.obs['vel_norms'] = numpy.linalg.norm(adata.layers['velocity'], axis=1)
        adata.obs['whole_data_centroid_sample_errors'] = np.zeros(len(adata))
        adata.obs['is_whole_centroid_neighbor'] = np.zeros(len(adata))
        
        adata.uns['neighbor_jacobs'] = [[] for _ in range(len(adata))]
        
        '''
        analyze each cluster
        '''
        whole_data_label = -1 # denotes using all data as a cluster
        label_set.add(whole_data_label)  # -1 denote for whole dataset
        raw_labels = labels
        labels = raw_labels.values
        
        for label in label_set:
            # skip if using only whole dataset
            if only_whole_data and label != whole_data_label:
                continue
            
            print('label:', label)
            if label == whole_data_label:
                # indices = np.full(len(labels), True)
                indices = np.arange(len(labels))
            else:
                # indices = labels == label
                indices = np.argwhere(labels == label)
                indices = indices.reshape([len(indices)])
                # print('indices shape:', indices.shape)
                
            # print('#samples in this cluster:', np.sum(indices))
            print('#samples in this cluster:', len(indices))
            sample_num = len(indices)
            
            # choose: PCA reduced by sklearn or reduced by packages?
            pca_model = None
            if use_pca:
                num_pc = max(1, sample_num//10)
                print('using %d principle components' % num_pc)
                raw_cluster_count_matrix = count_matrix[indices, ...]
                pca_model = PCA(
                    n_components=num_pc,
                    random_state=7).fit(raw_cluster_count_matrix)

                ratio = sum(pca_model.explained_variance_ratio_)
                print('total cluster pca explained variance ratio:', ratio)

                pca_count_matrix = pca_model.transform(raw_cluster_count_matrix)                
                label_count_matrix = pca_count_matrix
                label_velocities = pca_model.transform(
                    velocities[indices, ...])
            else:
                label_count_matrix = adata.X.todense()[indices, ...]
                label_velocities = velocities[indices, ...]


            '''
            Neighbor MAR part
            '''
            mses, r2s, max_eigenval_reals, feature_norms, bias_norms, jacobs = neighbor_MAR(label_count_matrix, label_velocities, neighbor_num=config.MAR_neighbor_num, pca_model=pca_model)
            cluster_centroid_sample_mses, is_centroid_neighbor_indicator, closest_sample_id_in_indices = centroid_neighbor_MAR(label_count_matrix, label_velocities, neighbor_num=config.MAR_neighbor_num, pca_model=pca_model)

            analyze_group_jacobian(adata[indices], jacobs, pca_model=pca_model, group_name='cluster'+str(label))
            eigen_stable_subset = np.array(max_eigenval_reals) < 0
            analyze_group_jacobian(adata[indices[eigen_stable_subset]],
                                   np.array(jacobs)[eigen_stable_subset, ...],
                                   pca_model=pca_model,
                                   group_name='cluster'+str(label) + 'eigenStable')
            # dont let whole data cluster overwrites everything
            if label != whole_data_label:
                adata.obs['clusters_neighbor_MAR_mse'][indices] = mses
                adata.obs['clusters_neighbor_MAR_r2'][indices] = r2s

                adata.obs['clusters_centroid_neighborhood_sample_mses'][indices[is_centroid_neighbor_indicator]] = cluster_centroid_sample_mses

                # debug check
                # print(indices)
                # print(adata.obs['clusters_centroid_neighborhood_sample_mses'][indices[is_centroid_neighbor_indicator]].shape)
                # print(is_centroid_neighbor_indicator.shape, cluster_centroid_sample_mses.shape)
                # print(cluster_centroid_sample_mses)
                # print(is_centroid_neighbor_indicator)                
                adata.obs['is_centroid_neighbor_indicator'][indices[is_centroid_neighbor_indicator]] = label

            else:
                adata.obs['whole_neighbor_MAR_mse'][indices] = mses
                adata.obs['whole_neighbor_MAR_r2'][indices] = r2s
                adata.obs['whole_neighbor_max_eigenVal_real'][indices] = max_eigenval_reals
                adata.obs['whole_neighbor_avg_vel_norms'] = feature_norms
                adata.obs['whole_neighbor_bias_norms'] = bias_norms
                adata.obs['whole_data_centroid_sample_errors'][indices[is_centroid_neighbor_indicator]] = cluster_centroid_sample_mses
                adata.obs['is_whole_centroid_neighbor'][indices[is_centroid_neighbor_indicator]] = 1
                adata.obs['is_whole_centroid_neighbor'][indices[closest_sample_id_in_indices]] = 3

                for i in range(len(indices)):
                    index = indices[i]
                    adata.uns['neighbor_jacobs'][index] = jacobs[i]

                
                
            '''
            whole cluster MAR
            '''
            model = Lasso(alpha=config.lasso_alpha).fit(label_count_matrix, label_velocities)
            predicted_velocities = model.predict(label_count_matrix)
            diff = predicted_velocities - label_velocities
            diff = np.sum(diff**2, axis=-1)
            errors[indices] = diff
            analyze_specific_cluster(adata,
                                     indices,
                                     predicted_velocities,
                                     diff,
                                     model,
                                     label,
                                     pca_model,
                                     emt_genes)

            r2_score = model.score(label_count_matrix, label_velocities)
            # explained_variance = sklearn.metrics\
            #                             .explained_variance_score(label_velocities,
            #                                                       predicted_velocities)
            mse = sklearn.metrics\
                         .mean_squared_error(label_velocities,
                                             predicted_velocities)
            print('label:%d, r^2 score: %.5f, mse:%.5f'
                  % (label, r2_score, mse))
            # scv.pl.scatter(adata, color=[ 'root_cells', 'end_points', 'errors', 'kmeans_labels'], save='error_root_end_points.png')
            if label == whole_data_label:
                adata.obs['whole_data_MAR_squared_error'] = errors
                pass
            else:
                adata.obs['cluster_squared_error'] = errors
                
        if not only_whole_data:
            scv.pl.scatter(
                adata,
                color=[
                    'root_cells',
                    'end_points',
                    'cluster_squared_error',
                    'whole_data_MAR_squared_error',
                    'Clusters'],
                save='error_root_end_points.png',
                show=False)
            scv.pl.velocity_embedding_stream(
                adata,
                color=[
                    'whole_neighbor_MAR_r2',
                    'whole_neighbor_MAR_mse',
                    'whole_neighbor_max_eigenVal_real',
                    'whole_data_MAR_squared_error',
                    'clusters_neighbor_MAR_r2',
                    'clusters_neighbor_MAR_mse',
                    'clusters_centroid_neighborhood_sample_mses',
                    'is_centroid_neighbor_indicator',
                    'Clusters',
                    'vel_norms'],
                colorbar=True,
                save='neighbor_MAR_stats.png',
                show=False)

        else:
            scv.pl.scatter(
                adata,
                color=[
                    'root_cells',
                    'end_points',
                    'whole_data_MAR_squared_error',
                    'Clusters'],
                save='error_root_end_points.png',
                show=False)
            scv.pl.velocity_embedding_stream(adata,
                                             color=[
                                                 'whole_neighbor_MAR_r2',
                                                 'whole_neighbor_MAR_mse',
                                                 'whole_neighbor_bias_norms',
                                                 'whole_neighbor_avg_vel_norms',
                                                 'whole_neighbor_max_eigenVal_real',
                                                 'whole_data_MAR_squared_error',
                                                 'Time',
                                                 'Clusters',
                                                 'vel_norms'],
                                             save='neighbor_MAR_stats.png',
                                             show=False)
        scv.pl.scatter(
            adata,
            color=[
                'root_cells',
                'end_points',
                'whole_data_centroid_sample_errors',
                'is_whole_centroid_neighbor',
                'vel_norms'],
            save='artificial_center_MAR.png',
            show=False)
        adata.obs.to_csv('./figures/adata_obs.csv')
        
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

    main_MAR(only_whole_data=config.only_whole_data, use_pca=config.use_pca)
    # main_graphlasso()


if __name__ == '__main__':
    main()
