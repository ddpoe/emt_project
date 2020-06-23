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
    
def main():
    loom_data_path = '../data/a549_tgfb1.loom'
    meta_path = '../data/a549_tgfb1_meta.csv'
    adata=dyn.read_loom(loom_data_path)
    meta=pd.read_csv(meta_path)
    print('data read complete', flush=True)
    cells=meta['Unnamed: 0'].to_numpy()
    
    # print('flag1')
    # treatment=np.array([[meta['Time'][np.squeeze(np.argwhere(cells==cell))]][0] for cell in adata.obs_names])
    # print('flag2')
    # adata.obs['treatment']=treatment
    # print('flag3')
    
    count_matrix=adata.X.todense()[:2, ...]
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    # scv.pl.velocity_embedding_stream(adata)
    velocities=adata.layers['velocity']
    count_matrix=adata.X.A
    run_graphlasso(count_matrix, prefix='count_matrix')
    run_graphlasso(velocities, prefix='velocity')

if __name__ == '__main__':
    main()
