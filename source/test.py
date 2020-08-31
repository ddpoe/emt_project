from utils import *
from fokker_planck_analysis import *


def test_process_kazu_data():
    adata = scv.read_loom(config.kazu_loom_data_path)
    # add concentration information
    adata = process_kazu_loom_data(
        adata,
        config.kazu_cbc_gbc_mapping_path,
        config.kazu_gbc_info_path)
    adata = process_kazu_loom_data(
        adata,
        config.kazu_cbc_gbc_mapping_path,
        config.kazu_gbc_info_path)


def test_dynamo():
    adata = scv.read_loom(config.a549_loom_data_path)
    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata)

    dyn.tl.reduceDimension(adata)
    dyn.tl.cell_velocities(adata)
    dyn.tl.cell_velocities(adata, basis='pca')
    # adata.obsm['velocity_pca'].shape
    dyn.tl.cell_wise_confidence(adata)
    dyn.vf.VectorField(adata)
    # dyn.vf.VectorField(adata, basis='pca')
    vectorfield = adata.uns['VecFld']['VecFld']
    dyn.pl.topography(adata)
    dyn.pl.save_fig(path='./results/test', ext='png')


def test_fokker():
    use_real_data = True
    use_real_data = False
    if use_real_data:
        meta = pd.read_csv(config.a549_meta_path)
        adata = scv.read_loom(config.a549_loom_data_path)
        adata = filter_a549_MET_samples(
            adata, meta, include_a549_days=config.include_a549_days)
        scv.pp.filter_and_normalize(adata, n_top_genes=config.n_top_genes)
        scv.pp.moments(adata)
        scv.tl.velocity(adata, perc=config.perc)
        X = adata.X.todense()
        V = adata.layers['velocity']
    else:
        X = np.identity(10)
        V = np.identity(10)
    fp_analyze(X, V)


if __name__ == '__main__':
    # test_process_kazu_data()
    test_fokker()
