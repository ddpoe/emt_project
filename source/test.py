from utils import *


def test_process_kazu_data():
    adata = scv.read_loom(config.kazu_loom_data_path)
    # add concentration information
    adata = process_kazu_loom_data(adata, config.kazu_cbc_gbc_mapping_path, config.kazu_gbc_info_path)
    adata = process_kazu_loom_data(adata, config.kazu_cbc_gbc_mapping_path, config.kazu_gbc_info_path)


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
    
if __name__ == '__main__':
    test_process_kazu_data()
    
