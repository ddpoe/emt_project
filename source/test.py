from utils import *


def test_process_kazu_data():
    adata = scv.read_loom(config.kazu_loom_data_path)
    # add concentration information
    adata = process_kazu_loom_data(adata, config.kazu_cbc_gbc_mapping_path, config.kazu_gbc_info_path)
    adata = process_kazu_loom_data(adata, config.kazu_cbc_gbc_mapping_path, config.kazu_gbc_info_path)


if __name__ == '__main__':
    test_process_kazu_data()
    
