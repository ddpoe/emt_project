import scipy
import scipy.io
import pandas as pd
import os
import numpy as np
import anndata

def test():
    dir_path = '../data/MCF10A_japan_data/'

    # load GSE data
    mat = scipy.io.mmread(os.path.join(dir_path, 'matrix.mtx'))
    barcodes = pd.read_csv(os.path.join(dir_path, 'barcodes.tsv'), header=None, delimiter='\t')
    cbc_gbc_frame = pd.read_csv(os.path.join(dir_path, 'CBC_GBC_summary.txt'), delimiter='\t')
    genes = pd.read_csv(os.path.join(dir_path, 'genes.tsv'), header=None, delimiter='\t')
    
    gse_cbc_codes = list(barcodes[0])
    gse_cbc_codes = [s[:-2] for s in gse_cbc_codes] # get rid of '-1'
    gse_cbc2gbc = {}
    for i in range(len(cbc_gbc_frame)):
        # print(i)
        cbc = cbc_gbc_frame['CBC'][i]
        gbc = cbc_gbc_frame['GBC'][i]
        gse_cbc2gbc[cbc] = gbc
    print('matched #cell barcodes between CBC_GBC_summary and GSE data:', len(set(gse_cbc2gbc.keys()) & set(gse_cbc_codes)))

    
    # check common barcodes between our processed dataset from raw and GSE data
    ann_obj = anndata.read_h5ad('../data/MCF10A/MCF10A/data/scVelo/MCF10A_filtered.h5ad')
    cbc_codes = ann_obj.obs_names
    print('matched #cell barcodes between GSE and our processed data:', len(set(cbc_codes) & set(gse_cbc_codes)))

    
if __name__ == '__main__':
    test()
