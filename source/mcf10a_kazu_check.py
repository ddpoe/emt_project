import scipy
import scipy.io
import pandas as pd
import os
import numpy as np
import anndata

def test():
    dir_path = '../data/MCF10A_exp1/'
    # load GSE data
    mat = scipy.io.mmread(os.path.join(dir_path, 'matrix.mtx'))
    barcodes = pd.read_csv(os.path.join(dir_path, 'barcodes.tsv'), header=None, delimiter='\t')
    cbc_gbc_frame = pd.read_csv(os.path.join(dir_path, 'CBC_GBC_summary.txt'), delimiter='\t')
    genes = pd.read_csv(os.path.join(dir_path, 'genes.tsv'), header=None, delimiter='\t')
    
    gse_cbc_codes = list(barcodes[0])
    gse_cbc_codes = [s[:-2] for s in gse_cbc_codes] # get rid of '-1'
    kazu_cbc2gbc = {}
    for i in range(len(cbc_gbc_frame)):
        # print(i)
        cbc = cbc_gbc_frame['CBC'][i]
        gbc = cbc_gbc_frame['GBC'][i]
        kazu_cbc2gbc[cbc] = gbc
        
    print('matched #cell barcodes between CBC_GBC_summary and GSE data:', len(set(kazu_cbc2gbc.keys()) & set(gse_cbc_codes)))
    print(list(kazu_cbc2gbc.keys())[:5])
    print('------')
    print(gse_cbc_codes[:5])
    
    # check common barcodes between our processed dataset from raw and GSE data
    ann_obj = anndata.read_h5ad(os.path.join(dir_path, 'MCF10A_filtered.h5ad'))
    cbc_codes = ann_obj.obs_names
    print('matched #cell barcodes between GSE and our processed data:', len(set(cbc_codes) & set(gse_cbc_codes)))
    
    old_cbc_gbc_frame = pd.read_csv(os.path.join(dir_path, 'CBC_GBC_summary.txt'), delimiter='\t')
    old_cbcs = list(old_cbc_gbc_frame['CBC'])
    print(old_cbcs[:5])
    print('matched #cell barcodes between two versions:', len(set(old_cbcs) & set(kazu_cbc2gbc.keys())))

def write_barcodes():
    dir_path = '../data/MCF10A_exp1/'
    kazu_onedrive_exp1_barcodes = pd.read_csv(os.path.join(dir_path, 'barcodes.tsv'), header=None, delimiter='\t')[0]
    kazu_onedrive_exp1_barcodes = [s[:-2] for s in kazu_onedrive_exp1_barcodes]
    with open(os.path.join(dir_path, 'kalisto_whitelist_barcodes.txt'), 'w+') as f:
        for barcode in kazu_onedrive_exp1_barcodes:
            f.write(barcode + '\n')

            
def get_mapping(frame):
    res = {}
    for i in range(len(frame)):
        cbc = frame['CBC'][i]
        gbc = frame['GBC'][i]
        res[cbc] = gbc
    return res


def test_kazu_new():
    dir_path = '../data/MCF10A_exp1/'
    exp1_cbc_gbc_frame = pd.read_csv(os.path.join(dir_path, 'CBC_GBC_summary.txt'), delimiter='\t')
    exp2_cbc_gbc_frame = pd.read_csv(os.path.join(dir_path, 'CBC_GBC.txt'), delimiter='\t')
    kazu_onedrive_exp1_cbc2gbc_mapping = get_mapping(exp1_cbc_gbc_frame)
    kazu_onedrive_exp1_cbcs = list(kazu_onedrive_exp1_cbc2gbc_mapping.keys())
    kazu_onedrive_exp2_cbc2gbc_mapping = get_mapping(exp2_cbc_gbc_frame)
    kazu_onedrive_exp2_cbcs = list(kazu_onedrive_exp2_cbc2gbc_mapping.keys())

    ann_obj = anndata.read_loom(os.path.join(dir_path, 'possorted_genome_bam_RIG79.loom'))
    adata_cbc_codes = [x[x.find(':')+1:] for x in list(ann_obj.obs_names)]
    # ann_obj = anndata.read_h5ad(os.path.join(dir_path, 'adata.h5ad'))
    # adata_cbc_codes = list(ann_obj.obs_names)


    kazu_onedrive_exp1_barcodes = pd.read_csv(os.path.join(dir_path, 'barcodes.tsv'), header=None, delimiter='\t')[0]
    kazu_onedrive_exp1_barcodes = [s[:-2] for s in kazu_onedrive_exp1_barcodes]
    whitelist_barcodes = pd.read_csv(os.path.join(dir_path, '10xv2_whitelist.txt'), header=None, delimiter='\t')[0]
    print('sample outputs:')
    print(sorted(list(adata_cbc_codes))[:10])
    print(sorted(kazu_onedrive_exp1_cbcs)[:10])
    print(sorted(kazu_onedrive_exp1_barcodes)[:10])
    
    print('total obs (cbc) in ann data:', len(adata_cbc_codes))
    print('len of barcode onedrive file:', len(kazu_onedrive_exp1_barcodes))
    print('len of mapping onedrive file:', len(kazu_onedrive_exp1_cbcs))
    print('len of whitelist:', len(whitelist_barcodes))
    print('matched #cell barcodes between adata and kazu onedrive exp1:', len(set(adata_cbc_codes) & set(kazu_onedrive_exp1_cbcs)))
    print('matched #cell barcodes between adata and kazu onedrive exp2:', len(set(adata_cbc_codes) & set(kazu_onedrive_exp2_cbcs)))
    print('matched #cell barcodes between kazu onedrive barcodes and adata:', len(set(adata_cbc_codes) & set(kazu_onedrive_exp1_barcodes)))
    print('matched #cell barcodes between kazu onedrive barcodes and kazu onedrive CBC_GBC mapping:', len(set(kazu_onedrive_exp1_cbcs) & set(kazu_onedrive_exp1_barcodes)))
    
    print('matched #cell barcodes between whitelist barcodes and annData:', len(set(whitelist_barcodes) & set(adata_cbc_codes)))    
    print('matched #cell barcodes between whitelist barcodes and kazu onedrive barcodes:', len(set(whitelist_barcodes) & set(kazu_onedrive_exp1_barcodes)))
    print('matched #cell barcodes between whitelist barcodes and kazu onedrive CBC_GBC mapping:', len(set(kazu_onedrive_exp1_cbcs) & set(whitelist_barcodes)))

    
if __name__ == '__main__':
    # test()
    test_kazu_new()
    # write_barcodes()
