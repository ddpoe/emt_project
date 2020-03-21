#!/home/dap182/anaconda3/envs/kb/bin/python

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import sys
import anndata
from anndata import AnnData

if len(sys.argv)==1:
    print()
    print(sys.argv[0],"<filename> <gene list> <output file>")
    print()
    exit()

filename=sys.argv[1]

output_path=sys.argv[3]

print()
print('Beginning to load in the data')
print()
if filename.endswith('.h5ad'):
    adata=anndata.read_h5ad(filename)
elif filename.endswith('.loom'):
    adata=anndata.read_loom(filename)
else:
    print('file type not currently supported')
    exit()

gene_list=np.loadtxt(sys.argv[2],dtype=str)

gene_idxs=np.hstack([np.where(adata.var_names==i)[0] for i in gene_list])

adata=adata[:,gene_idxs]

print()
print('Data subsetted and beginning to save file')
print()
if output_path.endswith('.h5ad'):
    AnnData.write_h5ad(adata,output_path)
elif output_path.endswith('.loom'):
    AnnData.write_loom(adata,output_path)
else:
    print('output file type not currently supported')
    exit()
