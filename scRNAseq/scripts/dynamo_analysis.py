#!/home/dap182/anaconda3/envs/kb/bin/python

from zipfile import ZipFile
import sys
import numpy as np
import matplotlib.pyplot as plt 
import dynamo as dyn
from anndata import AnnData

if len(sys.argv)==1:
    print()
    print(sys.argv[0],"<filename>")
    print()
    exit()

filename=sys.argv[1]

dataset_name=filename.split("/")[-1].split(".")[0]
#file_extension=filename.split(".")[-1]
file_extension='h5ad'

processed_output_path=dataset_name+'_processed.'+file_extension
P_output_path=dataset_name+'_P.'+'npy'
idx_output_path=dataset_name+'_idx.'+'npy'

print()
print('Loading Data')
print()
if filename.endswith('.h5ad'):
    adata=dyn.read_h5ad(filename)
elif filename.endswith('.loom'):
    adata=dyn.read_loom(filename)
else:
    print('file type not currently supported')

print()
print('Initial Preprocessing')
print()
dyn.pp.recipe_monocle(adata, n_top_genes=2000)

print()
print('Computing Dynamics')
print()
dyn.tl.dynamics(adata, mode='moment')

print()
print('Reducing Dimension')
print()
dyn.tl.reduceDimension(adata)

print()
print('Computing Cell Velocity')
print()
dyn.tl.cell_velocities(adata)

kmc=adata.uns['kmc']

print()
print('Saving KMC Probability Matrix')
print()
np.savetxt(P_output_path,kmc.P.todense())

print()
print('Saving KMC Index\'s')
print()
np.save(idx_output_path,kmc.Idx)

print()
print('Saving Processed Data')
print()
if file_extension=='h5ad':
    AnnData.write(adata,processed_output_path)
elif file_extension=='loom':
    AnnData.write(adata,processed_output_path)


    
