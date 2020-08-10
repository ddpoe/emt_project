import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--only-whole-data', type=bool, help='do MAR in whole data only, or in each cluster as well?')
parser.add_argument('--use-pca', type=bool, help='use pca in each cluster or not')
parser.add_argument('--include-a549-days', help='a list of days', nargs='*', default='all')
parser.add_argument('--use-emt-gene-filter', type=bool, help='')
parser.add_argument('--selected-genes-jacobian', type=str, help='', default=[], nargs='*')
parser.add_argument('--MAR-neighbor-num', type=int, help='')
parser.add_argument('--use-dataset', type=str, help='', choices=['a549', 'pancreas', 'cook'])
parser.add_argument('--lasso-alpha', type=float, help='', default=0.0001)
args = parser.parse_args()
print('arguments:', args)
only_whole_data = args.only_whole_data
use_pca = args.use_pca
use_emt_gene_filter = args.use_emt_gene_filter
include_a549_days = args.include_a549_days
selected_genes_jacobian = args.selected_genes_jacobian
MAR_neighbor_num = args.MAR_neighbor_num
use_dataset = args.use_dataset
lasso_alpha = args.lasso_alpha

# only_whole_data = False # 
# use_pca = True
# use_emt_gene_filter = False
# day0_only = ['0d', '8h', '1d', '3d', '7d'] # which time subset do we want to use? set to 'all' to include all day data
# day0_only = 'all' # which time subset do we want to use? set to 'all' to include all day data
# selected_genes_jacobian = ['FN1', 'SNAI2', 'VIM', 'GEM']
# selected_genes_jacobian = ['FN1', 'PMEPA1', 'SETBP1', 'PLA2G4A', 'TM4SF20' , 'SMIM14']
two_gene_graph_dir = './figures/two_gene_vector_field'
# MAR_neighbor_num = 40
# use_pancreas_data = False
# use_dataset = 'a549'
# use_dataset = 'pancreas'
# calculate_velocity_with_all_gene = True
