import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    '--only-whole-data',
    help='do MAR in whole data only, or in each cluster as well?',
    default=False,
    action='store_true')
parser.add_argument(
    '--use-pca',
    help='use pca in each cluster or not',
    default=False,
    action='store_true')
parser.add_argument(
    '--include-a549-days',
    help='a list of days',
    nargs='*',
    default='all')
parser.add_argument(
    '--use-emt-gene-filter',
    help='',
    default=False,
    action='store_true')
parser.add_argument(
    '--selected-genes-jacobian',
    type=str,
    help='',
    default=[],
    nargs='*')
parser.add_argument('--MAR-neighbor-num', type=int, help='')
parser.add_argument(
    '--use-dataset',
    type=str,
    help='',
    choices=[
        'a549',
        'pancreas',
        'cook',
         'kazu_mcf10a'])
parser.add_argument('--lasso-alpha', type=float, help='', default=0.0001)
parser.add_argument('--result-dir', type=str, help='', default='./results')
parser.add_argument(
    '--kazu-dosage-range',
    nargs=2,
    type=float,
    help='',
    default=[
        0,
         float('inf')])
parser.add_argument(
    '--mode',
    type=str,
    help='',
    default='analyze_MAR',
    choices=[
        'analyze_MAR',
        'analyze_PCA',
        'analyze_fokker_planck'])
parser.add_argument('--perc', type=int, help='', default=5)

args = parser.parse_args()
emt_gene_path = '/home/ke/emt_project/data/gene_lists/emt_genes_weikang.txt'
a549_loom_data_path = '/home/ke/emt_project/data/a549_tgfb1.loom'
a549_meta_path = '/home/ke/emt_project/data/a549_tgfb1_meta.csv'
kazu_loom_data_path = '/home/ke/emt_project/data/MCF10A_exp1/possorted_genome_bam_RIG79.loom'
kazu_cbc_gbc_mapping_path = '/home/ke/emt_project/data/MCF10A_exp1/CBC_GBC_summary.txt'
kazu_gbc_info_path = '/home/ke/emt_project/data/MCF10A_exp1/gbc_dosage_info.txt'

mode = args.mode
only_whole_data = args.only_whole_data
use_pca = args.use_pca
use_emt_gene_filter = args.use_emt_gene_filter
include_a549_days = args.include_a549_days
selected_genes_jacobian = args.selected_genes_jacobian
MAR_neighbor_num = args.MAR_neighbor_num
use_dataset = args.use_dataset
lasso_alpha = args.lasso_alpha
result_dir = args.result_dir
kazu_dosage_range = args.kazu_dosage_range
perc = args.perc
n_top_genes = 2000
fp_num_pc = 30


random_state = 7
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


def gen_config_folder_name():
    arg_strs = ['mode=' + mode.replace("_", '-'),
                'dataset=' + use_dataset,
                'use-emt-gene-filter=' + str(use_emt_gene_filter),
                'MAR-neighbor-num=' + str(MAR_neighbor_num),
                'lasso-alpha=' + str(lasso_alpha),
                'perc=' + str(perc)]

    if use_dataset == 'a549':
        arg_strs.append('include-a549-days=' + str(include_a549_days))
    elif use_dataset == 'kazu_mcf10a':
        arg_strs.append('kazu-dosage-range=' + str(kazu_dosage_range))
    name = '_'.join(arg_strs)
    return name
