from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from genes_ncbi_9606_proteincoding import GENEID2NT as GeneID2nt_homo
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
import collections as cx


obo_fname = download_go_basic_obo()

fin_gene2go = download_ncbi_associations()

obodag = GODag("go-basic.obo")

CDK1_gene_list=list(np.loadtxt('CDK1_top_effectors.txt',dtype=str))

#def load_data(directory):
    #F_adjusted=np

# Read NCBI's gene2go. Store annotations in a list of namedtuples (9606 is the tax ID for humans)
objanno = Gene2GoReader(fin_gene2go, taxids=[9606])

# Get namespace2association where:
#    namespace is:
#        BP: biological_process               
#        MF: molecular_function
#        CC: cellular_component
#    assocation is a dict:
#        key: NCBI GeneID
#        value: A set of GO IDs associated with that gene
ns2assoc = objanno.get_ns2assc()

for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated human  genes".format(NS=nspc, N=len(id2gos)))

print()
print(len(GeneID2nt_homo))
print()

goeaobj = GOEnrichmentStudyNS(
        GeneID2nt_homo.keys(), # List of mouse protein-coding genes
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method

symbols=np.zeros(len(GeneID2nt_homo.keys()),dtype='U100')
geneids=np.zeros(len(GeneID2nt_homo.keys()),dtype=int)

#Creating a lookup table to convert the gene symbols to the gene ids needed for the gene enrichment analysis
for idx,key in enumerate(GeneID2nt_homo.keys()):
    symbols[idx]=GeneID2nt_homo[key].Symbol
    geneids[idx]=GeneID2nt_homo[key].GeneID

boolean_symbol=np.isin(symbols,CDK1_gene_list)

matches_idx=np.where(boolean_symbol)[0]

geneids_matches=list(geneids[matches_idx])

goea_quiet_all = goeaobj.run_study(geneids_matches, prt=None)
goea_quiet_sig = [r for r in goea_quiet_all if r.p_fdr_bh < 0.05]

print('{N} of {M:,} results were significant'.format(
        N=len(goea_quiet_sig),
        M=len(goea_quiet_all)
    ))

print('Significant results: {E} enriched, {P} purified'.format(
    E=sum(1 for r in goea_quiet_sig if r.enrichment=='e'),
    P=sum(1 for r in goea_quiet_sig if r.enrichment=='p')))

ctr = cx.Counter([r.NS for r in goea_quiet_sig])
print('Significant results[{TOTAL}] = {BP} BP + {MF} MF + {CC} CC'.format(
    TOTAL=len(goea_quiet_sig),
    BP=ctr['BP'],  # biological_process
    MF=ctr['MF'],  # molecular_function
    CC=ctr['CC'])) # cellular_component

#goeaobj.wr_xlsx("CDK1_test.xlsx", goea_quiet_sig)
goeaobj.wr_txt("CDK1_test.txt", goea_quiet_sig)


goid_subset = [
    'GO:0003723', # MF D04 RNA binding (32 genes)
    'GO:0044822', # MF D05 poly(A) RNA binding (86 genes)
    'GO:0003729', # MF D06 mRNA binding (11 genes)
    'GO:0019843', # MF D05 rRNA binding (6 genes)
    'GO:0003746', # MF D06 translation elongation factor activity (5 genes)
]
plot_gos("nbt3102_MF_RNA_genecnt.png", 
    goid_subset, # Source GO ids
    obodag, 
    goea_results=goea_quiet_all) # Use pvals for coloring
