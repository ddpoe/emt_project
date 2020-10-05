The code in this folder was modified from the following tutorials from the goatools github. One is the processing of the genes and go catergories so they are converted to a dictonary that can be used with the program (I think. I did not look into what their processing script does). The other one takes this in with some other information that can be downloaded from the script and gives the go categories for the gene list given. I also think I downloaded another one from the github as well. 

The tutorial uses mouse genes as an example so I changed the tax id from 10090 to 9606, which is the tax id for humans

The only script that needs to be dealt with is the gene_enrichment.py. You pass in a gene list and you will need to convert this list to the gene ids that coorespond to the gene symbols. There is a section of this code that does this. I load in a gene list into a variable called CDK1_genes or something like that. If you change this variable throughout the script if should work fine. 

Tutorials used:
https://github.com/tanghaibao/goatools/blob/main/notebooks/backround_genes_ncbi.ipynb
https://github.com/tanghaibao/goatools/blob/main/notebooks/goea_nbt3102.ipynb
