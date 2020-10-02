import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
%matplotlib pyplot

F_adjusted=np.load('../fokker-plank-results/a549_0day/F_gene_adjusted.npy',allow_pickle=True)
genes=np.loadtxt('../fokker-plank-results/a549_0day/fp_gene_list.txt',dtype=str)



#Gives a 1 x gene_num array with the value at the desired quantile range. I will use these values in inequality statements
#to trim the F_adjusted matrix 
thresholds_upper=np.quantile(F_adjusted,0.9,axis=1)
thresholds_lower=np.quantile(F_adjusted,0.1,axis=1)

trimmed_matrix=np.copy(F_adjusted)

#This trims the matrix using the thresholds given above
for i in range(F_adjusted.shape[0]):
    mask = (F_adjusted[i] < thresholds_lower[i]) | (F_adjusted[i] > thresholds_upper[i]) 
    trimmed_matrix[i]=trimmed_matrix[i]*mask
    
print("The number of average sig genes per genes are",np.count_nonzero(trimmed_matrix)/F_adjusted.shape[0])
print("Percentage of nonzero elements in trimmed F: ",np.count_nonzero(trimmed_matrix)/F_adjusted.shape[0]**2*100,'%')
print("Percentage of nonzero elements in adjusted F: ",np.count_nonzero(F_adjusted)/F_adjusted.shape[0]**2*100,'%')

#example of how I got the CDK1  "significant" genes
#CDK1_arg=np.where(genes=='CDK1')

#CDK1_top_regulators=genes[np.argwhere(trimmed_matrix[:,CDK1_arg[0]]!=0).squeeze()]

#np.savetxt('CDK1_top_effectors.txt',CDK1_top_regulators,fmt='%s')
