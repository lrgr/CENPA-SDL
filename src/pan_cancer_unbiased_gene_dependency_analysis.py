#This script does the following:
# (1) Find the cell lines with high expression level of gene of interest (GOI). 
# Here "high expression level" is defined as the top 1/3-quantile with regards to RPKM value)
# (2) For each gene in DepMap gene dependency data, perform wilcoxon ranksum test 
# to see if cell lines with high GOI expression level has higher or worse dependency level than the remaining cell lines.

import pandas as pd
import numpy as np
import scipy as sp
import argparse
import sys
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument('-g','--gene',type = str,required = True)
parser.add_argument('-rs','--rpkm_subset',type = str,required = True)
parser.add_argument('-gds','--gene_dependency_subset',type = str,required = True)
parser.add_argument('-o','--output',type = str, required=True)
args = parser.parse_args(sys.argv[1:])

#load datasets
df_rpkm = pd.read_csv(args.rpkm_subset,sep = '\t',header=0, index_col=0 )
df_gd = pd.read_csv(args.gene_dependency_subset,sep = '\t', header = 0, index_col =0)

# (1) Find the cell lines with high expression level of gene of interest (GOI). 
# find cells with high GOI expresison level (defined as top 1/3 quartile) and the remaining cells (the remaining 2/3-quantile)
cutoff = round((1/3)*len(df_rpkm))
df_high_GOI_expression = pd.DataFrame(df_rpkm.iloc[0:cutoff])
high_GOI_cells = list(df_high_GOI_expression.index)
non_high_GOI_cells = set(df_rpkm.index)-set(high_GOI_cells)

print('Among {} cell lines, {} have high {} expression level'.format(len(df_rpkm),len(high_GOI_cells),args.gene))

# (2) For each gene in DepMap gene dependency data, perform wilcoxon ranksum test 
# to see if cell lines with high GOI expression level has higher or worse dependency level than the remaining cell lines.
# create empty dict for information needed
stats = dict()
pvals  = dict()
high_means = dict()
high_stds = dict()
remain_means = dict()
remain_stds = dict()

# iterate over each of the 17634 genes in the DepMap dataset and perform wilcoxon ranksum test between gene dependency score 
# of cell lines with high expression level of GOI and the remaining cell lines. The results are saved as dict entries.
for gene in list(df_gd.columns):
    df_g = pd.DataFrame(df_gd[gene])
    high_GOI_gen_dep = pd.DataFrame(df_g.loc[high_GOI_cells])
    non_high_GOI_gen_dep = pd.DataFrame(df_g.loc[non_high_GOI_cells])
    high_means[gene] = np.mean(high_GOI_gen_dep[gene])
    high_stds[gene] = np.std(high_GOI_gen_dep[gene])
    remain_means[gene] = np.mean(non_high_GOI_gen_dep[gene])
    remain_stds[gene]= np.std(non_high_GOI_gen_dep[gene])
    if high_means[gene]> remain_means[gene]:
        mwu_alt = 'greater'
    else:
        mwu_alt = 'less'
    mwu = mannwhitneyu(high_GOI_gen_dep[gene],non_high_GOI_gen_dep[gene], use_continuity = False, alternative =mwu_alt)
    stats[gene] = mwu[0]
    pvals[gene] = mwu[1]
   
# combine the results into a data frame
results = [stats,pvals,high_means,high_stds,remain_means,remain_stds]
df_results = pd.DataFrame(results).T
df_results.columns = ['mwu_stat','mwu_pval','high_mean','high_std','remain_mean','remain_std']

# select significant data based on pval and saved it to another df
df_results_sig = pd.DataFrame(df_results[df_results['mwu_pval']<=0.05])

# add a column to show whether the difference in gene dependency scores indicates proliferative advantange or disadvantage
# higher scores = proliferative advantage, lower scores = disadvantage
df_results_sig['adv_or_disadv'] = df_results_sig['high_mean']>df_results_sig['remain_mean']
df_results_sig =df_results_sig.replace({True:'adv',False:'disadv'})

# add a column to show the absolute difference between gene dependency scores of cells with high GOI expression level and the remaining cells
df_results_sig['diff'] = abs(df_results_sig['high_mean']-df_results_sig['remain_mean'])

#sort results by absolute difference
df_results_sig = df_results_sig.sort_values(by = ['diff'],ascending= False)

#save results to tsv
df_results_sig.to_csv(args.output,sep= '\t')