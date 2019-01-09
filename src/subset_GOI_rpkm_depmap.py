#This script does the following: 
# (1) Find the cell lines that are covered by both CCLE RNAseq data (RPKM) and DepMap gene dependecy data
# (2) Save the RPKM data and gene dependency data of cell lines found in both datasets to seperate tsv files for downstream analysis
# (3) Plot the expression distribution of gene of interests in cell lines to show if the selected cell lines has a skewed distribution in gene expression

import pandas as pd
import numpy as np
import scipy as sp
import argparse
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument('-g','--gene',type = str,required = True)
parser.add_argument('-r','--rpkm',type = str,required = True)
parser.add_argument('-gd','--gene_dependency_data',type = str, required=True)
parser.add_argument('-md','--metadata',type = str, required= True)
parser.add_argument('-or','--output_rpkm',type = str, required=True)
parser.add_argument('-odm','--output_depmap',type = str, required=True)
parser.add_argument('-p','--plot',type = str, required=True)
args = parser.parse_args(sys.argv[1:])

# (1) Find the cell lines that are covered by both CCLE RNAseq data (RPKM) and DepMap gene dependecy data
# load cell line metadata 
df_md = pd.read_csv(args.metadata,sep = ',')

# load gene dependency data
df_gd = pd.read_csv(args.gene_dependency_data, sep=',')

# load sorted rpkm of gene of interest from output dir
df_rpkm = pd.read_csv(args.rpkm,sep = '\t',header = 0)
df_rpkm = df_rpkm.drop(columns=['Unnamed: 0'])

# generate dict of DepMap ID and CCLE name
dict_dm_ccle= df_md.set_index('DepMap_ID')['CCLE_Name'].to_dict()

# replace DepMap ID of gene dependecy data with CCLE name 
df_gd = df_gd.replace({'DepMap_ID':dict_dm_ccle})

# clean DepMap data header as the header contains both the gene name and the numerous gene ID in parenthesis
new_col_name = ['CCLE_name']
old_col_name = list(df_gd.columns)
old_col_name.remove('DepMap_ID')
for gene in old_col_name:
    new_gene_name = gene.split(' ')[0]
    new_col_name.append(new_gene_name)
df_gd.columns = new_col_name

# generate list of cell lines found in both gene dependency data and the RPKM data
overlapping_cells = list(set(df_rpkm['Cell Line']).intersection(set(df_gd['CCLE_name'])))

# (2) Save the RPKM data and gene dependency data of cell lines found in both datasets to seperate tsv files for downstream analysis
# subset the rpkm data so that it only contains RPKM data of cell lines found in both RPKM and gene dependency data
df_rpkm = df_rpkm.set_index('Cell Line')
df_rpkm_subet = pd.DataFrame(df_rpkm.loc[overlapping_cells])
df_rpkm_subet = df_rpkm_subet.sort_values(by = ['RPKM'],ascending=False)
df_rpkm_subet.to_csv(args.output_rpkm,sep = '\t')

print('RPKM subset completed')

#subset the depmap data so that it only contains depmap data of cell lines found in both RPKM and gene dependency data
df_gd = df_gd.set_index('CCLE_name')
df_gd_subset = df_gd.loc[overlapping_cells]
df_gd_subset.to_csv(args.output_depmap,sep = '\t')

print('DepMap subset completed')

# (3) Plot the expression distribution of gene of interests in cell lines to show if the selected cell lines has a skewed distribution in gene expression
#combine the RPKM into the same dataframe for plotting
ccle_rpkm = df_rpkm.reset_index()
ccle_rpkm['Data'] = ['CCLE' for cell in range(len(ccle_rpkm))]

subset_rpkm = df_rpkm_subet.reset_index()
subset_rpkm['Data'] = ['CCLE overlapping with DepMap' for cell in range(len(subset_rpkm))]

## find cells in CCLE RPKM data that do not overlap with cells in Dep Map data
subset_status = [cell in overlapping_cells for cell in list(ccle_rpkm['Cell Line']) ]
ccle_rpkm['subset status'] = subset_status
non_subset_rpkm = pd.DataFrame(ccle_rpkm[ccle_rpkm['subset status']==False])
non_subset_rpkm['Data'] = ['Non-overlapping' for cell in range(len(non_subset_rpkm))]
ccle_rpkm = ccle_rpkm.drop(columns = ['subset status'])
non_subset_rpkm = non_subset_rpkm.drop(columns = ['subset status'])
frames = [ccle_rpkm,subset_rpkm,non_subset_rpkm]
plot_ready= pd.concat(frames)

print('plot ready dataframe completed')

# compute wilcoxon ranksum p value, and mean and std of each category of the data
mwu = mannwhitneyu(non_subset_rpkm['RPKM'],subset_rpkm['RPKM'],use_continuity = False, alternative ='greater')
ccle_mean = str(round(np.mean(ccle_rpkm['RPKM']),3))
ccle_std = str(round(np.std(ccle_rpkm['RPKM']),3))
non_subset_mean = str(round(np.mean(non_subset_rpkm['RPKM']),3))
non_subset_std = str(round(np.std(non_subset_rpkm['RPKM']),3))
subset_mean = str(round(np.mean(subset_rpkm['RPKM']),4))
subset_std = str(round(np.std(subset_rpkm['RPKM']),3))
mwu_pvalue = str(round(mwu[1],3))

print('wilcoxon ranksum test completed')

#plot stripplot and boxplot of the distribution
sns.set(style="whitegrid", color_codes=True,font_scale=2.5)
strip_color = {'CCLE':'#0404B4','CCLE overlapping with DepMap':'#B40404','Non-overlapping':'#00802B'}
box_color = {'CCLE':'#819FF7','CCLE overlapping with DepMap':'#F78181','Non-overlapping':'#80FFAA'}
fig, ax = plt.subplots()
fig.set_size_inches(20, 20)
sns.stripplot(x='Data', y = 'RPKM',data = plot_ready, 
              jitter= True, size = 6, 
              palette = strip_color).set_title('Comparing {} Expression Distribution'.format(args.gene),fontsize = 40,y=1.03)
sns.boxplot(x='Data', y = 'RPKM',
            data = plot_ready,
            width = 0.35, 
            palette = box_color )

ax.set_xlabel('Data',fontsize=35)
ax.set_ylabel('RPKM',fontsize=35)
ax.set_ylim(-5,50)

#add wilcoxon ranksum test results, mean, and std of each category of data to the plot
plt.text(0,45.5,'mean ='+ccle_mean,horizontalalignment='center',fontsize=24)
plt.text(0,44.2,'std ='+ccle_std,horizontalalignment='center',fontsize=24)
plt.text(1,45.5,'mean ='+subset_mean,horizontalalignment='center',fontsize=24)
plt.text(1,44.2,'std = '+subset_std,horizontalalignment='center',fontsize=24)
plt.text(2,45.5,'mean ='+non_subset_mean,horizontalalignment='center',fontsize=24)
plt.text(2,44.2,'std = '+non_subset_std,horizontalalignment='center',fontsize=24)
plt.text(1,48.5,'Wilcoxon Rank Sum between overlapping and non-overlapping p-value = '+mwu_pvalue,horizontalalignment='center',fontsize=22)

# save image
plt.savefig(args.plot,
            bbox_inches='tight', 
            pad_inches=0.5)

print('RPKM distribution image saved')

