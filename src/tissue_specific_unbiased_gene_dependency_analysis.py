# This script performs the following tasks:
# (1) Count the number of cell lines in each tissue type in the truncated gene dependency (GD) data 
# (the truncated GD data has only data of cel lines found in both CCLE RNAseq data and GE data )
# (2) Find the top 3 tissue types with the most cell lines in the data
# (3) For each of the three tissue types, perform Wilcoxon rank sum test to compare the CENPA expression levels 
# between cell lines that belong to the tissue type and those that do not
# (4) Plot the CENPA expression distribution between cell lines belonging to the tissue of interest and those that don't
# (5) For each of the three tissue types, perform tissue-specific wilcoxon ranksum test for difference in gene dependency scores 
# between high CENPA expressing cell lines and the remaining cell lines


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
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-g','--gene',type = str,required = True)
parser.add_argument('-rs','--rpkm_subset',type = str,required = True)
parser.add_argument('-gds','--gene_dependency_subset',type = str,required = True)
parser.add_argument('-od','--output_dir',type = str,required = True)

args = parser.parse_args(sys.argv[1:])

#load datasets
df_rpkm = pd.read_csv(args.rpkm_subset,sep = '\t',header=0, index_col=0 )
df_gd = pd.read_csv(args.gene_dependency_subset,sep = '\t', header = 0, index_col =0)

#(1) Count the number of cell lines in each tissue type in the truncated (GD) data 
Tissue = [cell.split('_',1)[1] for cell in list(df_gd.index)]
df_tissue = pd.DataFrame(dict(Counter(Tissue)),index=['number']).T
df_tissue = df_tissue.sort_values(by = ['number'],ascending= False)
print(df_tissue[0:10])

# (2) Find the top 3 tissue types with the most cell lines in the data
t = [(' ').join(df_tissue.index[i].split('_')).capitalize() for i in range(0,3)]
print('The top 3 tissue types with the most cell lines: (1) {} (2) {} (3) {}'.format(*t))
TOI = [df_tissue.index[i] for i in range(0,3)]

#(3) For each of the three tissue types, perform Wilcoxon rank sum test to compare the CENPA expression levels 
# between cell lines that belong to the tissue type and those that do not
def split_df_by_tissue_type(df,toi):
    toi_lines = [cell for cell in list(df.index) if cell.split('_',1)[1]==toi]
    non_toi_lines = list(set(df.index)-set(toi_lines))
    df_toi = pd.DataFrame(df.loc[toi_lines])
    df_non_toi = pd.DataFrame(df.loc[non_toi_lines])
    return (df_toi,df_non_toi,toi_lines,non_toi_lines)

def wilcoxon_across_tissue_rpkm(df,toi):
    split = split_df_by_tissue_type(df,toi)
    df_toi = split[0]
    df_non_toi = split[1]
    mean_toi = np.mean(df_toi['RPKM'])
    std_toi = np.std(df_toi['RPKM'])
    mean_non_toi = np.mean(df_non_toi['RPKM'])
    std_non_toi = np.std(df_non_toi['RPKM'])
    mean_all = np.mean(df['RPKM'])
    std_all = np.std(df['RPKM'])
    if mean_toi > mean_non_toi:
        mwu_alt = 'greater'
    else:
        mwu_alt = 'less'
    mwu = mannwhitneyu(df_toi['RPKM'],df_non_toi['RPKM'],use_continuity = False, alternative =mwu_alt)
    return (df_toi,df_non_toi,mwu[0],mwu[1],mean_all,std_all,mean_toi,std_toi,mean_non_toi,std_non_toi)

# store the results as dict of tuples
tissue_results = dict()
for tissue in TOI:
    tissue_results[tissue] = wilcoxon_across_tissue_rpkm(df_rpkm,tissue)

print('Wilcoxom Rank Sum Tests performed to examine respectively if {} expression levels in cell lines that belong to (1) {} or (2) {} or (3) {} are significantly different from those that do not'.format(args.gene,*t))

# (4) Plot the CENPA expression distribution between cell lines belonging to the tissue of interest and those that don't
def plot_data_prep(df,results,toi):
    # clean tissure name
    toi_cap = (' ').join(toi.split('_')).capitalize()
    toi_lower = (' ').join(toi.split('_')).lower()
    colname = 'Data'
    # assign data type 
    df_all = df.reset_index()
    df_all[colname] = ['All' for cell in range(len(df_all))]
    df_toi = results[0].reset_index()
    df_toi [colname] = ['{}'.format(toi_cap) for cell in range(len(df_toi))]
    df_non_toi= results[1].reset_index()
    df_non_toi[colname] = ['Non-{}'.format(toi_lower) for cell in range(len(df_non_toi))]
    #combine the above df into one df
    frames = [df_all,df_toi,df_non_toi]
    plot_ready= pd.concat(frames)
    return(plot_ready)

def round_3(x):
    return str(round(x,3))

def plot(df,results,toi):
    # generate strings of data type
    toi_cap = (' ').join(toi.split('_')).capitalize()
    toi_lower = (' ').join(toi.split('_')).lower()
    type_all = 'All'
    type_t = toi_cap
    type_nt = 'Non-{}'.format(toi_lower) 

    # title and filename
    title = '{} vs. non-{} {} expression distribution.png'.format(toi_cap,toi_lower,args.gene)
    file_name= '_'.join(title.split(' '))

    # plot canvas setup
    sns.set(style="whitegrid", color_codes=True,font_scale=2.5)
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20)

    # palette
    strip_color = {type_all:'#0404B4',type_t:'#B40404',type_nt:'#00802B'}
    box_color = {type_all:'#819FF7',type_t:'#F78181',type_nt:'#80FFAA'}

    #plot
    sns.stripplot(x='Data', y = 'RPKM',data = df, 
                  jitter= True, 
                  size = 6,
                  palette = strip_color).set_title(title,fontsize = 35,y=1.03)

    sns.boxplot(x='Data', y = 'RPKM',
                data = df,
                width = 0.35,
                palette = box_color)
    
    # label and text
    ax.set_xlabel('Data',fontsize=30)
    ax.set_ylabel('RPKM',fontsize=30)
    ax.set_ylim(-5,50)
    # (df_toi,df_non_toi,mwu[0],mwu[1],mean_all,std_all,mean_toi,std_toi,mean_non_toi,std_non_toi)
    plt.text(0,45.5,'mean = '+round_3(results[4]),horizontalalignment='center',fontsize=24)
    plt.text(0,44.2,'std = '+round_3(results[5]),horizontalalignment='center',fontsize=24)
    plt.text(1,45.5,'mean = '+round_3(results[6]),horizontalalignment='center',fontsize=24)
    plt.text(1,44.2,'std = '+round_3(results[7]),horizontalalignment='center',fontsize=24)
    plt.text(2,45.5,'mean = '+round_3(results[8]),horizontalalignment='center',fontsize=24)
    plt.text(2,44.2,'std = '+round_3(results[9]),horizontalalignment='center',fontsize=24)
    plt.text(1,48.5,'Wilcoxon Rank Sum between selectd and non-selected p-value = '+round_3(results[3]),horizontalalignment='center',fontsize=22)
    plt.savefig(args.output_dir+'{}'.format(file_name),
            bbox_inches='tight', 
            pad_inches=0.5)

for tissue in TOI:
    results = tissue_results[tissue]
    plot_ready = plot_data_prep(df_rpkm,results,tissue)
    plot(plot_ready,results,tissue)
    t = (' ').join(tissue.split('_')).lower()
    print ('plotted {} expression level distribution of {}'.format(args.gene,t))

# (5) For each of the three tissue types, perform tissue-specific wilcoxon ranksum test for difference in gene dependency scores 
# between high CENPA expressing cell lines and the remaining cell lines

def wilcoxon_within_tissue_gd(df,high_GOI_cells,non_high_GOI_cells):
    stats = dict()
    pvals  = dict()
    high_means = dict()
    high_stds = dict()
    remain_means = dict()
    remain_stds = dict()
    for gene in list(df.columns):
        df_g = pd.DataFrame(df[gene])
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
    results = [stats,pvals,high_means,high_stds,remain_means,remain_stds]
    df_results = pd.DataFrame(results).T
    df_results.columns = ['mwu_stat','mwu_pval','high_mean','high_std','remain_mean','remain_std']
    return(df_results)

def results_processing(df):
    # select statistically significant results by filtering with pval <= 0.05
    df_results_sig = pd.DataFrame(df[df['mwu_pval']<=0.05])

    # add a column to show whether the difference in gene dependency scores indicates proliferative advantange or disadvantage
    # higher scores = proliferative advantage, lower scores = disadvantage
    df_results_sig['adv_or_disadv'] = df_results_sig['high_mean']>df_results_sig['remain_mean']
    df_results_sig =df_results_sig.replace({True:'adv',False:'disadv'})

    # add a column to show the absolute difference between gene dependency scores of cells with high GOI expression level and the remaining cells
    df_results_sig['diff'] = abs(df_results_sig['high_mean']-df_results_sig['remain_mean'])

    #sort results by absolute difference
    df_results_sig = df_results_sig.sort_values(by = ['diff'],ascending= False)
    return(df_results_sig)

for tissue in TOI:
    results = tissue_results[tissue]
    rpkm_toi = results[0]
    cutoff = round((1/3)*len(rpkm_toi))
    high_GOI_cells = list(pd.DataFrame(rpkm_toi.iloc[0:cutoff]).index)
    non_high_GOI_cells = list(set(rpkm_toi.index)-set(high_GOI_cells))
    wilcoxon_results = wilcoxon_within_tissue_gd(df_gd,high_GOI_cells,non_high_GOI_cells)
    processed_results = results_processing(wilcoxon_results)
    processed_results.to_csv(args.output_dir+'Wilcoxon_rank_sum_gd_{}_{}.tsv'.format(tissue.args.gene),sep = '\t')
    print('Wilcoxon rank sum test for gene dependency scores of {} cell lines completed'.format(tissue))
