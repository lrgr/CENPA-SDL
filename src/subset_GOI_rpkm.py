import pandas as pd
import numpy as np
import scipy as sp
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-r','--rpkm',type = str,required = True)
parser.add_argument('-gd','--gene_dependency_data',type = str, required=True)
parser.add_argument('-md','--metadata',type = str, required= True)
parser.add_argument('-o','--output',type = str, required=True)
args = parser.parse_args(sys.argv[1:])

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
gd_cells = list(set(df_rpkm['Cell Line']).intersection(set(df_gd['CCLE_name'])))

# subset the rpkm data so that it only contains RPKM data of cell lines found in both RPKM and gene dependency data
df_rpkm = df_rpkm.set_index('Cell Line')
df_rpkm_subet = df_rpkm.loc[gd_cells]
df_rpkm_subet = df_rpkm_subet.sort_values(by = ['RPKM'],ascending=False)
df_rpkm_subet.to_csv(args.output,sep = '\t')

