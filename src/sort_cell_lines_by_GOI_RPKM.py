import pandas as pd
import numpy as np
import scipy as sp
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-g','--gene',type = str,required = True)
parser.add_argument('-e','--expression_data',type = str, required = True) 
parser.add_argument('-o','--output',type = str, required=True)
args = parser.parse_args(sys.argv[1:])

# load CCLE RNA seq data (RPKM)
# skipped first two rows as they contains information of dataset version which is not needed for the analysis
# header name are cell line names
df = pd.read_csv(args.expression_data, sep = '\t', header=0, skiprows = range(0,2)) 

###SORT BY EXPRESSION VALUE OF GENE OF INTEREST"
# extract expression values of the gene of interest and store it as a seperate dataframe GOI
# each row is a cell line and the RPKM value of the cell line 
GOI = df[df['Description']==args.gene].T

# clean GOI by removing unwanted header info that comes from df transposition and renaming columns
GOI.reset_index(inplace=True)
GOI.columns= ['Cell Line','RPKM']
GOI = GOI.drop(GOI.index[0:2])

# sort by RPKM value in descending order
GOI_sorted = GOI.sort_values(by=['RPKM'],ascending=False) 
GOI_sorted.reset_index(inplace=True,drop=True) #reset index

# save GOI_sorted to csv
GOI_sorted.to_csv(args.output,sep = '\t')



