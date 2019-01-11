# About
This project uses DepMap datasets and tumor data to identify SDL partners of CENPA and is conducted in collaboration with the Basrai Group at National Cancer Institute. 

# Usage
This repo requires `conda`. Install and activate the `conda` environment with:
```
conda env create -f environment.yml
conda activate CENPA-SDL
```
Download the raw data and generate output with:
```
snakemake all
```
The outputs include both the pan-cancer and tissue-specific analysis

# Data
The data used for this analysis can be found in `data/` which will contain:

- `CCLE_RPKM.gct`: RNAseq gene expression data for 1019 cancer cell lines (in RPKM), Cancer Cell Line Encyclopedia, the Broad Institute, released on January 2nd, 2019

- `DepMap_CRISPR.csv`: Genetic Dependency, CRISPR, Public 18Q4, CERES, the DepMap project, the Broad Institute, released in November 2018

- `DepMap_RNAi.csv`: Genetic Dependency, Combined RNAi, the DepMap project, the Broad Institute, released in November 2018
- `DepMap_cell_line_metadata.csv`: Cell line metadata about cell lines in the 18Q4 release, including mapping between DepMap ID and CCLE names, the DepMap project, the Broad Institute
- `TCGA_data_ISLE.RData`: Processed TCGA data from ISLE with information of 8479 samples. Information needed includes RNAseq for of 19001 genes and survival information (survival time in days + survival outcome in terms of alive or dead) of each sample.
