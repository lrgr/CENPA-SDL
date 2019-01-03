from os.path import join
SRC_DIR = 'src'
DATA_DIR = 'data'
DEPMAP_CRISPR_URL = 'https://ndownloader.figshare.com/files/13396070'
DEPMAP_RNAI_URL = 'https://ndownloader.figshare.com/files/11489681'
TCGA_URL = 'ftp://ftp.umiacs.umd.edu/pub/jooslee/prob.TCGA.RData'
CCLE_RPKM_URL = 'https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz'
DEPMAP_CRISPR_DATA = join(DATA_DIR,'DepMap_CRISPR.csv')
DEPMAP_RNAI_DATA = join(DATA_DIR,'DepMap_RNAi.csv')
TCGA_DATA = join(DATA_DIR,'TCGA_data_ISLE.RData')
CCLE_RPKM_DATA_GZ = join(DATA_DIR,'CCLE_RPKM.gct.gz')
CCLE_RPKM_DATA = join(DATA_DIR,'CCLE_RPKM.gct')

rule all: 
    input:
        DEPMAP_CRISPR_DATA,
        DEPMAP_RNAI_DATA,
        TCGA_DATA

rule download_datasets:
    params:
        crispr_url = DEPMAP_CRISPR_URL,
        rnai_url = DEPMAP_RNAI_URL,
        tcga_url = TCGA_URL,
        ccle_rpkm_url = CCLE_RPKM_URL
    output:
        crispr_data = DEPMAP_CRISPR_DATA,
        rnai_data = DEPMAP_RNAI_DATA,
        tcga_data = TCGA_DATA,
        ccle_rpkm_data_gz = CCLE_RPKM_DATA_GZ
        ccle_rpkm_data = CCLE_RPKM_DATA
    shell:
        '''
        wget -O {output.crispr_data} {params.crispr_url}
        wget -O {output.rnai_data} {params.rnai_url} 
        wget -O {output.tcga_data} {params.tcga_url}
        wget -O {output.ccle_rpkm_data_gz} {params.ccle_rpkm_url}
        gzip -d {output.ccle_rpkm_data_gz}
        '''

