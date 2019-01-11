from os.path import join

# Directories and file names
SRC_DIR = 'src'
DATA_DIR = 'data'
OUTPUT_DIR = 'output'
DEPMAP_CRISPR_DATA = join(DATA_DIR,'DepMap_CRISPR.csv')
DEPMAP_RNAI_DATA = join(DATA_DIR,'DepMap_RNAi.csv')
TCGA_DATA = join(DATA_DIR,'TCGA_data_ISLE.RData')
CCLE_RPKM_DATA_GZ = join(DATA_DIR,'CCLE_RPKM.gct.gz')
CCLE_RPKM_DATA = join(DATA_DIR,'CCLE_RPKM.gct')
DEPMAP_CELL_LINE_METADATA = join(DATA_DIR,'DepMap_cell_line_metadata.csv')

#Scripts
GOI_RPKM_SCRIPT = join(SRC_DIR,'sort_cell_lines_by_GOI_RPKM.py')
SUBSET_RPKM_SCRIPT = join(SRC_DIR,'subset_GOI_rpkm_depmap.py')
PAN_CANCER_UNBIASED_ANALYSIS = join(SRC_DIR,'pan_cancer_unbiased_gene_dependency_analysis.py.py')
TISSUE_SPECIFIC_UNBIASED_ANALYSIS= join(SRC_DIR,'tissue_specific_unbiased_gene_dependency_analysis.py')

# URL of files needed
DEPMAP_CRISPR_URL = 'https://ndownloader.figshare.com/files/13396070'
DEPMAP_RNAI_URL = 'https://ndownloader.figshare.com/files/11489681'
DEPMAP_CELL_LINE_METADATA_URL = 'https://depmap.org/portal/download/api/download/external?file_name=processed_portal_downloads%2Fdepmap-public-cell-line-metadata-183e.1%2FDepMap-2018q4-celllines.csv'
TCGA_URL = 'ftp://ftp.umiacs.umd.edu/pub/jooslee/prob.TCGA.RData'
CCLE_RPKM_URL = 'https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz'

# Configurations
config['GOI']=config.get('GOI','CENPA') #GOI=gene_of_interest
config['DepMap_data']=config.get('DepMap_data','CRISPR')

# Outputs
GOI_RPKM_SORTED = join(OUTPUT_DIR,'{}_RPKM_sorted.tsv'.format(config.get('GOI')))
GOI_RPKM_SUBSET = join(OUTPUT_DIR,'{}_RPKM_sorted_subset_depmap_{}.tsv'.format(config.get('GOI'),config.get('DepMap_data')))
DEPMAP_SUBSET = join(OUTPUT_DIR,'{}_depmap_{}_subset.tsv'.format(config.get('GOI'),config.get('DepMap_data')))
RPKM_DISTRIBUTION_PLOT = join(OUTPUT_DIR,'comparing_{}_expression_distribution.png'.format(config.get('GOI')))
WILCOXON_RANKSUM_RESULTS = join(OUTPUT_DIR,'wilcoxon_ranksum_results.tsv')
TISSUE_SPECFIC_DIR ='tissue_specific_unbiased_analysis'
TISSUE_SPECFIC_OUTPUT_DIR = join(OUTPUT_DIR,TISSUE_SPECFIC_DIR)

rule all: 
    input:
        DEPMAP_CRISPR_DATA,
        DEPMAP_RNAI_DATA,
        TCGA_DATA,
        CCLE_RPKM_DATA,
        DEPMAP_CELL_LINE_METADATA,
        GOI_RPKM_SORTED,
        GOI_RPKM_SUBSET,
        DEPMAP_SUBSET,
        RPKM_DISTRIBUTION_PLOT,
        WILCOXON_RANKSUM_RESULTS,
        TISSUE_SPECFIC_OUTPUT_DIR

rule tissue_specific_unbiased_gene_dependency_analysis:
    input:
        GOI_RPKM_SUBSET,
        DEPMAP_SUBSET
    output:
        TISSUE_SPECFIC_OUTPUT_DIR
    params:
        GOI = config.get('GOI')
    shell:
        '''
        python {TISSUE_SPECIFIC_UNBIASED_ANALYSIS} \
        -g {params.GOI} \
        -rs {GOI_RPKM_SUBSET} \
        -gds {DEPMAP_SUBSET} \
        -od {TISSUE_SPECFIC_OUTPUT_DIR}


rule pan_cancer_unbiased_gene_dependency_analysis:
    input:
        GOI_RPKM_SUBSET,
        DEPMAP_SUBSET
    output:
        WILCOXON_RANKSUM_RESULTS
    params:
        GOI = config.get('GOI')
    shell:
        '''
        python {PAN_CANCER_UNBIASED_ANALYSIS} \
        -g {params.GOI} \
        -rs {GOI_RPKM_SUBSET} \
        -gds {DEPMAP_SUBSET} \
        -o {WILCOXON_RANKSUM_RESULTS}
        '''

rule subset_GOI_rpkm_depmap:
    input:
        DEPMAP_CRISPR_DATA,
        DEPMAP_CELL_LINE_METADATA, 
        GOI_RPKM_SORTED
    output:
        GOI_RPKM_SUBSET,
        DEPMAP_SUBSET,
        RPKM_DISTRIBUTION_PLOT
    params:
        GOI = config.get('GOI')
    shell:
        '''
        python {SUBSET_RPKM_SCRIPT} \
        -g {params.GOI} \
        -r {GOI_RPKM_SORTED} \
        -gd {DEPMAP_CRISPR_DATA} \
        -md {DEPMAP_CELL_LINE_METADATA} \
        -or {GOI_RPKM_SUBSET} \
        -odm {DEPMAP_SUBSET} \
        -p {RPKM_DISTRIBUTION_PLOT}
        '''
  

rule gen_cell_line_list_by_GOI_RPKM:
    input:
        CCLE_RPKM_DATA
    params:
        goi = config.get('GOI')
    output:
        GOI_RPKM_SORTED
    shell:
        'python {GOI_RPKM_SCRIPT} -g {params.goi} -e {input} -o {output}'

rule download_cell_line_metadata:
    params:
        url = DEPMAP_CELL_LINE_METADATA_URL
    output:
        DEPMAP_CELL_LINE_METADATA
    shell:
        'wget -O {output} {params.url}'


rule download_tcga_isle:
    params:
        url = TCGA_URL
    output: 
        TCGA_DATA
    shell:
        'wget -O {output} {params.url}'

rule download_depmap:
    params:
        crispr_url = DEPMAP_CRISPR_URL,
        rnai_url = DEPMAP_RNAI_URL
    output:
        crispr_data = DEPMAP_CRISPR_DATA,
        rnai_data = DEPMAP_RNAI_DATA
    shell:
        '''
        wget -O {output.crispr_data} {params.crispr_url}
        wget -O {output.rnai_data} {params.rnai_url}
        '''

rule download_ccle_rpkm:
    params:
        url = CCLE_RPKM_URL
    output:
        CCLE_RPKM_DATA
    shell:
        '''
        wget -O {CCLE_RPKM_DATA_GZ} {params.url}
        gzip -d {CCLE_RPKM_DATA_GZ}
        '''


