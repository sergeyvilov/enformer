import numpy as np
import os
import pandas as pd

#liftover gene coordinates from GRCh38 to GRCh37

gene_list = '/s/project/mll/sergey/effect_prediction/promoter_mutations/biomart/genes_GRCh38.tsv'

progress_dir = '/s/project/mll/sergey/effect_prediction/promoter_mutations/genes_GRCh37/' #output dir

liftover_dir = '/s/project/mll/sergey/effect_prediction/tools/liftOver/'

rule all:
    input:
        progress_dir + 'genes_GRCh37.bed',


rule get_tss_bed:
    input:
        csv = gene_list
    output:
        bed = progress_dir + 'genes_GRCh38.bed6'
    shell:
        r'''
        cat {input.csv} |tail -n +2|awk 'BEGIN{{OFS="\t"}} {{print $5,$2-1,$3,$1}}' \
        |sort -V -k1,1 -k2,2|uniq|grep -E '^[0-9XYMT]+\s'| sed -e 's/^/chr/' -e 's/chrMT/chrM/'  > {output.bed}
        '''

rule liftover:
    input:
        bed = progress_dir + 'genes_GRCh38.bed6',
        chain_file = liftover_dir + 'hg38ToHg19.over.chain.gz' #chain file to convert positions see https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
    output:
        bed = progress_dir + 'genes_GRCh37.liftover.bed',
        umap = progress_dir + 'genes_GRCh37.umap'
    log:
        progress_dir + 'logs/liftover.log'
    shell:
        r'''
        {liftover_dir}/liftOver {input.bed}  {input.chain_file} {output.bed}  {output.umap}  > {log} 2>&1
        '''

rule change_chrom_names:
    input:
        bed = progress_dir + 'genes_GRCh37.liftover.bed',
    output:
        bed = progress_dir + 'genes_GRCh37.bed',
    shell:
        r'''
        cat {input.bed}|sed -e 's/^chr//' -e 's/^M/MT/'|grep -E '^[0-9XYMT]+\s' > {output.bed}
        '''
