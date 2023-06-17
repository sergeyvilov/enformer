import numpy as np
import os
import pandas as pd


progress_dir = '/s/project/mll/sergey/effect_prediction/outrider/aml_association/' #output dir

GRCh37_fa = '/s/project/mll/sergey/MLL_data/GRCh37.fa'

transcripts_canonocal_bed = '/s/project/mll/sergey/effect_prediction/outrider/gene_coords/transcripts_canonocal_GRCh37.bed'

ass_gene_list='/s/project/mll/sergey/effect_prediction/outrider/aml_association/ass_250.csv'

associated_genes=pd.read_csv(ass_gene_list).ass_gene.values


rule all:
    input:
        progress_dir + 'ass_genes.fa',


rule intersect_associated_genes:
    input:
        bed = transcripts_canonocal_bed
    output:
        bed = temp(progress_dir + 'associated_genes_GRCh37.tss.bed')
    run:
        df = pd.read_csv(input.bed, sep="\t", names=['chrom','tss_start','tss_end','geneID','score','strand'])
        df = df[df.geneID.isin(associated_genes)]
        df.to_csv(output.bed, sep="\t", header=None, index=None)

rule get_window_pos:
    input:
        bed = progress_dir + 'associated_genes_GRCh37.tss.bed',
        chrom_sizes = 'GRCh37.chrom_sizes',
    output:
        bed = progress_dir + 'associated_genes_GRCh37.bed',
        tmp = temp(progress_dir + 'associated_genes_GRCh37.bed.tmp')
    params:
        half_window = 196000
    shell:
        r'''
        bedtools flank -i {input.bed} -g {input.chrom_sizes} -s -b {params.half_window} > {output.tmp}
        cat {output.tmp}|awk '{{printf "%s%s",$0,(NR%2?"\t":RS)}}'|awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$9,$4,$5,$6}}'  > {output.bed}
        '''

rule extract_seq:
    input:
        bed = progress_dir + 'associated_genes_GRCh37.bed',
    output:
        fa = progress_dir + 'ass_genes.fa',
    shell:
        r'''
        bedtools getfasta -fi {GRCh37_fa} -bed {input.bed} -fo {output.fa} -nameOnly
        '''
