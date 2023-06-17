import numpy as np
import os
import pandas as pd


progress_dir = '/lustre/groups/epigenereg01/workspace/projects/vale/outrider/gene_coords/' #output dir

GRCh37_fa = '/lustre/groups/epigenereg01/workspace/projects/vale/calling_new/MLL/resources_GRCh37/GRCh37.fa'

liftover_dir = '/lustre/groups/epigenereg01/workspace/projects/vale/tools/liftOver/'

#outrider_results='/lustre/groups/epigenereg01/workspace/projects/vale/outrider/outrider_210223/tables/OUTRIDER_results_AML_panel.tsv'
#leukaemia_genes=pd.read_csv(outrider_results)['geneID'].str.replace('\..*','', regex=True).unique()

outrider_results='/lustre/groups/epigenereg01/workspace/projects/vale/outrider/outrider_210223/tables/120_random_insignif.csv'
insignificant_genes=pd.read_csv(outrider_results)['geneID'].str.replace('\..*','', regex=True).unique()


rule all:
    input:
        progress_dir + 'promoters_196000_left.bed',
        progress_dir + 'promoters_196000_symm.bed',
        progress_dir + 'insignificant_genes_GRCh37.fa',
        #progress_dir + 'aml_genes_GRCh37.fa',

rule get_tss_bed:
    input:
        csv = progress_dir + 'transcripts_canonocal_GRCh38.csv'
    output:
        bed = progress_dir + 'transcripts_canonocal_GRCh38.bed6'
    shell:
        r'''
        cat {input.csv} |tail -n +2|awk 'BEGIN{{FS=",";OFS="\t"}}{{if ($4==1) {{$4="+"}} else {{$4="-"}};print $8,$3-1,$3,$1,".",$4}}' \
        |sort -V -k1,1 -k2,2|uniq|grep -E '^[0-9XYMT]+\s'| sed -e 's/^/chr/' -e 's/chrMT/chrM/'  > {output.bed}
        '''

rule liftover:
    input:
        bed = progress_dir + 'transcripts_canonocal_GRCh38.bed6',
        chain_file = liftover_dir + 'chain_files/hg38ToHg19.over.chain.gz' #chain file to convert positions see https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
    output:
        bed = progress_dir + 'transcripts_canonocal_GRCh37.liftover.bed',
        umap = temp(progress_dir + 'transcripts_canonocal_GRCh37.umap')
    log:
        progress_dir + 'logs/liftover.log'
    shell:
        r'''
        {liftover_dir}/liftOver {input.bed}  {input.chain_file} {output.bed}  {output.umap}  > {log} 2>&1
        '''

rule change_chrom_names:
    input:
        bed = progress_dir + 'transcripts_canonocal_GRCh37.liftover.bed',
    output:
        bed = progress_dir + 'transcripts_canonocal_GRCh37.bed',
    shell:
        r'''
        cat {input.bed}|sed -e 's/^chr//' -e 's/^M/MT/'|grep -E '^[0-9XYMT]+\s' > {output.bed}
        '''

rule get_promoters_symm:
    input:
        bed = progress_dir + 'transcripts_canonocal_GRCh37.bed',
        chrom_sizes = 'GRCh37.chrom_sizes',
    output:
        bed = progress_dir + 'promoters_196000_symm.bed',
    params:
        workdir = progress_dir
    shell:
        r'''
        for promoter_size in 2000 5000 10000 25000 50000 196000;do
            bedtools flank -i {input.bed} -g {input.chrom_sizes} -s -l ${{promoter_size}} -r ${{promoter_size}} |\
            awk '{{printf "%s%s",$0,(NR%2?"\t":RS)}}'|awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$9,$4,$5,$6}}'  > {params.workdir}promoters_${{promoter_size}}_symm.bed
        done
        '''

rule get_promoters_left:
    input:
        bed = progress_dir + 'transcripts_canonocal_GRCh37.bed',
        chrom_sizes = 'GRCh37.chrom_sizes',
    output:
        bed = progress_dir + 'promoters_196000_left.bed',
    params:
        workdir = progress_dir
    shell:
        r'''
        for promoter_size in 2000 5000 10000 25000 50000 196000;do
            bedtools flank -i {input.bed} -g {input.chrom_sizes} -s -l ${{promoter_size}} -r 0 > {params.workdir}promoters_${{promoter_size}}_left.bed
        done
        '''


rule intersect_insignificant_genes:
    input:
        bed = progress_dir + 'transcripts_canonocal_GRCh37.bed'
    output:
        bed = progress_dir + 'insignificant_genes_GRCh37.tss.bed'
    run:
        df = pd.read_csv(input.bed, sep="\t", names=['chrom','tss_start','tss_end','geneID','score','strand'])
        df = df[df.geneID.isin(insignificant_genes)]
        df.to_csv(output.bed, sep="\t", header=None, index=None)

rule get_window_pos:
    input:
        bed = progress_dir + 'insignificant_genes_GRCh37.tss.bed',
        chrom_sizes = 'GRCh37.chrom_sizes',
    output:
        bed = progress_dir + 'insignificant_genes_GRCh37.bed',
        tmp = temp(progress_dir + 'insignificant_genes_GRCh37.bed.tmp')
    params:
        half_window = 196000
    shell:
        r'''
        bedtools flank -i {input.bed} -g {input.chrom_sizes} -s -b {params.half_window} > {output.tmp}
        cat {output.tmp}|awk '{{printf "%s%s",$0,(NR%2?"\t":RS)}}'|awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$9,$4,$5,$6}}'  > {output.bed}
        '''

rule extract_seq:
    input:
        bed = progress_dir + 'insignificant_genes_GRCh37.bed',
    output:
        fa = progress_dir + 'insignificant_genes_GRCh37.fa',
    shell:
        r'''
        bedtools getfasta -fi {GRCh37_fa} -bed {input.bed} -fo {output.fa} -nameOnly
        '''

# rule intersect_leukaemia_genes:
#     input:
#         bed = progress_dir + 'transcripts_canonocal_GRCh37.bed'
#     output:
#         bed = progress_dir + 'aml_genes_GRCh37.tss.bed'
#     run:
#         df = pd.read_csv(input.bed, sep="\t", names=['chrom','tss_start','tss_end','geneID','score','strand'])
#         df = df[df.geneID.isin(leukaemia_genes)]
#         df.to_csv(output.bed, sep="\t", header=None, index=None)
#
# rule get_window_pos:
#     input:
#         bed = progress_dir + 'aml_genes_GRCh37.tss.bed',
#         chrom_sizes = 'GRCh37.chrom_sizes',
#     output:
#         bed = progress_dir + 'aml_genes_GRCh37.bed',
#         tmp = temp(progress_dir + 'aml_genes_GRCh37.bed.tmp')
#     params:
#         half_window = 196000
#     shell:
#         r'''
#         bedtools flank -i {input.bed} -g {input.chrom_sizes} -s -b {params.half_window} > {output.tmp}
#         cat {output.tmp}|awk '{{printf "%s%s",$0,(NR%2?"\t":RS)}}'|awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$9,$4,$5,$6}}'  > {output.bed}
#         '''
#
# rule extract_seq:
#     input:
#         bed = progress_dir + 'aml_genes_GRCh37.bed',
#     output:
#         fa = progress_dir + 'aml_genes_GRCh37.fa',
#     shell:
#         r'''
#         bedtools getfasta -fi {GRCh37_fa} -bed {input.bed} -fo {output.fa} -nameOnly
#         '''
