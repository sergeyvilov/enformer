import pandas as pd
import numpy as np


counts_tsv = '/s/project/mll/sergey/MLL_data/mll5k_counts.tsv.gz'

gnomAD_vcf = '/s/project/mll/sergey/effect_prediction/tools/gnomAD/gnomAD.GRCh38.concat.GRCh37.vcf.gz'

progress_dir = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/coding_mutations/' #output dir

input_data_dir='/s/project/mll/rawdata/'

vcf_list='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_samples/analysed_vcfs.tsv'

genes_csv = '/s/project/mll/sergey/effect_prediction/tools/Ensembl/GRCh37/gene_transcript_pairs.tsv'

all_patients=pd.read_csv(vcf_list, names=['sample_ID', 'VCF_file'], sep=' ')['VCF_file'].str.replace('.vcf.gz','').unique()

rule all:
    input:
        expand(progress_dir + 'filter_counts/{patient}.format.tsv', patient=all_patients),


rule filter_calls:
    #take stop_gain, frameshift_variant and NMD_transcript_variant
    input:
        vcf = input_data_dir + 'Somatic_Analysis/{patient}.vcf.gz',
    output:
        vcf = progress_dir + 'passed/{patient}.vcf.gz',
        tbi = progress_dir + 'passed/{patient}.vcf.gz.tbi',
    shell:
        r'''
        bcftools view {input.vcf} --max-alleles 2 -v "snps,indels" -i 'FILTER="PASS" && (CSQT~"NMD_transcript_variant" || CSQT~"stop_gained" || CSQT~"frameshift_variant")'   -Oz -o {output.vcf}
        tabix -f {output.vcf}
        '''

rule annotate_counts:
    input:
        vcf = progress_dir + 'passed/{patient}.vcf.gz',
        tbi = progress_dir + 'passed/{patient}.vcf.gz.tbi',
        header = 'headers/MLL_counts_header.txt',
    output:
        vcf = progress_dir + 'MLL_counts/{patient}.vcf.gz',
        tbi = progress_dir + 'MLL_counts/{patient}.vcf.gz.tbi',
    params:
        counts_tsv = counts_tsv,
    shell:
        r'''
        bcftools annotate --threads 4 \
        -c 'CHROM,POS,REF,ALT,MLL_counts' \
        -x ^INFO/CSQT \
        -h {input.header}  \
        -a {params.counts_tsv} \
        {input.vcf} \
        -Oz -o {output.vcf}
        tabix -f {output.vcf}
        '''

rule annotate_gnomAD:
    #add allele frequency from gnomAD
    #needs a lightweight gnomAD version, where all variants are in a signle VCF file,
    #the INFO field of this file has only gnomAD AF, variants not passing gnomAD filters are removed
    input:
        vcf = progress_dir + 'MLL_counts/{patient}.vcf.gz',
        tbi = progress_dir + 'MLL_counts/{patient}.vcf.gz.tbi',
        header = 'headers/gnomAD_header.txt',
    output:
        vcf = progress_dir + 'gnomAD/{patient}.vcf.gz',
        tbi = progress_dir + 'gnomAD/{patient}.vcf.gz.tbi',
    params:
        gnomAD_vcf = gnomAD_vcf,
    shell:
        r'''
        bcftools annotate --threads 4 \
        -h {input.header} \
        -c 'ID,INFO/gnomAD_AF:=INFO/AF' \
        -a {params.gnomAD_vcf} \
        {input.vcf} \
        -Oz -o {output.vcf}
        tabix -f {output.vcf}
        '''

rule filter_counts:
    input:
        vcf = progress_dir + 'gnomAD/{patient}.vcf.gz',
        tbi = progress_dir + 'gnomAD/{patient}.vcf.gz.tbi',
    output:
        tsv = progress_dir + 'filter_counts/{patient}.tsv',
    shell:
        r'''
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQT\n' -i '(gnomAD_AF<=5e-4||gnomAD_AF=".")&&(MLL_counts[0]<18||MLL_counts[1]<27||MLL_counts=".")' {input.vcf} \
        |sed 's/$/\t'"{wildcards.patient}"'.vcf.gz/' > {output.tsv}
        '''

rule annotate_genes:
    input:
        tsv = progress_dir + 'filter_counts/{patient}.tsv',
        genes_csv = genes_csv
    output:
        tsv = progress_dir + 'filter_counts/{patient}.format.tsv',
    run:
        genes_df = pd.read_csv(input.genes_csv, usecols=[0,1], names=['gene_ENS', 'transcript_ENS'], sep='\t').set_index('transcript_ENS')
        with open(input.tsv, 'r') as inp:
            with open(output.tsv, 'w') as outp:
                for line in inp:
                    chrom,pos,id,ref,alt,csqt,vcf_name = line.split('\t')
                    for group in csqt.split(','):
                        _, gene_HGNC, transcript_ENS, mutations = group.split('|')
                        mutations = [mutations.find(mut_type)>-1 for  mut_type in ['frameshift_variant', 'stop_gained', 'NMD_transcript_variant']]
                        mutations = np.array(mutations)
                        if any(mutations) and transcript_ENS in genes_df.index:
                            gene_ENS = genes_df.loc[transcript_ENS].gene_ENS
                            mutations = '\t'.join(mutations.astype(int).astype(str))
                            outp.write(f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{mutations}\t{gene_ENS}\t{vcf_name}\n")
