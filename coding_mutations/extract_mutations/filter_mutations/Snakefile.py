import pandas as pd
import numpy as np

progress_dir = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/coding_mutations/MLL_filtering/' #output dir

input_data_dir='/s/project/mll/preprocess_202302/Somatic_Analysis_abSplice_gz/'

vcf_list='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_samples/analysed_vcfs.tsv'

genes_csv = '/s/project/mll/sergey/effect_prediction/tools/Ensembl/GRCh37/gene_transcript_pairs.tsv'

all_patients=pd.read_csv(vcf_list, names=['sample_ID', 'VCF_file'], sep='\t')['VCF_file'].str.replace('.vcf.gz','').unique()

rule all:
    input:
        expand(progress_dir + 'filter_counts/{patient}.format.tsv', patient=all_patients),


rule filter_calls:
    #take stop_gain, frameshift_variant and NMD_transcript_variant
    input:
        vcf = input_data_dir + '{patient}.abSplice.vcf.gz',
    output:
        tsv = progress_dir + 'filter_counts/{patient}.tsv',
    shell:
        r'''
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQT\n' \
        -i 'CSQT~"NMD_transcript_variant" || CSQT~"stop_gained" || CSQT~"frameshift_variant" || CSQT~"splice_"' \
        {input.vcf} |sed 's/$/\t'"{wildcards.patient}"'.vcf.gz/' > {output.tsv}
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
                        mutations = [mutations.find(mut_type)>-1 for  mut_type in ['frameshift_variant', 'stop_gained', 'NMD_transcript_variant', 'splice_donor_variant', 'splice_acceptor_variant','splice_region_variant']]
                        mutations = np.array(mutations)
                        if any(mutations) and transcript_ENS in genes_df.index:
                            gene_ENS = genes_df.loc[transcript_ENS].gene_ENS
                            mutations = '\t'.join(mutations.astype(int).astype(str))
                            outp.write(f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{mutations}\t{gene_ENS}\t{vcf_name}\n")
