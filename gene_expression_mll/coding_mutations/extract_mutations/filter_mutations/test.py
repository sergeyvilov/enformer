import pandas as pd
import numpy as np

genes_csv = '/s/project/mll/sergey/effect_prediction/promoter_mutations/biomart/all_transcripts_genes.tsv.gz'
genes_csv = '/s/project/mll/sergey/effect_prediction/tools/Ensembl/GRCh37/gene_transcript_pairs.tsv'

input = '/s/project/mll/sergey/effect_prediction/promoter_mutations/coding_mutations/filter_counts/comb-MaleCon_comb-MLL_19699_G1_P1.somatic.tsv'
output = 'test.tsv'

#genes_df = pd.read_csv(genes_csv, usecols=[0,1], skiprows=1, names=['gene_ENS', 'transcript_ENS'], sep='\t').set_index('transcript_ENS')
genes_df = pd.read_csv(genes_csv, usecols=[0,1], names=['gene_ENS', 'transcript_ENS'], sep='\t').set_index('transcript_ENS')
with open(input, 'r') as inp:
    with open(output, 'w') as outp:
        for line in inp:
            chrom,pos,id,ref,alt,csqt,vcf_name = line.split('\t')
            for group in csqt.split(','):
                _, gene_HGNC, transcript_ENS, mutations = group.split('|')
                mutations = [mutations.find(mut_type)>-1 for  mut_type in ['frameshift_variant', 'stop_gained', 'NMD_transcript_variant']]
                mutations = np.array(mutations)
                if not transcript_ENS in genes_df.index:
                    print(transcript_ENS)
                if any(mutations) and transcript_ENS in genes_df.index:
                    gene_ENS = genes_df.loc[transcript_ENS].gene_ENS
                    mutations = '\t'.join(mutations.astype(int).astype(str))
                    outp.write(f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{mutations}\t{gene_ENS}\t{vcf_name}\n")
