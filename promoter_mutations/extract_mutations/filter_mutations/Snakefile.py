import pandas as pd

promoter_suffix = '2000_symm'

promoters_bed = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/tables/promoters_{promoter_suffix}.bed'

progress_dir = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/{promoter_suffix}/MLL_filtering/' #output dir

input_data_dir='/s/project/mll/preprocess_202302/Somatic_Analysis_abSplice_gz/'

vcf_list='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_samples/analysed_vcfs.tsv'

all_patients=pd.read_csv(vcf_list, names=['sample_ID', 'VCF_file'], sep='\t')['VCF_file'].str.replace('.vcf.gz','').unique()

rule all:
    input:
        expand(progress_dir + 'genes/{patient}.tsv', patient=all_patients),


rule filter_calls:
    #use only mutations in promoter regions
    input:
        vcf = input_data_dir + '{patient}.abSplice.vcf.gz',
        bed = promoters_bed,
    output:
        vcf = progress_dir + 'passed/{patient}.vcf.gz',
        tbi = progress_dir + 'passed/{patient}.vcf.gz.tbi',
    shell:
        r'''
        bcftools view {input.vcf} -R {input.bed} -Oz -o {output.vcf}
        tabix -f {output.vcf}
        '''
        
rule annotate_genes:
    #add gene annotations using promoter bed file
    #each promoter mutations is annotated with the corresponding gene(s)
    input:
        vcf = progress_dir + 'passed/{patient}.vcf.gz',
        tbi = progress_dir + 'passed/{patient}.vcf.gz.tbi',
        header = 'headers/gene_header.txt',
        bed = promoters_bed,
    output:
        tsv = progress_dir + 'genes/{patient}.tsv',
    shell:
        r'''
        bcftools annotate --threads 4 \
        -h {input.header} \
        -c 'CHROM,FROM,TO,=gene' \
        -a {input.bed} \
        {input.vcf} \
        |grep -v '#'|cut -f1,2,3,4,5,8|sed 's/\t[^\t]*gene=/\t/'|sed 's/$/\t'"{wildcards.patient}"'.vcf.gz/' > {output.tsv}
        '''
