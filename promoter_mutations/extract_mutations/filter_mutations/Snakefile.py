import pandas as pd

promoter_suffix = '2000_left'

counts_tsv = '/s/project/mll/sergey/MLL_data/mll5k_counts.tsv.gz'

gnomAD_vcf = '/s/project/mll/sergey/effect_prediction/tools/gnomAD/gnomAD.GRCh38.concat.GRCh37.vcf.gz'

promoters_bed = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/tables/promoters_{promoter_suffix}.bed'

progress_dir = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/{promoter_suffix}/' #output dir

input_data_dir='/s/project/mll/rawdata/'

vcf_list='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_samples/analysed_vcfs.tsv'

all_patients=pd.read_csv(vcf_list, names=['sample_ID', 'VCF_file'], sep=' ')['VCF_file'].str.replace('.vcf.gz','').unique()

rule all:
    input:
        expand(progress_dir + 'gnomAD/{patient}.vcf.gz', patient=all_patients),
        expand(progress_dir + 'filter_gnomAD/{patient}.tsv', patient=all_patients),
        expand(progress_dir + 'filter_counts/{patient}.tsv', patient=all_patients),


rule filter_calls:
    #use only mutations in promoter regions
    input:
        vcf = input_data_dir + 'Somatic_Analysis/{patient}.vcf.gz',
        bed = promoters_bed,
    output:
        vcf = progress_dir + 'passed/{patient}.vcf.gz',
        tbi = progress_dir + 'passed/{patient}.vcf.gz.tbi',
    shell:
        r'''
        bcftools view {input.vcf} --max-alleles 2 -v "snps,indels" -i 'FILTER="PASS"' -R {input.bed} -Oz -o {output.vcf}
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
        -x INFO \
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

rule annotate_genes:
    #add gene annotations using promoter bed file
    #each promoter mutations is annotated with the corresponding gene(s)
    input:
        vcf = progress_dir + 'gnomAD/{patient}.vcf.gz',
        tbi = progress_dir + 'gnomAD/{patient}.vcf.gz.tbi',
        header = 'headers/gene_header.txt',
        bed = promoters_bed,
    output:
        vcf = progress_dir + 'genes/{patient}.vcf.gz',
        tbi = progress_dir + 'genes/{patient}.vcf.gz.tbi',
    shell:
        r'''
        bcftools annotate --threads 4 \
        -h {input.header} \
        -c 'CHROM,FROM,TO,=gene' \
        -a {input.bed} \
        {input.vcf} \
        -Oz -o {output.vcf}
        tabix -f {output.vcf}
        '''

rule filter_gnomAD:
    input:
        vcf = progress_dir + 'genes/{patient}.vcf.gz',
        tbi = progress_dir + 'genes/{patient}.vcf.gz.tbi',
    output:
        tsv = progress_dir + 'filter_gnomAD/{patient}.tsv',
    shell:
        r'''
        bcftools view -H -i '(gnomAD_AF<=5e-4||gnomAD_AF=".")' {input.vcf} \
        |cut -f1,2,3,4,5,8|sed 's/\t[^\t]*gene=/\t/'|sed 's/$/\t'"{wildcards.patient}"'.vcf.gz/' > {output.tsv}
        '''

rule filter_counts:
    input:
        vcf = progress_dir + 'genes/{patient}.vcf.gz',
        tbi = progress_dir + 'genes/{patient}.vcf.gz.tbi',
    output:
        tsv = progress_dir + 'filter_counts/{patient}.tsv',
    shell:
        r'''
        bcftools view -H -i '(gnomAD_AF<=5e-4||gnomAD_AF=".")&&(MLL_counts[0]<18||MLL_counts[1]<27||MLL_counts=".")' {input.vcf} \
        |cut -f1,2,3,4,5,8|sed 's/\t[^\t]*gene=/\t/'|sed 's/$/\t'"{wildcards.patient}"'.vcf.gz/' > {output.tsv}
        '''
