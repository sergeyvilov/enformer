import pandas as pd

promoter_suffix = '2000_symm'

counts_tsv = '/lustre/groups/epigenereg01/workspace/projects/vale/calling_new/MLL/reoccurence/mll5k_counts.tsv.gz'

gnomAD_vcf = '/lustre/groups/epigenereg01/workspace/projects/vale/tools/gnomAD/v3.1.1_GRCh38/af_only/GRCh37/gnomAD.GRCh38.concat.GRCh37.vcf.gz'

promoters_bed = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/tables/promoters_{promoter_suffix}.bed'

progress_dir = f'/s/project/mll/sergey/effect_prediction/promoter_mutations/{promoter_suffix}/' #output dir

input_data_dir='/s/project/mll/rawdata/Somatic_Analysis/'

vcf_list='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_samples/analysed_vcfs.tsv'

all_patients=pd.read_csv(vcf_list, names=['sample_ID', 'VCF_file'])['VCF_file'].str.replace('.vcf.gz','').unique()

rule all:
    input:
        expand(progress_dir + 'filtered/{patient}.vcf.gz', patient=all_patients),

rule filter_calls:
    input:
        vcf = input_data_dir + 'Somatic_Analysis/{patient}.vcf.gz',
    output:
        vcf = progress_dir + 'passed/{patient}.vcf.gz',
        tbi = progress_dir + 'passed/{patient}.vcf.gz.tbi',
    shell:
        r'''
        bcftools view {input.vcf} --max-alleles 2 -v "snps,indels" -i 'FILTER="PASS"' -Oz -o {output.vcf}
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
        workdir = progress_dir + 'gnomAD/'
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


rule apply_filteres:
    input:
        vcf = progress_dir + 'gnomAD/{patient}.vcf.gz',
        tbi = progress_dir + 'gnomAD/{patient}.vcf.gz.tbi',
    output:
        vcf = progress_dir + 'filtered/{patient}.vcf.gz',
        tbi = progress_dir + 'filtered/{patient}.vcf.gz.tbi',
    shell:
        r'''
        bcftools view -i '(gnomAD_AF<=5e-4||gnomAD_AF=".")&&(MLL_counts[0]<18||MLL_counts[1]<27||MLL_counts=".")' {input.vcf} -Oz -o {output.vcf}
        tabix -f {output.vcf}
        '''
