#!/bin/bash

promoter_length=$1
promoter_dir=$2
max_gnomAD_AF=$3

bed_file='/lustre/groups/epigenereg01/workspace/projects/vale/outrider/gene_coords/insignificant_genes_GRCh37.bed'
promoters="/lustre/groups/epigenereg01/workspace/projects/vale/outrider/gene_coords/promoters_${promoter_length}_${promoter_dir}.bed"
vcf_dir='/lustre/groups/epigenereg01/workspace/projects/vale/outrider/Strelka_variants/gnomAD'
output_dir="/lustre/groups/epigenereg01/workspace/projects/vale/outrider/promoter_variants/insignificant_genes/max_gnomAD_${max_gnomAD_AF}"

mkdir -p ${output_dir}

output_tsv="${output_dir}/promoters_${promoter_length}_${promoter_dir}.tsv"

tot_genes=$(wc -l $bed_file|cut -d" " -f1)
tot_vcfs=$(ls $vcf_dir/*.vcf.gz|wc -l)

n_gene=1

> $output_tsv

while  read -r _ region_start region_stop gene_name _ strand;do
  read -r chrom promoter_start promoter_stop _ <<< $(cat $promoters|grep $gene_name)
  echo "Extracting promoter variants for $gene_name::$chrom:$promoter_start-$promoter_stop ($n_gene/$tot_genes)"
  #n_vcf=1
  for vcf in $(ls $vcf_dir/*.vcf.gz);do
    vcf_name=$(basename $vcf)
    #echo "Processing $vcf_name ($n_vcf/$tot_vcfs)"
    bcftools query -i '(gnomAD_AF="."||gnomAD_AF<='$max_gnomAD_AF')&&(MLL_counts[0]<18||MLL_counts[1]<27||MLL_counts=".")' $vcf  -r "$chrom:$promoter_start-$promoter_stop" -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO\n"|sed "s/$/\t$vcf_name\t$gene_name\t$region_start\t$region_stop/" >> $output_tsv
  done
  n_gene=$((n_gene+1))
done  <$bed_file
