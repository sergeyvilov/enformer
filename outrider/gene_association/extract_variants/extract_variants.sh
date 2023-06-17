#!/bin/bash

promoter_length=5000
promoter_dir=symm
max_gnomAD_AF=5e-4

gene_list='/data/ouga/home/ag_gagneur/l_vilov/workspace/enformer/outrider/gene_association/AML_associated.tsv'

promoters="/s/project/mll/sergey/effect_prediction/outrider/gene_coords/promoters_${promoter_length}_${promoter_dir}.bed"

vcf_dir='/s/project/mll/sergey/effect_prediction/outrider/calling/gnomAD'

output_dir="/s/project/mll/sergey/effect_prediction/outrider/aml_association/promoter_variants/max_gnomAD_${max_gnomAD_AF}"

mkdir -p ${output_dir}

output_tsv="${output_dir}/promoters_${promoter_length}_${promoter_dir}.tsv"

tot_genes=$(wc -l $gene_list|cut -d" " -f1)

n_gene=1

> $output_tsv

while  read -r sample_name gene_name vcf_name _ ;do

  read -r chrom promoter_start promoter_stop _ <<< $(cat $promoters|grep $gene_name)
  
  echo "Extracting promoter variants for $vcf_name:$gene_name ($n_gene/$tot_genes)"
  
    vcf=$vcf_dir/$vcf_name
    
    #echo "Processing $vcf_name ($n_vcf/$tot_vcfs)"
    bcftools query -i '(gnomAD_AF="."||gnomAD_AF<='$max_gnomAD_AF')&&(MLL_counts[0]<18||MLL_counts[1]<27||MLL_counts=".")' $vcf  -r "$chrom:$promoter_start-$promoter_stop" -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO\n"|sed "s/$/\t$sample_name\t$gene_name/" >> $output_tsv
  
  n_gene=$((n_gene+1))
  
done  < <(tail -n+2 $gene_list)
