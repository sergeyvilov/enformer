#!/bin/bash

#combine mutations in a tsv table after filtering (filter_mutations)

promoter_suffix='2000_symm'

vcfdir="/s/project/mll/sergey/effect_prediction/promoter_mutations/${promoter_suffix}/filtered/"

output_tsv="/s/project/mll/sergey/effect_prediction/promoter_mutations/${promoter_suffix}/mutations.tsv"

for vcf in $(find "$vcfdir" -name '*.vcf.gz');do
    vcf_name=$(basename $vcf)
    bcftools view -H $vcf|cut -f1,2,3,4,5,8|sed 's/\t[^\t]*gene=/\t/'|sed 's/$/\t'"$vcf_name"'/'
done > $output_tsv