#!/bin/bash

#combine mutations in a tsv table after filtering (filter_mutations)

promoter_suffix='2000_symm' 

vcfdir="/s/project/mll/sergey/effect_prediction/promoter_mutations/${promoter_suffix}/filter_gnomAD/"

output_tsv="/s/project/mll/sergey/effect_prediction/promoter_mutations/${promoter_suffix}/mutations_filter_gnomAD.tsv"

cat $vcfdir/* > $output_tsv

vcfdir="/s/project/mll/sergey/effect_prediction/promoter_mutations/${promoter_suffix}/filter_counts/"

output_tsv="/s/project/mll/sergey/effect_prediction/promoter_mutations/${promoter_suffix}/mutations_filter_counts.tsv"

cat $vcfdir/* > $output_tsv