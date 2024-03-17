#!/bin/bash

#combine mutations in a tsv table after filtering (filter_mutations)


#vcfdir="/s/project/mll/sergey/effect_prediction/promoter_mutations/coding_mutations/filter_counts/"

vcfdir="/s/project/mll/sergey/effect_prediction/promoter_mutations/coding_mutations/MLL_filtering/filter_counts/"

cat $vcfdir/*.format.tsv > $vcfdir/mutations_filter_counts.tsv
