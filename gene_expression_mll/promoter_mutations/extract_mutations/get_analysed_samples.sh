#!/bin/bash

#get a list of unique samples (MLL_ID), used by outrider and activation tools

output_dir='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_samples/'

outrider_dir='/s/project/mll/sergey/effect_prediction/outrider/input_data/outrider'
activation_dir='/s/project/mll/sergey/effect_prediction/outrider/input_data/activation'

vcf_matching_jessy='/s/project/mll/sergey/effect_prediction/outrider/input_data/sample_annotation_for_drop.tsv'

mkdir -p $output_dir

cd $output_dir

echo 'Collecting OUTRIDER samples'

zless $outrider_dir/outrider_all.csv.gz|cut -d',' -f2|tail -n+2|sort|uniq > outrider_samples.csv

echo 'Collecting Activation samples'

zless $activation_dir/res_filter_out_all.csv.gz|cut -d',' -f2|tail -n+2|sort|uniq > activation_samples.csv

cat outrider_samples.csv activation_samples.csv|sort|uniq > analysed_samples.csv

rm outrider_samples.csv activation_samples.csv

cat $vcf_matching_jessy |tail -n+2|cut -f1,25|sed 's/"//g'|sort -k1,1 > vcf_matching.tsv

join analysed_samples.csv vcf_matching.tsv > analysed_vcfs.tsv

rm analysed_samples.csv




