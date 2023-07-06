#!/bin/bash

#get a list of unique genes, used by outrider and activation tools

output_dir='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_genes/'

outrider_dir='/s/project/mll/sergey/effect_prediction/outrider/input_data/outrider'
activation_dir='/s/project/mll/sergey/effect_prediction/outrider/input_data/activation'

output_dir='/s/project/mll/sergey/effect_prediction/promoter_mutations/analysed_genes/'

mkdir -p $output_dir

cd $output_dir

echo 'Collecting OUTRIDER genes'

zless $outrider_dir/outrider_all.csv.gz|cut -d',' -f1|tail -n+2|sort|uniq > outrider_genes.csv

echo 'Collecting Activation genes'

zless $activation_dir/res_filter_out_all.csv.gz|cut -d',' -f1|tail -n+2|sort|uniq > activation_genes.csv

cat outrider_genes.csv activation_genes.csv|sort|uniq > analysed_genes.csv

rm outrider_genes.csv activation_genes.csv


