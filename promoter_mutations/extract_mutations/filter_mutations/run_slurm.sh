#!/bin/bash
#exclude servers with no AVX support (otherwise calling takes ages)!
#mutect2 log should have the line "Using CPU-supported AVX-512 instructions"
#if the line is absent, the node should also be excluded as AVX works somehow more slowly on these machines

source /home/icb/sergey.vilov/.bashrc
conda activate vale-bio
rm -f slurm_logs/*
mkdir slurm_logs
snakemake -s Snakefile.py -k --restart-times 3 --rerun-incomplete --use-conda --latency-wait 180 --cluster-config cluster.yaml --cluster 'sbatch -p cpu_p \
--mem={cluster.mem} --time=2-00:00:00 --threads-per-core={cluster.threads_per_core} --nice=10000 -c {cluster.cores} -o slurm_logs/%j.out' -j 20
#-x ibis216-010-0[68-71] \
