#!/bin/bash

#SBATCH -J enformer
#SBATCH -p gpu_p
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --nice=10000
#SBATCH -c 2
#SBATCH -t 2-00:00:00
#SBATCH --mem=16G
#SBATCH -o logs/%a.o
#SBATCH -e logs/%a.e

source ~/.bashrc; conda activate enformer

c=0

for promoter_length in 2000 5000 10000 25000 50000;do

  for promoter_dir in left symm; do

    for max_gnomAD_AF in 5e-4 1;do

        if [ ${SLURM_ARRAY_TASK_ID} -eq $c ]; then

          python -u run_enformer.py $promoter_length $promoter_dir $max_gnomAD_AF

        fi

        c=$((c+1))

      done
  done
done
