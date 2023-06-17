#!/bin/bash


#SBATCH -J extract_variants
#SBATCH -p cpu_p
#SBATCH --nice=10000
#SBATCH -c 2
#SBATCH -t 2-00:00:00
#SBATCH --mem=2G
#SBATCH -o logs/%a.o
#SBATCH -e logs/%a.e



c=0

echo ${SLURM_ARRAY_TASK_ID}

for promoter_length in 2000 5000 10000 25000 50000;do

  for promoter_dir in left symm; do

    for max_gnomAD_AF in 5e-4 1;do

        if [ ${SLURM_ARRAY_TASK_ID} -eq $c ]; then

          ./extract_variants.sh $promoter_length $promoter_dir $max_gnomAD_AF

        fi

        c=$((c+1))

      done
  done
done
