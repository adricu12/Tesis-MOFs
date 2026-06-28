#!/bin/bash

#SBATCH --job-name=job-id
#SBATCH --nodes=int
#SBATCH --ntasks=int
#SBATCH --cpus-per-task=int
#SBATCH --time=hh:mm:ss
#SBATCH --mem=intG
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source ~/.bash_profile

module purge
module add [all pacakges needed]

conda activate host_guest

complexes_from_file HOST GUEST \
                    -nh 1 -nm 1 -nc 20000 -s SEED \
                    -r job-id_results \
                    --hdf5 \
                    >> job-id_results.out
mv job-id_results/complexes.h5 ./job-id.h5
rm -r job-id_results