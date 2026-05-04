#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=16
#SBATCH --mem=64gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/revbayes_rooted_non-clock.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/05_revbayes
#SBATCH -J rb_rooted_non-clock

# Bayesian rooted tree inference with RevBayes
# Usage: sbatch 05_revbayes/scripts/run_revbayes_rooted_non-clock.slurm


module load gnu/12
module load prebin/kry
module load openmpi
module load boost

NUM_RUNS=4
NUM_MCMC_ITER=1000000

mpirun -np ${SLURM_NTASKS} rb-mpi-dev scripts/mcmc_rooted_non-clock.Rev --args ${NUM_RUNS} ${NUM_MCMC_ITER}
