#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=8
#SBATCH --mem=32gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/revbayes_rooted_non-clock.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/05_revbayes
#SBATCH -J mito_revbayes_rooted_non-clock

# Bayesian rooted tree inference with RevBayes
# Usage: sbatch 05_revbayes/scripts/run_revbayes_rooted_non-clock.slurm


module load gnu/12
module load prebin/kry
module load openmpi
module load boost

NUM_RUNS=4
NUM_MCMC_ITER=100000

mpirun -np ${SLURM_NTASKS} rb-mpi-dev scripts/mcmc_rooted_non-clock.Rev --args ${NUM_RUNS} ${NUM_MCMC_ITER}
