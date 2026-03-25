#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=2
#SBATCH --mem=2gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/revbayes_timetree.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/05_revbayes
#SBATCH -J mito_revbayes_timetree

# Bayesian time tree inference with RevBayes
# Usage: sbatch 05_revbayes/scripts/run_revbayes_timetree.slurm


module load gnu/12
module load prebin/kry
module load openmpi
module load boost

NUM_RUNS=2
NUM_MCMC_ITER=100000

mpirun -np ${SLURM_NTASKS} rb-mpi-dev scripts/mcmc_timetree_mito.Rev ${NUM_RUNS} ${NUM_MCMC_ITER}
