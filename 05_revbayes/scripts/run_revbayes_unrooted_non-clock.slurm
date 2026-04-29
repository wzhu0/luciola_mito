#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=8
#SBATCH --mem=32gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/revbayes_unrooted_non-clock.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/05_revbayes
#SBATCH -J mito_revbayes_unrooted_non-clock

# Bayesian unrooted tree inference with RevBayes
# Usage: sbatch 05_revbayes/scripts/run_revbayes_unrooted_non-clock.slurm


module load gnu/12
module load prebin/kry
module load openmpi
module load boost

mpirun -np ${SLURM_NTASKS} rb-mpi-dev scripts/mcmc_unrooted_non-clock.Rev
