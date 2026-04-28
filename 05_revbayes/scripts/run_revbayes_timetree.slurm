#!/bin/bash
#SBATCH --partition=lemmium
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=8
#SBATCH --mem=32gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/revbayes_timetree.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/05_revbayes
#SBATCH -J mito_revbayes_timetree

# Bayesian time tree inference with RevBayes
# Usage: sbatch 05_revbayes/scripts/run_revbayes_timetree.slurm


module load gnu/12
module load prebin/lem
module load openmpi
module load boost

mpirun -np ${SLURM_NTASKS} rb-mpi-dev scripts/mcmc_timetree_mito.Rev
