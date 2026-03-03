#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/revbayes_timetree.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/04_revbayes
#SBATCH -J mito_revbayes_timetree

# Bayesian time tree inference with RevBayes
# GTR+Gamma+Inv, 28 partitions, UCLN relaxed clock, birth-death tree prior
# Submit as: sbatch 04_revbayes/scripts/run_revbayes_timetree.slurm


module load gnu/12
module load prebin/kry

mkdir -p output logs

echo "Starting RevBayes at $(date)"
echo "Working directory: $(pwd)"
echo "Nodes: $SLURM_NNODES, CPUs: $SLURM_CPUS_PER_TASK"

rb-dev scripts/mcmc_combined.Rev

echo "RevBayes finished at $(date)"
