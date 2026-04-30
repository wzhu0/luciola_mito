#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/convergence.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J convergence

module purge
module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate convergence

PREFIX=$1

Rscript 05_revbayes/scripts/assess_convergence.R ${PREFIX}
