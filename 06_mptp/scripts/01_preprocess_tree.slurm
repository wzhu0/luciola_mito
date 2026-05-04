#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/preprocess_mptp.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J preprocess_mptp

module purge
module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate convergence

Rscript 06_mptp/scripts/01_preprocess_tree.R 05_revbayes/output/Luciola_mito_unrooted_non-clock_MCMC150000_MCC.tree
