#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=2
#SBATCH --mem=64gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/sum_rb_trees.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/05_revbayes/output
#SBATCH -J sum_rb_trees

# RevBayes tree posterior summary
# Usage: sbatch sbatch 05_revbayes/scripts/sum_trees.slurm Luciola_mito_unrooted_non-clock_MCMC150000


module load gnu/12
module load prebin/kry
module load openmpi
module load boost

rb-dev /home/wzhu/luciola/mito/05_revbayes/scripts/sum_trees.Rev --args ${1}
