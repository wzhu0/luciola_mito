#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=2
#SBATCH --mem=32gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/revbayes_sum_trees.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito/05_revbayes/output
#SBATCH -J mito_revbayes_sum_tree

# RevBayes tree posterior summary
# Usage: sbatch 05_revbayes/output/sum_trees.Rev


module load gnu/12
module load prebin/kry
module load openmpi
module load boost

rb-dev /home/wzhu/luciola/mito/05_revbayes/scripts/sum_trees.Rev --args ${1}
