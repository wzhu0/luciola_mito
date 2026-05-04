#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/soft_convergence_%j.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J soft_convergence

# Assess MCMC convergence on a subset of ongoing RevBayes output.
#
# Useful when the required chain length for convergence is uncertain:
# run the MCMC for longer than needed, then call this script periodically
# with increasing NUM_MCMC_ITER to check convergence without interrupting
# the chain.
#
# Usage: sbatch soft_convergence.slurm <NUM_MCMC_ITER> <PREFIX>
#

module purge
module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate convergence

NUM_MCMC_ITER=$1
PREFIX=$2
OUTDIR="/home/wzhu/luciola/mito/05_revbayes/output/"
N_LINES=$(( NUM_MCMC_ITER / 10 + 2 ))

echo "${NUM_MCMC_ITER} MCMC iterations"

for run in 1 2 3 4; do
    for ext in log trees; do
        infile="${OUTDIR}${PREFIX}_run_${run}.${ext}"
        outfile="${OUTDIR}${PREFIX}_MCMC${NUM_MCMC_ITER}_run_${run}.${ext}"
        head -n "${N_LINES}" "${infile}" > "${outfile}"
    done
done

Rscript 05_revbayes/scripts/assess_convergence.R ${PREFIX}_MCMC${NUM_MCMC_ITER}

