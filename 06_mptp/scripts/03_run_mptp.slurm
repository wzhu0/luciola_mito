#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/mptp_%a.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mptp
#SBATCH --array=1-4

module load gnu/12
module load prebin/kry

DATA_DIR="06_mptp/data"
OUT_DIR="06_mptp/output"
BASE="Luciola_mito_unrooted_non-clock_MCMC150000_MCC_rooted"
MINBR=0.000253

mkdir -p ${OUT_DIR}

case ${SLURM_ARRAY_TASK_ID} in
    1)
        echo "ML - full tree"
        mptp --ml --multi \
            --minbr ${MINBR} \
            --tree_file ${DATA_DIR}/${BASE}.nwk \
            --output_file ${OUT_DIR}/${BASE}_ml
        ;;
    2)
        echo "ML - field-only tree"
        mptp --ml --multi \
            --minbr ${MINBR} \
            --tree_file ${DATA_DIR}/${BASE}_field_only.nwk \
            --output_file ${OUT_DIR}/${BASE}_field_only_ml
        ;;
    3)
        echo "MCMC - full tree"
        mptp \
            --mcmc 50000000 \
            --multi \
            --mcmc_sample 50000 \
            --mcmc_burnin 10000000 \
            --mcmc_runs 10 \
            --mcmc_startml \
            --minbr ${MINBR} \
            --tree_file ${DATA_DIR}/${BASE}.nwk \
            --output_file ${OUT_DIR}/${BASE}_mcmc
        ;;
    4)
        echo "MCMC - field-only tree"
        mptp \
            --mcmc 50000000 \
            --multi \
            --mcmc_sample 50000 \
            --mcmc_burnin 10000000 \
            --mcmc_runs 10 \
            --mcmc_startml \
            --minbr ${MINBR} \
            --tree_file ${DATA_DIR}/${BASE}_field_only.nwk \
            --output_file ${OUT_DIR}/${BASE}_field_only_mcmc
        ;;
esac

echo "Done: task ${SLURM_ARRAY_TASK_ID}"