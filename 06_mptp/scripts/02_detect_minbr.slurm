#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --output=/home/wzhu/luciola/mito/06_mptp/output/mptp_minbr.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mptp_minbr

module load gnu/12
module load prebin/kry

DATA_DIR="06_mptp/data"
OUT_DIR="06_mptp/output"
BASE="Luciola_mito_unrooted_non-clock_MCMC150000_MCC_rooted"

mkdir -p ${OUT_DIR}

echo "Minimum branch length for full tree..."
mptp --minbr_auto ${DATA_DIR}/concat_full.fasta \
    --tree_file ${DATA_DIR}/${BASE}.nwk \
    --output_file ${OUT_DIR}/minbr_full

echo ""
echo "Minimum branch length for field-only tree..."
mptp --minbr_auto ${DATA_DIR}/concat_field_only.fasta \
    --tree_file ${DATA_DIR}/${BASE}_field_only.nwk \
    --output_file ${OUT_DIR}/minbr_field_only

echo ""
echo "Done. Check the output above for --minbr values to use in 03_run_mptp.slurm"