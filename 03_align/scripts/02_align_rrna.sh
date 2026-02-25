#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=4
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/02_align_rrna.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_mafft

# Step 2b: Align rrnS and rrnL with MAFFT (auto strategy)
# Run after 01_extract.sh is complete
# Submit as: sbatch 03_align/scripts/02_align_rrna.sh

set -euo pipefail

module load gnu/12
module load prebin/kry
module load mafft/7.526

# ── Paths ─────────────────────────────────────────────────────────────────────
EXTRACT_DIR="03_align/01_extracted"
ALIGN_DIR="03_align/02_aligned"

mkdir -p "${ALIGN_DIR}" logs

# ── Align each rRNA gene ──────────────────────────────────────────────────────
for gene in rrnS rrnL; do
    IN="${EXTRACT_DIR}/${gene}.fasta"

    if [[ ! -f "$IN" ]]; then
        echo "ERROR: ${IN} not found. Run 01_extract.sh first."
        exit 1
    fi

    echo "MAFFT: aligning ${gene}"
    mafft --auto \
        --thread 4 \
        "$IN" \
        > "${ALIGN_DIR}/${gene}.fasta" \
        2> "${ALIGN_DIR}/${gene}_mafft.log"

    echo "Done: ${gene} -> ${ALIGN_DIR}/${gene}.fasta"
done
