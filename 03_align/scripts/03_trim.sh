#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/03_trim.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_trimal

# Step 3: Trim all 15 alignments with trimAl automated1
# Run after both 02_align_pcg.sh and 02_align_rrna.sh are complete
# Submit as: sbatch 03_align/scripts/03_trim.sh

set -euo pipefail

module load gnu/12
module load prebin/kry
module load trimal/1.5.1

# ── Paths ─────────────────────────────────────────────────────────────────────
ALIGN_DIR="03_align/02_aligned"
TRIM_DIR="03_align/03_trimmed"

mkdir -p "${TRIM_DIR}" logs

PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)

# ── Trim PCGs (input is _NT.fasta from MACSE) ─────────────────────────────────
for gene in "${PCG_GENES[@]}"; do
    IN="${ALIGN_DIR}/${gene}_NT.fasta"
    if [[ ! -f "$IN" ]]; then
        echo "ERROR: ${IN} not found. Run 02_align_pcg.sh first."
        exit 1
    fi
    echo "trimAl: ${gene}"
    trimal -in "$IN" -out "${TRIM_DIR}/${gene}.fasta" -automated1 \
        2> "${TRIM_DIR}/${gene}_trimal.log"
done

# ── Trim rRNAs ────────────────────────────────────────────────────────────────
for gene in "${RRNA_GENES[@]}"; do
    IN="${ALIGN_DIR}/${gene}.fasta"
    if [[ ! -f "$IN" ]]; then
        echo "ERROR: ${IN} not found. Run 02_align_rrna.sh first."
        exit 1
    fi
    echo "trimAl: ${gene}"
    trimal -in "$IN" -out "${TRIM_DIR}/${gene}.fasta" -automated1 \
        2> "${TRIM_DIR}/${gene}_trimal.log"
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Trimmed alignment summary ==="
printf "%-10s | %14s | %s\n" "Gene" "Trimmed length" "Sequences"
echo "-----------|----------------|----------"
for gene in "${PCG_GENES[@]}" "${RRNA_GENES[@]}"; do
    FA="${TRIM_DIR}/${gene}.fasta"
    if [[ -f "$FA" ]]; then
        NSEQ=$(grep -c "^>" "$FA")
        LEN=$(grep -v "^>" "$FA" | head -1 | tr -d '\n' | wc -c)
        printf "%-10s | %14s | %s\n" "$gene" "$LEN" "$NSEQ"
    else
        printf "%-10s | %14s | %s\n" "$gene" "MISSING" "0"
    fi
done

echo ""
echo "Done. Trimmed alignments in: ${TRIM_DIR}/"
echo "Next: concatenate and partition for RevBayes"
