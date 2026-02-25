#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/02_align_pcg_%a.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_macse
#SBATCH --array=0-12

# Step 2a: Align 13 PCGs with MACSE v2 (insect mito code)
# Each array job handles one gene in parallel
# Run after 01_extract.sh is complete
# Submit as: sbatch 03_align/scripts/02_align_pcg.sh

set -euo pipefail

module load gnu/12
module load prebin/kry
module load macse/2.07

# ── Paths ─────────────────────────────────────────────────────────────────────
EXTRACT_DIR="03_align/01_extracted"
ALIGN_DIR="03_align/02_aligned"

mkdir -p "${ALIGN_DIR}" logs

# ── Gene for this array task ──────────────────────────────────────────────────
PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
GENE="${PCG_GENES[$SLURM_ARRAY_TASK_ID]}"

echo "Array task ${SLURM_ARRAY_TASK_ID}: aligning ${GENE}"

IN="${EXTRACT_DIR}/${GENE}.fasta"

if [[ ! -f "$IN" ]]; then
    echo "ERROR: ${IN} not found. Run 01_extract.sh first."
    exit 1
fi

# ── MACSE alignment ───────────────────────────────────────────────────────────
# -gc_def 5 = insect mitochondrial genetic code
macse -prog alignSequences \
    -seq "$IN" \
    -gc_def 5 \
    -out_NT "${ALIGN_DIR}/${GENE}_NT.fasta" \
    -out_AA "${ALIGN_DIR}/${GENE}_AA.fasta" \
    2>&1 | tee "${ALIGN_DIR}/${GENE}_macse.log"

# Replace MACSE frameshift character '!' with gap '-' for downstream tools
sed -i 's/!/-/g' "${ALIGN_DIR}/${GENE}_NT.fasta"

echo "Done: ${GENE} -> ${ALIGN_DIR}/${GENE}_NT.fasta"
