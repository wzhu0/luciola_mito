#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --array=0-11
#SBATCH --output=/home/wzhu/luciola/mito/logs/mitgard_%a.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mitgard

# Run MITGARD mitogenome assembly on merged RNAseq samples
# Output: 01_assembly/output_rnaseq/<sample>/
# Usage: sbatch 01_assembly/scripts/01_run_mitgard.sh

set -euo pipefail

set +u
source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate mitgard
set -u

WORK_DIR="/home/wzhu/luciola/mito"
RAW_DIR="${WORK_DIR}/00_raw_reads"
OUT_DIR="${WORK_DIR}/01_assembly/output_rnaseq"
REFERENCE="${WORK_DIR}/utils/Luciola_lusitanica.fasta"

SAMPLES=(
    f_1  f_3  f_5  f_6  f_16  f_17
    m_7  m_9  m_10 m_12 m_14  m_15
)

key="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
sex="${key%%_*}"
id="${key#*_}"
sample="ItCol_${sex}_${id}"

R1="${RAW_DIR}/${sample}_1.fq.gz"
R2="${RAW_DIR}/${sample}_2.fq.gz"

echo "Job ${SLURM_ARRAY_TASK_ID}: assembling ${sample}"

# check input files
for f in "$R1" "$R2" "$REFERENCE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing $f"
        exit 1
    fi
done

# skip if already completed
if [[ -f "${OUT_DIR}/${sample}/${sample}.fasta" ]]; then
    echo "SKIPPING ${sample} — output already exists"
    exit 0
fi

mkdir -p "${OUT_DIR}/${sample}"
cd "${OUT_DIR}/${sample}"

MITGARD.py \
    -s "$sample" \
    -1 "$R1" \
    -2 "$R2" \
    -R "$REFERENCE" \
    -L R \
    -c 8 \
    -M 28G

echo "DONE: ${sample}"
