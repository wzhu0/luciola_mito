#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --array=0-11
#SBATCH --output=/home/wzhu/luciola/mito/logs/rnaseq_merge_%a.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J rnaseq_merge

# Merge body (A) and head (T) RNAseq libraries per individual
# Submit as: sbatch 01_assembly/scripts/00_merge_rnaseq.sh


set -euo pipefail

RAW_DIR="00_raw_reads"

SAMPLES=(
    f_1  f_3  f_5  f_6  f_16  f_17
    m_7  m_9  m_10 m_12 m_14  m_15
)

key="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
sex="${key%%_*}"
id="${key#*_}"

prefix="${sex}_LuLu_${id}"
outname="ItCol_${sex}_${id}"

R1_body="${RAW_DIR}/${prefix}A_1_val_1.fq.gz"
R2_body="${RAW_DIR}/${prefix}A_2_val_2.fq.gz"
R1_head="${RAW_DIR}/${prefix}T_1_val_1.fq.gz"
R2_head="${RAW_DIR}/${prefix}T_2_val_2.fq.gz"

echo "Job ${SLURM_ARRAY_TASK_ID}: merging ${outname}"

for f in "$R1_body" "$R2_body" "$R1_head" "$R2_head"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing $f"
        exit 1
    fi
done

OUT_R1="${RAW_DIR}/${outname}_1.fq.gz"
OUT_R2="${RAW_DIR}/${outname}_2.fq.gz"

if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
    echo "SKIPPING ${outname} — merged files already exist"
    exit 0
fi

cat "$R1_body" "$R1_head" > "$OUT_R1"
cat "$R2_body" "$R2_head" > "$OUT_R2"

echo "DONE: ${outname}"
