#!/bin/bash
# Submit one SLURM job per sample for mitogenome assembly with GetOrganelle
# Usage: bash 01_assembly/scripts/run_getOrganelle.sh

WORK_DIR="$HOME/luciola/mito"
READS_DIR="$WORK_DIR/00_raw_reads"
LOG_DIR="$WORK_DIR/logs"
OUT_DIR="$WORK_DIR/01_assembly/output"
SLURM_DIR="$WORK_DIR/01_assembly/slurm_scripts"
SAMPLE_LIST="$WORK_DIR/utils/sample_list.txt"

mkdir -p "$LOG_DIR" "$OUT_DIR" "$SLURM_DIR"

while IFS= read -r line || [[ -n "$line" ]]; do
    # skip empty lines
    [[ -z "${line// }" ]] && continue

    sample=$(basename "$line")

    # find reads
    r1=$(find "$READS_DIR" -maxdepth 1 -name "${sample}_1_val_1*.fq.gz" | head -1)
    r2=$(find "$READS_DIR" -maxdepth 1 -name "${sample}_2_val_2*.fq.gz" | head -1)

    # check R1 and R2 exist
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: reads not found for $sample, skipping"
        continue
    fi

    # skip if already completed (circular or scaffold)
    if grep -q "Thank you!" "$OUT_DIR/$sample/get_org.log.txt" 2>/dev/null; then
        echo "SKIPPING $sample — already completed"
        continue
    fi

    # skip if job is currently running or pending
    if squeue -u "$USER" -o "%o" -h | grep -qF "asm_${sample}.out"; then
        echo "SKIPPING $sample — job already in queue or running"
        continue
    fi

    slurm_script="$SLURM_DIR/asm_${sample}.slurm"

    cat > "$slurm_script" <<EOF
#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --time=3-00:00:00
#SBATCH --output=${LOG_DIR}/asm_${sample}.out
#SBATCH --qos=normal_prio
#SBATCH -D ${WORK_DIR}
#SBATCH -J asm_${sample}

module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate getorganelle

get_organelle_from_reads.py \\
    -1 ${r1} \\
    -2 ${r2} \\
    -R 10 \\
    -k 21,45,65,85,105 \\
    -F animal_mt \\
    -s ${WORK_DIR}/utils/Luciola_lusitanica.fasta \\
    --genes ${WORK_DIR}/utils/Luciola_lusitanica.fasta \\
    -o ${OUT_DIR}/${sample} \\
    -t 2
EOF

    echo "Submitting $sample ..."
    sbatch "$slurm_script"

done < "$SAMPLE_LIST"

echo "All jobs submitted."
