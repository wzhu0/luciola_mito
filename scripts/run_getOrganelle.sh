#!/bin/bash

WORK_DIR="$HOME/luciola/mito"
DATA_DIR="$WORK_DIR/data"
LOG_DIR="$WORK_DIR/logs"
OUT_DIR="$WORK_DIR/getorganelle_out"
SLURM_DIR="$WORK_DIR/slurm_scripts"
SAMPLE_LIST="$WORK_DIR/utils/sample_list.txt"

mkdir -p "$LOG_DIR" "$OUT_DIR" "$SLURM_DIR"

while IFS= read -r line || [[ -n "$line" ]]; do
    # skip empty lines
    [[ -z "${line// }" ]] && continue

    sample=$(basename "$line")

    # determine suffix: .fixed if 2020 or pooled sample, otherwise not
    if [[ "$line" == 2020/* || "$sample" == *pooled* ]]; then
        r1="$DATA_DIR/${sample}_1_val_1.fixed.fq.gz"
        r2="$DATA_DIR/${sample}_2_val_2.fixed.fq.gz"
    else
        r1="$DATA_DIR/${sample}_1_val_1.fq.gz"
        r2="$DATA_DIR/${sample}_2_val_2.fq.gz"
    fi

    # check R1 and R2 exist
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "WARNING: reads not found for $sample, skipping"
        continue
    fi

    # skip if already completed (circular or scaffold)
    if grep -q "Thank you!" "$OUT_DIR/$sample/get_org.log.txt" 2>/dev/null; then
        echo "SKIPPING $sample — already completed"
        continue
    fi

    # skip if job is currently running or pending in SLURM (match by output log path)
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
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=wenjie.zhu@lmu.de
#SBATCH --qos=normal_prio
#SBATCH -D ${WORK_DIR}
#SBATCH -J asm_${sample:0:12}

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
