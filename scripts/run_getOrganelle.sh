#!/bin/bash

DATA_DIR="$HOME/luciola/mito/data"
WORK_DIR="$HOME/luciola/mito"
LOG_DIR="$WORK_DIR/logs"
OUT_DIR="$WORK_DIR/getorganelle_out"
SLURM_DIR="$WORK_DIR/slurm_scripts"

mkdir -p "$LOG_DIR" "$OUT_DIR" "$SLURM_DIR"

# find all R1 files, derive sample name from each
find "$DATA_DIR" -maxdepth 1 -name "*_1_val_1*.fq.gz" | sort | while read -r r1; do
    
    # derive R2 and sample name
    r2="${r1/_1_val_1/_2_val_2}"
    filename=$(basename "$r1")
    sample="${filename%%_1_val_1*}"  # strip everything from _1_val_1 onward

    # check R2 exists
    if [[ ! -f "$r2" ]]; then
        echo "WARNING: R2 not found for $sample, skipping"
        continue
    fi

    # skip if output already exists
    if [[ -f "$OUT_DIR/$sample/get_org.log.txt" ]] && grep -q "Thanks for using GetOrganelle" "$OUT_DIR/$sample/get_org.log.txt" 2>/dev/null; then
        echo "SKIPPING $sample — already completed successfully"
        continue
    fi

    slurm_script="$SLURM_DIR/asm_${sample}.sh"

    cat > "$slurm_script" <<EOF
#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
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

done

echo "All jobs submitted."
