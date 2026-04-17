#!/bin/bash
# Mitogenome annotation with MITOS2
# Usage: bash 02_annotation/scripts/01_run_mitos2.sh

WORK_DIR="$(pwd)"
ASSEMBLY_DIR="$WORK_DIR/01_assembly/assemblies"
OUT_DIR="$WORK_DIR/02_annotation/output"
LOG_DIR="$WORK_DIR/logs"
SLURM_DIR="$WORK_DIR/02_annotation/slurm_scripts"
REFDIR="$WORK_DIR/utils/mitos2_db"
REFSEQ="refseq63m"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$SLURM_DIR"

for fasta in "$ASSEMBLY_DIR"/*.fasta; do
    sample=$(basename "$fasta" .fasta)

    # skip if already completed
    if [[ -d "$OUT_DIR/$sample" ]]; then
        echo "SKIPPING $sample — output already exists"
        continue
    fi

    # skip if job already running
    if squeue -u "$USER" -o "%o" -h | grep -qF "mitos2_${sample}.out"; then
        echo "SKIPPING $sample — job already in queue or running"
        continue
    fi

    slurm_script="$SLURM_DIR/mitos2_${sample}.slurm"

    cat > "$slurm_script" << EOF
#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=${LOG_DIR}/mitos2_${sample}.out
#SBATCH --qos=normal_prio
#SBATCH -D ${WORK_DIR}
#SBATCH -J mt2_${sample}

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate mitos2

mkdir -p ${OUT_DIR}/${sample}

runmitos \\
    -i ${fasta} \\
    -c 5 \\
    -o ${OUT_DIR}/${sample} \\
    -r ${REFSEQ} \\
    -R ${REFDIR} \\
    --best
EOF

    echo "Submitting $sample ..."
    sbatch "$slurm_script"

done

echo "All MITOS2 jobs submitted."
