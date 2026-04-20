#!/bin/bash
# Mitogenome annotation with MITOS2 for RNAseq (MITGARD) assemblies
# Usage: bash 02_annotation/scripts/01_run_mitos2_rnaseq.sh

WORK_DIR="$(pwd)"
ASSEMBLY_DIR="$WORK_DIR/01_assembly/assemblies"
OUT_DIR="$WORK_DIR/02_annotation/output"
LOG_DIR="$WORK_DIR/logs"
SLURM_DIR="$WORK_DIR/02_annotation/slurm_scripts"
REFDIR="$WORK_DIR/utils/mitos2_db"
REFSEQ="refseq63m"

SAMPLES=(
    ItCol_f_1  ItCol_f_3  ItCol_f_5  ItCol_f_6  ItCol_f_16  ItCol_f_17
    ItCol_m_7  ItCol_m_9  ItCol_m_10 ItCol_m_12 ItCol_m_14  ItCol_m_15
)

mkdir -p "$OUT_DIR" "$LOG_DIR" "$SLURM_DIR"

for sample in "${SAMPLES[@]}"; do
    fasta="$ASSEMBLY_DIR/${sample}.fasta"

    if [[ ! -f "$fasta" ]]; then
        echo "MISSING $sample — fasta not found, skipping"
        continue
    fi

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
