#!/bin/bash

WORK_DIR="$HOME/luciola/mito"
ASSEMBLY_DIR="$WORK_DIR/01_assembly/assemblies"
OUT_DIR="$WORK_DIR/02_annotation/output/mitofinder"
LOG_DIR="$WORK_DIR/logs"
SLURM_DIR="$WORK_DIR/02_annotation/slurm_scripts"
SAMPLE_LIST="$WORK_DIR/utils/sample_list.txt"
REF="$WORK_DIR/utils/Luciola_lusitanica.gb"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$SLURM_DIR"

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue

    sample=$(basename "$line")
    fasta="$ASSEMBLY_DIR/${sample}.fasta"

    if [[ ! -f "$fasta" ]]; then
        echo "WARNING: fasta not found for $sample, skipping"
        continue
    fi

    # skip if already completed
    if [[ -f "$OUT_DIR/$sample/${sample}_MitoFinder_arwen_Final_Results/${sample}_mtDNA_contig.gb" ]]; then
        echo "SKIPPING $sample — already completed"
        continue
    fi

    # skip if job already running
    if squeue -u "$USER" -o "%o" -h | grep -qF "mitofinder_${sample}.out"; then
        echo "SKIPPING $sample — job already in queue or running"
        continue
    fi

    slurm_script="$SLURM_DIR/mitofinder_${sample}.slurm"

    cat > "$slurm_script" <<EOF
#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=3-00:00:00
#SBATCH --output=${LOG_DIR}/mitofinder_${sample}.out
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=wenjie.zhu@lmu.de
#SBATCH --qos=normal_prio
#SBATCH -D ${OUT_DIR}
#SBATCH -J mf_${sample}

module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate mitofinder

mitofinder \\
    -j ${sample} \\
    -a ${fasta} \\
    -r ${REF} \\
    -o 5 \\
    -p 2 \\
    --adjust-direction \\
    -t arwen
EOF

    echo "Submitting $sample ..."
    sbatch "$slurm_script"

done < "$SAMPLE_LIST"

echo "All MitoFinder jobs submitted."
