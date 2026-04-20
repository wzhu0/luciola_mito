#!/bin/bash
# Collect MITGARD assemblies into 01_assembly/assemblies/
# Usage: bash 01_assembly/scripts/03_collect_rnaseq_assemblies.sh

OUT_DIR="01_assembly/output_rnaseq"
ASSEMBLY_DIR="01_assembly/assemblies"

SAMPLES=(
    ItCol_f_1  ItCol_f_3  ItCol_f_5  ItCol_f_6  ItCol_f_16  ItCol_f_17
    ItCol_m_7  ItCol_m_9  ItCol_m_10 ItCol_m_12 ItCol_m_14  ItCol_m_15
)

mkdir -p "$ASSEMBLY_DIR"

for sample in "${SAMPLES[@]}"; do
    src="${OUT_DIR}/${sample}/${sample}_mitogenome.fa"
    dst="${ASSEMBLY_DIR}/${sample}.fasta"

    if [[ ! -f "$src" ]]; then
        echo "MISSING: $src — skipping"
        continue
    fi

    if [[ -f "$dst" ]]; then
        echo "SKIPPING $sample — already in assemblies/"
        continue
    fi

    cp "$src" "$dst"
    echo "COLLECTED: $sample"
done

echo "Done."
