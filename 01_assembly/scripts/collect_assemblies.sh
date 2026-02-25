#!/bin/bash

WORK_DIR="$HOME/luciola/mito"
OUT_DIR="$WORK_DIR/01_assembly/output"
FASTA_DIR="$WORK_DIR/01_assembly/assemblies"
SAMPLE_LIST="$WORK_DIR/utils/sample_list.txt"

mkdir -p "$FASTA_DIR"

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue

    sample=$(basename "$line")
    fasta="$OUT_DIR/$sample/animal_mt.K105.complete.graph1.1.path_sequence.fasta"

    if [[ ! -f "$fasta" ]]; then
        echo "WARNING: fasta not found for $sample, skipping"
        continue
    fi

    cp "$fasta" "$FASTA_DIR/${sample}.fasta"
    echo "Collected $sample"

done < "$SAMPLE_LIST"

echo "Done. All fastas in $FASTA_DIR"

