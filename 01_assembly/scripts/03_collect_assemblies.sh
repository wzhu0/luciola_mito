#!/bin/bash
# Collect final assembly FASTAs into 01_assembly/assemblies/
# Usage: bash 01_assembly/scripts/03_collect_assemblies.sh

WORK_DIR="$HOME/luciola/mito"
OUT_DIR="$WORK_DIR/01_assembly/output"
FASTA_DIR="$WORK_DIR/01_assembly/assemblies"
SAMPLE_LIST="$WORK_DIR/utils/sample_list.txt"

mkdir -p "$FASTA_DIR"

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue

    sample=$(basename "$line")

    complete="$OUT_DIR/$sample/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
    scaffold="$OUT_DIR/$sample/animal_mt.K105.scaffolds.graph1.1.path_sequence.fasta"

    if [[ -f "$complete" ]]; then
        cp "$complete" "$FASTA_DIR/${sample}.fasta"
        echo "COPIED   $sample (complete)"
    elif [[ -f "$scaffold" ]]; then
        cp "$scaffold" "$FASTA_DIR/${sample}.fasta"
        echo "COPIED   $sample (scaffold)"
    else
        echo "MISSING  $sample — no assembly found"
    fi

done < "$SAMPLE_LIST"

echo "Done. Assemblies in $FASTA_DIR"
