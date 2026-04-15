#!/bin/bash
# Summarise GenBank assembly quality for public mitogenomes
# Usage: bash 01_assembly/scripts/02b_check_genbank_assembly.sh

WORK_DIR="$(pwd)"
ASSEMBLIES_DIR="$WORK_DIR/01_assembly/assemblies"
SAMPLE_LIST="$WORK_DIR/utils/sample_gb.txt"

printf "%-45s %-10s %-10s\n" "SAMPLE" "LENGTH" "NSEQS"
printf "%-45s %-10s %-10s\n" "------" "------" "-----"

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue

    sample="$line"
    fasta="$ASSEMBLIES_DIR/${sample}.fasta"

    if [[ ! -f "$fasta" ]]; then
        printf "%-45s %-10s %-10s\n" "$sample" "NOT FOUND" "N/A"
        continue
    fi

    nseqs=$(grep -c ">" "$fasta")
    length=$(grep -v ">" "$fasta" | tr -d '\n' | wc -c)

    printf "%-45s %-10s %-10s\n" "$sample" "$length" "$nseqs"

done < "$SAMPLE_LIST"
