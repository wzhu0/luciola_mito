#!/bin/bash
# Append GenBank gene-level sequences to per-gene extracted FASTAs
# Usage: bash 03_align/scripts/02_add_gb_samples.sh

GB_DIR="03_align/00_gb_genes"
EXTRACT_DIR="03_align/01_extracted"

for fasta in "${GB_DIR}"/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    echo "Processing ${sample} ..."

    gene=""
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "${line// }" ]] && continue

        if [[ "$line" == ">"* ]]; then
            gene="${line#>}"
        else
            out="${EXTRACT_DIR}/${gene}.fasta"
            if [[ ! -f "$out" ]]; then
                echo "  WARNING: ${out} not found — skipping gene ${gene}"
                continue
            fi
            # check if sample already present to avoid duplicates
            if grep -q "^>${sample}$" "$out"; then
                echo "  SKIPPING ${gene} — ${sample} already present"
                continue
            fi
            echo ">${sample}" >> "$out"
            echo "$line"     >> "$out"
            echo "  Appended ${gene}"
        fi
    done < "$fasta"

done

echo "Done"
