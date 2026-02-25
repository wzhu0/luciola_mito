#!/bin/bash

WORK_DIR="$HOME/luciola/mito"
OUT_DIR="$WORK_DIR/getorganelle_out"
SAMPLE_LIST="$WORK_DIR/utils/sample_list.txt"

# print header
printf "%-45s %-12s %-10s %-10s\n" "SAMPLE" "STATUS" "LENGTH" "COVERAGE"
printf "%-45s %-12s %-10s %-10s\n" "------" "------" "------" "--------"

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue

    sample=$(basename "$line")
    log="$OUT_DIR/$sample/get_org.log.txt"

    # check if completed
    if ! grep -q "Thank you!" "$log" 2>/dev/null; then
        printf "%-45s %-12s %-10s %-10s\n" "$sample" "INCOMPLETE" "N/A" "N/A"
        continue
    fi

    # status: circular or scaffold
    status=$(grep "Result status of animal_mt:" "$log" | sed 's/.*: //')

    # coverage
    coverage=$(grep "Average animal_mt base-coverage" "$log" | sed 's/.*= //')

    # length and number of paths
    fasta_complete="$OUT_DIR/$sample/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
    fasta_scaffold="$OUT_DIR/$sample/animal_mt.K105.scaffolds.graph1.1.path_sequence.fasta"

    if [[ -f "$fasta_complete" ]]; then
        fasta="$fasta_complete"
    elif [[ -f "$fasta_scaffold" ]]; then
        fasta="$fasta_scaffold"
    else
        fasta=""
    fi

    if [[ -n "$fasta" ]]; then
        length=$(grep -v ">" "$fasta" | tr -d '\n' | wc -c)
    else
        length="N/A"
    fi

    # check for multiple paths
    npaths=$(grep "Writing PATH" "$log" | wc -l)
    if [[ "$npaths" -gt 1 ]]; then
        status="${status} (${npaths} paths)"
    fi

    printf "%-45s %-12s %-10s %-10s\n" "$sample" "$status" "$length" "$coverage"

done < "$SAMPLE_LIST"
