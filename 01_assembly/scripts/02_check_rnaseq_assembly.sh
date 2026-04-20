#!/bin/bash
# Summarise MITGARD assembly results for all RNAseq samples
# Usage: bash 01_assembly/scripts/02_check_rnaseq_assembly.sh

OUT_DIR="01_assembly/output_rnaseq"

SAMPLES=(
    ItCol_f_1  ItCol_f_3  ItCol_f_5  ItCol_f_6  ItCol_f_16  ItCol_f_17
    ItCol_m_7  ItCol_m_9  ItCol_m_10 ItCol_m_12 ItCol_m_14  ItCol_m_15
)

printf "%-20s %-12s %-10s\n" "SAMPLE" "STATUS" "LENGTH"
printf "%-20s %-12s %-10s\n" "------" "------" "------"

for sample in "${SAMPLES[@]}"; do
    fa="${OUT_DIR}/${sample}/${sample}_mitogenome.fa"
    log_pattern="logs/mitgard_*.out"

    # find the log file for this sample
    log=$(grep -l "assembling ${sample}" logs/mitgard_*.out 2>/dev/null | head -1)

    if [[ ! -f "$fa" ]]; then
        printf "%-20s %-12s %-10s\n" "$sample" "MISSING" "N/A"
        continue
    fi

    # check for DONE in log
    if [[ -n "$log" ]] && grep -q "^DONE: ${sample}" "$log"; then
        status="COMPLETE"
    else
        status="INCOMPLETE"
    fi

    nseq=$(grep -c "^>" "$fa")
    length=$(grep -v "^>" "$fa" | tr -d '\n' | wc -c)

    # flag if length is outside expected range
    if [[ "$length" -lt 15400 || "$length" -gt 18500 ]]; then
        status="${status}(!)"
    fi

    # flag multiple sequences
    if [[ "$nseq" -gt 1 ]]; then
        status="${status} (${nseq} seqs)"
    fi

    printf "%-20s %-12s %-10s\n" "$sample" "$status" "$length"
done

echo ""
echo "NOTE: Traceback errors in logs (collections.Hashable) are a known"
echo "Python 3.10+ compatibility issue in MITGARD and are non-fatal."
echo "All samples with STATUS=COMPLETE and LENGTH 15400-18500 bp are good to proceed."
