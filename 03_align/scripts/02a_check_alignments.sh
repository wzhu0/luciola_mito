#!/bin/bash
# Sanity check for gene alignments: sequence counts and alignment lengths
# Usage: bash 03_align/scripts/02a_check_alignments.sh

ALIGN_DIR="03_align/02_aligned"
SAMPLE_LIST="utils/sample_list.txt"

TOTAL=$(awk 'NF' "$SAMPLE_LIST" | wc -l)
PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)

printf "%-10s | %6s | %10s\n" "GENE" "NSEQ" "LENGTH(bp)"
printf "%-10s | %6s | %10s\n" "----" "-----" "----------"

ALL_OK=true

for gene in "${PCG_GENES[@]}"; do
    fa="${ALIGN_DIR}/${gene}_NT.fasta"
    if [[ ! -f "$fa" ]]; then
        printf "%-10s | %6s | %10s\n" "$gene" "MISSING" "-"
        ALL_OK=false; continue
    fi
    nseq=$(grep -c "^>" "$fa")
    len=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' "$fa" | sort -u | tr '\n' ',' | sed 's/,$//')
    [[ "$nseq" -ne "$TOTAL" ]] && ALL_OK=false && flag=" <--" || flag=""
    printf "%-10s | %6s | %10s%s\n" "$gene" "${nseq}/${TOTAL}" "$len" "$flag"
done

for gene in "${RRNA_GENES[@]}"; do
    fa="${ALIGN_DIR}/${gene}.fasta"
    if [[ ! -f "$fa" ]]; then
        printf "%-10s | %6s | %10s\n" "$gene" "MISSING" "-"
        ALL_OK=false; continue
    fi
    nseq=$(grep -c "^>" "$fa")
    len=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' "$fa" | sort -u | tr '\n' ',' | sed 's/,$//')
    [[ "$nseq" -ne "$TOTAL" ]] && ALL_OK=false && flag=" <--" || flag=""
    printf "%-10s | %6s | %10s%s\n" "$gene" "${nseq}/${TOTAL}" "$len" "$flag"
done

echo ""
$ALL_OK && echo "All alignments OK." || echo "WARNING: issues found. Check flagged genes above."
