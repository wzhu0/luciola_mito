#!/bin/bash
# Sanity check for extracted gene sequences before alignment
# Verifies sequence counts and lengths for all 15 genes
# MOD3 FAILs are informational only — MACSE handles reading frame correction
# Only NSEQ mismatches and missing files require action before proceeding
# Usage: bash 03_align/scripts/03_check_extracted.sh
# Run from: ~/luciola/mito/

EXTRACT_DIR="03_align/01_extracted"
TOTAL=$(ls "${EXTRACT_DIR}"/*.fasta 2>/dev/null | head -1 | xargs grep -c "^>" 2>/dev/null || echo "?")
TOTAL=$(grep -c "^>" "${EXTRACT_DIR}/nad3.fasta" 2>/dev/null || echo "?")

PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)

# Expected minimum lengths for rRNAs
RRNA_MIN_rrnL=1251
RRNA_MIN_rrnS=747

printf "%-10s | %6s | %-4s | %s\n" "GENE" "NSEQ" "MOD3" "LENGTH(bp)"
printf "%-10s | %6s | %-4s | %s\n" "----" "-----" "----" "----------"

HAS_ERRORS=false

for gene in "${PCG_GENES[@]}"; do
    fa="${EXTRACT_DIR}/${gene}.fasta"
    if [[ ! -f "$fa" ]]; then
        printf "%-10s | %6s | %-4s | %s\n" "$gene" "MISSING" "-" "-"
        HAS_ERRORS=true; continue
    fi

    nseq=$(grep -c "^>" "$fa")
    [[ "$nseq" -ne "$TOTAL" ]] && HAS_ERRORS=true && count_flag=" <--" || count_flag=""

    lens=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq) print length(seq)}' "$fa" \
           | sort | uniq -c | sort -rn)
    summary=$(echo "$lens" | awk '{printf "%s(%s) ", $2, $1}')

    mod3_bad=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq) print length(seq)}' "$fa" \
               | awk '$1%3!=0' | sort -u | tr '\n' ',' | sed 's/,$//')
    [[ -n "$mod3_bad" ]] && mod3_flag="FAIL" || mod3_flag="OK"

    printf "%-10s | %6s | %-4s | %s%s\n" "$gene" "${nseq}/${TOTAL}" "$mod3_flag" "$summary" "$count_flag"
    [[ -n "$mod3_bad" ]] && echo "  INFO: ${gene} has lengths not divisible by 3: ${mod3_bad} (handled by MACSE)"
done

for gene in "${RRNA_GENES[@]}"; do
    fa="${EXTRACT_DIR}/${gene}.fasta"
    if [[ ! -f "$fa" ]]; then
        printf "%-10s | %6s | %-4s | %s\n" "$gene" "MISSING" "-" "-"
        HAS_ERRORS=true; continue
    fi

    nseq=$(grep -c "^>" "$fa")
    [[ "$nseq" -ne "$TOTAL" ]] && HAS_ERRORS=true && count_flag=" <--" || count_flag=""

    [[ "$gene" == "rrnL" ]] && MIN=$RRNA_MIN_rrnL || MIN=$RRNA_MIN_rrnS
    truncated=$(awk -v min="$MIN" \
        '/^>/{if(name && seq && length(seq)<min) print name"\t"length(seq); name=substr($0,2); seq=""} \
         !/^>/{seq=seq$0} \
         END{if(name && seq && length(seq)<min) print name"\t"length(seq)}' "$fa")

    lens=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq) print length(seq)}' "$fa" \
           | sort | uniq -c | sort -rn)
    summary=$(echo "$lens" | awk '{printf "%s(%s) ", $2, $1}')

    printf "%-10s | %6s | %-4s | %s%s\n" "$gene" "${nseq}/${TOTAL}" "n/a" "$summary" "$count_flag"
    if [[ -n "$truncated" ]]; then
        HAS_ERRORS=true
        echo "  WARNING: truncated ${gene} sequences (< ${MIN} bp):"
        echo "$truncated" | awk '{print "    "$0}'
    fi
done

echo ""
if $HAS_ERRORS; then
    echo "WARNING: NSEQ mismatches or missing files found above — fix before proceeding."
    echo "NOTE: MOD3 FAILs are expected and handled by MACSE, not an error."
else
    echo "All extracted sequences present. MOD3 FAILs (if any) will be handled by MACSE."
fi
