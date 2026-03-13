#!/bin/bash
# Sanity check for extracted gene sequences before alignment
# Verifies sequence counts and lengths for all 15 genes
# Usage: bash 03_align/scripts/01a_check_extracted.sh

EXTRACT_DIR="03_align/01_extracted"
SAMPLE_LIST="utils/sample_list.txt"

TOTAL=$(awk 'NF' "$SAMPLE_LIST" | wc -l)
PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)

# Expected lengths for rRNAs — warn if any sequence deviates
RRNA_MIN_rrnL=1251
RRNA_MIN_rrnS=747

printf "%-10s | %6s | %s\n" "GENE" "NSEQ" "LENGTH(bp)"
printf "%-10s | %6s | %s\n" "----" "-----" "----------"

ALL_OK=true

for gene in "${PCG_GENES[@]}"; do
    fa="${EXTRACT_DIR}/${gene}.fasta"
    if [[ ! -f "$fa" ]]; then
        printf "%-10s | %6s | %s\n" "$gene" "MISSING" "-"
        ALL_OK=false; continue
    fi
    nseq=$(grep -c "^>" "$fa")
    lens=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' "$fa" | sort | uniq -c | sort -rn)
    summary=$(echo "$lens" | awk '{printf "%s(%s) ", $2, $1}')
    [[ "$nseq" -ne "$TOTAL" ]] && ALL_OK=false && flag=" <--" || flag=""
    printf "%-10s | %6s | %s%s\n" "$gene" "${nseq}/${TOTAL}" "$summary" "$flag"
done

for gene in "${RRNA_GENES[@]}"; do
    fa="${EXTRACT_DIR}/${gene}.fasta"
    if [[ ! -f "$fa" ]]; then
        printf "%-10s | %6s | %s\n" "$gene" "MISSING" "-"
        ALL_OK=false; continue
    fi
    nseq=$(grep -c "^>" "$fa")
    [[ "$nseq" -ne "$TOTAL" ]] && ALL_OK=false && flag=" <--" || flag=""

    # Check for truncated sequences
    if [[ "$gene" == "rrnL" ]]; then MIN=$RRNA_MIN_rrnL; else MIN=$RRNA_MIN_rrnS; fi
    truncated=$(awk -v min="$MIN" '/^>/{name=substr($0,2)} !/^>/{seq=seq$0}
        /^>/{if(seq && length(seq)<min) print name"\t"length(seq); seq=""}
        END{if(seq && length(seq)<min) print name"\t"length(seq)}' "$fa")

    lens=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' "$fa" | sort | uniq -c | sort -rn)
    summary=$(echo "$lens" | awk '{printf "%s(%s) ", $2, $1}')
    printf "%-10s | %6s | %s%s\n" "$gene" "${nseq}/${TOTAL}" "$summary" "$flag"

    if [[ -n "$truncated" ]]; then
        ALL_OK=false
        echo "  WARNING: truncated ${gene} sequences (< ${MIN} bp):"
        echo "$truncated" | awk '{print "    "$0}'
    fi
done

echo ""
$ALL_OK && echo "All extracted sequences OK." || echo "WARNING: issues found. Check above before proceeding to alignment."
