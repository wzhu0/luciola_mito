#!/bin/bash
# Mitogenome assembly curation QC
# Usage: bash mito_curation_qc.sh
# Run from: ~/luciola/mito/

SAMPLE_LIST="utils/sample_list.txt"
ASSEMBLY_DIR="01_assembly/assemblies"
ANNOTATION_DIR="02_annotation/output/mitos2"
OUT="02_annotation/curation_qc.tsv"

# Expected values
MIN_LEN=15500
MAX_LEN=18500
EXPECTED_PCG=13
EXPECTED_RRNA=2
EXPECTED_TRNA=22

# Header
echo -e "SAMPLE\tASM_LEN\tLEN_OK\tN_CONTIGS\tCONTIG_OK\tPCG_COUNT\tPCG_OK\tRRNA_COUNT\tRRNA_OK\tTRNA_COUNT\tTRNA_OK\tDUPLICATED_GENES\tISSUES" > "$OUT"

while read -r sample; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue

    ISSUES=""

    # ── 1. Find assembly FASTA ──────────────────────────────────────────────
    # sample_list.txt may have year/ prefix (e.g. 2022/AlVe-...), strip it
    sample_id="${sample##*/}"

    # Try {sample_id}.fasta first, then {sample_id}.fa
    if [[ -f "${ASSEMBLY_DIR}/${sample_id}.fasta" ]]; then
        ASM_FA="${ASSEMBLY_DIR}/${sample_id}.fasta"
    elif [[ -f "${ASSEMBLY_DIR}/${sample_id}.fa" ]]; then
        ASM_FA="${ASSEMBLY_DIR}/${sample_id}.fa"
    else
        ASM_FA=$(find "${ASSEMBLY_DIR}" -name "${sample_id}*" 2>/dev/null | head -1)
    fi

    if [[ -z "$ASM_FA" ]]; then
        echo -e "${sample}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNO_ASSEMBLY_FOUND" >> "$OUT"
        continue
    fi

    # ── 2. Genome length and contig count ───────────────────────────────────
    N_CONTIGS=$(grep -c "^>" "$ASM_FA")
    ASM_LEN=$(grep -v "^>" "$ASM_FA" | tr -d '\n' | wc -c)

    CONTIG_OK="YES"
    LEN_OK="YES"
    [[ "$N_CONTIGS" -ne 1 ]] && CONTIG_OK="NO" && ISSUES="${ISSUES}MULTIPLE_CONTIGS(${N_CONTIGS});"
    [[ "$ASM_LEN" -lt "$MIN_LEN" ]] && LEN_OK="NO" && ISSUES="${ISSUES}TOO_SHORT(${ASM_LEN});"
    [[ "$ASM_LEN" -gt "$MAX_LEN" ]] && LEN_OK="NO" && ISSUES="${ISSUES}TOO_LONG(${ASM_LEN});"

    # ── 3. MITOS2 GFF annotation ─────────────────────────────────────────────
    GFF="${ANNOTATION_DIR}/${sample_id}/result.gff"
    if [[ ! -f "$GFF" ]]; then
        echo -e "${sample}\t${ASM_LEN}\t${LEN_OK}\t${N_CONTIGS}\t${CONTIG_OK}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNO_GFF_FOUND" >> "$OUT"
        continue
    fi

    # Count features using Name= field on 'gene' feature lines (MITOS2 format)
    PCG_NAMES="nad1|nad2|nad3|nad4|nad4l|nad5|nad6|cox1|cox2|cox3|atp6|atp8|cob"
    RRNA_NAMES="rrnS|rrnL|rrn12|rrn16"

    N_PCG=$(awk '$3=="gene"' "$GFF" | grep -iE "Name=(${PCG_NAMES});" | wc -l)
    N_RRNA=$(awk '$3=="ncRNA_gene"' "$GFF" | grep -iE "Name=(rrnS|rrnL);" | wc -l)
    N_TRNA=$(awk '$3=="ncRNA_gene"' "$GFF" | grep -i "Name=trn" | wc -l)

    PCG_OK="YES";  RRNA_OK="YES"; TRNA_OK="YES"
    [[ "$N_PCG"  -ne "$EXPECTED_PCG"  ]] && PCG_OK="NO"  && ISSUES="${ISSUES}PCG=${N_PCG};"
    [[ "$N_RRNA" -ne "$EXPECTED_RRNA" ]] && RRNA_OK="NO" && ISSUES="${ISSUES}RRNA=${N_RRNA};"
    [[ "$N_TRNA" -ne "$EXPECTED_TRNA" ]] && TRNA_OK="NO" && ISSUES="${ISSUES}TRNA=${N_TRNA};"

    # ── 4. Check for duplicated gene annotations ─────────────────────────────
    DUPS=$(cat <(awk '$3=="gene"' "$GFF" | grep -iE "Name=(${PCG_NAMES});") \
               <(awk '$3=="ncRNA_gene"' "$GFF" | grep -iE "Name=(rrnS|rrnL);") \
           | grep -oP "Name=\K[^;]+" \
           | sort | uniq -d | tr '\n' ',' | sed 's/,$//')
    [[ -n "$DUPS" ]] && ISSUES="${ISSUES}DUPLICATED:${DUPS};"

    # ── 5. Summary ───────────────────────────────────────────────────────────
    [[ -z "$ISSUES" ]] && ISSUES="OK"

    echo -e "${sample}\t${ASM_LEN}\t${LEN_OK}\t${N_CONTIGS}\t${CONTIG_OK}\t${N_PCG}\t${PCG_OK}\t${N_RRNA}\t${RRNA_OK}\t${N_TRNA}\t${TRNA_OK}\t${DUPS:-none}\t${ISSUES}" >> "$OUT"

done < "$SAMPLE_LIST"

echo ""
echo "QC complete: $OUT"
echo ""
echo "=== SUMMARY ==="
echo "Total samples:    $(grep -v "^SAMPLE" "$OUT" | wc -l)"
echo "Clean (OK):       $(grep -v "^SAMPLE" "$OUT" | grep -c "	OK$")"
echo "With issues:      $(grep -v "^SAMPLE" "$OUT" | grep -vc "	OK$")"
echo ""
echo "=== SAMPLES WITH ISSUES ==="
grep -v "^SAMPLE" "$OUT" | grep -v "	OK$" | awk -F'\t' '{print $1"\t"$NF}'
