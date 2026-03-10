#!/bin/bash
# Annotation QC: check gene counts and assembly length for all samples
# Usage: bash 02_annotation/scripts/02_mito_curation_qc.sh
# Run from: ~/luciola/mito/

SAMPLE_LIST="utils/sample_list.txt"
ASSEMBLY_DIR="01_assembly/assemblies"
ANNOTATION_DIR="02_annotation/output"
OUT="02_annotation/curation_qc.tsv"

MIN_LEN=15400
MAX_LEN=18500
EXPECTED_PCG=13
EXPECTED_RRNA=2
EXPECTED_TRNA=22
PCG_NAMES="nad1|nad2|nad3|nad4|nad4l|nad5|nad6|cox1|cox2|cox3|atp6|atp8|cob"

echo -e "SAMPLE\tASM_LEN\tLEN_OK\tN_CONTIGS\tCONTIG_OK\tPCG_COUNT\tPCG_OK\tRRNA_COUNT\tRRNA_OK\tTRNA_COUNT\tTRNA_OK\tDUPLICATED_GENES\tISSUES" > "$OUT"

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue
    sample=$(basename "$line")
    ISSUES=""

    fasta="${ASSEMBLY_DIR}/${sample}.fasta"
    if [[ ! -f "$fasta" ]]; then
        echo -e "${sample}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNO_ASSEMBLY" >> "$OUT"; continue
    fi

    N_CONTIGS=$(grep -c "^>" "$fasta")
    ASM_LEN=$(grep -v "^>" "$fasta" | tr -d '\n' | wc -c)
    LEN_OK="YES"; CONTIG_OK="YES"
    [[ "$N_CONTIGS" -ne 1 ]]      && CONTIG_OK="NO" && ISSUES="${ISSUES}MULTIPLE_CONTIGS(${N_CONTIGS});"
    [[ "$ASM_LEN" -lt "$MIN_LEN" ]] && LEN_OK="NO"   && ISSUES="${ISSUES}TOO_SHORT(${ASM_LEN});"
    [[ "$ASM_LEN" -gt "$MAX_LEN" ]] && LEN_OK="NO"   && ISSUES="${ISSUES}TOO_LONG(${ASM_LEN});"

    gff="${ANNOTATION_DIR}/${sample}/result.gff"
    if [[ ! -f "$gff" ]]; then
        echo -e "${sample}\t${ASM_LEN}\t${LEN_OK}\t${N_CONTIGS}\t${CONTIG_OK}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNO_GFF" >> "$OUT"; continue
    fi

    N_PCG=$(awk  '$3=="gene"'      "$gff" | grep -iE "Name=(${PCG_NAMES});" | wc -l)
    N_RRNA=$(awk '$3=="ncRNA_gene"' "$gff" | grep -iE "Name=(rrnS|rrnL);"   | wc -l)
    N_TRNA=$(awk '$3=="ncRNA_gene"' "$gff" | grep -i  "Name=trn"             | wc -l)
    PCG_OK="YES"; RRNA_OK="YES"; TRNA_OK="YES"
    [[ "$N_PCG"  -ne "$EXPECTED_PCG"  ]] && PCG_OK="NO"  && ISSUES="${ISSUES}PCG=${N_PCG};"
    [[ "$N_RRNA" -ne "$EXPECTED_RRNA" ]] && RRNA_OK="NO" && ISSUES="${ISSUES}RRNA=${N_RRNA};"
    [[ "$N_TRNA" -ne "$EXPECTED_TRNA" ]] && TRNA_OK="NO" && ISSUES="${ISSUES}TRNA=${N_TRNA};"

    DUPS=$(cat <(awk '$3=="gene"'      "$gff" | grep -iE "Name=(${PCG_NAMES});") \
               <(awk '$3=="ncRNA_gene"' "$gff" | grep -iE "Name=(rrnS|rrnL);") \
           | grep -oP "Name=\K[^;]+" | sort | uniq -d | tr '\n' ',' | sed 's/,$//')
    [[ -n "$DUPS" ]] && ISSUES="${ISSUES}DUPLICATED:${DUPS};"

    [[ -z "$ISSUES" ]] && ISSUES="OK"
    echo -e "${sample}\t${ASM_LEN}\t${LEN_OK}\t${N_CONTIGS}\t${CONTIG_OK}\t${N_PCG}\t${PCG_OK}\t${N_RRNA}\t${RRNA_OK}\t${N_TRNA}\t${TRNA_OK}\t${DUPS:-none}\t${ISSUES}" >> "$OUT"

done < "$SAMPLE_LIST"

echo "QC complete: $OUT"
echo ""
echo "Total:      $(grep -vc "^SAMPLE" "$OUT")"
echo "Clean:      $(grep -vc "^SAMPLE" "$OUT" | xargs -I{} grep -c "	OK$" "$OUT")"
echo ""
echo "=== Samples with issues ==="
grep -v "^SAMPLE" "$OUT" | grep -v "	OK$" | awk -F'\t' '{print $1"\t"$NF}'
