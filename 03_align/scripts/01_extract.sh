#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/01_extract.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_extract

# Step 1: Extract 13 PCGs + rrnS + rrnL from each assembly using MITOS2 GFF
# Handles genes that wrap around the origin of circular mitogenomes
# Run from: ~/luciola/mito/
# Submit as: sbatch 03_align/scripts/01_extract.sh

set -euo pipefail

module load gnu/12
module load prebin/kry
module load bedtools/2.29.1
module load samtools/1.2.0

# ── Paths ─────────────────────────────────────────────────────────────────────
SAMPLE_LIST="utils/sample_list.txt"
ASSEMBLY_DIR="01_assembly/assemblies"
MITOS2_DIR="02_annotation/output/mitos2"
EXTRACT_DIR="03_align/01_extracted"
TMP_DIR="tmp"

mkdir -p "${EXTRACT_DIR}" logs "${TMP_DIR}"

# ── Gene lists ────────────────────────────────────────────────────────────────
PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)
ALL_GENES=("${PCG_GENES[@]}" "${RRNA_GENES[@]}")

# Initialise empty output files
for gene in "${ALL_GENES[@]}"; do
    > "${EXTRACT_DIR}/${gene}.fasta"
done

# ── Extract ───────────────────────────────────────────────────────────────────
while read -r sample; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    sample_id="${sample##*/}"

    ASM_FA="${ASSEMBLY_DIR}/${sample_id}.fasta"
    GFF="${MITOS2_DIR}/${sample_id}/result.gff"

    if [[ ! -f "$ASM_FA" || ! -f "$GFF" ]]; then
        echo "WARNING: missing files for ${sample_id}, skipping"
        continue
    fi

    # Rename contig to sample_id in temp files
    # Use awk for GFF to avoid regex issues with special chars in contig names
    CONTIG=$(head -1 "$ASM_FA" | sed 's/^>//')
    TMP_FA=$(mktemp ~/luciola/mito/${TMP_DIR}/mito_XXXXXX.fasta)
    TMP_GFF=$(mktemp ~/luciola/mito/${TMP_DIR}/mito_XXXXXX.gff)

    sed "s/^>${CONTIG}/>${sample_id}/" "$ASM_FA" > "$TMP_FA"
    awk -v old="$CONTIG" -v new="$sample_id" \
        'BEGIN{OFS="\t"} $1==old{$1=new} {print}' "$GFF" > "$TMP_GFF"

    samtools faidx "$TMP_FA"

    # Get genome length from fai
    GENOME_LEN=$(awk -v chrom="$sample_id" '$1==chrom{print $2}' "${TMP_FA}.fai")

    # Extract each gene
    for gene in "${ALL_GENES[@]}"; do
        # Feature type differs between PCGs and rRNAs
        if [[ " ${PCG_GENES[*]} " == *" ${gene} "* ]]; then
            FEAT="gene"
        else
            FEAT="ncRNA_gene"
        fi

        GENE_LINE=$(awk -v feat="$FEAT" -v gname="$gene" \
            '$3==feat && $9 ~ "Name="gname";"' "$TMP_GFF" | head -1)

        if [[ -z "$GENE_LINE" ]]; then
            echo "WARNING: ${gene} not found in ${sample_id}"
            continue
        fi

        CHROM=$(echo  "$GENE_LINE" | awk '{print $1}')
        START=$(echo  "$GENE_LINE" | awk '{print $4-1}')  # convert to 0-based
        END=$(echo    "$GENE_LINE" | awk '{print $5}')
        STRAND=$(echo "$GENE_LINE" | awk '{print $7}')

        if [[ "$END" -le "$GENOME_LEN" ]]; then
            # ── Normal case: gene fits within genome bounds ───────────────────
            echo -e "${CHROM}\t${START}\t${END}\t${sample_id}\t0\t${STRAND}" \
                | bedtools getfasta -fi "$TMP_FA" -bed - -s -name \
                | sed "s/^>.*$/>${sample_id}/" \
                >> "${EXTRACT_DIR}/${gene}.fasta"
        else
            # ── Wraparound case: gene crosses the origin of circular genome ───
            echo "  NOTE: ${gene} in ${sample_id} wraps around origin (${START}-${END}, genome=${GENOME_LEN})"
            OVERFLOW=$(( END - GENOME_LEN ))
            TMP_P1=$(mktemp ~/luciola/mito/${TMP_DIR}/piece1_XXXXXX.fa)
            TMP_P2=$(mktemp ~/luciola/mito/${TMP_DIR}/piece2_XXXXXX.fa)

            # Extract both pieces on + strand first, then RC combined if needed
            echo -e "${CHROM}\t${START}\t${GENOME_LEN}\t${sample_id}\t0\t+" \
                | bedtools getfasta -fi "$TMP_FA" -bed - -name \
                | grep -v "^>" > "$TMP_P1"
            echo -e "${CHROM}\t0\t${OVERFLOW}\t${sample_id}\t0\t+" \
                | bedtools getfasta -fi "$TMP_FA" -bed - -name \
                | grep -v "^>" > "$TMP_P2"

            # Concatenate pieces into one sequence
            COMBINED=$(cat "$TMP_P1" "$TMP_P2" | tr -d '\n')

            if [[ "$STRAND" == "-" ]]; then
                # Reverse complement
                COMBINED=$(echo "$COMBINED" | tr 'ACGTacgt' 'TGCAtgca' | rev)
            fi

            echo ">${sample_id}" >> "${EXTRACT_DIR}/${gene}.fasta"
            echo "$COMBINED"     >> "${EXTRACT_DIR}/${gene}.fasta"

            rm -f "$TMP_P1" "$TMP_P2"
        fi
    done

    rm -f "$TMP_FA" "$TMP_GFF" "${TMP_FA}.fai"
    echo "Extracted: ${sample_id}"

done < "$SAMPLE_LIST"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Extraction summary ==="
ALL_OK=true
for gene in "${ALL_GENES[@]}"; do
    N=$(grep -c "^>" "${EXTRACT_DIR}/${gene}.fasta" 2>/dev/null || echo 0)
    STATUS=""
    [[ "$N" -ne 32 ]] && STATUS=" <-- MISSING $(( 32 - N ))" && ALL_OK=false
    echo "  ${gene}: ${N}/32 sequences${STATUS}"
done
echo ""
if $ALL_OK; then
    echo "All genes extracted for all 32 samples."
else
    echo "WARNING: some genes are missing sequences. Check warnings above."
fi
echo ""
echo "Done. Sequences in: ${EXTRACT_DIR}/"
echo "Next: sbatch 03_align/scripts/02_align_pcg.sh"
echo "      sbatch 03_align/scripts/02_align_rrna.sh"
