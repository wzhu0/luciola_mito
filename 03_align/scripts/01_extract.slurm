#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/03_align_extract.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_extract

# Step 1: Extract 13 PCGs + rrnS + rrnL from each assembly using MITOS2 GFF
# Handles genes that wrap around the origin of circular mitogenomes
# rrnL and rrnS boundaries are recomputed from flanking tRNAs (trnL1, trnV)
# to correct frequent MITOS2 5'-end truncation of rRNA genes
# orientations are handled based on + and - strand 
# Submit as: sbatch 03_align/scripts/01_extract.slurm

set -euo pipefail

module load gnu/12
module load prebin/kry
module load bedtools/2.29.1
module load samtools/1.2.0

SAMPLE_LIST="utils/sample_list.txt"
ASSEMBLY_DIR="01_assembly/assemblies"
MITOS2_DIR="02_annotation/output"
EXTRACT_DIR="03_align/01_extracted"
TMP_DIR="tmp"

mkdir -p "${EXTRACT_DIR}" logs "${TMP_DIR}"

# Gene lists 
PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)
ALL_GENES=("${PCG_GENES[@]}" "${RRNA_GENES[@]}")

# Initialise empty output files
for gene in "${ALL_GENES[@]}"; do
    > "${EXTRACT_DIR}/${gene}.fasta"
done

# Extract
# Helper: extract a sequence region, handling wraparound and strand 
# Usage: extract_region CHROM START END STRAND GENOME_LEN TMP_FA OUT_FA SAMPLE
extract_region() {
    local chrom=$1 start=$2 end=$3 strand=$4 genome_len=$5
    local tmp_fa=$6 out_fa=$7 sample=$8

    if [[ "$start" -lt 0 ]]; then
        # 5' end wraps before origin (negative start)
        local actual_start=$(( genome_len + start ))
        local p1=$(mktemp ${TMP_DIR}/piece1_XXXXXX.fa)
        local p2=$(mktemp ${TMP_DIR}/piece2_XXXXXX.fa)
        echo -e "${chrom}\t${actual_start}\t${genome_len}\t${sample}\t0\t+" \
            | bedtools getfasta -fi "$tmp_fa" -bed - -name | grep -v "^>" > "$p1"
        echo -e "${chrom}\t0\t${end}\t${sample}\t0\t+" \
            | bedtools getfasta -fi "$tmp_fa" -bed - -name | grep -v "^>" > "$p2"
        local seq=$(cat "$p1" "$p2" | tr -d '\n')
        [[ "$strand" == "-" ]] && seq=$(echo "$seq" | tr 'ACGTacgt' 'TGCAtgca' | rev)
        echo ">${sample}" >> "$out_fa"; echo "$seq" >> "$out_fa"
        rm -f "$p1" "$p2"
    elif [[ "$end" -gt "$genome_len" ]]; then
        # 3' end wraps past genome end
        local overflow=$(( end - genome_len ))
        local p1=$(mktemp ${TMP_DIR}/piece1_XXXXXX.fa)
        local p2=$(mktemp ${TMP_DIR}/piece2_XXXXXX.fa)
        echo -e "${chrom}\t${start}\t${genome_len}\t${sample}\t0\t+" \
            | bedtools getfasta -fi "$tmp_fa" -bed - -name | grep -v "^>" > "$p1"
        echo -e "${chrom}\t0\t${overflow}\t${sample}\t0\t+" \
            | bedtools getfasta -fi "$tmp_fa" -bed - -name | grep -v "^>" > "$p2"
        local seq=$(cat "$p1" "$p2" | tr -d '\n')
        [[ "$strand" == "-" ]] && seq=$(echo "$seq" | tr 'ACGTacgt' 'TGCAtgca' | rev)
        echo ">${sample}" >> "$out_fa"; echo "$seq" >> "$out_fa"
        rm -f "$p1" "$p2"
    else
        # Normal case
        echo -e "${chrom}\t${start}\t${end}\t${sample}\t0\t${strand}" \
            | bedtools getfasta -fi "$tmp_fa" -bed - -s -name \
            | sed "s/^>.*$/>${sample}/" >> "$out_fa"
    fi
}

# Main extraction loop
while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue
    sample=$(basename "$line")

    ASM_FA="${ASSEMBLY_DIR}/${sample}.fasta"
    GFF="${MITOS2_DIR}/${sample}/result.gff"

    if [[ ! -f "$ASM_FA" || ! -f "$GFF" ]]; then
        echo "WARNING: missing files for ${sample}, skipping"
        continue
    fi

    CONTIG=$(head -1 "$ASM_FA" | sed 's/^>//')
    TMP_FA=$(mktemp ${TMP_DIR}/mito_XXXXXX.fasta)
    TMP_GFF=$(mktemp ${TMP_DIR}/mito_XXXXXX.gff)

    awk -v old="$CONTIG" -v new="$sample" \
    'NR==1 && $0==">"old{print ">"new; next} {print}' "$ASM_FA" > "$TMP_FA"
    awk -v old="$CONTIG" -v new="$sample" \
        'BEGIN{OFS="\t"} $1==old{$1=new} {print}' "$GFF" > "$TMP_GFF"

    samtools faidx "$TMP_FA"

    # Get genome length from fai
    GENOME_LEN=$(awk -v chrom="$sample" '$1==chrom{print $2}' "${TMP_FA}.fai")

    # Extract each gene
    for gene in "${ALL_GENES[@]}"; do
        # Feature type differs between PCGs and rRNAs
        [[ " ${PCG_GENES[*]} " == *" ${gene} "* ]] && FEAT="gene" || FEAT="ncRNA_gene"

        GENE_LINE=$(awk -v feat="$FEAT" -v gname="$gene" \
            '$3==feat && $9 ~ "Name="gname";"' "$TMP_GFF" | head -1)

        if [[ -z "$GENE_LINE" ]]; then
            echo "WARNING: ${gene} not found in ${sample}"
            continue
        fi

        CHROM=$(echo  "$GENE_LINE" | awk '{print $1}')
        STRAND=$(echo "$GENE_LINE" | awk '{print $7}')

        # rrnL and rrnS: recompute boundaries from flanking tRNAs
        # Handles both + and - strand orientations
        # + strand gene order: rrnS → trnV → rrnL → trnL1
        #   rrnL 3' end: trnL1_start + 19  (20bp overlap)
        #   rrnS 3' end: trnV_start - 1    (adjacent, no overlap)
        # - strand gene order (coords): trnL1 → rrnL → trnV → rrnS
        #   rrnL 3' end: trnL1_end - 19
        #   rrnS 3' end: trnV_end + 1        
        if [[ "$gene" == "rrnL" ]]; then
            TRNL1_LINE=$(awk '$3=="ncRNA_gene" && $9 ~ "Name=trnL1;"' "$TMP_GFF" | head -1)
            if [[ -z "$TRNL1_LINE" ]]; then
                echo "WARNING: trnL1 not found in ${sample}, falling back to GFF coords for rrnL"
                START=$(echo "$GENE_LINE" | awk '{print $4-1}')
                END=$(echo   "$GENE_LINE" | awk '{print $5}')
            elif [[ "$STRAND" == "+" ]]; then
                TRNL1_START=$(echo "$TRNL1_LINE" | awk '{print $4}')
                END=$(( TRNL1_START + 19 ))
                START=$(( END - 1251 ))  # 0-based
            else
                TRNL1_END=$(echo "$TRNL1_LINE" | awk '{print $5}')
                START=$(( TRNL1_END - 19 - 1 ))  # 0-based
                END=$(( START + 1251 ))
            fi

        elif [[ "$gene" == "rrnS" ]]; then
            TRNV_LINE=$(awk '$3=="ncRNA_gene" && $9 ~ "Name=trnV;"' "$TMP_GFF" | head -1)
            if [[ -z "$TRNV_LINE" ]]; then
                echo "WARNING: trnV not found in ${sample}, falling back to GFF coords for rrnS"
                START=$(echo "$GENE_LINE" | awk '{print $4-1}')
                END=$(echo   "$GENE_LINE" | awk '{print $5}')
            elif [[ "$STRAND" == "+" ]]; then
                TRNV_START=$(echo "$TRNV_LINE" | awk '{print $4}')
                END=$(( TRNV_START - 1 ))
                START=$(( END - 747 ))  # 0-based
            else
                TRNV_END=$(echo "$TRNV_LINE" | awk '{print $5}')
                START=$(( TRNV_END ))  # 0-based, adjacent
                END=$(( START + 747 ))
            fi

        else
            START=$(echo "$GENE_LINE" | awk '{print $4-1}')  # convert to 0-based
            END=$(echo   "$GENE_LINE" | awk '{print $5}')
        fi

        # Log wraparound cases
        if [[ "$START" -lt 0 || "$END" -gt "$GENOME_LEN" ]]; then
            echo "  NOTE: ${gene} in ${sample} wraps around origin (start=${START}, end=${END}, genome=${GENOME_LEN})"
        fi

        extract_region "$CHROM" "$START" "$END" "$STRAND" "$GENOME_LEN" \
            "$TMP_FA" "${EXTRACT_DIR}/${gene}.fasta" "$sample"
    done

    rm -f "$TMP_FA" "$TMP_GFF" "${TMP_FA}.fai"
    echo "DONE     ${sample}"

done < "$SAMPLE_LIST"

echo ""
echo "=== Extraction summary ==="
TOTAL=$(awk 'NF' "$SAMPLE_LIST" | wc -l)
ALL_OK=true
for gene in "${ALL_GENES[@]}"; do
    N=$(grep -c "^>" "${EXTRACT_DIR}/${gene}.fasta" 2>/dev/null || echo 0)
    STATUS=""
    [[ "$N" -ne "$TOTAL" ]] && STATUS=" <-- MISSING $(( TOTAL - N ))" && ALL_OK=false
    echo "  ${gene}: ${N}/${TOTAL}${STATUS}"
done
echo ""
$ALL_OK && echo "All genes extracted for all samples." \
        || echo "WARNING: some genes missing. Check warnings above."
