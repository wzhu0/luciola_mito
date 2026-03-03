#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/reextract_rrnL.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_reextract_rrnL

# Re-extract rrnL for 7 samples where MITOS2 annotation missed the 5' end
# true_end   = trnL1_start + 19  (20 bp rrnL/trnL1 overlap)
# true_start = true_end - 1251 + 1
# Submit as: sbatch 03_align/scripts/reextract_rrnL.sh

module load gnu/12
module load prebin/kry
module load samtools/1.2.0

ASSEMBLY_DIR="/home/wzhu/luciola/mito/01_assembly/output"
EXTRACTED="/home/wzhu/luciola/mito/03_align/01_extracted"
TMP_DIR="/home/wzhu/luciola/mito/tmp"
TARGET_LEN=1251
OVERLAP=19

mkdir -p "$TMP_DIR" /home/wzhu/luciola/mito/logs

OUT="${EXTRACTED}/rrnL_corrected.fasta"
> "$OUT"

echo "Re-extracting rrnL (target: ${TARGET_LEN} bp)..."

# ── AlVe-LuNov-34-m-2022-2 ───────────────────────────────────────────────────
sample="AlVe-LuNov-34-m-2022-2"
FA="${ASSEMBLY_DIR}/${sample}/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
contig="327+(circular)"
true_end=$(( 12207 + OVERLAP ))
true_start=$(( true_end - TARGET_LEN + 1 ))
echo "  ${sample}: ${true_start}-${true_end}"
samtools faidx "$FA" "${contig}:${true_start}-${true_end}" \
    | awk -v s="$sample" '/^>/{print ">"s} !/^>/{print}' >> "$OUT"

# ── AlVe-LuNov-36-m-2022-4 ───────────────────────────────────────────────────
sample="AlVe-LuNov-36-m-2022-4"
FA="${ASSEMBLY_DIR}/${sample}/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
contig="6282+(circular)"
true_end=$(( 9093 + OVERLAP ))
true_start=$(( true_end - TARGET_LEN + 1 ))
echo "  ${sample}: ${true_start}-${true_end}"
samtools faidx "$FA" "${contig}:${true_start}-${true_end}" \
    | awk -v s="$sample" '/^>/{print ">"s} !/^>/{print}' >> "$OUT"

# ── AlVe-LuNov-38-m-2022-6 ───────────────────────────────────────────────────
sample="AlVe-LuNov-38-m-2022-6"
FA="${ASSEMBLY_DIR}/${sample}/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
contig="373+(circular)"
true_end=$(( 9093 + OVERLAP ))
true_start=$(( true_end - TARGET_LEN + 1 ))
echo "  ${sample}: ${true_start}-${true_end}"
samtools faidx "$FA" "${contig}:${true_start}-${true_end}" \
    | awk -v s="$sample" '/^>/{print ">"s} !/^>/{print}' >> "$OUT"

# ── MoDo-LuNov-31-m-2022-1 ───────────────────────────────────────────────────
sample="MoDo-LuNov-31-m-2022-1"
FA="${ASSEMBLY_DIR}/${sample}/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
contig="712-(circular)"
true_end=$(( 11954 + OVERLAP ))
true_start=$(( true_end - TARGET_LEN + 1 ))
echo "  ${sample}: ${true_start}-${true_end}"
samtools faidx "$FA" "${contig}:${true_start}-${true_end}" \
    | awk -v s="$sample" '/^>/{print ">"s} !/^>/{print}' >> "$OUT"

# ── MoDo-LuNov-32-m-2022-2 ───────────────────────────────────────────────────
sample="MoDo-LuNov-32-m-2022-2"
FA="${ASSEMBLY_DIR}/${sample}/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
contig="207-(circular)"
true_end=$(( 5541 + OVERLAP ))
true_start=$(( true_end - TARGET_LEN + 1 ))
echo "  ${sample}: ${true_start}-${true_end}"
samtools faidx "$FA" "${contig}:${true_start}-${true_end}" \
    | awk -v s="$sample" '/^>/{print ">"s} !/^>/{print}' >> "$OUT"

# ── ItPie-LuLus-20-m-2022-5 (wraps origin) ───────────────────────────────────
sample="ItPie-LuLus-20-m-2022-5"
FA="${ASSEMBLY_DIR}/${sample}/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
contig="189-(circular)"
genome_len=16470
true_end=$(( 975 + OVERLAP ))
true_start=$(( true_end - TARGET_LEN + 1 ))
# true_start is negative, so wraps: actual_start = genome_len + true_start
actual_start=$(( genome_len + true_start ))
echo "  ${sample}: wraps origin ${actual_start}-${genome_len} + 1-${true_end}"
TMP1="${TMP_DIR}/itpie_piece1.fasta"
TMP2="${TMP_DIR}/itpie_piece2.fasta"
samtools faidx "$FA" "${contig}:${actual_start}-${genome_len}" | grep -v "^>" | tr -d '\n' > "$TMP1"
samtools faidx "$FA" "${contig}:1-${true_end}"                | grep -v "^>" | tr -d '\n' > "$TMP2"
echo ">${sample}" >> "$OUT"
cat "$TMP1" "$TMP2" >> "$OUT"
echo "" >> "$OUT"
rm -f "$TMP1" "$TMP2"

# ── PoAlf-LuLus-65-m-2022-11 ─────────────────────────────────────────────────
sample="PoAlf-LuLus-65-m-2022-11"
FA="${ASSEMBLY_DIR}/${sample}/animal_mt.K105.complete.graph1.1.path_sequence.fasta"
contig="543_539_545_127_15539_1615-(circular)"
true_end=$(( 3082 + OVERLAP ))
true_start=$(( true_end - TARGET_LEN + 1 ))
echo "  ${sample}: ${true_start}-${true_end}"
samtools faidx "$FA" "${contig}:${true_start}-${true_end}" \
    | awk -v s="$sample" '/^>/{print ">"s} !/^>/{print}' >> "$OUT"

# ── Verify lengths ─────────────────────────────────────────────────────────────
echo ""
echo "=== Extracted rrnL lengths ==="
awk '/^>/{if(name) print name"\t"length(seq); name=substr($0,2); seq=""}
     !/^>/{seq=seq $0}
     END{if(name) print name"\t"length(seq)}' "$OUT"
echo ""
echo "All should be ${TARGET_LEN} bp"
echo "Output: ${OUT}"
