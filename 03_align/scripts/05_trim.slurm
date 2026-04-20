#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/05_trim.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J trim

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate trimal

ALIGN_DIR="03_align/02_aligned"
TRIM_DIR="03_align/03_trimmed"
UTILS_DIR="03_align/utils"

mkdir -p "${TRIM_DIR}" "${UTILS_DIR}" logs

PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)

# helper: get alignment length from fasta (wrap-aware)
aln_len() {
    local fa=$1
    local nseq=$(grep -c "^>" "$fa")
    local total=$(grep -v "^>" "$fa" | tr -d '\n' | wc -c)
    echo $((total / nseq))
}

# helper: count ? per sequence (wrap-aware, uses first sequence)
count_masked() {
    awk '/^>/{if(seq!="")exit} /^[^>]/{seq=seq$0} END{gsub(/[^?]/,"",seq); print length(seq)}' "$1"
}

# PCGs: mask removed columns as ? rather than removing them
for gene in "${PCG_GENES[@]}"; do
    IN="${ALIGN_DIR}/${gene}_NT.fasta"
    COLS="${UTILS_DIR}/${gene}_colnumbering.txt"
    OUT="${TRIM_DIR}/${gene}.fasta"

    echo "Masking: ${gene}"

    trimal -in "$IN" -out /dev/null -gt 0.3 -colnumbering 2>/dev/null | \
        grep "^#ColumnsMap" > "$COLS"

    python3 - "$IN" "$COLS" "$OUT" << 'PYEOF'
import sys

fasta_in  = sys.argv[1]
cols_file = sys.argv[2]
fasta_out = sys.argv[3]

with open(cols_file) as f:
    content = f.read()
kept = set(
    int(x.strip())
    for x in content.replace('#ColumnsMap', '').split(',')
    if x.strip().lstrip('-').isdigit()
)

seqs = {}
order = []
with open(fasta_in) as f:
    name = None
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            name = line
            order.append(name)
            seqs[name] = []
        else:
            seqs[name].extend(list(line))

aln_len = len(seqs[order[0]])

with open(fasta_out, 'w') as out:
    for name in order:
        seq = seqs[name]
        masked = ''.join(seq[i] if i in kept else '?' for i in range(aln_len))
        out.write(f"{name}\n")
        for i in range(0, len(masked), 60):
            out.write(masked[i:i+60] + '\n')

n_masked = aln_len - len(kept)
print(f"  {aln_len} cols total, {len(kept)} kept, {n_masked} masked as ?")
PYEOF
done

# rRNAs: standard trimAl
for gene in "${RRNA_GENES[@]}"; do
    IN="${ALIGN_DIR}/${gene}.fasta"
    echo "trimAl: ${gene}"
    trimal -in "$IN" -out "${TRIM_DIR}/${gene}.fasta" -gt 0.3
done

# Summary
echo ""
echo "=== Trimmed alignment summary ==="
printf "%-10s | %10s | %10s | %10s | %s\n" "Gene" "Pre" "Post" "?-masked" "MOD3"
echo "-----------|------------|------------|------------|------"

for gene in "${PCG_GENES[@]}" "${RRNA_GENES[@]}"; do
    if [[ " ${RRNA_GENES[*]} " == *" ${gene} "* ]]; then
        PRE_FA="${ALIGN_DIR}/${gene}.fasta"
    else
        PRE_FA="${ALIGN_DIR}/${gene}_NT.fasta"
    fi
    POST_FA="${TRIM_DIR}/${gene}.fasta"

    PRE=$(aln_len "$PRE_FA")
    POST=$(aln_len "$POST_FA")
    MASKED=$(count_masked "$POST_FA")

    if [[ " ${PCG_GENES[*]} " == *" ${gene} "* ]]; then
        [[ $((POST % 3)) -ne 0 ]] && MOD3="FAIL <--" || MOD3="OK"
    else
        MOD3="n/a"
    fi

    printf "%-10s | %10s | %10s | %10s | %s\n" "$gene" "$PRE" "$POST" "$MASKED" "$MOD3"
done

echo ""
echo "Done. Trimmed alignments in: ${TRIM_DIR}/"
echo "Next: concatenate and partition for RevBayes"

