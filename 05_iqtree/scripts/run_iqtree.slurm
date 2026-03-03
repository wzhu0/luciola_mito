#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --output=/home/wzhu/luciola/mito/logs/05_iqtree.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_iqtree

# Step 5b: Concatenate trimmed alignments and run IQ-TREE
# Produces an unrooted ML tree for quick topology check
# Run after 03_trim.sh is complete
# Submit as: sbatch 05_iqtree/scripts/05_concat_iqtree.sh

module load gnu/12
module load prebin/kry

TRIM_DIR="03_align/03_trimmed"
CONCAT_DIR="05_iqtree"

mkdir -p "${CONCAT_DIR}" logs

PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)
ALL_GENES=("${PCG_GENES[@]}" "${RRNA_GENES[@]}")

# ── Step 1: Concatenate all 15 alignments into one FASTA ─────────────────────
echo "Concatenating alignments..."

# Get sample list from first gene
mapfile -t SAMPLES < <(grep "^>" "${TRIM_DIR}/nad1.fasta" | sed 's/^>//')

# Build a python script to do the concatenation reliably
python3 << EOF
import os

trim_dir = "${TRIM_DIR}"
out_dir  = "${CONCAT_DIR}"
genes    = "${ALL_GENES[*]}".split()

# Read all alignments into dict: gene -> {sample: seq}
aln = {}
for gene in genes:
    aln[gene] = {}
    fa = os.path.join(trim_dir, gene + ".fasta")
    name = None
    with open(fa) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:]
                aln[gene][name] = ""
            else:
                aln[gene][name] += line

# Get sample order from first gene
samples = list(aln[genes[0]].keys())

# Write concatenated FASTA
with open(os.path.join(out_dir, "concat.fasta"), "w") as out:
    for sample in samples:
        seq = "".join(aln[gene].get(sample, "") for gene in genes)
        out.write(">" + sample + "\n")
        out.write(seq + "\n")

# Report
first_seq = "".join(aln[gene].get(samples[0], "") for gene in genes)
print(f"Concatenated: {len(samples)} taxa, {len(first_seq)} sites")

# Write per-gene lengths for partition file
with open(os.path.join(out_dir, "gene_lengths.txt"), "w") as out:
    for gene in genes:
        seqs = list(aln[gene].values())
        out.write(f"{gene}\t{len(seqs[0])}\n")
EOF

# ── Step 2: Write partition file for IQ-TREE ─────────────────────────────────
echo "Writing partition file..."
PART_FILE="${CONCAT_DIR}/partitions.txt"
> "$PART_FILE"

POS=1
while IFS=$'\t' read -r gene len; do
    END=$(( POS + len - 1 ))

    # PCGs: split by codon position (pos1+2 combined, pos3 separate)
    if [[ " ${PCG_GENES[*]} " == *" ${gene} "* ]]; then
        echo "DNA, ${gene}_pos12 = ${POS}-${END}\2, $((POS+1))-${END}\2" >> "$PART_FILE"
        echo "DNA, ${gene}_pos3  = $((POS+2))-${END}\2"                  >> "$PART_FILE"
    else
        echo "DNA, ${gene} = ${POS}-${END}" >> "$PART_FILE"
    fi

    POS=$(( END + 1 ))
done < "${CONCAT_DIR}/gene_lengths.txt"

echo "Partition file written: ${PART_FILE}"
echo "Total concatenated length: $((POS - 1)) bp"

# ── Step 3: Run IQ-TREE ───────────────────────────────────────────────────────
echo ""
echo "Running IQ-TREE..."

iqtree2 \
    -s  "${CONCAT_DIR}/concat.fasta" \
    -p  "${PART_FILE}" \
    -m  GTR+G+I \
    -B  1000 \
    --redo \
    -T  ${SLURM_CPUS_PER_TASK} \
    --prefix "${CONCAT_DIR}/luciola_mito"

echo ""
echo "IQ-TREE complete. Results in: ${CONCAT_DIR}/"
echo "Tree file: ${CONCAT_DIR}/luciola_mito.treefile"