#!/bin/bash
#SBATCH --partition=lemmium
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/04_iqtree.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_iqtree

module load gnu/12
module load prebin/lem

TRIM_DIR="03_align/03_trimmed"
DATA_DIR="04_iqtree/data"
OUT_DIR="04_iqtree/output"

mkdir -p "${DATA_DIR}" "${OUT_DIR}" logs

PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)
ALL_GENES=("${PCG_GENES[@]}" "${RRNA_GENES[@]}")

# Step 1: Concatenate all 15 alignments into one FASTA
python3 << EOF
import os

trim_dir = "${TRIM_DIR}"
data_dir = "${DATA_DIR}"
genes    = "${ALL_GENES[*]}".split()

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

samples = list(aln[genes[0]].keys())

with open(os.path.join(data_dir, "concat.fasta"), "w") as out:
    for sample in samples:
        seq = "".join(aln[gene].get(sample, "") for gene in genes)
        out.write(">" + sample + "\n")
        out.write(seq + "\n")

first_seq = "".join(aln[gene].get(samples[0], "") for gene in genes)
print(f"Concatenated: {len(samples)} taxa, {len(first_seq)} sites")

with open(os.path.join(data_dir, "gene_lengths.txt"), "w") as out:
    for gene in genes:
        seqs = list(aln[gene].values())
        out.write(f"{gene}\t{len(seqs[0])}\n")
EOF

# Step 2: Write partition file matching RevBayes 4-partition scheme
# PCG codon pos 1, 2, 3 pooled across all PCGs; rRNAs concatenated
PART_FILE="${DATA_DIR}/partitions.txt"
> "$PART_FILE"

POS1_SITES=""
POS2_SITES=""
POS3_SITES=""
RRNA_SITES=""

POS=1
while IFS=$'\t' read -r gene len; do
    END=$(( POS + len - 1 ))

    if [[ " ${PCG_GENES[*]} " == *" ${gene} "* ]]; then
        POS1_SITES="${POS1_SITES:+${POS1_SITES}, }${POS}-${END}\3"
        POS2_SITES="${POS2_SITES:+${POS2_SITES}, }$((POS+1))-${END}\3"
        POS3_SITES="${POS3_SITES:+${POS3_SITES}, }$((POS+2))-${END}\3"
    else
        RRNA_SITES="${RRNA_SITES:+${RRNA_SITES}, }${POS}-${END}"
    fi

    POS=$(( END + 1 ))
done < "${DATA_DIR}/gene_lengths.txt"

echo "DNA, PCG_pos1 = ${POS1_SITES}" >> "$PART_FILE"
echo "DNA, PCG_pos2 = ${POS2_SITES}" >> "$PART_FILE"
echo "DNA, PCG_pos3 = ${POS3_SITES}" >> "$PART_FILE"
echo "DNA, rRNA     = ${RRNA_SITES}" >> "$PART_FILE"

echo "Partition file written: ${PART_FILE}"
cat "$PART_FILE"

# Step 3: Run IQ-TREE
# - GTR+G4+I to match RevBayes model
# - 1000 ultrafast bootstrap replicates
# - 10 independent tree search runs
# - model applied per partition
iqtree2 \
    -s  "${DATA_DIR}/concat.fasta" \
    -p  "${PART_FILE}" \
    -m  GTR+G4+I \
    -B  1000 \
    --runs 10 \
    --redo \
    -T  ${SLURM_CPUS_PER_TASK} \
    --prefix "${OUT_DIR}/luciola_mito"

