#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
##SBATCH --time=3-00:00:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/04_iqtree.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_iqtree

# Step 4: Concatenate trimmed alignments and run IQ-TREE
# Produces an unrooted ML tree with UFBoot support values
# Run after 03_align/scripts/03_trim.sh is complete
# Usage: sbatch 04_iqtree/scripts/run_iqtree.slurm

module load gnu/12
module load prebin/kry

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
with open(os.path.join(data_dir, "concat.fasta"), "w") as out:
    for sample in samples:
        seq = "".join(aln[gene].get(sample, "") for gene in genes)
        out.write(">" + sample + "\n")
        out.write(seq + "\n")

# Report
first_seq = "".join(aln[gene].get(samples[0], "") for gene in genes)
print(f"Concatenated: {len(samples)} taxa, {len(first_seq)} sites")

# Write per-gene lengths for partition file
with open(os.path.join(data_dir, "gene_lengths.txt"), "w") as out:
    for gene in genes:
        seqs = list(aln[gene].values())
        out.write(f"{gene}\t{len(seqs[0])}\n")
EOF

# Step 2: Write partition file for IQ-TREE 
PART_FILE="${DATA_DIR}/partitions.txt"
> "$PART_FILE"

POS=1
while IFS=$'\t' read -r gene len; do
    END=$(( POS + len - 1 ))

    # PCGs: three separate codon position partitions
    if [[ " ${PCG_GENES[*]} " == *" ${gene} "* ]]; then
        echo "DNA, ${gene}_pos1 = ${POS}-${END}\3"         >> "$PART_FILE"
        echo "DNA, ${gene}_pos2 = $((POS+1))-${END}\3"     >> "$PART_FILE"
        echo "DNA, ${gene}_pos3 = $((POS+2))-${END}\3"     >> "$PART_FILE"
    else
        echo "DNA, ${gene} = ${POS}-${END}"                 >> "$PART_FILE"
    fi

    POS=$(( END + 1 ))
done < "${DATA_DIR}/gene_lengths.txt"


# Step 3: Run IQ-TREE
iqtree2 \
    -s  "${DATA_DIR}/concat.fasta" \
    -p  "${PART_FILE}" \
    -m  MFP \
    -B  1000 \
    --redo \
    -T  ${SLURM_CPUS_PER_TASK} \
    --prefix "${OUT_DIR}/luciola_mito"
