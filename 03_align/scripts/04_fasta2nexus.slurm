#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/03_align_fasta2nexus.out
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J mito_fasta2nex

set -euo pipefail

module load gnu/12
module load prebin/kry

TRIM_DIR="03_align/03_trimmed"
NEX_DIR="05_revbayes/data"

mkdir -p "${NEX_DIR}" logs

PCG_GENES=(nad1 nad2 nad3 nad4 nad4l nad5 nad6 cox1 cox2 cox3 atp6 atp8 cob)
RRNA_GENES=(rrnS rrnL)
ALL_GENES=("${PCG_GENES[@]}" "${RRNA_GENES[@]}")

for gene in "${ALL_GENES[@]}"; do
    IN="${TRIM_DIR}/${gene}.fasta"
    OUT="${NEX_DIR}/${gene}.nex"

    if [[ ! -f "$IN" ]]; then
        echo "ERROR: ${IN} not found."
        exit 1
    fi

    NTAX=$(grep -c "^>" "$IN")

    # Get alignment length: read lines after first header until next header
    # Using python to avoid awk double-print issue
    NCHAR=$(python3 -c "
seq = ''
found = False
with open('${IN}') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if found: break
            found = True
        elif found:
            seq += line
print(len(seq))
")

    echo "Converting ${gene}: ${NTAX} taxa, ${NCHAR} chars"

    echo "#NEXUS"                               > "$OUT"
    echo ""                                    >> "$OUT"
    echo "Begin data;"                         >> "$OUT"
    echo "    Dimensions ntax=${NTAX} nchar=${NCHAR};" >> "$OUT"
    echo "    Format datatype=DNA missing=? gap=-;" >> "$OUT"
    echo "    Matrix"                          >> "$OUT"

    awk '
        /^>/ {
            if (name != "") print name "\t" seq
            name = substr($0, 2); seq = ""
        }
        !/^>/ { seq = seq $0 }
        END   { if (name != "") print name "\t" seq }
    ' "$IN" >> "$OUT"

    printf "    ;\nEnd;\n" >> "$OUT"

    echo "  -> ${OUT}"
done

echo ""
echo "=== NEXUS conversion summary ==="
for gene in "${ALL_GENES[@]}"; do
    NEX="${NEX_DIR}/${gene}.nex"
    [[ -f "$NEX" ]] && echo "  ${gene}.nex: OK" || echo "  ${gene}.nex: MISSING"
done
echo "Done. NEXUS files in: ${NEX_DIR}/"
