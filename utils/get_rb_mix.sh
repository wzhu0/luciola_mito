#!/usr/bin/env bash

PREFIX=$1
NUM_REPLICATES=$2
SUFFIX=$3

input_files=()
for i in $(seq 1 "$NUM_REPLICATES"); do
    input_files+=("${PREFIX}_run_${i}.${SUFFIX}")
done

for f in "${input_files[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: input file not found: $f" >&2
        exit 1
    fi
done

outfile="${PREFIX}.${SUFFIX}"
ncols=$(head -n 1 "${input_files[0]}" | awk '{print NF}')

# Write header: insert Replicate_ID as second column
head -n 1 "${input_files[0]}" | awk '{
    print $1 "\tReplicate_ID\t" substr($0, index($0, $2))
}' > "$outfile"

# Strip headers from all replicates into temp files
tmpdir=$(mktemp -d)
for i in $(seq 0 $(( NUM_REPLICATES - 1 ))); do
    tail -n +2 "${input_files[$i]}" > "${tmpdir}/rep_${i}.tmp"
done

# Interleave: paste all rep tmp files line by line, then split fields back out
paste $(for i in $(seq 0 $(( NUM_REPLICATES - 1 ))); do echo "${tmpdir}/rep_${i}.tmp"; done) \
| awk -v nrep="$NUM_REPLICATES" -v ncols="$ncols" 'BEGIN{OFS="\t"; gidx=0}{
    for (r = 0; r < nrep; r++) {
        start = r * ncols + 1      # 1-based field index for this replicate
        row = gidx "\t" r          # global index + replicate ID
        for (j = start + 1; j <= start + ncols - 1; j++) {
            row = row OFS $j       # skip field `start` (original Iteration column)
        }
        print row
        gidx++
    }
}' >> "$outfile"

rm -rf "$tmpdir"

