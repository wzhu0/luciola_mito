#!/bin/bash
# Copy reads from hpc_exchange to "$HOME/luciola/mito/00_raw_reads"
# Usage: bash utils/cp_reads.sh

SOURCE_DIR="/mnt/netvolumes/srva229/bayes/hpc_exchange/shared/data_luciola"
DEST_DIR="$HOME/luciola/mito/00_raw_reads"
SAMPLE_LIST="$HOME/luciola/mito/utils/sample_list.txt"

mkdir -p "$DEST_DIR"

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line// }" ]] && continue

    year_dir="${SOURCE_DIR}/$(dirname "$line")"
    sample=$(basename "$line")

    r1=$(find "$year_dir" -maxdepth 1 -name "${sample}_1_val_1*.fq.gz" | head -1)
    r2=$(find "$year_dir" -maxdepth 1 -name "${sample}_2_val_2*.fq.gz" | head -1)

    # Check if the raw read files exist in hpc_exchange
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "MISSING  $sample"; continue
    fi

    dest_r1="$DEST_DIR/$(basename "$r1")"
    dest_r2="$DEST_DIR/$(basename "$r2")"

    # Skip if both exist in the destination dir
    if [[ -f "$dest_r1" && -f "$dest_r2" ]]; then
        echo "SKIP     $sample"; continue
    fi

    cp "$r1" "$dest_r1" && cp "$r2" "$dest_r2" \
        && echo "COPIED   $sample" \
        || { echo "ERROR    $sample"; exit 1; }

done < "$SAMPLE_LIST"

echo "Done."