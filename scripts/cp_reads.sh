SOURCE_DIR="/mnt/netvolumes/srva229/bayes/hpc_exchange/shared/data_luciola/2025"
DEST_DIR="$HOME/luciola/mito/data"
SAMPLE_LIST="$HOME/luciola/mito/utils/sample_list.txt"

while IFS= read -r line || [[ -n "$line" ]]; do
    # skip empty lines
    [[ -z "${line// }" ]] && continue
    
    sample_path="$SOURCE_DIR/${line#./}"
    sample_name=$(basename "$sample_path")
    subdir=$(dirname "$sample_path")
    
    # find R1 and R2 — with or without .fixed
    r1=$(find "$subdir" -maxdepth 1 -name "${sample_name}_1_val_1*.fq.gz" | head -1)
    r2=$(find "$subdir" -maxdepth 1 -name "${sample_name}_2_val_2*.fq.gz" | head -1)
    
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: reads not found for $sample_name, skipping"
        continue
    fi
    
    echo "Copying $sample_name ..."
    cp "$r1" "$r2" "$DEST_DIR/"

done < "$SAMPLE_LIST"

echo "Done."