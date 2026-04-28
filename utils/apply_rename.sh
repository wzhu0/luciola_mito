#!/bin/bash
# Usage: bash utils/apply_rename.sh 05_revbayes/output/Luciola_mito_timetree_MCC.tree

TARGET="$1"
MAP="utils/Luciola_rename_map.tsv"

if [[ -z "$TARGET" ]]; then
    echo "Usage: $0 <target_file>"
    exit 1
fi

while IFS=$'\t' read -r original replacement; do
    original="${original//$'\r'/}"
    replacement="${replacement//$'\r'/}"
    perl -pi -e "s/\Q${original}\E/${replacement}/g" "$TARGET"
done < "$MAP"

echo "Done."