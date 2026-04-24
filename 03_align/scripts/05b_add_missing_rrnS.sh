#!/bin/bash
# Add SwLa-LuSpp-22-f-2020_pooled_101 as all-? to trimmed rrnS alignment

TRIM_DIR="03_align/03_trimmed"
TARGET="SwLa-LuSpp-22-f-2020_pooled_101"
FA="${TRIM_DIR}/rrnS.fasta"

# check if already present
if grep -q "^>${TARGET}$" "$FA"; then
    echo "${TARGET} already in rrnS alignment, skipping"
    exit 0
fi

# calculate alignment length from first sequence (wrap-aware)
ALN_LEN=$(awk '/^>/{if(seq!="")exit} /^[^>]/{seq=seq$0} END{print length(seq)}' "$FA")
echo "rrnS alignment length: ${ALN_LEN}"

# generate ? string of that length
QSEQ=$(printf '%0.s?' $(seq 1 $ALN_LEN))

# append to fasta, wrapping at 60 chars
echo ">${TARGET}" >> "$FA"
echo "$QSEQ" | fold -w 60 >> "$FA"

echo "Added ${TARGET} as ${ALN_LEN} ? to ${FA}"

