#!/usr/bin/env bash

set -u

THREADS=8

FAILED_LOG="antismash_failed.txt"
DONE_LOG="antismash_done.txt"

rm -f "$FAILED_LOG" "$DONE_LOG"

find . -name "*.clean.gbk" -print0 | while IFS= read -r -d '' gbk
do
    base=$(basename "$gbk" .clean.gbk)
    outdir="${base}_antismash"

    echo "========================================"
    echo "Running antiSMASH for: $gbk"
    echo "Output directory: $outdir"
    echo "========================================"

    if [ -d "$outdir" ]; then
        echo "Removing existing output: $outdir"
        rm -rf "$outdir"
    fi

    antismash \
        -c "$THREADS" \
        --output-dir "$outdir" \
        --output-basename "$base" \
        "$gbk"

    status=$?

    if [ "$status" -eq 0 ]; then
        echo "$gbk" >> "$DONE_LOG"
        echo "Finished: $gbk"
    else
        echo "$gbk" >> "$FAILED_LOG"
        echo "FAILED: $gbk"
    fi

done

echo "========================================"
echo "Batch antiSMASH finished."
echo "Done samples: $DONE_LOG"
echo "Failed samples: $FAILED_LOG"
echo "========================================"
