#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
import sys

root = Path(".")

input_files = sorted(root.rglob("*.fixed.gbk"))

if not input_files:
    print("No *.fixed.gbk files found.")
    sys.exit(1)

summary_file = "removed_no_gene_records.tsv"

with open(summary_file, "w") as summary:
    summary.write("file\tkept_records\tremoved_records\tremoved_record_ids\n")

    for gbk in input_files:
        out_gbk = gbk.with_name(gbk.stem.replace(".fixed", "") + ".clean.gbk")

        kept = []
        removed = []

        print(f"Cleaning: {gbk}")

        for record in SeqIO.parse(str(gbk), "genbank"):
            cds_count = sum(1 for f in record.features if f.type == "CDS")

            if cds_count > 0:
                kept.append(record)
            else:
                removed.append(record.id)

        if kept:
            SeqIO.write(kept, str(out_gbk), "genbank")
            print(f"  Kept records: {len(kept)}")
            print(f"  Removed records without CDS: {len(removed)}")
            print(f"  Output: {out_gbk}")
        else:
            print(f"  WARNING: no records with CDS found in {gbk}, no clean file written")

        summary.write(
            f"{gbk}\t{len(kept)}\t{len(removed)}\t{','.join(removed)}\n"
        )

print(f"Done. Summary written to: {summary_file}")
