#!/usr/bin/env python3

from pathlib import Path
import re
import sys

root = Path(".")

def safe_prefix(name):
    name = re.sub(r"[^A-Za-z0-9_]", "_", name)
    if not name:
        name = "SEQ"
    return name[:10]

def count_origin_sequence(record_lines):
    """
    Count actual nucleotide/amino acid letters in the ORIGIN section.
    This avoids trusting the broken LOCUS line.
    """
    in_origin = False
    length = 0

    for line in record_lines:
        if line.startswith("ORIGIN"):
            in_origin = True
            continue

        if in_origin:
            if line.startswith("//"):
                break

            # Remove numbers, spaces, punctuation; count only letters
            seq = re.sub(r"[^A-Za-z]", "", line)
            length += len(seq)

    return length

def get_date_from_locus(line):
    m = re.search(r"\d{2}-[A-Z]{3}-\d{4}", line)
    if m:
        return m.group(0)
    return "01-JAN-2000"

def get_topology_from_locus(line):
    if "circular" in line.lower():
        return "circular"
    return "linear"

def fix_file(gbk_path):
    prefix = safe_prefix(gbk_path.stem)
    out_path = gbk_path.with_name(gbk_path.stem + ".fixed.gbk")

    with open(gbk_path, "r", errors="replace") as f:
        lines = f.readlines()

    records = []
    current = []

    for line in lines:
        if line.startswith("LOCUS") and current:
            records.append(current)
            current = []

        current.append(line)

        if line.startswith("//"):
            records.append(current)
            current = []

    if current:
        records.append(current)

    fixed_records = []
    record_number = 0

    for record in records:
        locus_indices = [i for i, line in enumerate(record) if line.startswith("LOCUS")]

        if not locus_indices:
            fixed_records.extend(record)
            continue

        record_number += 1
        locus_index = locus_indices[0]
        old_locus = record[locus_index]

        seq_len = count_origin_sequence(record)

        if seq_len == 0:
            print(f"WARNING: {gbk_path} record {record_number} has sequence length 0", file=sys.stderr)

        date = get_date_from_locus(old_locus)
        topology = get_topology_from_locus(old_locus)

        # LOCUS ID must be short. Max 16 characters is safest.
        new_id = f"{prefix}_{record_number:05d}"

        new_locus = f"LOCUS       {new_id:<16} {seq_len:>11} bp    DNA     {topology:<8} BCT {date}\n"

        record[locus_index] = new_locus
        fixed_records.extend(record)

    with open(out_path, "w") as out:
        out.writelines(fixed_records)

    print(f"Fixed: {gbk_path} -> {out_path}")

gbk_files = sorted(
    p for p in root.rglob("*.gbk")
    if not p.name.endswith(".fixed.gbk")
)

if not gbk_files:
    print("No .gbk files found.")
    sys.exit(0)

for gbk in gbk_files:
    fix_file(gbk)

print("Done.")
