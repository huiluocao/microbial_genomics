#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict


BIOSAMPLE_RE = re.compile(r"^SAM[A-Z]*[0-9]+$")


def read_tsv(path):
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    return fieldnames, rows


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Wrote: {path}")
    print(f"Rows:  {len(rows)}")


def pick_col(cols, candidates, required=False, label="column"):
    for c in candidates:
        if c in cols:
            return c

    if required:
        raise ValueError(
            f"Could not find required {label}. Tried: {', '.join(candidates)}"
        )

    return None


def split_values(value):
    if value is None:
        return []

    value = str(value).strip()

    if not value:
        return []

    return [x.strip() for x in value.split(";") if x.strip()]


def looks_like_biosample(value):
    if value is None:
        return False

    value = str(value).strip()

    if not value:
        return False

    return bool(BIOSAMPLE_RE.match(value))


def unique_join(values):
    seen = set()
    out = []

    for value in values:
        for part in split_values(value):
            if part not in seen:
                seen.add(part)
                out.append(part)

    return ";".join(out)


def count_joined(value):
    return len(split_values(value))


def add_value(d, key, value):
    if not key:
        return

    for v in split_values(value):
        d[key].append(v)


def get_first_value(row, col):
    if not col:
        return ""
    return row.get(col, "").strip()


def sample_related_cols(cols):
    """
    Columns likely to contain sample/BioSample accessions.

    Important for your files:
      sample_accession can be SAMD00770387
      secondary_sample_accession can be DRS378499
    """
    preferred = [
        "biosample_accession",
        "bio_sample_accession",
        "sample_accession",
        "secondary_sample_accession",
        "secondary_accession",
        "accession",
    ]

    out = []

    for c in preferred:
        if c in cols and c not in out:
            out.append(c)

    for c in cols:
        lc = c.lower()
        if "sample" in lc and "accession" in lc and c not in out:
            out.append(c)

    return out


def accession_like_cols(cols):
    """
    Wider set of accession columns, used only for BioSample detection/mapping.
    """
    preferred = [
        "biosample_accession",
        "bio_sample_accession",
        "sample_accession",
        "secondary_sample_accession",
        "secondary_accession",
        "accession",
        "run_accession",
        "experiment_accession",
        "study_accession",
        "secondary_study_accession",
        "assembly_accession",
        "analysis_accession",
    ]

    out = []

    for c in preferred:
        if c in cols and c not in out:
            out.append(c)

    for c in cols:
        if "accession" in c.lower() and c not in out:
            out.append(c)

    return out


def detect_biosample_from_row(row, cols, sample_to_biosample=None):
    """
    Detect BioSample accession from a metadata row.

    Correctly handles:
      sample_accession = SAMD00770387
      secondary_sample_accession = DRS378499

    Strategy:
      1. Look in sample-related columns for SAMN/SAMEA/SAMD accession.
      2. Look in all accession-like columns for SAMN/SAMEA/SAMD accession.
      3. If no SAM accession is present, map ERS/SRS/DRS/etc to BioSample
         using sample_to_biosample.
    """
    if sample_to_biosample is None:
        sample_to_biosample = {}

    # First, check sample-related columns.
    for c in sample_related_cols(cols):
        for v in split_values(row.get(c, "")):
            if looks_like_biosample(v):
                return v

    # Then, check all accession-like columns.
    for c in accession_like_cols(cols):
        for v in split_values(row.get(c, "")):
            if looks_like_biosample(v):
                return v

    # Finally, try mapping non-SAM sample IDs to BioSample.
    for c in sample_related_cols(cols):
        for v in split_values(row.get(c, "")):
            if v in sample_to_biosample:
                return sample_to_biosample[v]

    for c in accession_like_cols(cols):
        for v in split_values(row.get(c, "")):
            if v in sample_to_biosample:
                return sample_to_biosample[v]

    return ""


def get_sample_related_values(row, cols):
    values = []

    for c in sample_related_cols(cols):
        for v in split_values(row.get(c, "")):
            values.append(v)

    return values


def print_columns(label, cols):
    print(f"\n{label} columns:")
    for i, c in enumerate(cols, start=1):
        print(f"  {i}\t{c}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create compact ENA metadata summary keyed by BioSample accession. "
            "Handles cases where sample_accession itself is SAMN/SAMEA/SAMD."
        )
    )

    parser.add_argument("--runs", required=True, help="ENA read_run TSV file")
    parser.add_argument("--assemblies", required=True, help="ENA assembly TSV file")
    parser.add_argument("--samples", required=True, help="ENA sample TSV file")
    parser.add_argument("--outprefix", required=True, help="Output prefix")

    args = parser.parse_args()

    print("Reading input files...")

    run_cols, runs = read_tsv(args.runs)
    asm_cols, assemblies = read_tsv(args.assemblies)
    sample_cols, samples = read_tsv(args.samples)

    print(f"Runs:       {len(runs)}")
    print(f"Assemblies: {len(assemblies)}")
    print(f"Samples:    {len(samples)}")

    print_columns("Run file", run_cols)
    print_columns("Assembly file", asm_cols)
    print_columns("Sample file", sample_cols)

    # ------------------------------------------------------------------
    # Detect useful run columns
    # ------------------------------------------------------------------

    run_accession_col = pick_col(
        run_cols,
        ["run_accession", "accession"],
        required=True,
        label="run accession column",
    )

    run_experiment_col = pick_col(
        run_cols,
        ["experiment_accession"],
        required=False,
    )

    run_study_col = pick_col(
        run_cols,
        ["study_accession", "secondary_study_accession"],
        required=False,
    )

    run_tax_id_col = pick_col(
        run_cols,
        ["tax_id", "taxon_id"],
        required=False,
    )

    run_tax_name_col = pick_col(
        run_cols,
        ["scientific_name", "taxonomic_name", "organism", "species"],
        required=False,
    )

    run_platform_col = pick_col(
        run_cols,
        ["instrument_platform"],
        required=False,
    )

    run_model_col = pick_col(
        run_cols,
        ["instrument_model"],
        required=False,
    )

    run_library_strategy_col = pick_col(
        run_cols,
        ["library_strategy"],
        required=False,
    )

    run_library_layout_col = pick_col(
        run_cols,
        ["library_layout"],
        required=False,
    )

    run_fastq_col = pick_col(
        run_cols,
        ["fastq_ftp"],
        required=False,
    )

    # ------------------------------------------------------------------
    # Detect useful assembly columns
    # ------------------------------------------------------------------

    asm_accession_col = pick_col(
        asm_cols,
        ["assembly_accession", "accession", "analysis_accession"],
        required=True,
        label="assembly accession column",
    )

    asm_study_col = pick_col(
        asm_cols,
        ["study_accession", "secondary_study_accession"],
        required=False,
    )

    asm_tax_id_col = pick_col(
        asm_cols,
        ["tax_id", "taxon_id"],
        required=False,
    )

    asm_tax_name_col = pick_col(
        asm_cols,
        ["scientific_name", "taxonomic_name", "organism", "species"],
        required=False,
    )

    asm_name_col = pick_col(
        asm_cols,
        ["assembly_name", "name"],
        required=False,
    )

    asm_level_col = pick_col(
        asm_cols,
        ["assembly_level"],
        required=False,
    )

    asm_description_col = pick_col(
        asm_cols,
        ["description"],
        required=False,
    )

    # ------------------------------------------------------------------
    # Detect useful sample columns
    # ------------------------------------------------------------------

    sample_accession_col = pick_col(
        sample_cols,
        ["sample_accession", "accession"],
        required=True,
        label="sample accession column",
    )

    sample_study_col = pick_col(
        sample_cols,
        ["study_accession", "secondary_study_accession"],
        required=False,
    )

    sample_tax_id_col = pick_col(
        sample_cols,
        ["tax_id", "taxon_id"],
        required=False,
    )

    sample_tax_name_col = pick_col(
        sample_cols,
        ["scientific_name", "taxonomic_name", "organism", "species", "description"],
        required=False,
    )

    print("\nDetected main columns:")
    print(f"  Run accession column:          {run_accession_col}")
    print(f"  Run experiment column:         {run_experiment_col}")
    print(f"  Run study column:              {run_study_col}")
    print(f"  Run taxonomy ID column:        {run_tax_id_col}")
    print(f"  Run taxonomy name column:      {run_tax_name_col}")
    print(f"  Assembly accession column:     {asm_accession_col}")
    print(f"  Assembly study column:         {asm_study_col}")
    print(f"  Assembly taxonomy ID column:   {asm_tax_id_col}")
    print(f"  Assembly taxonomy name column: {asm_tax_name_col}")
    print(f"  Assembly name column:          {asm_name_col}")
    print(f"  Assembly level column:         {asm_level_col}")
    print(f"  Sample accession column:       {sample_accession_col}")
    print(f"  Sample taxonomy/name column:   {sample_tax_name_col}")

    # ------------------------------------------------------------------
    # Build sample accession -> BioSample map from sample file
    # ------------------------------------------------------------------

    sample_to_biosample = {}
    all_biosamples = set()
    samples_without_biosample = 0

    for row in samples:
        biosample_acc = detect_biosample_from_row(row, sample_cols)

        if not biosample_acc:
            samples_without_biosample += 1
            continue

        all_biosamples.add(biosample_acc)

        # Map all sample-related values in the row to this BioSample.
        # This allows DRS/ERS/SRS -> SAMD/SAMEA/SAMN if both occur.
        for v in get_sample_related_values(row, sample_cols):
            sample_to_biosample[v] = biosample_acc

    print(f"\nBioSamples detected from sample file:       {len(all_biosamples)}")
    print(f"Sample/SRA-to-BioSample mappings created:  {len(sample_to_biosample)}")
    print(f"Sample rows without BioSample detected:    {samples_without_biosample}")

    # ------------------------------------------------------------------
    # Aggregation dictionaries keyed by BioSample accession
    # ------------------------------------------------------------------

    biosample_to_sample_accessions = defaultdict(list)
    biosample_to_run_accessions = defaultdict(list)
    biosample_to_experiment_accessions = defaultdict(list)
    biosample_to_assembly_accessions = defaultdict(list)
    biosample_to_study_accessions = defaultdict(list)
    biosample_to_tax_ids = defaultdict(list)
    biosample_to_taxonomy_names = defaultdict(list)
    biosample_to_assembly_names = defaultdict(list)
    biosample_to_assembly_levels = defaultdict(list)
    biosample_to_assembly_descriptions = defaultdict(list)
    biosample_to_platforms = defaultdict(list)
    biosample_to_instrument_models = defaultdict(list)
    biosample_to_library_strategies = defaultdict(list)
    biosample_to_library_layouts = defaultdict(list)
    biosample_to_fastq_ftps = defaultdict(list)

    unmapped_runs = []
    unmapped_assemblies = []

    # ------------------------------------------------------------------
    # Add sample metadata
    # ------------------------------------------------------------------

    for row in samples:
        biosample_acc = detect_biosample_from_row(row, sample_cols, sample_to_biosample)

        if not biosample_acc:
            continue

        all_biosamples.add(biosample_acc)

        for v in get_sample_related_values(row, sample_cols):
            add_value(biosample_to_sample_accessions, biosample_acc, v)

        if sample_study_col:
            add_value(
                biosample_to_study_accessions,
                biosample_acc,
                row.get(sample_study_col, ""),
            )

        if sample_tax_id_col:
            add_value(
                biosample_to_tax_ids,
                biosample_acc,
                row.get(sample_tax_id_col, ""),
            )

        if sample_tax_name_col:
            add_value(
                biosample_to_taxonomy_names,
                biosample_acc,
                row.get(sample_tax_name_col, ""),
            )

    # ------------------------------------------------------------------
    # Add run metadata
    # ------------------------------------------------------------------

    runs_without_biosample = 0

    for row in runs:
        run_acc = get_first_value(row, run_accession_col)
        biosample_acc = detect_biosample_from_row(row, run_cols, sample_to_biosample)

        if not biosample_acc:
            runs_without_biosample += 1
            unmapped_runs.append(
                {
                    "run_accession": run_acc,
                    "sample_accession": row.get("sample_accession", ""),
                    "secondary_sample_accession": row.get("secondary_sample_accession", ""),
                    "study_accession": row.get("study_accession", ""),
                    "tax_id": row.get("tax_id", ""),
                    "scientific_name": row.get("scientific_name", ""),
                    "reason": "No SAMN/SAMEA/SAMD BioSample detected and no mapping found",
                }
            )
            continue

        all_biosamples.add(biosample_acc)

        add_value(biosample_to_run_accessions, biosample_acc, run_acc)

        for v in get_sample_related_values(row, run_cols):
            add_value(biosample_to_sample_accessions, biosample_acc, v)

        if run_experiment_col:
            add_value(
                biosample_to_experiment_accessions,
                biosample_acc,
                row.get(run_experiment_col, ""),
            )

        if run_study_col:
            add_value(
                biosample_to_study_accessions,
                biosample_acc,
                row.get(run_study_col, ""),
            )

        if run_tax_id_col:
            add_value(
                biosample_to_tax_ids,
                biosample_acc,
                row.get(run_tax_id_col, ""),
            )

        if run_tax_name_col:
            add_value(
                biosample_to_taxonomy_names,
                biosample_acc,
                row.get(run_tax_name_col, ""),
            )

        if run_platform_col:
            add_value(
                biosample_to_platforms,
                biosample_acc,
                row.get(run_platform_col, ""),
            )

        if run_model_col:
            add_value(
                biosample_to_instrument_models,
                biosample_acc,
                row.get(run_model_col, ""),
            )

        if run_library_strategy_col:
            add_value(
                biosample_to_library_strategies,
                biosample_acc,
                row.get(run_library_strategy_col, ""),
            )

        if run_library_layout_col:
            add_value(
                biosample_to_library_layouts,
                biosample_acc,
                row.get(run_library_layout_col, ""),
            )

        if run_fastq_col:
            add_value(
                biosample_to_fastq_ftps,
                biosample_acc,
                row.get(run_fastq_col, ""),
            )

    # ------------------------------------------------------------------
    # Add assembly metadata
    # ------------------------------------------------------------------

    assemblies_without_biosample = 0

    for row in assemblies:
        asm_acc = get_first_value(row, asm_accession_col)
        biosample_acc = detect_biosample_from_row(row, asm_cols, sample_to_biosample)

        if not biosample_acc:
            assemblies_without_biosample += 1
            unmapped_assemblies.append(
                {
                    "assembly_accession": asm_acc,
                    "sample_accession": row.get("sample_accession", ""),
                    "secondary_sample_accession": row.get("secondary_sample_accession", ""),
                    "study_accession": row.get("study_accession", ""),
                    "tax_id": row.get("tax_id", ""),
                    "scientific_name": row.get("scientific_name", ""),
                    "reason": "No SAMN/SAMEA/SAMD BioSample detected and no mapping found",
                }
            )
            continue

        all_biosamples.add(biosample_acc)

        add_value(biosample_to_assembly_accessions, biosample_acc, asm_acc)

        for v in get_sample_related_values(row, asm_cols):
            add_value(biosample_to_sample_accessions, biosample_acc, v)

        if asm_study_col:
            add_value(
                biosample_to_study_accessions,
                biosample_acc,
                row.get(asm_study_col, ""),
            )

        if asm_tax_id_col:
            add_value(
                biosample_to_tax_ids,
                biosample_acc,
                row.get(asm_tax_id_col, ""),
            )

        if asm_tax_name_col:
            add_value(
                biosample_to_taxonomy_names,
                biosample_acc,
                row.get(asm_tax_name_col, ""),
            )

        if asm_name_col:
            add_value(
                biosample_to_assembly_names,
                biosample_acc,
                row.get(asm_name_col, ""),
            )

        if asm_level_col:
            add_value(
                biosample_to_assembly_levels,
                biosample_acc,
                row.get(asm_level_col, ""),
            )

        if asm_description_col:
            add_value(
                biosample_to_assembly_descriptions,
                biosample_acc,
                row.get(asm_description_col, ""),
            )

    # ------------------------------------------------------------------
    # Build output rows
    # ------------------------------------------------------------------

    rows = []

    for biosample_acc in sorted(all_biosamples):
        sample_accessions = unique_join(
            biosample_to_sample_accessions.get(biosample_acc, [])
        )

        run_accessions = unique_join(
            biosample_to_run_accessions.get(biosample_acc, [])
        )

        experiment_accessions = unique_join(
            biosample_to_experiment_accessions.get(biosample_acc, [])
        )

        assembly_accessions = unique_join(
            biosample_to_assembly_accessions.get(biosample_acc, [])
        )

        study_accessions = unique_join(
            biosample_to_study_accessions.get(biosample_acc, [])
        )

        tax_ids = unique_join(
            biosample_to_tax_ids.get(biosample_acc, [])
        )

        taxonomy_names = unique_join(
            biosample_to_taxonomy_names.get(biosample_acc, [])
        )

        assembly_names = unique_join(
            biosample_to_assembly_names.get(biosample_acc, [])
        )

        assembly_levels = unique_join(
            biosample_to_assembly_levels.get(biosample_acc, [])
        )

        assembly_descriptions = unique_join(
            biosample_to_assembly_descriptions.get(biosample_acc, [])
        )

        instrument_platforms = unique_join(
            biosample_to_platforms.get(biosample_acc, [])
        )

        instrument_models = unique_join(
            biosample_to_instrument_models.get(biosample_acc, [])
        )

        library_strategies = unique_join(
            biosample_to_library_strategies.get(biosample_acc, [])
        )

        library_layouts = unique_join(
            biosample_to_library_layouts.get(biosample_acc, [])
        )

        fastq_ftps = unique_join(
            biosample_to_fastq_ftps.get(biosample_acc, [])
        )

        row = {
            "biosample_accession": biosample_acc,
            "sample_accessions": sample_accessions,
            "run_accessions": run_accessions,
            "experiment_accessions": experiment_accessions,
            "assembly_accessions": assembly_accessions,
            "study_accessions": study_accessions,
            "tax_ids": tax_ids,
            "taxonomy_names": taxonomy_names,
            "assembly_names": assembly_names,
            "assembly_levels": assembly_levels,
            "assembly_descriptions": assembly_descriptions,
            "instrument_platforms": instrument_platforms,
            "instrument_models": instrument_models,
            "library_strategies": library_strategies,
            "library_layouts": library_layouts,
            "fastq_ftps": fastq_ftps,
            "n_sample_accessions": count_joined(sample_accessions),
            "n_runs": count_joined(run_accessions),
            "n_experiments": count_joined(experiment_accessions),
            "n_assemblies": count_joined(assembly_accessions),
            "has_reads": "yes" if run_accessions else "no",
            "has_assembly": "yes" if assembly_accessions else "no",
            "has_reads_and_assembly": (
                "yes" if run_accessions and assembly_accessions else "no"
            ),
        }

        rows.append(row)

    # ------------------------------------------------------------------
    # Write main output
    # ------------------------------------------------------------------

    out_fields = [
        "biosample_accession",
        "sample_accessions",
        "run_accessions",
        "experiment_accessions",
        "assembly_accessions",
        "study_accessions",
        "tax_ids",
        "taxonomy_names",
        "assembly_names",
        "assembly_levels",
        "assembly_descriptions",
        "instrument_platforms",
        "instrument_models",
        "library_strategies",
        "library_layouts",
        "fastq_ftps",
        "n_sample_accessions",
        "n_runs",
        "n_experiments",
        "n_assemblies",
        "has_reads",
        "has_assembly",
        "has_reads_and_assembly",
    ]

    out_file = f"{args.outprefix}_biosample_summary.tsv"
    write_tsv(out_file, out_fields, rows)

    # ------------------------------------------------------------------
    # Write unmapped diagnostic files
    # ------------------------------------------------------------------

    unmapped_run_file = f"{args.outprefix}_unmapped_runs.tsv"
    unmapped_run_fields = [
        "run_accession",
        "sample_accession",
        "secondary_sample_accession",
        "study_accession",
        "tax_id",
        "scientific_name",
        "reason",
    ]

    write_tsv(unmapped_run_file, unmapped_run_fields, unmapped_runs)

    unmapped_asm_file = f"{args.outprefix}_unmapped_assemblies.tsv"
    unmapped_asm_fields = [
        "assembly_accession",
        "sample_accession",
        "secondary_sample_accession",
        "study_accession",
        "tax_id",
        "scientific_name",
        "reason",
    ]

    write_tsv(unmapped_asm_file, unmapped_asm_fields, unmapped_assemblies)

    # ------------------------------------------------------------------
    # Final report
    # ------------------------------------------------------------------

    biosamples_with_samples = sum(
        1 for r in rows if r["sample_accessions"]
    )

    biosamples_with_runs = sum(
        1 for r in rows if r["run_accessions"]
    )

    biosamples_with_assemblies = sum(
        1 for r in rows if r["assembly_accessions"]
    )

    biosamples_with_both = sum(
        1 for r in rows if r["run_accessions"] and r["assembly_accessions"]
    )

    print("\nSummary:")
    print(f"  BioSamples in output:                 {len(rows)}")
    print(f"  BioSamples with sample accessions:    {biosamples_with_samples}")
    print(f"  BioSamples with reads/runs:           {biosamples_with_runs}")
    print(f"  BioSamples with assemblies:           {biosamples_with_assemblies}")
    print(f"  BioSamples with both reads+assembly:  {biosamples_with_both}")
    print(f"  Samples without BioSample detected:   {samples_without_biosample}")
    print(f"  Runs without BioSample after mapping: {runs_without_biosample}")
    print(f"  Assemblies without BioSample mapping: {assemblies_without_biosample}")

    print("\nDone.")


if __name__ == "__main__":
    main()
