#!/usr/bin/env python3

import argparse
import pandas as pd


def read_tsv(path, name):
    """
    Read TSV safely as strings.
    """
    print(f"Reading {name}: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    df = df.fillna("")
    print(f"  {name}: {len(df)} rows, {len(df.columns)} columns")
    return df


def prefix_columns(df, prefix, keep):
    """
    Prefix columns except join keys.
    """
    rename = {}
    for c in df.columns:
        if c not in keep:
            rename[c] = f"{prefix}_{c}"
    return df.rename(columns=rename)


def write_table(df, path):
    df.to_csv(path, sep="\t", index=False)
    print(f"Wrote: {path}")
    print(f"  Rows: {len(df)}")
    print(f"  Columns: {len(df.columns)}")


def main():
    parser = argparse.ArgumentParser(
        description="Link ENA read_run, read_experiment, sample, assembly and study metadata."
    )

    parser.add_argument("--runs", required=True)
    parser.add_argument("--experiments", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--assemblies", required=True)
    parser.add_argument("--studies", required=False)
    parser.add_argument("--outprefix", default="ena_linked")

    args = parser.parse_args()

    runs = read_tsv(args.runs, "runs")
    experiments = read_tsv(args.experiments, "experiments")
    samples = read_tsv(args.samples, "samples")
    assemblies = read_tsv(args.assemblies, "assemblies")

    studies = None
    if args.studies:
        studies = read_tsv(args.studies, "studies")

    print("\nAvailable key columns:")
    for name, df in [
        ("runs", runs),
        ("experiments", experiments),
        ("samples", samples),
        ("assemblies", assemblies),
        ("studies", studies),
    ]:
        if df is None:
            continue
        keys = [
            c for c in [
                "run_accession",
                "experiment_accession",
                "sample_accession",
                "secondary_sample_accession",
                "study_accession",
                "secondary_study_accession",
                "assembly_accession",
            ]
            if c in df.columns
        ]
        print(f"  {name}: {keys}")

    # ------------------------------------------------------------
    # 1. Link runs to experiments by experiment_accession
    # ------------------------------------------------------------

    linked = runs.copy()

    if (
        "experiment_accession" in linked.columns
        and "experiment_accession" in experiments.columns
    ):
        exp_prefixed = prefix_columns(
            experiments,
            prefix="experiment",
            keep=["experiment_accession"],
        )

        linked = linked.merge(
            exp_prefixed,
            on="experiment_accession",
            how="left",
        )

        print("\nLinked runs to experiments using experiment_accession")
    else:
        print("\nWARNING: Cannot link runs to experiments; experiment_accession missing")

    # ------------------------------------------------------------
    # 2. Link to samples by sample_accession
    # ------------------------------------------------------------

    if (
        "sample_accession" in linked.columns
        and "sample_accession" in samples.columns
    ):
        sample_prefixed = prefix_columns(
            samples,
            prefix="sample",
            keep=["sample_accession"],
        )

        linked = linked.merge(
            sample_prefixed,
            on="sample_accession",
            how="left",
        )

        print("Linked runs/experiments to samples using sample_accession")
    else:
        print("WARNING: Cannot link to samples; sample_accession missing")

    # ------------------------------------------------------------
    # 3. Link to studies by study_accession
    # ------------------------------------------------------------

    if studies is not None:
        if (
            "study_accession" in linked.columns
            and "study_accession" in studies.columns
        ):
            study_prefixed = prefix_columns(
                studies,
                prefix="study",
                keep=["study_accession"],
            )

            linked = linked.merge(
                study_prefixed,
                on="study_accession",
                how="left",
            )

            print("Linked to studies using study_accession")
        else:
            print("WARNING: Cannot link to studies; study_accession missing")

    # ------------------------------------------------------------
    # 4. Save full run-centered table
    # ------------------------------------------------------------

    write_table(
        linked,
        f"{args.outprefix}_runs_experiments_samples_studies.tsv",
    )

    # ------------------------------------------------------------
    # 5. Build assembly-centered table linked to samples
    # ------------------------------------------------------------

    assemblies_linked = assemblies.copy()

    if (
        "sample_accession" in assemblies_linked.columns
        and "sample_accession" in samples.columns
    ):
        sample_prefixed = prefix_columns(
            samples,
            prefix="sample",
            keep=["sample_accession"],
        )

        assemblies_linked = assemblies_linked.merge(
            sample_prefixed,
            on="sample_accession",
            how="left",
        )

        print("Linked assemblies to samples using sample_accession")
    else:
        print("WARNING: Cannot link assemblies to samples; sample_accession missing")

    if studies is not None:
        if (
            "study_accession" in assemblies_linked.columns
            and "study_accession" in studies.columns
        ):
            study_prefixed = prefix_columns(
                studies,
                prefix="study",
                keep=["study_accession"],
            )

            assemblies_linked = assemblies_linked.merge(
                study_prefixed,
                on="study_accession",
                how="left",
            )

            print("Linked assemblies to studies using study_accession")
        else:
            print("WARNING: Cannot link assemblies to studies; study_accession missing")

    write_table(
        assemblies_linked,
        f"{args.outprefix}_assemblies_samples_studies.tsv",
    )

    # ------------------------------------------------------------
    # 6. Build run-to-assembly links using sample_accession
    # ------------------------------------------------------------

    if (
        "sample_accession" in runs.columns
        and "sample_accession" in assemblies.columns
    ):
        run_to_assembly = runs.merge(
            assemblies,
            on="sample_accession",
            how="left",
            suffixes=("_run", "_assembly"),
        )

        write_table(
            run_to_assembly,
            f"{args.outprefix}_run_to_assembly_by_sample.tsv",
        )

        matched = run_to_assembly[
            run_to_assembly.filter(regex="assembly_accession").columns[0]
        ].astype(str).ne("").sum() if any(
            c.startswith("assembly_accession") for c in run_to_assembly.columns
        ) else 0

        print(f"Run-to-assembly rows with assembly accession: {matched}")

    else:
        print("WARNING: Cannot create run-to-assembly table; sample_accession missing")

    # ------------------------------------------------------------
    # 7. Summary tables
    # ------------------------------------------------------------

    summary = []

    def add_summary(name, df):
        summary.append(
            {
                "table": name,
                "rows": len(df),
                "columns": len(df.columns),
            }
        )

    add_summary("runs", runs)
    add_summary("experiments", experiments)
    add_summary("samples", samples)
    add_summary("assemblies", assemblies)
    if studies is not None:
        add_summary("studies", studies)
    add_summary("linked_runs_experiments_samples_studies", linked)
    add_summary("linked_assemblies_samples_studies", assemblies_linked)

    summary_df = pd.DataFrame(summary)
    write_table(summary_df, f"{args.outprefix}_summary.tsv")

    print("\nDone.")


if __name__ == "__main__":
    main()
