#!/usr/bin/env python3

import argparse
import csv
import io
import json
import re
import time
from collections import Counter

import requests


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"


def read_bioprojects(path):
    ids = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            ids.append(line)

    return ids


def clean_text(x):
    if x is None:
        return ""

    return " ".join(str(x).replace("\n", " ").replace("\r", " ").split())


def simplify_value(value):
    if value is None:
        return ""

    if isinstance(value, str):
        return clean_text(value)

    if isinstance(value, (int, float)):
        return value

    return json.dumps(value, ensure_ascii=False)


def request_get(url, params=None, retries=5, sleep=1.0):
    last_error = ""

    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, params=params, timeout=90)

            if r.status_code == 200:
                return r

            last_error = f"HTTP {r.status_code}: {r.text[:300]}"

        except Exception as e:
            last_error = str(e)

        wait = sleep * attempt
        print(f"Request failed attempt {attempt}/{retries}: {last_error}")
        print(f"Retrying in {wait:.1f} seconds...")
        time.sleep(wait)

    raise RuntimeError(f"Request failed after {retries} tries: {url} {params}")


def eutils_params(args, extra):
    params = dict(extra)
    params["tool"] = "bioproject_info_with_datatype"
    params["email"] = args.email

    if args.api_key:
        params["api_key"] = args.api_key

    return params


def search_bioproject_uid(bioproject_id, args):
    """
    Convert BioProject accession, e.g. PRJNA123456, to BioProject UID.
    """

    search_terms = [
        f"{bioproject_id}[Project Accession]",
        f"{bioproject_id}[Accession]",
        f"{bioproject_id}[All Fields]",
    ]

    for term in search_terms:
        params = eutils_params(args, {
            "db": "bioproject",
            "term": term,
            "retmode": "json",
            "retmax": 10,
        })

        r = request_get(f"{EUTILS}/esearch.fcgi", params=params)
        data = r.json()

        ids = data.get("esearchresult", {}).get("idlist", [])

        if ids:
            return ids[0]

        time.sleep(args.sleep)

    return ""


def get_bioproject_summary(uid, args):
    """
    Get BioProject metadata from NCBI ESummary.
    """

    params = eutils_params(args, {
        "db": "bioproject",
        "id": uid,
        "retmode": "json",
    })

    r = request_get(f"{EUTILS}/esummary.fcgi", params=params)
    data = r.json()

    result = data.get("result", {})

    if uid not in result:
        return {}

    return result[uid]


def fetch_sra_runinfo_for_project(bioproject_id):
    """
    Get SRA RunInfo for a BioProject.

    Useful columns include:
    Run, BioProject, BioSample, LibraryStrategy, LibrarySource,
    LibrarySelection, Platform, ScientificName, TaxID, etc.
    """

    params = {
        "acc": bioproject_id
    }

    r = request_get(RUNINFO_URL, params=params)
    text = r.text.strip()

    if not text:
        return []

    lower = text.lower()

    if lower.startswith("error") or lower.startswith("no data"):
        return []

    reader = csv.DictReader(io.StringIO(text))
    rows = []

    for row in reader:
        if not any(row.values()):
            continue

        clean_row = {
            k: clean_text(v)
            for k, v in row.items()
            if k is not None
        }

        rows.append(clean_row)

    return rows


def normalize(x):
    """
    Normalize library strategy/source strings.

    Examples:
    RNA-Seq -> RNASEQ
    WGS -> WGS
    Bisulfite-Seq -> BISULFITESEQ
    """

    if x is None:
        return ""

    return re.sub(r"[^A-Za-z0-9]+", "", str(x).upper())


def classify_sra_run(row):
    """
    Classify one SRA run into a broad data type.
    """

    strategy = normalize(row.get("LibraryStrategy", ""))
    source = normalize(row.get("LibrarySource", ""))
    selection = normalize(row.get("LibrarySelection", ""))

    combined = " ".join([strategy, source, selection])

    # RNA / expression
    if (
        "RNA" in combined
        or "TRANSCRIPT" in combined
        or strategy in {
            "RNASEQ",
            "MIRNASEQ",
            "NCRNASEQ",
            "FLCDNA",
            "EST",
            "CAGE",
        }
    ):
        return "RNA-seq / expression"

    # Genome sequencing
    if strategy in {
        "WGS",
        "WGA",
        "CLONE",
        "POOLCLONE",
        "FINISHING",
        "SYNTHETICLONGREAD",
    }:
        return "Genome"

    # Metagenome
    if (
        "META" in combined
        or strategy in {
            "METAGENOMIC",
            "METATRANSCRIPTOMIC",
        }
    ):
        return "Metagenome"

    # Amplicon
    if strategy in {
        "AMPLICON",
        "TARGETEDLOCUS",
    }:
        return "Amplicon"

    # Epigenomics / regulation
    if strategy in {
        "CHIPSEQ",
        "ATACSEQ",
        "DNASEHYPERSENSITIVITY",
        "MNASESEQ",
        "BISULFITESEQ",
        "METHYLSEQ",
        "HIC",
        "CHIA-PET",
    }:
        return "Epigenomics / regulation"

    # Exome / capture / targeted DNA
    if strategy in {
        "WXS",
        "TARGETEDCAPTURE",
    }:
        return "Targeted / exome"

    # Tn-seq / insertion etc.
    if strategy in {
        "TNSEQ",
        "TRASEQ",
        "BARCODEDSELECTED",
    }:
        return "Functional genomics"

    if strategy:
        return "Other"

    return "Unknown"


def summarize_sra_datatypes(run_rows):
    """
    Summarize SRA runs for one BioProject.
    """

    if not run_rows:
        return {
            "total_sra_runs": 0,
            "unique_biosamples": 0,
            "data_type_primary": "Unknown",
            "data_type_all": "Unknown",
            "data_type_counts": "",
            "library_strategy_counts": "",
            "library_source_counts": "",
            "platform_counts": "",
            "is_mixed_data_type": "no",
        }

    datatype_counts = Counter()
    strategy_counts = Counter()
    source_counts = Counter()
    platform_counts = Counter()
    biosamples = set()

    for row in run_rows:
        datatype = classify_sra_run(row)
        datatype_counts[datatype] += 1

        strategy = row.get("LibraryStrategy", "")
        source = row.get("LibrarySource", "")
        platform = row.get("Platform", "")
        biosample = row.get("BioSample", "")

        if strategy:
            strategy_counts[strategy] += 1

        if source:
            source_counts[source] += 1

        if platform:
            platform_counts[platform] += 1

        if biosample:
            biosamples.add(biosample)

    primary_datatype = datatype_counts.most_common(1)[0][0]

    data_type_all = "; ".join(
        datatype
        for datatype, count in datatype_counts.most_common()
    )

    data_type_counts = "; ".join(
        f"{datatype}:{count}"
        for datatype, count in datatype_counts.most_common()
    )

    library_strategy_counts = "; ".join(
        f"{strategy}:{count}"
        for strategy, count in strategy_counts.most_common()
    )

    library_source_counts = "; ".join(
        f"{source}:{count}"
        for source, count in source_counts.most_common()
    )

    platform_counts = "; ".join(
        f"{platform}:{count}"
        for platform, count in platform_counts.most_common()
    )

    is_mixed = "yes" if len(datatype_counts) > 1 else "no"

    return {
        "total_sra_runs": len(run_rows),
        "unique_biosamples": len(biosamples),
        "data_type_primary": primary_datatype,
        "data_type_all": data_type_all,
        "data_type_counts": data_type_counts,
        "library_strategy_counts": library_strategy_counts,
        "library_source_counts": library_source_counts,
        "platform_counts": platform_counts,
        "is_mixed_data_type": is_mixed,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Get BioProject information and infer data type from linked SRA runs."
    )

    parser.add_argument(
        "--projects",
        required=True,
        help="Text file with one BioProject ID per line."
    )

    parser.add_argument(
        "--email",
        required=True,
        help="Your email address. Required by NCBI."
    )

    parser.add_argument(
        "--api-key",
        default="",
        help="Optional NCBI API key."
    )

    parser.add_argument(
        "--output",
        default="bioproject_info_with_datatype.csv",
        help="Output CSV file."
    )

    parser.add_argument(
        "--sleep",
        type=float,
        default=0.34,
        help="Delay between NCBI requests. Use 0.34 without API key, 0.11 with API key."
    )

    args = parser.parse_args()

    bioprojects = read_bioprojects(args.projects)

    if not bioprojects:
        raise SystemExit("No BioProject IDs found in input file.")

    rows = []

    for i, bioproject_id in enumerate(bioprojects, start=1):
        print(f"[{i}/{len(bioprojects)}] Processing {bioproject_id}")

        row = {
            "query_bioproject": bioproject_id,
            "found": "no",
            "bioproject_uid": "",
        }

        # BioProject summary
        uid = search_bioproject_uid(bioproject_id, args)

        if uid:
            print(f"  BioProject UID: {uid}")

            row["found"] = "yes"
            row["bioproject_uid"] = uid

            summary = get_bioproject_summary(uid, args)

            for key, value in summary.items():
                row[key] = simplify_value(value)

        else:
            print(f"  BioProject not found in NCBI BioProject database")

        time.sleep(args.sleep)

        # SRA RunInfo and data type summary
        print(f"  Fetching SRA RunInfo")

        run_rows = fetch_sra_runinfo_for_project(bioproject_id)
        sra_summary = summarize_sra_datatypes(run_rows)

        row.update(sra_summary)

        print(f"  SRA runs: {row['total_sra_runs']}")
        print(f"  Data type: {row['data_type_all']}")

        rows.append(row)

        time.sleep(args.sleep)

    # Preferred column order
    preferred_cols = [
        "query_bioproject",
        "found",
        "bioproject_uid",
        "project_acc",
        "project_id",
        "project_name",
        "project_title",
        "project_description",
        "organism_name",
        "organism_taxid",
        "submission",
        "last_update",
        "total_sra_runs",
        "unique_biosamples",
        "data_type_primary",
        "data_type_all",
        "data_type_counts",
        "is_mixed_data_type",
        "library_strategy_counts",
        "library_source_counts",
        "platform_counts",
    ]

    all_cols = set()

    for row in rows:
        all_cols.update(row.keys())

    fieldnames = []

    for col in preferred_cols:
        if col in all_cols and col not in fieldnames:
            fieldnames.append(col)

    for col in sorted(all_cols):
        if col not in fieldnames:
            fieldnames.append(col)

    with open(args.output, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print()
    print("Done.")
    print(f"Wrote: {args.output}")
    print(f"BioProjects processed: {len(rows)}")


if __name__ == "__main__":
    main()
