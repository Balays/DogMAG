#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


FASTA_SUFFIXES = (
    ".fa.gz",
    ".fna.gz",
    ".fasta.gz",
    ".fa",
    ".fna",
    ".fasta",
    ".fas",
    ".fsa",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a GTDB-Tk-style MAG taxonomy summary for a dRep95 "
            "dereplicated genome collection from an existing full GTDB-Tk "
            "summary table. This does not run GTDB-Tk."
        )
    )
    parser.add_argument(
        "--drep95-genome-dir",
        default="CanMAG_depletion_panels_20260721_indexed/dRep95/dereplicated_genomes",
        help="Directory containing dRep95 representative FASTAs.",
    )
    parser.add_argument(
        "--drep95-genome-info",
        default="CanMAG_depletion_panels_20260721_indexed/dRep95/data_tables/genomeInfo.csv",
        help="dRep genomeInfo.csv with completeness/contamination/length/N50.",
    )
    parser.add_argument(
        "--full-summary",
        default="work/dogfirst_mag_taxonomy_gtdbtk_20260719/mag_taxonomy_gtdbtk_summary.tsv",
        help="Full MAG GTDB-Tk merged summary TSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="work/dogmag_drep95_gtdbtk_from_full",
        help="Output directory.",
    )
    return parser.parse_args()


def is_fasta(path: Path) -> bool:
    return path.is_file() and any(path.name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def strip_fasta_suffixes(name: str) -> str:
    token = Path(name).name
    changed = True
    while changed:
        changed = False
        for suffix in FASTA_SUFFIXES:
            if token.endswith(suffix):
                token = token[: -len(suffix)]
                changed = True
                break
    return token


def clean_id(text: str) -> str:
    token = strip_fasta_suffixes(text)
    token = re.sub(r"[^A-Za-z0-9_]+", "_", token)
    token = re.sub(r"_+", "_", token).strip("_")
    return token


def bin_group(text: str) -> str:
    match = re.match(r"^(bin\d+)", clean_id(text))
    return match.group(1) if match else ""


def aliases_for_name(text: str) -> set[str]:
    raw = Path(text).name
    stripped = strip_fasta_suffixes(raw)
    aliases = {raw, stripped, clean_id(raw), clean_id(stripped)}

    # dRep95 names often look like bin10017__bin10017_mag_polished.fa.fa.
    # The full GTDB run usually contains only bin10017_mag_polished.
    if "__" in stripped:
        right = stripped.split("__", 1)[1]
        aliases.add(right)
        aliases.add(clean_id(right))

    group = bin_group(stripped)
    if group:
        aliases.add(group)

    return {alias for alias in aliases if alias}


def read_table(path: Path, delimiter: str) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Empty or malformed table: {path}")
        return list(reader), list(reader.fieldnames)


def read_drep_quality(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    rows, _ = read_table(path, ",")
    quality: dict[str, dict[str, str]] = {}
    for row in rows:
        genome = row.get("genome", "")
        for alias in aliases_for_name(genome):
            quality[alias] = row
    return quality


def add_full_summary_indexes(
    rows: list[dict[str, str]],
) -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]]]:
    exact: dict[str, dict[str, str]] = {}
    by_group: dict[str, list[dict[str, str]]] = {}

    for row in rows:
        values = [
            row.get("genome_id", ""),
            row.get("file_name", ""),
            row.get("genome_path", ""),
            row.get("source_basalt_quality_key", ""),
        ]
        for value in values:
            if not value:
                continue
            for alias in aliases_for_name(value):
                exact.setdefault(alias, row)

        group = bin_group(row.get("genome_id", "") or row.get("file_name", ""))
        if group:
            by_group.setdefault(group, []).append(row)

    return exact, by_group


def choose_match(
    fasta: Path,
    exact_index: dict[str, dict[str, str]],
    group_index: dict[str, list[dict[str, str]]],
) -> tuple[dict[str, str] | None, str, str]:
    for alias in aliases_for_name(fasta.name):
        row = exact_index.get(alias)
        if row:
            return row, "exact_full_gtdb_match", row.get("genome_id", "")

    group = bin_group(fasta.name)
    if group and group in group_index:
        rows = sorted(group_index[group], key=lambda row: row.get("genome_id", ""))
        return rows[0], "same_bin_group_inferred", rows[0].get("genome_id", "")

    return None, "missing_from_full_gtdb", ""


def main() -> None:
    args = parse_args()
    genome_dir = Path(args.drep95_genome_dir)
    genome_info = Path(args.drep95_genome_info)
    full_summary = Path(args.full_summary)
    output_dir = Path(args.output_dir)

    if not genome_dir.is_dir():
        raise FileNotFoundError(f"dRep95 genome directory not found: {genome_dir}")
    if not full_summary.is_file():
        raise FileNotFoundError(f"Full GTDB-Tk summary not found: {full_summary}")

    full_rows, full_fields = read_table(full_summary, "\t")
    exact_index, group_index = add_full_summary_indexes(full_rows)
    quality = read_drep_quality(genome_info)

    fastas = sorted(path.resolve() for path in genome_dir.iterdir() if is_fasta(path))
    output_rows: list[dict[str, str]] = []
    match_rows: list[dict[str, str]] = []
    counts: Counter[str] = Counter()

    output_fields = full_fields

    for fasta in fastas:
        matched, source, source_genome_id = choose_match(fasta, exact_index, group_index)
        counts[source] += 1

        genome_id = clean_id(fasta.name)
        q = {}
        for alias in aliases_for_name(fasta.name):
            if alias in quality:
                q = quality[alias]
                break

        if matched:
            out = {field: matched.get(field, "") for field in full_fields}
        else:
            out = {field: "" for field in full_fields}

        out["genome_id"] = genome_id
        out["genome_path"] = str(fasta)
        out["file_name"] = fasta.name
        out["genome_size"] = q.get("length", out.get("genome_size", ""))
        out["completeness"] = q.get("completeness", out.get("completeness", ""))
        out["contamination"] = q.get("contamination", out.get("contamination", ""))
        out["n50"] = q.get("N50", out.get("n50", ""))
        output_rows.append(out)

        match_rows.append(
            {
                "drep95_genome_id": genome_id,
                "drep95_file": fasta.name,
                "bin_group": bin_group(fasta.name),
                "taxonomy_source": source,
                "taxonomy_source_genome_id": source_genome_id,
                "classification": out.get("classification", ""),
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "mag_taxonomy_gtdbtk_summary.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    match_path = output_dir / "drep95_gtdbtk_match_report.tsv"
    with match_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "drep95_genome_id",
                "drep95_file",
                "bin_group",
                "taxonomy_source",
                "taxonomy_source_genome_id",
                "classification",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(match_rows)

    summary_path = output_dir / "summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"drep95_fastas\t{len(fastas)}\n")
        for key in ("exact_full_gtdb_match", "same_bin_group_inferred", "missing_from_full_gtdb"):
            handle.write(f"{key}\t{counts[key]}\n")
        handle.write(f"full_summary_rows\t{len(full_rows)}\n")
        handle.write(f"full_summary\t{full_summary}\n")
        handle.write(f"drep95_genome_dir\t{genome_dir}\n")
        handle.write(f"drep95_genome_info\t{genome_info}\n")

    print(f"dRep95 FASTAs: {len(fastas)}")
    print(f"Exact full GTDB matches: {counts['exact_full_gtdb_match']}")
    print(f"Same-bin inferred matches: {counts['same_bin_group_inferred']}")
    print(f"Missing from full GTDB: {counts['missing_from_full_gtdb']}")
    print(f"Wrote {output_path}")
    print(f"Wrote {match_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
