#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge GTDB-Tk bacterial/archaeal summary files with MAG metadata."
    )
    parser.add_argument(
        "--metadata",
        default="work/mag_taxonomy/mag_taxonomy_input_metadata.tsv",
        help="Metadata TSV from prepare_mag_taxonomy_inputs.py.",
    )
    parser.add_argument(
        "--gtdbtk-dir",
        default="work/mag_taxonomy/gtdbtk_out",
        help="GTDB-Tk output directory.",
    )
    parser.add_argument(
        "--prefix",
        default="canmag",
        help="GTDB-Tk output prefix. Default: canmag",
    )
    parser.add_argument(
        "--output",
        default="work/mag_taxonomy/mag_taxonomy_gtdbtk_summary.tsv",
        help="Merged output TSV.",
    )
    return parser.parse_args()


def read_metadata(path: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Empty or malformed metadata TSV: {path}")
        rows = {row["genome_id"]: row for row in reader if row.get("genome_id")}
        return rows, list(reader.fieldnames)


def read_gtdb_summary(path: Path, domain: str) -> tuple[list[dict[str, str]], list[str]]:
    if not path.exists():
        return [], []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return [], []
        rows = []
        for row in reader:
            row["gtdb_domain_summary"] = domain
            rows.append(row)
        return rows, ["gtdb_domain_summary"] + list(reader.fieldnames)


def main() -> None:
    args = parse_args()
    metadata, metadata_fields = read_metadata(Path(args.metadata))
    gtdb_dir = Path(args.gtdbtk_dir)
    rows: list[dict[str, str]] = []
    gtdb_fields: list[str] = []

    for domain, suffix in (("bacteria", "bac120"), ("archaea", "ar53")):
        summary_path = gtdb_dir / f"{args.prefix}.{suffix}.summary.tsv"
        domain_rows, domain_fields = read_gtdb_summary(summary_path, domain)
        rows.extend(domain_rows)
        for field in domain_fields:
            if field not in gtdb_fields:
                gtdb_fields.append(field)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = metadata_fields + [field for field in gtdb_fields if field not in {"user_genome"}]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: r.get("user_genome", "")):
            genome_id = row.get("user_genome", "")
            out = dict(metadata.get(genome_id, {"genome_id": genome_id}))
            for field in gtdb_fields:
                if field != "user_genome":
                    out[field] = row.get(field, "")
            writer.writerow(out)

    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
