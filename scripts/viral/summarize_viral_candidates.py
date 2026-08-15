#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge geNomad viral summaries with the CanMAG viral contig context table "
            "and, optionally, CheckV quality summaries."
        )
    )
    parser.add_argument(
        "--context",
        required=True,
        help="viral_contig_context.tsv from prepare_viral_contig_inputs.py.",
    )
    parser.add_argument(
        "--genomad-virus-summary",
        required=True,
        help="geNomad *_virus_summary.tsv file.",
    )
    parser.add_argument(
        "--checkv-quality",
        default="",
        help="Optional CheckV quality_summary.tsv.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Merged output TSV.",
    )
    return parser.parse_args()


def read_table(path: Path, key: str) -> tuple[dict[str, dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Empty or malformed TSV: {path}")
        rows = {row[key]: row for row in reader if row.get(key)}
        return rows, list(reader.fieldnames)


def prefixed_fields(fields: list[str], prefix: str, exclude: set[str]) -> list[str]:
    return [f"{prefix}{field}" for field in fields if field not in exclude]


def get_context_row(context: dict[str, dict[str, str]], seq_name: str) -> dict[str, str]:
    row = context.get(seq_name)
    if row:
        return row
    if "|provirus_" in seq_name:
        parent = seq_name.split("|provirus_", 1)[0]
        return context.get(parent, {})
    return {}


def main() -> None:
    args = parse_args()
    context, context_fields = read_table(Path(args.context), "viral_contig_id")
    genomad, genomad_fields = read_table(Path(args.genomad_virus_summary), "seq_name")

    checkv: dict[str, dict[str, str]] = {}
    checkv_fields: list[str] = []
    if args.checkv_quality:
        checkv, checkv_fields = read_table(Path(args.checkv_quality), "contig_id")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = (
        ["seq_name"]
        + prefixed_fields(context_fields, "context_", {"viral_contig_id"})
        + prefixed_fields(genomad_fields, "genomad_", {"seq_name"})
        + prefixed_fields(checkv_fields, "checkv_", {"contig_id"})
    )

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for seq_name in sorted(genomad):
            row: dict[str, str] = {"seq_name": seq_name}
            context_row = get_context_row(context, seq_name)
            genomad_row = genomad[seq_name]
            checkv_row = checkv.get(seq_name, {})
            for field in context_fields:
                if field != "viral_contig_id":
                    row[f"context_{field}"] = context_row.get(field, "")
            for field in genomad_fields:
                if field != "seq_name":
                    row[f"genomad_{field}"] = genomad_row.get(field, "")
            for field in checkv_fields:
                if field != "contig_id":
                    row[f"checkv_{field}"] = checkv_row.get(field, "")
            writer.writerow(row)


if __name__ == "__main__":
    main()
