#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import Iterator, Optional, TextIO


DEFAULT_QUALITIES = ("Complete", "High-quality", "Medium-quality")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract a filtered viral candidate table, ID list, and FASTA."
    )
    parser.add_argument(
        "--candidates",
        required=True,
        help="Merged candidate TSV from summarize_viral_candidates.py.",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Source FASTA, usually geNomad *_virus.fna.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix for .tsv, .ids.txt, and .fna files.",
    )
    parser.add_argument(
        "--qualities",
        default=",".join(DEFAULT_QUALITIES),
        help=(
            "Comma-separated CheckV qualities to keep. "
            f"Default: {','.join(DEFAULT_QUALITIES)}"
        ),
    )
    parser.add_argument(
        "--min-completeness",
        type=float,
        default=0.0,
        help="Optional minimum CheckV completeness. Default: 0",
    )
    parser.add_argument(
        "--max-contamination",
        type=float,
        default=100.0,
        help="Optional maximum CheckV contamination. Default: 100",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: Optional[str] = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def as_float(text: str, fallback: float) -> float:
    try:
        if text in {"", "NA", "nan"}:
            return fallback
        return float(text)
    except ValueError:
        return fallback


def write_fasta_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for index in range(0, len(sequence), 80):
        handle.write(sequence[index : index + 80] + "\n")


def main() -> None:
    args = parse_args()
    keep_qualities = {item.strip() for item in args.qualities.split(",") if item.strip()}
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    kept_rows: list[dict[str, str]] = []
    with Path(args.candidates).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Empty or malformed TSV: {args.candidates}")
        fieldnames = list(reader.fieldnames)
        for row in reader:
            quality = row.get("checkv_checkv_quality", "")
            completeness = as_float(row.get("checkv_completeness", ""), 0.0)
            contamination = as_float(row.get("checkv_contamination", ""), 100.0)
            if quality not in keep_qualities:
                continue
            if completeness < args.min_completeness:
                continue
            if contamination > args.max_contamination:
                continue
            kept_rows.append(row)

    kept_ids = {row["seq_name"] for row in kept_rows}

    table_path = output_prefix.with_suffix(".tsv")
    ids_path = output_prefix.with_suffix(".ids.txt")
    fasta_path = output_prefix.with_suffix(".fna")

    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept_rows)

    with ids_path.open("w", encoding="utf-8", newline="\n") as handle:
        for seq_id in sorted(kept_ids):
            handle.write(seq_id + "\n")

    written = 0
    with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in iter_fasta(Path(args.fasta)):
            if header in kept_ids:
                write_fasta_record(handle, header, sequence)
                written += 1

    missing = len(kept_ids) - written
    print(f"Kept {len(kept_rows)} candidates; wrote {written} FASTA records to {fasta_path}.")
    if missing:
        print(f"Warning: {missing} kept IDs were not found in the FASTA.")


if __name__ == "__main__":
    main()
