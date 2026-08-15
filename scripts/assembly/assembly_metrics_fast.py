#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


VALID_SUFFIXES = (
    ".fasta",
    ".fa",
    ".fna",
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute assembly metrics for FASTA files in a directory as quickly as "
            "possible using one worker process per file."
        )
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing FASTA files",
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output TSV path, or - for stdout. Default: -",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Scan subdirectories recursively",
    )
    parser.add_argument(
        "--min-contig-length",
        type=int,
        default=0,
        help="Only include contigs at or above this length. Default: 0",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="Number of parallel worker processes. Default: CPU count",
    )
    return parser.parse_args()


def strip_fasta_suffix(name: str) -> str:
    for suffix in VALID_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def discover_files(input_dir: Path, recursive: bool) -> list[Path]:
    iterator = input_dir.rglob("*") if recursive else input_dir.iterdir()
    files = [
        path
        for path in iterator
        if path.is_file() and any(path.name.endswith(suffix) for suffix in VALID_SUFFIXES)
    ]
    return sorted(files)


def calc_n_stat(lengths: list[int], fraction: float) -> int:
    if not lengths:
        return 0
    threshold = sum(lengths) * fraction
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= threshold:
            return length
    return 0


def calc_median(lengths: list[int]) -> float:
    if not lengths:
        return 0.0
    values = sorted(lengths)
    mid = len(values) // 2
    if len(values) % 2 == 1:
        return float(values[mid])
    return (values[mid - 1] + values[mid]) / 2.0


def open_binary(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def scan_fasta(path: Path, min_contig_length: int) -> dict[str, str]:
    lengths: list[int] = []
    current_len = 0
    gc_count = 0
    n_count = 0
    total_bp_all = 0
    contig_count_all = 0
    has_record = False

    with open_binary(path) as handle:
        for raw_line in handle:
            if not raw_line:
                continue
            if raw_line.startswith(b">"):
                if has_record:
                    contig_count_all += 1
                    total_bp_all += current_len
                    if current_len >= min_contig_length:
                        lengths.append(current_len)
                current_len = 0
                has_record = True
                continue

            line = raw_line.strip()
            if not line:
                continue
            if not has_record:
                raise ValueError("Found sequence content before first FASTA header")

            current_len += len(line)
            gc_count += sum(base in b"GgCc" for base in line)
            n_count += sum(base in b"Nn" for base in line)

    if has_record:
        contig_count_all += 1
        total_bp_all += current_len
        if current_len >= min_contig_length:
            lengths.append(current_len)

    if not has_record:
        raise ValueError("No FASTA records found")

    total_bp = sum(lengths)
    gc_pct = (gc_count / total_bp_all * 100.0) if total_bp_all else 0.0
    n_pct = (n_count / total_bp_all * 100.0) if total_bp_all else 0.0

    return {
        "sample_id": strip_fasta_suffix(path.name),
        "file_name": path.name,
        "file_path": str(path.resolve()),
        "compression": "gzip" if path.suffix == ".gz" else "plain",
        "contig_count": str(len(lengths)),
        "contig_count_all": str(contig_count_all),
        "total_bp": str(total_bp),
        "total_bp_all": str(total_bp_all),
        "longest_contig": str(max(lengths, default=0)),
        "shortest_contig": str(min(lengths, default=0)),
        "mean_contig_length": f"{(total_bp / len(lengths)):.2f}" if lengths else "0.00",
        "median_contig_length": f"{calc_median(lengths):.2f}",
        "n50": str(calc_n_stat(lengths, 0.5)),
        "n90": str(calc_n_stat(lengths, 0.9)),
        "gc_pct": f"{gc_pct:.4f}",
        "n_pct": f"{n_pct:.4f}",
        "contigs_ge_1000": str(sum(1 for value in lengths if value >= 1000)),
        "contigs_ge_2500": str(sum(1 for value in lengths if value >= 2500)),
        "contigs_ge_5000": str(sum(1 for value in lengths if value >= 5000)),
        "min_contig_length_filter": str(min_contig_length),
        "status": "ok",
        "error": "",
    }


def scan_fasta_safe(args: tuple[Path, int]) -> dict[str, str]:
    path, min_contig_length = args
    try:
        return scan_fasta(path, min_contig_length)
    except Exception as exc:
        return {
            "sample_id": strip_fasta_suffix(path.name),
            "file_name": path.name,
            "file_path": str(path.resolve()),
            "compression": "gzip" if path.suffix == ".gz" else "plain",
            "contig_count": "0",
            "contig_count_all": "0",
            "total_bp": "0",
            "total_bp_all": "0",
            "longest_contig": "0",
            "shortest_contig": "0",
            "mean_contig_length": "0.00",
            "median_contig_length": "0",
            "n50": "0",
            "n90": "0",
            "gc_pct": "0.0000",
            "n_pct": "0.0000",
            "contigs_ge_1000": "0",
            "contigs_ge_2500": "0",
            "contigs_ge_5000": "0",
            "min_contig_length_filter": str(min_contig_length),
            "status": "error",
            "error": str(exc),
        }


def write_rows(output_path: str, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "sample_id",
        "file_name",
        "file_path",
        "compression",
        "contig_count",
        "contig_count_all",
        "total_bp",
        "total_bp_all",
        "longest_contig",
        "shortest_contig",
        "mean_contig_length",
        "median_contig_length",
        "n50",
        "n90",
        "gc_pct",
        "n_pct",
        "contigs_ge_1000",
        "contigs_ge_2500",
        "contigs_ge_5000",
        "min_contig_length_filter",
        "status",
        "error",
    ]

    if output_path == "-":
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        return

    output_file = Path(output_path).resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        print(f"Missing input directory: {input_dir}", file=sys.stderr)
        return 1

    files = discover_files(input_dir, args.recursive)
    if not files:
        print(f"No FASTA files found in: {input_dir}", file=sys.stderr)
        return 1

    worker_count = max(1, args.workers)
    job_args = [(path, args.min_contig_length) for path in files]

    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        rows = list(executor.map(scan_fasta_safe, job_args, chunksize=1))

    rows.sort(key=lambda row: row["file_name"])
    write_rows(args.output, rows)

    ok_count = sum(1 for row in rows if row["status"] == "ok")
    error_count = len(rows) - ok_count
    print(
        f"Processed {len(rows)} files with {worker_count} workers; ok={ok_count}; error={error_count}",
        file=sys.stderr,
    )
    return 0 if error_count == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
