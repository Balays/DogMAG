#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


VALID_SUFFIXES = (
    ".fasta",
    ".fa",
    ".fna",
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
)


@dataclass
class AssemblySource:
    source_path: Path
    collection: str
    platform: str
    assembler: str


@dataclass(frozen=True)
class AssemblyTask:
    source: AssemblySource
    assembly_id: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize collected final assemblies into a BASALT-ready working "
            "directory with a manifest, stats table, filtered FASTAs, and logs."
        )
    )
    parser.add_argument(
        "--input-root",
        default="final_assemblies",
        help="Root directory containing assembly collections. Default: final_assemblies",
    )
    parser.add_argument(
        "--output-root",
        default="work",
        help="Output root for normalized assemblies and reports. Default: work",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1500,
        help="Minimum contig length for the filtered FASTA. Default: 1500",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing normalized outputs",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned actions without writing outputs",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="Number of assemblies to process in parallel. Default: available CPUs",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def strip_fasta_suffix(name: str) -> str:
    for suffix in VALID_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def sanitize_token(text: str) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    token = re.sub(r"_+", "_", token)
    return token.lower() or "assembly"


def discover_sources(input_root: Path) -> list[AssemblySource]:
    collections = {
        "LRS_flye": ("LRS", "flye"),
        "SRS_megahit": ("SRS", "megahit"),
    }
    sources: list[AssemblySource] = []

    if input_root.name in collections:
        platform, assembler = collections[input_root.name]
        for child in sorted(input_root.iterdir()):
            if child.is_file() and child.name.endswith(VALID_SUFFIXES):
                sources.append(
                    AssemblySource(
                        source_path=child,
                        collection=input_root.name,
                        platform=platform,
                        assembler=assembler,
                    )
                )
        return sources

    for collection, (platform, assembler) in collections.items():
        collection_dir = input_root / collection
        if not collection_dir.is_dir():
            continue
        for child in sorted(collection_dir.iterdir()):
            if child.is_file() and child.name.endswith(VALID_SUFFIXES):
                sources.append(
                    AssemblySource(
                        source_path=child,
                        collection=collection,
                        platform=platform,
                        assembler=assembler,
                    )
                )
    return sources


def make_unique_assembly_id(base_id: str, used: set[str]) -> str:
    if base_id not in used:
        used.add(base_id)
        return base_id
    index = 2
    while True:
        candidate = f"{base_id}__{index}"
        if candidate not in used:
            used.add(candidate)
            return candidate
        index += 1


def open_text_auto(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def iter_fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    seq_parts: list[str] = []
    seen_record = False

    with open_text_auto(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(seq_parts)
                    if not sequence:
                        raise ValueError(f"Empty sequence for FASTA record {header!r}")
                    yield header, sequence
                header = line[1:].strip()
                if not header:
                    raise ValueError("Encountered FASTA header without an identifier")
                seq_parts = []
                seen_record = True
                continue
            if header is None:
                raise ValueError("Found sequence content before the first FASTA header")
            seq_parts.append(re.sub(r"\s+", "", line))

    if header is not None:
        sequence = "".join(seq_parts)
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record {header!r}")
        yield header, sequence
    elif not seen_record:
        raise ValueError("No FASTA records found")


def write_record(handle: TextIO, header: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def calc_n50(lengths: Iterable[int]) -> int:
    values = sorted((value for value in lengths if value > 0), reverse=True)
    if not values:
        return 0
    total = sum(values)
    threshold = total / 2
    running = 0
    for value in values:
        running += value
        if running >= threshold:
            return value
    return 0


def scan_fasta_lengths(path: Path) -> list[int]:
    return [len(sequence) for _, sequence in iter_fasta_records(path)]


def log_message(handle: TextIO, message: str) -> None:
    print(message)
    handle.write(message + "\n")
    handle.flush()


def process_source(
    source: AssemblySource,
    output_root: Path,
    min_length: int,
    assembly_id: str,
    force: bool,
    dry_run: bool,
) -> tuple[dict[str, str], dict[str, str], list[str]]:
    sample_id = strip_fasta_suffix(source.source_path.name)
    assembly_dir = output_root / "assemblies" / assembly_id
    ensure_dir(assembly_dir)

    normalized_path = assembly_dir / "assembly.fasta"
    filtered_path = assembly_dir / "assembly.filtered.fa"

    if force:
        if normalized_path.exists():
            normalized_path.unlink()
        if filtered_path.exists():
            filtered_path.unlink()

    notes: list[str] = [f"collection={source.collection}"]
    compression = "gzip" if source.source_path.suffix == ".gz" else "plain"
    status = "ok"
    log_lines: list[str] = []

    if normalized_path.exists() and filtered_path.exists() and not force:
        notes.append("Existing normalized outputs kept")
        status = "skipped_existing"
        try:
            contig_lengths = scan_fasta_lengths(normalized_path)
            filtered_lengths = scan_fasta_lengths(filtered_path)
        except Exception as exc:
            status = "failed_processing"
            notes.append(f"Could not read existing outputs: {exc}")
            contig_lengths = []
            filtered_lengths = []
        else:
            log_lines.append(
                f"[SKIP] {source.source_path} -> {normalized_path} (assembly_id={assembly_id}, status={status})",
            )
            manifest_row = {
                "assembly_id": assembly_id,
                "sample_id": sample_id,
                "platform": source.platform,
                "assembler": source.assembler,
                "source_path": str(source.source_path),
                "normalized_fasta_path": str(normalized_path),
                "filtered_fasta_path": str(filtered_path),
                "compression": compression,
                "status": status,
                "notes": "; ".join(notes),
            }
            stats_row = {
                "assembly_id": assembly_id,
                "sample_id": sample_id,
                "platform": source.platform,
                "assembler": source.assembler,
                "source_path": str(source.source_path),
                "normalized_fasta_path": str(normalized_path),
                "filtered_fasta_path": str(filtered_path),
                "contig_count": str(len(contig_lengths)),
                "total_bp": str(sum(contig_lengths)),
                "n50": str(calc_n50(contig_lengths)),
                "longest_contig": str(max(contig_lengths, default=0)),
                "contigs_ge_1000": str(sum(1 for value in contig_lengths if value >= 1000)),
                "contigs_ge_1500": str(sum(1 for value in contig_lengths if value >= 1500)),
                "contigs_ge_2000": str(sum(1 for value in contig_lengths if value >= 2000)),
                "filtered_contig_count": str(len(filtered_lengths)),
                "filtered_total_bp": str(sum(filtered_lengths)),
                "filtered_n50": str(calc_n50(filtered_lengths)),
                "filtered_longest_contig": str(max(filtered_lengths, default=0)),
                "status": status,
                "notes": "; ".join(notes),
            }
            return manifest_row, stats_row, log_lines

    if dry_run:
        log_lines.append(
            f"[DRY-RUN] {source.source_path} -> {normalized_path} and {filtered_path} (assembly_id={assembly_id})",
        )
        manifest_row = {
            "assembly_id": assembly_id,
            "sample_id": sample_id,
            "platform": source.platform,
            "assembler": source.assembler,
            "source_path": str(source.source_path),
            "normalized_fasta_path": str(normalized_path),
            "filtered_fasta_path": str(filtered_path),
            "compression": compression,
            "status": "dry_run",
            "notes": "; ".join(notes),
        }
        stats_row = {
            "assembly_id": assembly_id,
            "sample_id": sample_id,
            "platform": source.platform,
            "assembler": source.assembler,
            "source_path": str(source.source_path),
            "normalized_fasta_path": str(normalized_path),
            "filtered_fasta_path": str(filtered_path),
            "contig_count": "",
            "total_bp": "",
            "n50": "",
            "longest_contig": "",
            "contigs_ge_1000": "",
            "contigs_ge_1500": "",
            "contigs_ge_2000": "",
            "filtered_contig_count": "",
            "filtered_total_bp": "",
            "filtered_n50": "",
            "filtered_longest_contig": "",
            "status": "dry_run",
            "notes": "; ".join(notes),
        }
        return manifest_row, stats_row, log_lines

    contig_lengths: list[int] = []
    filtered_lengths: list[int] = []

    try:
        log_lines.append(
            f"[RUN] {source.source_path} -> {normalized_path} and {filtered_path} (assembly_id={assembly_id})",
        )
        with normalized_path.open("w", encoding="utf-8") as normalized_handle, filtered_path.open(
            "w", encoding="utf-8"
        ) as filtered_handle:
            for header, sequence in iter_fasta_records(source.source_path):
                length = len(sequence)
                contig_lengths.append(length)
                write_record(normalized_handle, header, sequence)
                if length >= min_length:
                    filtered_lengths.append(length)
                    write_record(filtered_handle, header, sequence)
    except Exception as exc:
        status = "failed_processing"
        notes.append(str(exc))
        log_lines.append(f"[ERROR] {source.source_path}: {exc}")
        if normalized_path.exists():
            normalized_path.unlink(missing_ok=True)
        if filtered_path.exists():
            filtered_path.unlink(missing_ok=True)
        manifest_row = {
            "assembly_id": assembly_id,
            "sample_id": sample_id,
            "platform": source.platform,
            "assembler": source.assembler,
            "source_path": str(source.source_path),
            "normalized_fasta_path": str(normalized_path),
            "filtered_fasta_path": str(filtered_path),
            "compression": compression,
            "status": status,
            "notes": "; ".join(notes),
        }
        stats_row = {
            "assembly_id": assembly_id,
            "sample_id": sample_id,
            "platform": source.platform,
            "assembler": source.assembler,
            "source_path": str(source.source_path),
            "normalized_fasta_path": str(normalized_path),
            "filtered_fasta_path": str(filtered_path),
            "contig_count": "0",
            "total_bp": "0",
            "n50": "0",
            "longest_contig": "0",
            "contigs_ge_1000": "0",
            "contigs_ge_1500": "0",
            "contigs_ge_2000": "0",
            "filtered_contig_count": "0",
            "filtered_total_bp": "0",
            "filtered_n50": "0",
            "filtered_longest_contig": "0",
            "status": status,
            "notes": "; ".join(notes),
        }
        return manifest_row, stats_row, log_lines

    if not contig_lengths:
        status = "failed_processing"
        notes.append("No contigs were written")
    elif not filtered_lengths:
        status = "filtered_empty"
        notes.append(f"No contigs passed min_length={min_length}")

    manifest_row = {
        "assembly_id": assembly_id,
        "sample_id": sample_id,
        "platform": source.platform,
        "assembler": source.assembler,
        "source_path": str(source.source_path),
        "normalized_fasta_path": str(normalized_path),
        "filtered_fasta_path": str(filtered_path),
        "compression": compression,
        "status": status,
        "notes": "; ".join(notes),
    }

    stats_row = {
        "assembly_id": assembly_id,
        "sample_id": sample_id,
        "platform": source.platform,
        "assembler": source.assembler,
        "source_path": str(source.source_path),
        "normalized_fasta_path": str(normalized_path),
        "filtered_fasta_path": str(filtered_path),
        "contig_count": str(len(contig_lengths)),
        "total_bp": str(sum(contig_lengths)),
        "n50": str(calc_n50(contig_lengths)),
        "longest_contig": str(max(contig_lengths, default=0)),
        "contigs_ge_1000": str(sum(1 for value in contig_lengths if value >= 1000)),
        "contigs_ge_1500": str(sum(1 for value in contig_lengths if value >= 1500)),
        "contigs_ge_2000": str(sum(1 for value in contig_lengths if value >= 2000)),
        "filtered_contig_count": str(len(filtered_lengths)),
        "filtered_total_bp": str(sum(filtered_lengths)),
        "filtered_n50": str(calc_n50(filtered_lengths)),
        "filtered_longest_contig": str(max(filtered_lengths, default=0)),
        "status": status,
        "notes": "; ".join(notes),
    }
    log_lines.append(
        f"[OK] {source.source_path} (assembly_id={assembly_id}, status={status})",
    )
    return manifest_row, stats_row, log_lines


def build_tasks(sources: list[AssemblySource]) -> list[AssemblyTask]:
    tasks: list[AssemblyTask] = []
    used_ids: set[str] = set()
    for source in sources:
        sample_id = strip_fasta_suffix(source.source_path.name)
        base_id = f"{sanitize_token(source.platform)}_{sanitize_token(source.assembler)}__{sanitize_token(sample_id)}"
        assembly_id = make_unique_assembly_id(base_id, used_ids)
        tasks.append(AssemblyTask(source=source, assembly_id=assembly_id))
    return tasks


def run_task(
    task: AssemblyTask,
    output_root: Path,
    min_length: int,
    force: bool,
    dry_run: bool,
) -> tuple[dict[str, str], dict[str, str], list[str]]:
    return process_source(
        source=task.source,
        output_root=output_root,
        min_length=min_length,
        assembly_id=task.assembly_id,
        force=force,
        dry_run=dry_run,
    )


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        print("--jobs must be at least 1", file=sys.stderr)
        return 1

    input_root = Path(args.input_root).resolve()
    output_root = Path(args.output_root).resolve()

    if not input_root.is_dir():
        print(f"Missing input root: {input_root}", file=sys.stderr)
        return 1

    ensure_dir(output_root / "assemblies")
    ensure_dir(output_root / "qc")
    ensure_dir(output_root / "logs")

    log_path = output_root / "logs" / "normalization.log"
    manifest_path = output_root / "assembly_manifest.tsv"
    stats_path = output_root / "qc" / "assembly_stats.tsv"

    sources = discover_sources(input_root)
    if not sources:
        print(f"No assembly files found under: {input_root}", file=sys.stderr)
        return 1
    tasks = build_tasks(sources)

    manifest_rows: list[dict[str, str]] = []
    stats_rows: list[dict[str, str]] = []

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(
            f"normalize_final_assemblies.py\tinput_root={input_root}\toutput_root={output_root}\tmin_length={args.min_length}\tdry_run={args.dry_run}\tjobs={args.jobs}\n"
        )
        if args.jobs == 1:
            for manifest_row, stats_row, log_lines in (
                run_task(
                    task=task,
                    output_root=output_root,
                    min_length=args.min_length,
                    force=args.force,
                    dry_run=args.dry_run,
                )
                for task in tasks
            ):
                for line in log_lines:
                    log_message(log_handle, line)
                manifest_rows.append(manifest_row)
                stats_rows.append(stats_row)
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
                results = executor.map(
                    run_task,
                    tasks,
                    [output_root] * len(tasks),
                    [args.min_length] * len(tasks),
                    [args.force] * len(tasks),
                    [args.dry_run] * len(tasks),
                )
                for manifest_row, stats_row, log_lines in results:
                    for line in log_lines:
                        log_message(log_handle, line)
                    manifest_rows.append(manifest_row)
                    stats_rows.append(stats_row)

    manifest_fields = [
        "assembly_id",
        "sample_id",
        "platform",
        "assembler",
        "source_path",
        "normalized_fasta_path",
        "filtered_fasta_path",
        "compression",
        "status",
        "notes",
    ]
    stats_fields = [
        "assembly_id",
        "sample_id",
        "platform",
        "assembler",
        "source_path",
        "normalized_fasta_path",
        "filtered_fasta_path",
        "contig_count",
        "total_bp",
        "n50",
        "longest_contig",
        "contigs_ge_1000",
        "contigs_ge_1500",
        "contigs_ge_2000",
        "filtered_contig_count",
        "filtered_total_bp",
        "filtered_n50",
        "filtered_longest_contig",
        "status",
        "notes",
    ]

    write_tsv(manifest_path, manifest_fields, manifest_rows)
    write_tsv(stats_path, stats_fields, stats_rows)

    ok_count = sum(1 for row in manifest_rows if row["status"] == "ok")
    warning_count = sum(1 for row in manifest_rows if row["status"] == "filtered_empty")
    failed_count = sum(1 for row in manifest_rows if row["status"] == "failed_processing")
    dry_run_count = sum(1 for row in manifest_rows if row["status"] == "dry_run")
    skipped_count = sum(1 for row in manifest_rows if row["status"] == "skipped_existing")

    print(f"Processed assemblies: {len(manifest_rows)}")
    print(f"OK: {ok_count}")
    print(f"Filtered empty: {warning_count}")
    print(f"Failed: {failed_count}")
    print(f"Skipped existing: {skipped_count}")
    print(f"Dry-run only: {dry_run_count}")
    print(f"Manifest: {manifest_path}")
    print(f"Stats: {stats_path}")
    print(f"Log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
