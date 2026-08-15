#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, TextIO


FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas", ".fsa", ".fasta.gz", ".fa.gz", ".fna.gz")


@dataclass(frozen=True)
class FastaRecord:
    header: str
    description: str
    sequence: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare viral-screening FASTAs and a contig context table from assemblies "
            "plus a BASALT contig assignment mask."
        )
    )
    parser.add_argument(
        "--assembly",
        action="append",
        default=[],
        help="Assembly FASTA path. Repeat for multiple assemblies.",
    )
    parser.add_argument(
        "--assembly-dir",
        action="append",
        default=[],
        help="Directory containing assembly FASTAs. Repeat for multiple directories.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively scan --assembly-dir directories.",
    )
    parser.add_argument(
        "--basalt-assignments",
        default="work/viral_contigs/basalt_contig_assignments.tsv",
        help="BASALT assignment TSV from make_basalt_contig_mask.py.",
    )
    parser.add_argument(
        "--output-dir",
        default="work/viral_contigs",
        help="Output directory. Default: work/viral_contigs",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1500,
        help="Minimum contig length to emit. Default: 1500",
    )
    parser.add_argument(
        "--prefix-headers",
        action="store_true",
        help="Prefix emitted FASTA IDs with assembly_id__ to avoid cross-assembly collisions.",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def is_fasta(path: Path) -> bool:
    return path.is_file() and any(path.name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def strip_fasta_suffix(name: str) -> str:
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def sanitize_token(text: str) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    token = re.sub(r"_+", "_", token)
    return token.lower() or "assembly"


def iter_fasta(path: Path) -> Iterator[FastaRecord]:
    header: Optional[str] = None
    description = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header=header, description=description, sequence="".join(chunks))
                text = line[1:]
                parts = text.split(maxsplit=1)
                header = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield FastaRecord(header=header, description=description, sequence="".join(chunks))


def discover_assemblies(args: argparse.Namespace) -> list[Path]:
    paths = [Path(path) for path in args.assembly]
    for raw_dir in args.assembly_dir:
        folder = Path(raw_dir)
        iterator = folder.rglob("*") if args.recursive else folder.iterdir()
        paths.extend(path for path in iterator if is_fasta(path))
    unique = sorted({path.resolve() for path in paths if is_fasta(path)}, key=str)
    if not unique:
        raise FileNotFoundError("No assembly FASTA files were found.")
    return unique


def load_basalt_assignments(path: Path) -> dict[str, list[dict[str, str]]]:
    if not path.exists():
        return {}
    assignments: dict[str, list[dict[str, str]]] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            contig_id = row.get("contig_id", "")
            if contig_id:
                assignments.setdefault(contig_id, []).append(row)
    return assignments


def load_basalt_assignments_by_md5(path: Path) -> dict[str, list[dict[str, str]]]:
    if not path.exists():
        return {}
    assignments: dict[str, list[dict[str, str]]] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seq_md5 = row.get("sequence_md5", "")
            if seq_md5:
                assignments.setdefault(seq_md5, []).append(row)
    return assignments


def write_fasta_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for index in range(0, len(sequence), 80):
        handle.write(sequence[index : index + 80] + "\n")


def summarize_assignment(rows: list[dict[str, str]]) -> dict[str, str]:
    if not rows:
        return {
            "basalt_status": "unbinned",
            "basalt_bins": "",
            "basalt_runs": "",
            "basalt_assemblies": "",
            "max_completeness": "",
            "min_contamination": "",
        }
    completeness_values = []
    contamination_values = []
    for row in rows:
        try:
            completeness_values.append(float(row.get("completeness", "")))
        except ValueError:
            pass
        try:
            contamination_values.append(float(row.get("contamination", "")))
        except ValueError:
            pass
    return {
        "basalt_status": "binned",
        "basalt_bins": ";".join(sorted({row.get("bin_id", "") for row in rows if row.get("bin_id", "")})),
        "basalt_runs": ";".join(sorted({row.get("basalt_run", "") for row in rows if row.get("basalt_run", "")})),
        "basalt_assemblies": ";".join(sorted({row.get("assembly_id", "") for row in rows if row.get("assembly_id", "")})),
        "max_completeness": f"{max(completeness_values):.3f}" if completeness_values else "",
        "min_contamination": f"{min(contamination_values):.3f}" if contamination_values else "",
    }


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    assignments = load_basalt_assignments(Path(args.basalt_assignments))
    md5_assignments = load_basalt_assignments_by_md5(Path(args.basalt_assignments))
    assemblies = discover_assemblies(args)

    all_fasta = output_dir / "viral_screen_all_contigs.fna"
    unbinned_fasta = output_dir / "viral_screen_unbinned_contigs.fna"
    binned_fasta = output_dir / "viral_screen_binned_contigs_for_prophages.fna"
    context_tsv = output_dir / "viral_contig_context.tsv"

    fieldnames = [
        "viral_contig_id",
        "original_contig_id",
        "assembly_id",
        "assembly_path",
        "length",
        "sequence_md5",
        "basalt_match_method",
        "basalt_status",
        "basalt_bins",
        "basalt_runs",
        "basalt_assemblies",
        "max_completeness",
        "min_contamination",
    ]

    with all_fasta.open("w", encoding="utf-8", newline="\n") as all_handle, unbinned_fasta.open(
        "w", encoding="utf-8", newline="\n"
    ) as unbinned_handle, binned_fasta.open("w", encoding="utf-8", newline="\n") as binned_handle, context_tsv.open(
        "w", encoding="utf-8", newline=""
    ) as context_handle:
        writer = csv.DictWriter(context_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for assembly_path in assemblies:
            assembly_id = sanitize_token(strip_fasta_suffix(assembly_path.name))
            for record in iter_fasta(assembly_path):
                length = len(record.sequence)
                if length < args.min_length:
                    continue
                viral_id = f"{assembly_id}__{record.header}" if args.prefix_headers else record.header
                seq_md5 = hashlib.md5(record.sequence.upper().encode("ascii")).hexdigest()
                assignment_rows = assignments.get(record.header, [])
                match_method = "contig_id" if assignment_rows else ""
                if not assignment_rows:
                    assignment_rows = md5_assignments.get(seq_md5, [])
                    match_method = "sequence_md5" if assignment_rows else ""
                summary = summarize_assignment(assignment_rows)
                write_fasta_record(all_handle, viral_id, record.sequence)
                if summary["basalt_status"] == "binned":
                    write_fasta_record(binned_handle, viral_id, record.sequence)
                else:
                    write_fasta_record(unbinned_handle, viral_id, record.sequence)
                writer.writerow(
                    {
                        "viral_contig_id": viral_id,
                        "original_contig_id": record.header,
                        "assembly_id": assembly_id,
                        "assembly_path": str(assembly_path),
                        "length": length,
                        "sequence_md5": seq_md5,
                        "basalt_match_method": match_method,
                        **summary,
                    }
                )


if __name__ == "__main__":
    main()
