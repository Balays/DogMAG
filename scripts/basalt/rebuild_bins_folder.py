#!/usr/bin/env python3
"""Rebuild BASALT Bins_folder.txt from the actual *_genomes folders.

This is intended for partially recovered BASALT workspaces where the saved
Bins_folder.txt has become misaligned with the assembly order after cleanup or
recovery passes. The script scans the workspace for top-level ``*_genomes``
directories that still contain bin FASTA files, groups them by assembly, and
rewrites ``Bins_folder.txt`` in the workspace's real assembly order.
"""

from __future__ import annotations

import argparse
import re
import shlex
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".fsa"}
GENOMES_SUFFIX = "_genomes"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rebuild BASALT Bins_folder.txt from existing *_genomes folders."
    )
    parser.add_argument("--workspace", required=True, help="Path to the BASALT workspace.")
    parser.add_argument(
        "--backup-suffix",
        default=".bak_rebuild",
        help="Suffix used when backing up an existing Bins_folder.txt.",
    )
    return parser.parse_args()


def log(msg: str) -> None:
    print(msg, flush=True)


def has_bin_fastas(folder: Path) -> bool:
    for child in folder.iterdir():
        if child.is_file() and child.suffix.lower() in FASTA_SUFFIXES:
            return True
    return False


def folder_prefix(folder_name: str) -> str:
    pattern = re.compile(
        r"_(?:0\.\d+_maxbin2|[0-9]+_(?:metabat|concoct)|100_semibin|1_SingleContig)_genomes$"
    )
    stripped = pattern.sub("", folder_name)
    return stripped if stripped != folder_name else folder_name


def assembly_key(name: str) -> str:
    return re.sub(r"^\d+_", "", name)


def folder_sort_key(folder_name: str) -> Tuple[int, float, str]:
    match = re.search(r"_(0\.\d+)_maxbin2_genomes$", folder_name)
    if match:
        return (0, float(match.group(1)), folder_name)

    match = re.search(r"_(\d+)_metabat_genomes$", folder_name)
    if match:
        return (1, float(match.group(1)), folder_name)

    match = re.search(r"_(\d+)_concoct_genomes$", folder_name)
    if match:
        return (2, float(match.group(1)), folder_name)

    if folder_name.endswith("_100_semibin_genomes"):
        return (3, 100.0, folder_name)

    if folder_name.endswith("_1_SingleContig_genomes"):
        return (4, 1.0, folder_name)

    return (9, 0.0, folder_name)


def parse_run_basalt_assemblies(run_script: Path) -> List[str]:
    text = run_script.read_text(encoding="utf-8", errors="replace")
    text = text.replace("\\\r\n", " ").replace("\\\n", " ")
    text = text.replace("\r", "")
    tokens = shlex.split(text)
    for idx, token in enumerate(tokens[:-1]):
        if token == "-a":
            return [item for item in tokens[idx + 1].split(",") if item]
    return []


def fallback_assembly_entries(workspace: Path) -> List[str]:
    assemblies_dir = workspace / "assemblies"
    if not assemblies_dir.exists():
        return []
    entries = sorted(
        p.name
        for p in assemblies_dir.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES
    )
    return [f"assemblies/{name}" for name in entries]


def assembly_entries(workspace: Path) -> List[str]:
    run_script = workspace / "run_basalt.sh"
    if run_script.exists():
        parsed = parse_run_basalt_assemblies(run_script)
        if parsed:
            return parsed
    return fallback_assembly_entries(workspace)


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    if not workspace.exists():
        raise FileNotFoundError(f"Workspace does not exist: {workspace}")

    entries = assembly_entries(workspace)
    if not entries:
        raise RuntimeError("Could not determine assembly order from run_basalt.sh or assemblies/")

    grouped: Dict[str, List[str]] = {}
    for folder in workspace.iterdir():
        if not folder.is_dir() or not folder.name.endswith(GENOMES_SUFFIX):
            continue
        if not has_bin_fastas(folder):
            continue
        grouped.setdefault(assembly_key(folder_prefix(folder.name)), []).append(folder.name)

    rows: List[Tuple[str, List[str]]] = []
    missing: List[str] = []
    for entry in entries:
        base = assembly_key(Path(entry).name)
        folders = sorted(grouped.get(base, []), key=folder_sort_key)
        rows.append((entry, folders))
        if not folders:
            missing.append(entry)

    bins_folder = workspace / "Bins_folder.txt"
    if bins_folder.exists():
        backup = bins_folder.with_name(bins_folder.name + args.backup_suffix)
        if not backup.exists():
            backup.write_text(bins_folder.read_text(encoding="utf-8"), encoding="utf-8")
            log(f"Backed up existing {bins_folder.name} to {backup.name}")

    with bins_folder.open("w", encoding="utf-8", newline="\n") as handle:
        for entry, folders in rows:
            handle.write(f"{entry}\t{folders!r}\n")

    kept_count = sum(1 for _, folders in rows if folders)
    log(f"Rebuilt {bins_folder.name} with {len(rows)} assembly row(s); {kept_count} row(s) have bin folders.")
    if missing:
        missing_path = workspace / "rebuild_bins_folder_missing.tsv"
        with missing_path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("assembly\n")
            for entry in missing:
                handle.write(f"{entry}\n")
        log(f"Wrote missing-assembly report to {missing_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
