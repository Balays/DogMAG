#!/usr/bin/env python3
"""Refresh CheckM2 outputs for existing BASALT ``*_genomes`` folders.

Use this when binner FASTA folders exist but their paired ``*_checkm`` folders
are missing ``quality_report.tsv``. The helper does not rerun binners or change
``Bins_folder.txt``; it only rebuilds stale CheckM2 output folders in place.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterator, Optional


FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".fsa")
CHECKM2_QUALITY_REPORT_HEADER = (
    "Name\tCompleteness\tContamination\tCompleteness_Model_Used\t"
    "Translation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\t"
    "Genome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\t"
    "Max_Contig_Length\tAdditional_Notes\n"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Refresh CheckM2 quality reports for BASALT binner folders."
    )
    parser.add_argument("--workspace", required=True, help="BASALT workspace path.")
    parser.add_argument("--threads", type=int, default=20, help="CheckM2 threads.")
    parser.add_argument(
        "--checkm2db",
        help="Path to uniref100.KO.1.dmnd or its containing CheckM2 directory.",
    )
    parser.add_argument(
        "--assemblies-prefix",
        help=(
            "Optional prefix filter, for example '1_' or "
            "'5_lrs_flye__b_barcode37.fa'."
        ),
    )
    parser.add_argument("--force", action="store_true", help="Rebuild even if report exists.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions only.")
    return parser.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def iter_fasta_files(folder: Path) -> Iterator[Path]:
    for child in folder.iterdir():
        if child.is_file() and child.suffix.lower() in FASTA_SUFFIXES:
            yield child


def detect_extension(folder: Path) -> Optional[str]:
    counts: dict[str, int] = {}
    for fasta in iter_fasta_files(folder):
        suffix = fasta.suffix.lower().lstrip(".")
        counts[suffix] = counts.get(suffix, 0) + 1
    if not counts:
        return None
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def infer_extension_from_name(folder: Path) -> str:
    name = folder.name.lower()
    if "maxbin2" in name or "concoct" in name:
        return "fasta"
    return "fa"


def normalize_checkm2db(path_value: Optional[str]) -> Optional[str]:
    path_value = path_value or os.environ.get("CHECKM2DB")
    if not path_value:
        return None

    path = Path(path_value)
    if path.is_dir() and (path / "uniref100.KO.1.dmnd").is_file():
        return str(path / "uniref100.KO.1.dmnd")
    return str(path)


def run_checkm2(
    workspace: Path,
    genomes_dir: Path,
    checkm_dir: Path,
    extension: str,
    threads: int,
    checkm2db: Optional[str],
    force: bool,
    dry_run: bool,
) -> bool:
    report = checkm_dir / "quality_report.tsv"
    if report.exists() and not force:
        log(f"[SKIP] {genomes_dir.name}: quality_report.tsv already exists")
        return False

    if checkm_dir.exists():
        log(f"[INFO] Remove stale {checkm_dir.name}")
        if not dry_run:
            shutil.rmtree(checkm_dir)

    if not any(iter_fasta_files(genomes_dir)):
        log(f"[EMPTY] {genomes_dir.name}: writing header-only quality_report.tsv")
        if not dry_run:
            checkm_dir.mkdir(parents=True, exist_ok=True)
            report.write_text(CHECKM2_QUALITY_REPORT_HEADER)
        return True

    cmd = [
        "checkm2",
        "predict",
        "-t",
        str(threads),
        "-i",
        genomes_dir.name,
        "-x",
        extension,
        "-o",
        checkm_dir.name,
    ]
    log("[CMD] " + " ".join(cmd))
    if dry_run:
        return True

    env = os.environ.copy()
    if checkm2db:
        env["CHECKM2DB"] = checkm2db
    subprocess.run(cmd, cwd=str(workspace), check=True, env=env)
    if not report.exists():
        raise RuntimeError(f"CheckM2 did not create {report}")
    return True


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    if not workspace.is_dir():
        raise NotADirectoryError(workspace)

    checkm2db = normalize_checkm2db(args.checkm2db)
    if checkm2db:
        log(f"[INFO] CHECKM2DB={checkm2db}")
    else:
        log("[WARN] CHECKM2DB was not provided; relying on CheckM2 global config")

    processed = 0
    skipped = 0
    for genomes_dir in sorted(workspace.glob("*_genomes")):
        if not genomes_dir.is_dir():
            continue
        if args.assemblies_prefix and not genomes_dir.name.startswith(args.assemblies_prefix):
            continue
        extension = detect_extension(genomes_dir) or infer_extension_from_name(genomes_dir)
        checkm_dir = workspace / genomes_dir.name.replace("_genomes", "_checkm")
        changed = run_checkm2(
            workspace,
            genomes_dir,
            checkm_dir,
            extension,
            args.threads,
            checkm2db,
            args.force,
            args.dry_run,
        )
        if changed:
            processed += 1
        else:
            skipped += 1

    log(f"[DONE] processed={processed} skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
