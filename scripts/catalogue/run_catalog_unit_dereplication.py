#!/usr/bin/env python3
"""Create Waltham-like and SGB-like dRep catalogue units from BASALT bins.

This is an internal comparison helper. It does not create official
Huttenhower/PhyloPhlAn SGBs; the 95% ANI output is deliberately labelled
"SGB-like" to keep the comparison honest.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
QUALITY_REPORT_NAMES = (
    "Best_binset_quality_report.tsv",
    "reselected_fasta_manifest.tsv",
    "OLC_quality_report.tsv",
    "combined_quality_report.tsv",
    "drep_genomeInfo.csv",
    "genomeInfo.csv",
)
PROJECT_LABEL_ALIASES = {
    "toti": "Dog_M0",
    "dog_m0": "Dog_M0",
    "dogm0": "Dog_M0",
}


@dataclass
class QualityRecord:
    genome_key: str
    bin_id: str
    completeness: float | None
    contamination: float | None
    genome_size: int | None
    n50: float | None
    project: str
    source_bin_id: str
    inferred_binner: str


@dataclass
class DrepSummary:
    representative_count: int | None
    cluster_count: int | None
    membership_rows: list[dict[str, str]]
    representative_genomes: set[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build quality-filtered CanMAG bin manifests and optionally run dRep "
            "at 99% ANI (strain-level representatives) and 95% ANI "
            "(SGB-like species-level clusters)."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="BASALT run directory, BestBinset directory, or FASTA directory.",
    )
    parser.add_argument(
        "--quality-report",
        default="",
        help=(
            "Optional quality table. Supports BASALT Best_binset_quality_report.tsv, "
            "combined_quality_report.tsv, or dRep genomeInfo CSV."
        ),
    )
    parser.add_argument("--out", required=True, help="Output directory.")
    parser.add_argument("--threads", type=int, default=58, help="dRep threads.")
    parser.add_argument(
        "--project-label",
        default="",
        help="Fallback project/source label when it cannot be inferred from names.",
    )
    parser.add_argument(
        "--drep-executable",
        default="dRep",
        help="dRep executable name or path. Default: dRep.",
    )
    parser.add_argument(
        "--ani99-drep-dir",
        default="",
        help="Optional existing 99%% ANI dRep output directory to summarize.",
    )
    parser.add_argument(
        "--ani95-drep-dir",
        default="",
        help="Optional existing 95%% ANI dRep output directory to summarize.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Stage quality-passing FASTAs and execute both dRep jobs.",
    )
    parser.add_argument(
        "--summarize-only",
        action="store_true",
        help="Do not stage or run dRep; parse existing dRep outputs under --out.",
    )
    return parser.parse_args()


def has_fasta_files(path: Path) -> bool:
    if not path.is_dir():
        return False
    return any(child.is_file() and is_fasta_name(child.name) for child in path.iterdir())


def is_fasta_name(name: str) -> bool:
    lower = name.lower()
    return any(lower.endswith(suffix) for suffix in FASTA_SUFFIXES)


def strip_fasta_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def fasta_stage_name(path: Path, used: set[str]) -> str:
    base = sanitize_filename(strip_fasta_suffix(path.name))
    candidate = f"{base}.fa"
    if candidate not in used:
        used.add(candidate)
        return candidate
    counter = 2
    while True:
        candidate = f"{base}__dup{counter}.fa"
        if candidate not in used:
            used.add(candidate)
            return candidate
        counter += 1


def sanitize_filename(value: str) -> str:
    safe_chars = []
    for char in value:
        if char.isalnum() or char in "._+-=":
            safe_chars.append(char)
        else:
            safe_chars.append("_")
    return "".join(safe_chars).strip("._") or "bin"


def discover_fasta_dir(input_path: Path) -> Path:
    if input_path.is_file():
        raise ValueError(f"--input must be a directory, not a file: {input_path}")
    if not input_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    candidates: list[tuple[int, Path]] = []
    if has_fasta_files(input_path):
        candidates.append((0, input_path))

    for priority, name in enumerate(("BestBinset", "Final_bestbinset", "dereplicated_genomes"), 1):
        child = input_path / name
        if has_fasta_files(child):
            candidates.append((priority, child))

    for child in input_path.iterdir():
        if child.is_dir() and child.name.endswith("BestBinsSet") and has_fasta_files(child):
            candidates.append((10, child))

    if not candidates:
        for root, dirs, _files in bounded_walk(input_path, max_depth=4):
            for dirname in dirs:
                candidate = Path(root) / dirname
                if (
                    dirname in {"BestBinset", "Final_bestbinset", "dereplicated_genomes"}
                    or dirname.endswith("BestBinsSet")
                ) and has_fasta_files(candidate):
                    candidates.append((20 + len(candidate.relative_to(input_path).parts), candidate))

    if not candidates:
        raise FileNotFoundError(
            f"Could not find a FASTA bin directory under {input_path}. "
            "Expected FASTAs directly under --input or in a BASALT BestBinset directory."
        )

    candidates.sort(key=lambda item: (item[0], len(item[1].parts), str(item[1])))
    return candidates[0][1].resolve()


def bounded_walk(root: Path, max_depth: int) -> Iterable[tuple[str, list[str], list[str]]]:
    root = root.resolve()
    for current, dirs, files in os.walk(root):
        current_path = Path(current)
        depth = len(current_path.relative_to(root).parts)
        if depth >= max_depth:
            dirs[:] = []
        yield current, dirs, files


def discover_quality_report(input_path: Path, fasta_dir: Path, provided: str) -> Path:
    if provided:
        path = Path(provided)
        if not path.is_absolute():
            path = Path.cwd() / path
        if not path.is_file():
            raise FileNotFoundError(f"Provided quality report does not exist: {path}")
        return path.resolve()

    direct_candidates: list[Path] = []
    for parent in (fasta_dir, input_path):
        for name in QUALITY_REPORT_NAMES:
            candidate = parent / name
            if candidate.is_file():
                direct_candidates.append(candidate)

    if direct_candidates:
        direct_candidates.sort(key=quality_report_priority)
        return direct_candidates[0].resolve()

    recursive_candidates: list[Path] = []
    for root, _dirs, files in bounded_walk(input_path, max_depth=4):
        for filename in files:
            if filename in QUALITY_REPORT_NAMES:
                recursive_candidates.append(Path(root) / filename)

    if not recursive_candidates:
        raise FileNotFoundError(
            f"Could not find a quality report under {input_path}. "
            f"Expected one of: {', '.join(QUALITY_REPORT_NAMES)}"
        )

    recursive_candidates.sort(key=quality_report_priority)
    return recursive_candidates[0].resolve()


def resolve_optional_path(value: str) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = Path.cwd() / path
    return path.resolve()


def quality_report_priority(path: Path) -> tuple[int, int, str]:
    name_priority = {name: index for index, name in enumerate(QUALITY_REPORT_NAMES)}
    return (name_priority.get(path.name, 99), len(path.parts), str(path))


def detect_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    first_line = sample.splitlines()[0] if sample.splitlines() else ""
    return "\t" if first_line.count("\t") >= first_line.count(",") else ","


def normalize_header(header: str) -> str:
    return header.strip().strip('"').strip().lower().replace(" ", "_").replace("-", "_")


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = str(value).strip().strip('"')
    if not value or value.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value: str | None) -> int | None:
    parsed = parse_float(value)
    if parsed is None:
        return None
    return int(round(parsed))


def first_present(row: dict[str, str], keys: Iterable[str]) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip().strip('"')
    return ""


def load_quality_records(path: Path) -> list[QualityRecord]:
    delimiter = detect_delimiter(path)
    records: list[QualityRecord] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Quality report has no header: {path}")

        for raw_row in reader:
            row = {normalize_header(k): (v or "") for k, v in raw_row.items() if k is not None}
            genome_key = first_present(
                row,
                (
                    "copied_file",
                    "selected_fasta",
                    "selected_file",
                    "genome",
                    "bin_id",
                    "candidate_id",
                    "name",
                    "bin",
                    "binid",
                ),
            )
            bin_id = first_present(
                row,
                (
                    "source_bin_id",
                    "bin_id",
                    "candidate_id",
                    "selected_fasta",
                    "name",
                    "genome",
                    "copied_file",
                    "bin",
                ),
            )
            if not genome_key:
                continue

            records.append(
                QualityRecord(
                    genome_key=Path(genome_key).name,
                    bin_id=bin_id or Path(genome_key).name,
                    completeness=parse_float(first_present(row, ("completeness", "comp"))),
                    contamination=parse_float(first_present(row, ("contamination", "contam", "con"))),
                    genome_size=parse_int(
                        first_present(row, ("genome_size", "genome_size_bp", "size_(bp)", "size", "length"))
                    ),
                    n50=parse_float(first_present(row, ("n50", "n_50"))),
                    project=first_present(row, ("project", "project_subset", "source_project")),
                    source_bin_id=first_present(row, ("source_bin_id", "bin_id")),
                    inferred_binner=first_present(row, ("inferred_binner", "binner", "binner_family")),
                )
            )
    if not records:
        raise ValueError(f"No quality records were parsed from {path}")
    return records


def build_quality_index(records: list[QualityRecord]) -> dict[str, QualityRecord]:
    index: dict[str, QualityRecord] = {}
    for record in records:
        for key in quality_keys(record.genome_key):
            index.setdefault(key, record)
        for key in quality_keys(record.bin_id):
            index.setdefault(key, record)
        if record.source_bin_id:
            for key in quality_keys(record.source_bin_id):
                index.setdefault(key, record)
    return index


def quality_keys(value: str) -> set[str]:
    value = Path(value.strip().strip('"')).name
    if not value:
        return set()
    keys = {value, strip_fasta_suffix(value)}
    if not is_fasta_name(value):
        keys.update({f"{value}.fa", f"{value}.fna", f"{value}.fasta"})
    if "__" in value:
        keys.add(value.split("__", 1)[1])
        keys.add(strip_fasta_suffix(value.split("__", 1)[1]))
    return {key for key in keys if key}


def discover_fastas(fasta_dir: Path) -> list[Path]:
    fastas = sorted(path.resolve() for path in fasta_dir.iterdir() if path.is_file() and is_fasta_name(path.name))
    if not fastas:
        raise FileNotFoundError(f"No FASTA bins found in {fasta_dir}")
    return fastas


def classify_quality(completeness: float | None, contamination: float | None) -> str:
    if completeness is None or contamination is None:
        return "missing_quality"
    if completeness >= 90.0 and contamination <= 5.0:
        return "high"
    if completeness >= 50.0 and contamination <= 10.0:
        return "medium_only"
    return "low_or_other"


def is_quality_passing(quality_class: str) -> bool:
    return quality_class in {"high", "medium_only"}


def parse_bin_id(bin_id: str) -> dict[str, str]:
    name = Path(bin_id).name
    project = ""
    if "__" in name:
        maybe_project, remainder = name.split("__", 1)
        if maybe_project:
            project = maybe_project
            name = remainder

    source_assembly = ""
    binner_label = ""
    for marker in (".fa_", ".fasta_", ".fna_"):
        if marker in name:
            left, right = name.split(marker, 1)
            source_assembly = left + marker[:-1]
            binner_label = strip_fasta_suffix(right)
            break
    if not source_assembly:
        source_assembly = strip_fasta_suffix(name)
        binner_label = ""

    bin_number = ""
    if "." in binner_label:
        maybe_label, maybe_number = binner_label.rsplit(".", 1)
        if maybe_number.isdigit() or maybe_number == "99999":
            binner_label = maybe_label
            bin_number = maybe_number

    source_index = ""
    source_label = source_assembly
    if "_" in source_assembly:
        prefix, remainder = source_assembly.split("_", 1)
        if prefix.isdigit():
            source_index = prefix
            source_label = remainder

    label_lower = binner_label.lower()
    binner_family = "unknown"
    if "maxbin2" in label_lower:
        binner_family = "maxbin2"
    elif "metabat" in label_lower:
        binner_family = "metabat"
    elif "concoct" in label_lower:
        binner_family = "concoct"
    elif "semibin" in label_lower:
        binner_family = "semibin"
    elif "singlecontig" in label_lower or "single_contig" in label_lower:
        binner_family = "single_contig"

    return {
        "project_from_name": project,
        "source_assembly": source_assembly,
        "source_index": source_index,
        "source_label": source_label,
        "binner_label": binner_label,
        "binner_family": binner_family,
        "bin_number": bin_number,
    }


def fasta_stats(path: Path) -> tuple[int, int]:
    lengths: list[int] = []
    current = 0
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                    current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    total = sum(lengths)
    if not lengths:
        return 0, 0
    half = total / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return total, length
    return total, 0


def build_manifest_rows(
    fastas: list[Path],
    quality_index: dict[str, QualityRecord],
    fallback_project: str,
) -> list[dict[str, str]]:
    stage_names: set[str] = set()
    rows: list[dict[str, str]] = []
    for fasta in fastas:
        record = None
        for key in quality_keys(fasta.name):
            record = quality_index.get(key)
            if record:
                break

        bin_id = record.bin_id if record else strip_fasta_suffix(fasta.name)
        parsed = parse_bin_id(bin_id or fasta.name)
        project = (
            (record.project if record else "")
            or parsed["project_from_name"]
            or infer_project_from_filename(fasta.name)
            or fallback_project
        )
        project = normalize_project_label(project)
        binner_family = (record.inferred_binner if record and record.inferred_binner else parsed["binner_family"])
        completeness = record.completeness if record else None
        contamination = record.contamination if record else None
        genome_size = record.genome_size if record else None
        n50 = record.n50 if record else None
        if genome_size is None or n50 is None:
            computed_size, computed_n50 = fasta_stats(fasta)
            genome_size = genome_size if genome_size is not None else computed_size
            n50 = n50 if n50 is not None else computed_n50

        quality_class = classify_quality(completeness, contamination)
        stage_name = fasta_stage_name(fasta, stage_names) if is_quality_passing(quality_class) else ""
        rows.append(
            {
                "input_path": str(fasta),
                "original_file": fasta.name,
                "staged_file": stage_name,
                "bin_id": bin_id,
                "source_project": project,
                "source_assembly": parsed["source_assembly"],
                "source_index": parsed["source_index"],
                "source_label": parsed["source_label"],
                "binner_family": binner_family,
                "binner_label": parsed["binner_label"],
                "bin_number": parsed["bin_number"],
                "completeness": format_optional_float(completeness),
                "contamination": format_optional_float(contamination),
                "genome_size": str(genome_size) if genome_size is not None else "",
                "n50": format_optional_float(n50),
                "quality_class": quality_class,
                "passes_quality_filter": "yes" if is_quality_passing(quality_class) else "no",
            }
        )
    return rows


def infer_project_from_filename(name: str) -> str:
    if "__" in name:
        return name.split("__", 1)[0]
    return ""


def normalize_project_label(value: str) -> str:
    cleaned = value.strip()
    if not cleaned:
        return cleaned
    return PROJECT_LABEL_ALIASES.get(cleaned.lower(), cleaned)


def format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.6g}"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_metric_tsv(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def stage_quality_bins(rows: list[dict[str, str]], stage_dir: Path) -> list[Path]:
    stage_dir.mkdir(parents=True, exist_ok=True)
    staged_paths: list[Path] = []
    for row in rows:
        if row["passes_quality_filter"] != "yes":
            continue
        source = Path(row["input_path"])
        target = stage_dir / row["staged_file"]
        if target.exists():
            if target.stat().st_size == source.stat().st_size and not source.name.lower().endswith(".gz"):
                staged_paths.append(target)
                continue
            raise FileExistsError(f"Staged file already exists with different content: {target}")
        if source.name.lower().endswith(".gz"):
            with gzip.open(source, "rb") as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
        else:
            try:
                os.link(source, target)
            except OSError:
                shutil.copy2(source, target)
        staged_paths.append(target)
    return staged_paths


def write_drep_genome_info(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["genome", "completeness", "contamination", "length", "N50"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            if row["passes_quality_filter"] != "yes":
                continue
            writer.writerow(
                {
                    "genome": row["staged_file"],
                    "completeness": row["completeness"],
                    "contamination": row["contamination"],
                    "length": row["genome_size"],
                    "N50": row["n50"],
                }
            )


def build_drep_command(
    executable: str,
    work_dir: Path,
    stage_dir: Path,
    genome_info: Path,
    threads: int,
    secondary_ani: float,
) -> str:
    genome_glob = f"{shell_quote(str(stage_dir))}/*.fa"
    return (
        f'{shell_quote(executable)} dereplicate {shell_quote(str(work_dir))} '
        f'-g {genome_glob} '
        f'-p {threads} -comp 50 -con 10 -pa 0.90 -sa {secondary_ani:.2f} '
        f'--genomeInfo {shell_quote(str(genome_info))}'
    )


def shell_quote(value: str) -> str:
    if not value:
        return "''"
    if all(char.isalnum() or char in "/._:+-=" for char in value):
        return value
    return "'" + value.replace("'", "'\"'\"'") + "'"


def write_run_scripts(out_dir: Path, commands: dict[str, str]) -> None:
    for label, command in commands.items():
        script = out_dir / f"run_drep_{label}.sh"
        script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n" + command + "\n", encoding="utf-8")


def run_drep_job(
    executable: str,
    work_dir: Path,
    staged_fastas: list[Path],
    genome_info: Path,
    threads: int,
    secondary_ani: float,
) -> None:
    if work_dir.exists():
        raise FileExistsError(
            f"dRep output already exists: {work_dir}. "
            "Use --summarize-only to parse existing results, or choose a new --out directory."
        )
    command = [
        executable,
        "dereplicate",
        str(work_dir),
        "-g",
        *[str(path) for path in staged_fastas],
        "-p",
        str(threads),
        "-comp",
        "50",
        "-con",
        "10",
        "-pa",
        "0.90",
        "-sa",
        f"{secondary_ani:.2f}",
        "--genomeInfo",
        str(genome_info),
    ]
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Could not find dRep executable '{executable}'. "
            "Install/activate dRep or run the generated run_drep_*.sh scripts in the correct environment."
        ) from exc
    except OSError as exc:
        raise RuntimeError(
            f"Failed to launch dRep, possibly because the expanded FASTA argument list is too long: {exc}. "
            "Try running the generated run_drep_*.sh scripts from bash/WSL."
        ) from exc


def parse_drep_output(drep_dir: Path) -> DrepSummary:
    cdb = drep_dir / "data_tables" / "Cdb.csv"
    wdb = drep_dir / "data_tables" / "Wdb.csv"
    derep_dir = drep_dir / "dereplicated_genomes"

    winner_by_cluster: dict[str, str] = {}
    score_by_genome: dict[str, str] = {}
    if wdb.is_file():
        with wdb.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                genome = row.get("genome", "")
                cluster = row.get("cluster", "")
                if cluster and genome:
                    winner_by_cluster[cluster] = genome
                    score_by_genome[genome] = row.get("score", "")

    membership_rows: list[dict[str, str]] = []
    clusters: set[str] = set()
    if cdb.is_file():
        with cdb.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                genome = row.get("genome", "")
                cluster = row.get("secondary_cluster", "")
                if not genome or not cluster:
                    continue
                clusters.add(cluster)
                representative = winner_by_cluster.get(cluster, "")
                membership_rows.append(
                    {
                        "cluster_id": cluster,
                        "genome": genome,
                        "representative": representative,
                        "is_representative": "yes" if representative == genome else "no",
                        "representative_score": score_by_genome.get(representative, ""),
                        "primary_cluster": row.get("primary_cluster", ""),
                        "threshold": row.get("threshold", ""),
                    }
                )

    representative_genomes = {Path(genome).name for genome in winner_by_cluster.values()}
    representative_count = len(representative_genomes) if representative_genomes else None
    if representative_count is None and derep_dir.is_dir():
        representative_genomes = {path.name for path in derep_dir.iterdir() if path.is_file() and is_fasta_name(path.name)}
        representative_count = len(representative_genomes)

    cluster_count = len(clusters) if clusters else representative_count
    return DrepSummary(representative_count, cluster_count, membership_rows, representative_genomes)


def summarize_counts(rows: list[dict[str, str]]) -> dict[str, int]:
    total = len(rows)
    high = sum(1 for row in rows if row["quality_class"] == "high")
    medium = sum(1 for row in rows if row["quality_class"] == "medium_only")
    missing = sum(1 for row in rows if row["quality_class"] == "missing_quality")
    low = total - high - medium - missing
    return {
        "total_input_bins": total,
        "high_quality_bins": high,
        "medium_only_bins": medium,
        "high_plus_medium_bins": high + medium,
        "low_or_other_bins": low,
        "missing_quality_bins": missing,
    }


def write_group_counts(path: Path, rows: list[dict[str, str]], ani99: DrepSummary, ani95: DrepSummary) -> None:
    fieldnames = [
        "group_type",
        "group",
        "input_bins",
        "high",
        "medium_only",
        "high_plus_medium",
        "low_or_other",
        "ani99_representatives",
        "ani95_sgb_like_clusters",
    ]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        for group_type, key in (("source_project", "source_project"), ("binner_family", "binner_family")):
            group = row.get(key, "") or "unknown"
            grouped.setdefault((group_type, group), []).append(row)

    out_rows: list[dict[str, str]] = []
    for (group_type, group), group_rows in sorted(grouped.items()):
        counts = summarize_counts(group_rows)
        out_rows.append(
            {
                "group_type": group_type,
                "group": group,
                "input_bins": str(counts["total_input_bins"]),
                "high": str(counts["high_quality_bins"]),
                "medium_only": str(counts["medium_only_bins"]),
                "high_plus_medium": str(counts["high_plus_medium_bins"]),
                "low_or_other": str(counts["low_or_other_bins"]),
                "ani99_representatives": count_representatives_for_group(ani99, group_rows),
                "ani95_sgb_like_clusters": count_representatives_for_group(ani95, group_rows),
            }
        )
    write_tsv(path, out_rows, fieldnames)


def count_representatives_for_group(summary: DrepSummary, group_rows: list[dict[str, str]]) -> str:
    if summary.representative_count is None:
        return "NA"
    group_stage_names = {Path(row["staged_file"]).name for row in group_rows if row.get("staged_file")}
    return str(len(summary.representative_genomes.intersection(group_stage_names)))


def write_comparison_context(path: Path, canmag_counts: dict[str, str]) -> None:
    rows = [
        {
            "concept": "Total generated/categorized MAGs",
            "branck": "61515",
            "waltham": "5753",
            "canmag": canmag_counts.get("total_input_bins", ""),
            "canmag_source": "Input bin count",
        },
        {
            "concept": "Quality-passing MAGs",
            "branck": "28981 high/medium",
            "waltham": "quality-filtered before 99% dRep",
            "canmag": canmag_counts.get("high_plus_medium_bins", ""),
            "canmag_source": "manifest_quality_pass_bins.tsv",
        },
        {
            "concept": "Species-level/SGB-like units",
            "branck": "2320 official SGBs",
            "waltham": "not main unit",
            "canmag": canmag_counts.get("ani95_sgb_like_clusters", ""),
            "canmag_source": "drep_ani95 winners; SGB-like, not official SGBs",
        },
        {
            "concept": "Strain-level representatives",
            "branck": "not main unit",
            "waltham": "1031",
            "canmag": canmag_counts.get("ani99_strain_level_representatives", ""),
            "canmag_source": "drep_ani99 winners",
        },
    ]
    write_tsv(path, rows, ["concept", "branck", "waltham", "canmag", "canmag_source"])


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_dir = discover_fasta_dir(input_path)
    quality_report = discover_quality_report(input_path, fasta_dir, args.quality_report)
    records = load_quality_records(quality_report)
    quality_index = build_quality_index(records)
    fastas = discover_fastas(fasta_dir)
    fallback_project = args.project_label or input_path.name
    manifest_rows = build_manifest_rows(fastas, quality_index, fallback_project)
    quality_rows = [row for row in manifest_rows if row["passes_quality_filter"] == "yes"]

    fieldnames = [
        "input_path",
        "original_file",
        "staged_file",
        "bin_id",
        "source_project",
        "source_assembly",
        "source_index",
        "source_label",
        "binner_family",
        "binner_label",
        "bin_number",
        "completeness",
        "contamination",
        "genome_size",
        "n50",
        "quality_class",
        "passes_quality_filter",
    ]
    write_tsv(out_dir / "manifest_all_bins.tsv", manifest_rows, fieldnames)
    write_tsv(out_dir / "manifest_quality_pass_bins.tsv", quality_rows, fieldnames)

    genome_info = out_dir / "drep_genomeInfo.csv"
    write_drep_genome_info(genome_info, quality_rows)

    stage_dir = out_dir / "quality_pass_bins"
    drep_ani99_dir = resolve_optional_path(args.ani99_drep_dir) or out_dir / "drep_ani99"
    drep_ani95_dir = resolve_optional_path(args.ani95_drep_dir) or out_dir / "drep_ani95"
    commands = {
        "ani99": build_drep_command(args.drep_executable, drep_ani99_dir, stage_dir, genome_info, args.threads, 0.99),
        "ani95": build_drep_command(args.drep_executable, drep_ani95_dir, stage_dir, genome_info, args.threads, 0.95),
    }
    write_run_scripts(out_dir, commands)

    staged_fastas: list[Path] = []
    if args.run and not args.summarize_only:
        staged_fastas = stage_quality_bins(quality_rows, stage_dir)
        run_drep_job(args.drep_executable, drep_ani99_dir, staged_fastas, genome_info, args.threads, 0.99)
        run_drep_job(args.drep_executable, drep_ani95_dir, staged_fastas, genome_info, args.threads, 0.95)

    ani99 = parse_drep_output(drep_ani99_dir)
    ani95 = parse_drep_output(drep_ani95_dir)
    write_tsv(
        out_dir / "ani99_cluster_membership.tsv",
        ani99.membership_rows,
        [
            "cluster_id",
            "genome",
            "representative",
            "is_representative",
            "representative_score",
            "primary_cluster",
            "threshold",
        ],
    )
    write_tsv(
        out_dir / "ani95_cluster_membership.tsv",
        ani95.membership_rows,
        [
            "cluster_id",
            "genome",
            "representative",
            "is_representative",
            "representative_score",
            "primary_cluster",
            "threshold",
        ],
    )

    counts = summarize_counts(manifest_rows)
    ani99_count = ani99.representative_count
    ani95_count = ani95.representative_count
    canmag_counts = {
        **{key: str(value) for key, value in counts.items()},
        "ani99_strain_level_representatives": str(ani99_count) if ani99_count is not None else "NA",
        "ani95_sgb_like_clusters": str(ani95_count) if ani95_count is not None else "NA",
    }

    write_metric_tsv(
        out_dir / "quality_summary.tsv",
        [
            ("input_path", str(input_path)),
            ("fasta_dir", str(fasta_dir)),
            ("quality_report", str(quality_report)),
            ("drep_ani99_dir", str(drep_ani99_dir)),
            ("drep_ani95_dir", str(drep_ani95_dir)),
            ("total_input_bins", str(counts["total_input_bins"])),
            ("high_quality_bins", str(counts["high_quality_bins"])),
            ("medium_only_bins", str(counts["medium_only_bins"])),
            ("high_plus_medium_bins", str(counts["high_plus_medium_bins"])),
            ("low_or_other_bins", str(counts["low_or_other_bins"])),
            ("missing_quality_bins", str(counts["missing_quality_bins"])),
            ("run_requested", "yes" if args.run else "no"),
            ("summarize_only", "yes" if args.summarize_only else "no"),
        ],
    )
    write_metric_tsv(
        out_dir / "catalog_unit_summary.tsv",
        [
            ("total_input_bins", str(counts["total_input_bins"])),
            ("quality_passing_bins", str(counts["high_plus_medium_bins"])),
            ("ani99_strain_level_representatives", str(ani99_count) if ani99_count is not None else "NA"),
            ("ani95_sgb_like_species_level_clusters", str(ani95_count) if ani95_count is not None else "NA"),
            ("ani99_label", "Waltham-like 99% ANI strain-level representatives"),
            ("ani95_label", "SGB-like 95% ANI species-level clusters; not official Huttenhower SGBs"),
        ],
    )
    write_group_counts(out_dir / "quality_counts_by_group.tsv", manifest_rows, ani99, ani95)
    write_comparison_context(out_dir / "comparison_context.tsv", canmag_counts)

    print(f"Input FASTA directory: {fasta_dir}")
    print(f"Quality report: {quality_report}")
    print(f"Total input bins: {counts['total_input_bins']}")
    print(f"Quality-passing bins: {counts['high_plus_medium_bins']}")
    print(f"99% ANI representatives: {ani99_count if ani99_count is not None else 'NA'}")
    print(f"95% ANI SGB-like clusters: {ani95_count if ani95_count is not None else 'NA'}")
    if not args.run:
        print("dRep was not run. Use --run or execute the generated run_drep_*.sh scripts.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
