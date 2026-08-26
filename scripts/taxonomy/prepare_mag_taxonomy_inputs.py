#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".fsa", ".fa.gz", ".fasta.gz", ".fna.gz")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare GTDB-Tk batch input and metadata for bacterial/archaeal MAG "
            "taxonomy from dereplicated genome FASTAs."
        )
    )
    parser.add_argument(
        "--genome-dir",
        default="final_drep_all_current_20260504/dereplicated_genomes",
        help="Directory containing dereplicated MAG FASTAs.",
    )
    parser.add_argument(
        "--quality-report",
        default="final_drep_input_all_current_bestbins_20260504/combined_quality_report.tsv",
        help="Combined BASALT quality report TSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="work/mag_taxonomy",
        help="Output directory. Default: work/mag_taxonomy",
    )
    parser.add_argument(
        "--min-completeness",
        type=float,
        default=50.0,
        help="Minimum completeness to include. Default: 50",
    )
    parser.add_argument(
        "--max-contamination",
        type=float,
        default=10.0,
        help="Maximum contamination to include. Default: 10",
    )
    parser.add_argument(
        "--include-low-quality",
        action="store_true",
        help="Include all genomes even if they fail the completeness/contamination cutoffs.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Threads to write into the GTDB-Tk command template. Default: 16",
    )
    parser.add_argument(
        "--path-prefix-map",
        action="append",
        default=[],
        metavar="FROM=TO",
        help=(
            "Rewrite output paths in batch/metadata files. Useful when preparing on "
            "another execution environment, e.g. /data/DogMAG=/work/DogMAG."
        ),
    )
    return parser.parse_args()


def is_fasta(path: Path) -> bool:
    return path.is_file() and any(path.name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def strip_fasta_suffix(name: str) -> str:
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def safe_genome_id(path: Path) -> str:
    token = strip_fasta_suffix(path.name)
    token = re.sub(r"[^A-Za-z0-9_]+", "_", token)
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "genome"


def as_float(text: str) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def load_quality(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return {}
        key = "copied_file" if "copied_file" in reader.fieldnames else reader.fieldnames[0]
        return {row[key]: row for row in reader if row.get(key)}


def parse_prefix_maps(raw_maps: list[str]) -> list[tuple[str, str]]:
    maps: list[tuple[str, str]] = []
    for raw in raw_maps:
        if "=" not in raw:
            raise ValueError(f"Invalid --path-prefix-map value: {raw}")
        left, right = raw.split("=", 1)
        maps.append((left.rstrip("\\/"), right.rstrip("\\/")))
    return maps


def display_path(path: Path, prefix_maps: list[tuple[str, str]]) -> str:
    text = str(path).replace("\\", "/")
    for left, right in prefix_maps:
        left_norm = left.replace("\\", "/")
        if text.lower().startswith(left_norm.lower()):
            return right + text[len(left_norm) :]
    return text


def main() -> None:
    args = parse_args()
    genome_dir = Path(args.genome_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix_maps = parse_prefix_maps(args.path_prefix_map)
    quality = load_quality(Path(args.quality_report))
    genomes = sorted(path.resolve() for path in genome_dir.iterdir() if is_fasta(path))

    batch_path = output_dir / "gtdbtk_batchfile.tsv"
    metadata_path = output_dir / "mag_taxonomy_input_metadata.tsv"
    command_path = output_dir / "run_gtdbtk_classify_wf.sh"
    excluded_path = output_dir / "excluded_genomes.tsv"

    metadata_fields = [
        "genome_id",
        "genome_path",
        "file_name",
        "project",
        "status",
        "original_bin_id",
        "genome_size",
        "completeness",
        "contamination",
        "n50",
    ]
    excluded_fields = metadata_fields + ["exclude_reason"]

    with batch_path.open("w", encoding="utf-8", newline="\n") as batch_handle, metadata_path.open(
        "w", encoding="utf-8", newline=""
    ) as meta_handle, excluded_path.open("w", encoding="utf-8", newline="") as excluded_handle:
        meta_writer = csv.DictWriter(meta_handle, fieldnames=metadata_fields, delimiter="\t")
        excluded_writer = csv.DictWriter(excluded_handle, fieldnames=excluded_fields, delimiter="\t")
        meta_writer.writeheader()
        excluded_writer.writeheader()

        used_ids: set[str] = set()
        for genome_path in genomes:
            base_id = safe_genome_id(genome_path)
            genome_id = base_id
            counter = 2
            while genome_id in used_ids:
                genome_id = f"{base_id}_{counter}"
                counter += 1
            used_ids.add(genome_id)

            q = quality.get(genome_path.name, {})
            completeness = as_float(q.get("Completeness", ""))
            contamination = as_float(q.get("Contamination", ""))
            row = {
                "genome_id": genome_id,
                "genome_path": display_path(genome_path, prefix_maps),
                "file_name": genome_path.name,
                "project": q.get("project", ""),
                "status": q.get("status", ""),
                "original_bin_id": q.get("original_bin_id", ""),
                "genome_size": q.get("Genome_size", ""),
                "completeness": q.get("Completeness", ""),
                "contamination": q.get("Contamination", ""),
                "n50": q.get("N50", ""),
            }

            reasons: list[str] = []
            if completeness is None:
                reasons.append("missing_completeness")
            elif completeness < args.min_completeness:
                reasons.append("low_completeness")
            if contamination is None:
                reasons.append("missing_contamination")
            elif contamination > args.max_contamination:
                reasons.append("high_contamination")

            if reasons and not args.include_low_quality:
                excluded_writer.writerow({**row, "exclude_reason": ";".join(reasons)})
                continue

            batch_handle.write(f"{display_path(genome_path, prefix_maps)}\t{genome_id}\n")
            meta_writer.writerow(row)

    with command_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")
        handle.write('GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:?set GTDBTK_DATA_PATH to the GTDB-Tk reference data directory}"\n')
        handle.write("export GTDBTK_DATA_PATH\n\n")
        handle.write("gtdbtk classify_wf \\\n")
        handle.write(f"  --batchfile {display_path(batch_path, prefix_maps)} \\\n")
        handle.write(f"  --out_dir {display_path(output_dir / 'gtdbtk_out', prefix_maps)} \\\n")
        handle.write(f"  --cpus {args.threads} \\\n")
        handle.write(f"  --pplacer_cpus {args.threads} \\\n")
        handle.write("  --prefix canmag \\\n")
        handle.write("  --force\n")

    print(f"Wrote {batch_path}")
    print(f"Wrote {metadata_path}")
    print(f"Wrote {command_path}")


if __name__ == "__main__":
    main()
