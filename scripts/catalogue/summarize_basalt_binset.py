#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Iterable


REPORT_NAMES = (
    "Best_binset_quality_report.tsv",
    "OLC_quality_report.tsv",
)

FASTA_SUFFIXES = (
    ".fa",
    ".fna",
    ".fasta",
    ".fa.gz",
    ".fna.gz",
    ".fasta.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize BASALT bin quality reports. The input can be a BASALT run "
            "directory, a binset directory, or a quality report TSV file."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="BASALT run directory, binset directory, or quality report TSV",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help=(
            "Directory for summary TSVs. Default: <resolved_binset_dir>/summary. "
            "Use - to skip writing files."
        ),
    )
    parser.add_argument(
        "--hq-completeness",
        type=float,
        default=90.0,
        help="HQ completeness cutoff. Default: 90",
    )
    parser.add_argument(
        "--hq-contamination",
        type=float,
        default=5.0,
        help="HQ contamination cutoff (strictly less than). Default: 5",
    )
    parser.add_argument(
        "--mq-completeness",
        type=float,
        default=50.0,
        help="MQ completeness cutoff. Default: 50",
    )
    parser.add_argument(
        "--mq-contamination",
        type=float,
        default=10.0,
        help="MQ contamination cutoff (strictly less than). Default: 10",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve inputs and planned outputs without writing summary files",
    )
    return parser.parse_args()


def path_depth(path: Path) -> int:
    return len(path.parts)


def resolve_quality_report(input_path: Path) -> tuple[Path, Path]:
    if input_path.is_file():
        return input_path.resolve(), input_path.resolve().parent

    if not input_path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    candidates: list[Path] = []
    for report_name in REPORT_NAMES:
        candidates.extend(input_path.rglob(report_name))

    candidates = [candidate.resolve() for candidate in candidates if candidate.is_file()]
    if not candidates:
        raise FileNotFoundError(
            f"Could not find any BASALT quality report under {input_path}. "
            f"Expected one of: {', '.join(REPORT_NAMES)}"
        )

    name_priority = {name: index for index, name in enumerate(REPORT_NAMES)}
    candidates.sort(
        key=lambda candidate: (
            name_priority.get(candidate.name, len(REPORT_NAMES)),
            path_depth(candidate.relative_to(input_path.resolve())),
            -candidate.stat().st_mtime,
            str(candidate),
        )
    )
    report_path = candidates[0]
    return report_path, report_path.parent


def infer_source_parts(source_assembly: str) -> tuple[str, str]:
    if "_" not in source_assembly:
        return "", source_assembly
    prefix, remainder = source_assembly.split("_", 1)
    return (prefix, remainder) if prefix.isdigit() else ("", source_assembly)


def parse_bin_id(bin_id: str) -> dict[str, str]:
    source_assembly = ""
    binner_label = ""
    bin_number = ""

    for suffix in (".fa_", ".fasta_", ".fna_"):
        if suffix in bin_id:
            left, right = bin_id.split(suffix, 1)
            source_assembly = left + suffix[:-1]
            binner_label = right
            break

    if not source_assembly:
        source_assembly = bin_id
        binner_label = ""

    if "." in binner_label:
        maybe_label, maybe_number = binner_label.rsplit(".", 1)
        if maybe_number.isdigit():
            binner_label = maybe_label
            bin_number = maybe_number

    source_index, source_label = infer_source_parts(source_assembly)
    binner_family = "unknown"
    label_lower = binner_label.lower()
    if "maxbin2" in label_lower:
        binner_family = "maxbin2"
    elif "metabat" in label_lower:
        binner_family = "metabat"
    elif "concoct" in label_lower:
        binner_family = "concoct"
    elif "semibin" in label_lower:
        binner_family = "semibin"
    elif "singlecontig" in label_lower:
        binner_family = "single_contig"

    return {
        "source_assembly": source_assembly,
        "source_index": source_index,
        "source_label": source_label,
        "binner_label": binner_label,
        "binner_family": binner_family,
        "bin_number": bin_number,
    }


def median_or_zero(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.median(values) if values else 0.0


def mean_or_zero(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.mean(values) if values else 0.0


def classify_bin(
    completeness: float,
    contamination: float,
    hq_completeness: float,
    hq_contamination: float,
    mq_completeness: float,
    mq_contamination: float,
) -> str:
    if completeness >= hq_completeness and contamination < hq_contamination:
        return "HQ"
    if completeness >= mq_completeness and contamination < mq_contamination:
        return "MQ"
    return "LQ"


def load_rows(
    report_path: Path,
    hq_completeness: float,
    hq_contamination: float,
    mq_completeness: float,
    mq_contamination: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with report_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Bin_ID", "Genome_size", "Completeness", "Contamination", "N50"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Quality report is missing required columns: {', '.join(sorted(missing))}"
            )

        for row in reader:
            parsed = parse_bin_id(row["Bin_ID"])
            completeness = float(row["Completeness"])
            contamination = float(row["Contamination"])
            genome_size = int(float(row["Genome_size"]))
            n50 = float(row["N50"])
            quality_class = classify_bin(
                completeness,
                contamination,
                hq_completeness,
                hq_contamination,
                mq_completeness,
                mq_contamination,
            )
            rows.append(
                {
                    "bin_id": row["Bin_ID"],
                    "genome_size": genome_size,
                    "completeness": completeness,
                    "contamination": contamination,
                    "n50": n50,
                    "quality_class": quality_class,
                    **parsed,
                }
            )
    return rows


def summarize_group(rows: list[dict[str, object]]) -> dict[str, object]:
    return {
        "total_bins": len(rows),
        "hq_bins": sum(1 for row in rows if row["quality_class"] == "HQ"),
        "mq_bins": sum(1 for row in rows if row["quality_class"] in {"HQ", "MQ"}),
        "lq_bins": sum(1 for row in rows if row["quality_class"] == "LQ"),
        "median_genome_size": int(round(median_or_zero(row["genome_size"] for row in rows))),
        "median_completeness": round(median_or_zero(row["completeness"] for row in rows), 2),
        "median_contamination": round(median_or_zero(row["contamination"] for row in rows), 2),
        "median_n50": int(round(median_or_zero(row["n50"] for row in rows))),
        "mean_genome_size": round(mean_or_zero(row["genome_size"] for row in rows), 2),
        "mean_completeness": round(mean_or_zero(row["completeness"] for row in rows), 2),
        "mean_contamination": round(mean_or_zero(row["contamination"] for row in rows), 2),
        "mean_n50": round(mean_or_zero(row["n50"] for row in rows), 2),
    }


def discover_bin_fastas(binset_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in binset_dir.iterdir()
        if path.is_file() and any(path.name.endswith(suffix) for suffix in FASTA_SUFFIXES)
    )


def build_overall_summary(
    rows: list[dict[str, object]],
    report_path: Path,
    binset_dir: Path,
    fasta_count: int,
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    summary = summarize_group(rows)
    return [
        {"metric": "input_path", "value": str(Path(args.input).resolve())},
        {"metric": "binset_dir", "value": str(binset_dir)},
        {"metric": "quality_report", "value": str(report_path)},
        {"metric": "bin_fasta_count", "value": fasta_count},
        {"metric": "report_row_count", "value": len(rows)},
        {"metric": "total_bins", "value": summary["total_bins"]},
        {"metric": "hq_bins", "value": summary["hq_bins"]},
        {"metric": "mq_bins", "value": summary["mq_bins"]},
        {"metric": "lq_bins", "value": summary["lq_bins"]},
        {
            "metric": "hq_fraction",
            "value": f"{(summary['hq_bins'] / len(rows)):.4f}" if rows else "0.0000",
        },
        {
            "metric": "mq_fraction",
            "value": f"{(summary['mq_bins'] / len(rows)):.4f}" if rows else "0.0000",
        },
        {"metric": "median_genome_size", "value": summary["median_genome_size"]},
        {"metric": "median_completeness", "value": summary["median_completeness"]},
        {"metric": "median_contamination", "value": summary["median_contamination"]},
        {"metric": "median_n50", "value": summary["median_n50"]},
        {"metric": "mean_genome_size", "value": summary["mean_genome_size"]},
        {"metric": "mean_completeness", "value": summary["mean_completeness"]},
        {"metric": "mean_contamination", "value": summary["mean_contamination"]},
        {"metric": "mean_n50", "value": summary["mean_n50"]},
        {
            "metric": "hq_definition",
            "value": f"completeness>={args.hq_completeness} and contamination<{args.hq_contamination}",
        },
        {
            "metric": "mq_definition",
            "value": f"completeness>={args.mq_completeness} and contamination<{args.mq_contamination}",
        },
    ]


def build_group_table(
    rows: list[dict[str, object]],
    group_key: str,
    extra_columns: list[str],
) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[group_key])].append(row)

    output_rows: list[dict[str, object]] = []
    for key in sorted(grouped):
        group_rows = grouped[key]
        summary = summarize_group(group_rows)
        out_row = {group_key: key}
        for column in extra_columns:
            out_row[column] = group_rows[0][column]
        out_row.update(summary)
        output_rows.append(out_row)
    return output_rows


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser()
    report_path, binset_dir = resolve_quality_report(input_path)
    fasta_files = discover_bin_fastas(binset_dir)
    rows = load_rows(
        report_path,
        hq_completeness=args.hq_completeness,
        hq_contamination=args.hq_contamination,
        mq_completeness=args.mq_completeness,
        mq_contamination=args.mq_contamination,
    )

    default_output_dir = binset_dir / "summary"
    write_outputs = args.output_dir != "-"
    output_dir = Path(args.output_dir).expanduser() if args.output_dir else default_output_dir

    overall_rows = build_overall_summary(rows, report_path, binset_dir, len(fasta_files), args)
    per_sample_rows = build_group_table(
        rows,
        group_key="source_assembly",
        extra_columns=["source_index", "source_label"],
    )
    per_binner_rows = build_group_table(
        rows,
        group_key="binner_family",
        extra_columns=[],
    )
    detailed_rows = [
        {
            "bin_id": row["bin_id"],
            "source_assembly": row["source_assembly"],
            "source_index": row["source_index"],
            "source_label": row["source_label"],
            "binner_family": row["binner_family"],
            "binner_label": row["binner_label"],
            "bin_number": row["bin_number"],
            "genome_size": row["genome_size"],
            "completeness": f"{row['completeness']:.2f}",
            "contamination": f"{row['contamination']:.2f}",
            "n50": int(round(float(row["n50"]))),
            "quality_class": row["quality_class"],
        }
        for row in rows
    ]

    print(f"Input:\t{input_path.resolve()}")
    print(f"Binset dir:\t{binset_dir}")
    print(f"Quality report:\t{report_path}")
    print(f"Bin FASTA files:\t{len(fasta_files)}")
    print(f"Report rows:\t{len(rows)}")
    print(
        "HQ definition:\t"
        f"completeness>={args.hq_completeness} and contamination<{args.hq_contamination}"
    )
    print(
        "MQ definition:\t"
        f"completeness>={args.mq_completeness} and contamination<{args.mq_contamination}"
    )
    overall = summarize_group(rows)
    print(f"Total bins:\t{overall['total_bins']}")
    print(f"HQ bins:\t{overall['hq_bins']}")
    print(f"MQ bins:\t{overall['mq_bins']}")
    print(f"Median N50:\t{overall['median_n50']}")

    if args.dry_run:
        if write_outputs:
            print(f"[DRY-RUN] Would write summary tables to:\t{output_dir}")
        else:
            print("[DRY-RUN] Output writing disabled")
        return 0

    if write_outputs:
        write_tsv(
            output_dir / "overall_summary.tsv",
            overall_rows,
            fieldnames=["metric", "value"],
        )
        write_tsv(
            output_dir / "per_sample_summary.tsv",
            per_sample_rows,
            fieldnames=[
                "source_assembly",
                "source_index",
                "source_label",
                "total_bins",
                "hq_bins",
                "mq_bins",
                "lq_bins",
                "median_genome_size",
                "median_completeness",
                "median_contamination",
                "median_n50",
                "mean_genome_size",
                "mean_completeness",
                "mean_contamination",
                "mean_n50",
            ],
        )
        write_tsv(
            output_dir / "per_binner_summary.tsv",
            per_binner_rows,
            fieldnames=[
                "binner_family",
                "total_bins",
                "hq_bins",
                "mq_bins",
                "lq_bins",
                "median_genome_size",
                "median_completeness",
                "median_contamination",
                "median_n50",
                "mean_genome_size",
                "mean_completeness",
                "mean_contamination",
                "mean_n50",
            ],
        )
        write_tsv(
            output_dir / "bins_detailed.tsv",
            detailed_rows,
            fieldnames=[
                "bin_id",
                "source_assembly",
                "source_index",
                "source_label",
                "binner_family",
                "binner_label",
                "bin_number",
                "genome_size",
                "completeness",
                "contamination",
                "n50",
                "quality_class",
            ],
        )
        print(f"Summary tables written to:\t{output_dir}")
    else:
        print("Output writing skipped")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
