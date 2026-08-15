#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


READ_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
)

SKIP_DIR_NAMES = {
    ".fastq_read_stats",
    "fastqc",
    "multiqc_data",
}

ASSEMBLY_COLLECTIONS = {
    "LRS_flye": ("LRS", "flye"),
    "SRS_megahit": ("SRS", "megahit"),
}


@dataclass
class FastqRecord:
    read_id: str
    basename: str
    source_path: Path
    relative_path: str
    top_group: str
    read_type: str
    read_layout: str
    mate_label: str
    token: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a FASTQ inventory and propose links between normalized assemblies "
            "and raw read files."
        )
    )
    parser.add_argument(
        "--assembly-manifest",
        default=None,
        help="Normalized assembly manifest TSV. Defaults to work/assembly_manifest.tsv if present.",
    )
    parser.add_argument(
        "--assembly-root",
        default=None,
        help=(
            "Assembly root directory to scan directly, for example final_assemblies or "
            "final_assemblies/LRS_flye. Used when no manifest is supplied."
        ),
    )
    parser.add_argument(
        "--fastq-root",
        default="fastq",
        help="FASTQ root directory to scan. Default: fastq",
    )
    parser.add_argument(
        "--extra-fastq-list",
        action="append",
        default=[],
        help=(
            "Plain-text manifest of additional FASTQ paths to include even when they "
            "are not physically present under --fastq-root. May be given multiple times."
        ),
    )
    parser.add_argument(
        "--output-root",
        default="work",
        help="Output root for manifests and logs. Default: work",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned outputs without writing files",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    parser.add_argument(
        "--refine",
        action="store_true",
        help=(
            "Refine an existing assembly_read_links.tsv by replacing only rows that "
            "are currently unmatched or review, while keeping matched rows unchanged"
        ),
    )
    parser.add_argument(
        "--link-mode",
        choices=("all", "best"),
        default="all",
        help=(
            "For assemblies with multiple plausible FASTQ links, keep all equal-top "
            "candidates or only the single best row. Default: all"
        ),
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="Number of assembly chunks to score in parallel. Default: available CPUs",
    )
    return parser.parse_args()


def log_message(handle: TextIO, message: str) -> None:
    print(message)
    handle.write(message + "\n")
    handle.flush()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def strip_read_suffix(name: str) -> str:
    for suffix in READ_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def sanitize_token(text: str) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    token = re.sub(r"_+", "_", token)
    return token.lower()


def paired_group_key(text: str) -> str:
    text = strip_read_suffix(text)
    text = re.sub(r"(_R1_filt|_R2_filt|_R1|_R2|\.R1|\.R2)$", "", text, flags=re.IGNORECASE)
    return sanitize_token(text)


def detect_read_layout(base: str) -> tuple[str, str]:
    if re.search(r"(?:^|[_\-\.])(r?1)(?:$|[_\-\.])", base, flags=re.IGNORECASE):
        return "paired", "R1"
    if re.search(r"(?:^|[_\-\.])(r?2)(?:$|[_\-\.])", base, flags=re.IGNORECASE):
        return "paired", "R2"
    return "single", ""


def infer_read_type(path: Path, base: str) -> str:
    joined = f"{path.as_posix()} {base}".lower()
    read_layout, _ = detect_read_layout(base)

    if read_layout == "paired":
        return "short"

    short_context_markers = (
        "short_reads",
        "short_reads_pairs",
        "dmd_wgs_merged",
        "toti_srs_",
        "novaseq",
        "miseq",
        "illumina",
        "fastp",
        "merged",
    )
    if any(marker in joined for marker in short_context_markers):
        return "short"

    if (
        "long_reads" in joined
        or "ont" in joined
        or "barcode" in joined
        or "zymo_hmw" in joined
        or re.search(r"(?:^|[_\-.])bc\d+(?:$|[_\-.])", joined)
        or re.search(r"(?:^|[_\-.])mn(?:$|[_\-.])", joined)
        or "pass" in joined
    ):
        return "long"
    return "short"


def should_skip_path(path: Path) -> bool:
    if any(part in SKIP_DIR_NAMES for part in path.parts):
        return True
    if any(part.startswith(".") and part != "." for part in path.parts):
        return True
    lowered = path.name.lower()
    return lowered.endswith((".html", ".zip", ".txt", ".xlsx", ".json", ".tsv", ".tsv.gz", ".jpg", ".png"))


def make_fastq_record(
    source_path: Path,
    *,
    basename_override: str | None = None,
    relative_path: str,
    top_group: str,
    seen_ids: set[str],
) -> FastqRecord:
    base = basename_override if basename_override is not None else strip_read_suffix(source_path.name)
    layout, mate = detect_read_layout(base)
    token = sanitize_token(base)
    read_id = token
    counter = 2
    while read_id in seen_ids:
        read_id = f"{token}__{counter}"
        counter += 1
    seen_ids.add(read_id)
    return FastqRecord(
        read_id=read_id,
        basename=base,
        source_path=source_path,
        relative_path=relative_path,
        top_group=top_group,
        read_type=infer_read_type(source_path, base),
        read_layout=layout,
        mate_label=mate,
        token=token,
    )


def discover_fastqs(fastq_root: Path) -> list[FastqRecord]:
    records: list[FastqRecord] = []
    seen_ids: set[str] = set()
    seen_sources: set[Path] = set()

    for path in sorted(fastq_root.rglob("*")):
        if not path.is_file():
            continue
        if should_skip_path(path):
            continue
        if not path.name.endswith(READ_SUFFIXES):
            continue
        resolved_source = path.resolve()
        if resolved_source in seen_sources:
            continue
        seen_sources.add(resolved_source)

        relative = path.relative_to(fastq_root).as_posix()
        top_group = relative.split("/", 1)[0]
        records.append(
            make_fastq_record(
                path,
                basename_override=strip_read_suffix(path.name),
                relative_path=relative,
                top_group=top_group,
                seen_ids=seen_ids,
            )
        )
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def discover_fastqs_from_list(
    manifest_path: Path,
    *,
    seen_ids: set[str],
    seen_sources: set[str],
) -> list[FastqRecord]:
    records: list[FastqRecord] = []
    pattern = re.compile(r"([A-Za-z0-9._-]+\.(?:fastq|fq)(?:\.gz)?)\s+->\s+(.+)$")

    with manifest_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if "->" not in line:
                continue
            match = pattern.search(line)
            if not match:
                continue

            link_name = match.group(1).strip()
            source_str = match.group(2).strip()
            source_key = source_str.lower()
            if source_key in seen_sources:
                continue
            seen_sources.add(source_key)

            source_path = Path(source_str)
            relative_path = f"__external__/{manifest_path.name}/{link_name}"
            records.append(
                make_fastq_record(
                    source_path,
                    basename_override=strip_read_suffix(link_name),
                    relative_path=relative_path,
                    top_group="external_manifest",
                    seen_ids=seen_ids,
                )
            )

    return records


def strip_fasta_suffix(name: str) -> str:
    for suffix in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def make_unique_assembly_id(base_id: str, used: set[str]) -> str:
    if base_id not in used:
        used.add(base_id)
        return base_id
    counter = 2
    while True:
        candidate = f"{base_id}__{counter}"
        if candidate not in used:
            used.add(candidate)
            return candidate
        counter += 1


def discover_assembly_rows(assembly_root: Path) -> list[dict[str, str]]:
    roots: list[tuple[Path, str, str, str]] = []

    if assembly_root.name in ASSEMBLY_COLLECTIONS:
        platform, assembler = ASSEMBLY_COLLECTIONS[assembly_root.name]
        roots.append((assembly_root, assembly_root.name, platform, assembler))
    else:
        for collection, (platform, assembler) in ASSEMBLY_COLLECTIONS.items():
            collection_dir = assembly_root / collection
            if collection_dir.is_dir():
                roots.append((collection_dir, collection, platform, assembler))

    if not roots:
        return []

    rows: list[dict[str, str]] = []
    used_ids: set[str] = set()
    valid_suffixes = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")

    for collection_dir, collection, platform, assembler in roots:
        for path in sorted(collection_dir.iterdir()):
            if not path.is_file() or not path.name.endswith(valid_suffixes):
                continue
            sample_id = strip_fasta_suffix(path.name)
            base_id = f"{sanitize_token(platform)}_{sanitize_token(assembler)}__{sanitize_token(sample_id)}"
            assembly_id = make_unique_assembly_id(base_id, used_ids)
            rows.append(
                {
                    "assembly_id": assembly_id,
                    "sample_id": sample_id,
                    "platform": platform,
                    "assembler": assembler,
                    "status": "ok",
                    "source_path": str(path.resolve()),
                    "notes": f"collection={collection}",
                }
            )
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def merge_refined_links(
    existing_rows: list[dict[str, str]],
    new_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    existing_by_assembly: dict[str, list[dict[str, str]]] = {}
    for row in existing_rows:
        existing_by_assembly.setdefault(row["assembly_id"], []).append(row)

    new_by_assembly: dict[str, list[dict[str, str]]] = {}
    for row in new_rows:
        new_by_assembly.setdefault(row["assembly_id"], []).append(row)

    merged: list[dict[str, str]] = []
    assembly_ids = []
    seen: set[str] = set()
    for row in existing_rows + new_rows:
        assembly_id = row["assembly_id"]
        if assembly_id not in seen:
            seen.add(assembly_id)
            assembly_ids.append(assembly_id)

    for assembly_id in assembly_ids:
        existing = existing_by_assembly.get(assembly_id, [])
        new = new_by_assembly.get(assembly_id, [])

        if existing and any(row.get("match_status") == "matched" for row in existing):
            merged.extend(existing)
        elif new:
            merged.extend(new)
        else:
            merged.extend(existing)

    return merged


def date_token(text: str) -> str:
    match = re.search(r"(20\d{6}|20\d{4})", text)
    return match.group(1).lower() if match else ""


def barcode_token(text: str) -> str:
    match = re.search(r"(?:barcode|bc)(\d+)", text, flags=re.IGNORECASE)
    return f"barcode{match.group(1).zfill(2)}".lower() if match else ""


def kennel_series_token(text: str) -> str:
    match = re.search(r"\b([ab]\d+)_barcode\d+\b", text, flags=re.IGNORECASE)
    return match.group(1).lower() if match else ""


def is_all_kennel_barcode_name(text: str) -> bool:
    return bool(re.search(r"\b[ab]\d+_barcode\d+\b", text, flags=re.IGNORECASE))


def score_candidate(assembly_row: dict[str, str], record: FastqRecord) -> tuple[int, list[str]]:
    sample_id = assembly_row["sample_id"]
    platform = assembly_row["platform"]
    sample_token = sanitize_token(sample_id)
    base_token = record.token
    score = 0
    reasons: list[str] = []

    if platform == "SRS" and record.read_type != "short":
        return -999, ["platform_mismatch"]
    if platform == "LRS" and record.read_type != "long":
        return -999, ["platform_mismatch"]

    if sample_token == base_token:
        score += 120
        reasons.append("exact_token_match")

    if sample_token in base_token or base_token in sample_token:
        score += 60
        reasons.append("substring_match")

    sample_date = date_token(sample_id)
    read_date = date_token(record.basename)
    if sample_date and read_date and sample_date == read_date:
        score += 40
        reasons.append("date_match")

    sample_barcode = barcode_token(sample_id)
    read_barcode = barcode_token(record.basename)
    if sample_barcode and read_barcode and sample_barcode == read_barcode:
        score += 80
        reasons.append("barcode_match")

    sample_series = kennel_series_token(sample_id)
    read_series = kennel_series_token(record.basename)
    if sample_series and read_series and sample_series == read_series:
        score += 60
        reasons.append("series_match")

    lower_sample = sample_id.lower()
    lower_base = record.basename.lower()

    if platform == "LRS":
        sample_is_all_kennel = is_all_kennel_barcode_name(sample_id)
        read_is_all_kennel = is_all_kennel_barcode_name(record.basename)
        read_is_serteperti = "serteperti" in lower_base
        sample_is_serteperti = "serteperti" in lower_sample

        if sample_is_all_kennel:
            if not (read_is_all_kennel or read_is_serteperti):
                return -999, ["wrong_project_context"]
            if sample_barcode and read_barcode and sample_barcode != read_barcode:
                return -999, ["barcode_mismatch"]
            if read_is_all_kennel:
                score += 25
                reasons.append("all_kennel_context")
                if record.relative_path.startswith("All_Kennel_ONT_WGS/"):
                    score += 15
                    reasons.append("local_placeholder_preferred")
            if read_is_serteperti and read_barcode:
                score += 35
                reasons.append("barcode_project_context")
        else:
            if read_is_serteperti:
                score += 20
                reasons.append("serteperti_context")
            if "concatenated" in lower_base:
                score += 10
                reasons.append("concatenated_context")
            if sample_is_serteperti and sample_barcode and read_barcode and sample_barcode == read_barcode:
                score += 80
                reasons.append("barcode_match")

    if platform == "SRS":
        if "merged" in lower_sample and "merged" in lower_base:
            score += 15
            reasons.append("merged_context")
        if lower_sample.startswith("toti_") and "toti" in lower_base:
            score += 20
            reasons.append("toti_context")
        if lower_sample.startswith("dmd_") and lower_base.startswith("dmd_"):
            score += 20
            reasons.append("dmd_context")
        sample_tokens = {token for token in sample_token.split("_") if token}
        base_tokens = {token for token in paired_group_key(record.basename).split("_") if token}
        overlap = len(sample_tokens & base_tokens)
        if overlap >= 3:
            score += overlap * 20
            reasons.append("token_overlap")
        elif overlap >= 2:
            score += overlap * 10
            reasons.append("token_overlap")
        if record.read_layout == "paired":
            score += 40
            reasons.append("paired_member")

    return score, reasons


def classify_match(score: int, candidate_count: int) -> tuple[str, str]:
    if score < 0 or candidate_count == 0:
        return "unmatched", "no_usable_candidates"
    if score >= 120 and candidate_count == 1:
        return "matched", "high"
    if score >= 90:
        return "matched", "medium"
    if score >= 60:
        return "review", "low"
    return "unmatched", "no_confident_match"


def choose_links(
    assembly_rows: list[dict[str, str]],
    fastq_records: list[FastqRecord],
    *,
    link_mode: str,
) -> list[dict[str, str]]:
    links: list[dict[str, str]] = []

    for row in assembly_rows:
        if row.get("status") not in {"ok", "skipped_existing"}:
            links.append(
                {
                    "assembly_id": row["assembly_id"],
                    "sample_id": row["sample_id"],
                    "platform": row["platform"],
                    "assembler": row["assembler"],
                    "match_status": "skipped_assembly",
                    "confidence": "na",
                    "read_type": "",
                    "read_layout": "",
                    "read_id": "",
                    "read_path": "",
                    "candidate_count": "0",
                    "score": "",
                    "match_reason": f"assembly_status={row.get('status', '')}",
                }
            )
            continue

        scored: list[tuple[int, list[str], FastqRecord]] = []
        for record in fastq_records:
            score, reasons = score_candidate(row, record)
            if score > 0:
                scored.append((score, reasons, record))

        scored.sort(key=lambda item: (-item[0], item[2].relative_path))
        candidate_count = len(scored)
        top_score = scored[0][0] if scored else -1
        match_status, confidence = classify_match(top_score, candidate_count)

        if match_status == "unmatched":
            links.append(
                {
                    "assembly_id": row["assembly_id"],
                    "sample_id": row["sample_id"],
                    "platform": row["platform"],
                    "assembler": row["assembler"],
                    "match_status": match_status,
                    "confidence": confidence,
                    "read_type": "",
                    "read_layout": "",
                    "read_id": "",
                    "read_path": "",
                    "candidate_count": str(candidate_count),
                    "score": str(top_score if top_score >= 0 else ""),
                    "match_reason": ";".join(scored[0][1]) if scored else "none",
                }
            )
            continue

        if row["platform"] == "SRS":
            paired_groups: dict[str, list[tuple[int, list[str], FastqRecord]]] = {}
            for item in scored:
                record = item[2]
                if record.read_layout == "paired":
                    paired_groups.setdefault(paired_group_key(record.basename), []).append(item)

            complete_groups: list[tuple[int, list[tuple[int, list[str], FastqRecord]]]] = []
            for _, group_items in paired_groups.items():
                mate_labels = {item[2].mate_label for item in group_items}
                if {"R1", "R2"}.issubset(mate_labels):
                    group_best = max(item[0] for item in group_items) + 100
                    complete_groups.append((group_best, group_items))

            if complete_groups:
                complete_groups.sort(key=lambda item: -item[0])
                best_group_score = complete_groups[0][0]
                if link_mode == "best":
                    chosen_groups = [complete_groups[0][1]]
                else:
                    chosen_groups = [
                        group_items
                        for group_score, group_items in complete_groups
                        if group_score == best_group_score
                    ]
                for group_items in chosen_groups:
                    picked: dict[str, tuple[int, list[str], FastqRecord]] = {}
                    for item in sorted(group_items, key=lambda value: (-value[0], value[2].relative_path)):
                        mate = item[2].mate_label
                        if mate not in picked:
                            picked[mate] = item
                    for mate in ("R1", "R2"):
                        if mate not in picked:
                            continue
                        best = picked[mate]
                        links.append(
                            {
                                "assembly_id": row["assembly_id"],
                                "sample_id": row["sample_id"],
                                "platform": row["platform"],
                                "assembler": row["assembler"],
                                "match_status": match_status,
                                "confidence": confidence,
                                "read_type": best[2].read_type,
                                "read_layout": best[2].read_layout,
                                "read_id": best[2].read_id,
                                "read_path": str(best[2].source_path),
                                "candidate_count": str(candidate_count),
                                "score": str(best[0]),
                                "match_reason": ";".join(best[1] + ["paired_group_selected"]),
                            }
                        )
            else:
                best = scored[0]
                links.append(
                    {
                        "assembly_id": row["assembly_id"],
                        "sample_id": row["sample_id"],
                        "platform": row["platform"],
                        "assembler": row["assembler"],
                        "match_status": match_status,
                        "confidence": confidence,
                        "read_type": best[2].read_type,
                        "read_layout": best[2].read_layout,
                        "read_id": best[2].read_id,
                        "read_path": str(best[2].source_path),
                        "candidate_count": str(candidate_count),
                        "score": str(best[0]),
                        "match_reason": ";".join(best[1]),
                    }
                )
        elif link_mode == "best":
            best = scored[0]
            links.append(
                {
                    "assembly_id": row["assembly_id"],
                    "sample_id": row["sample_id"],
                    "platform": row["platform"],
                    "assembler": row["assembler"],
                    "match_status": match_status,
                    "confidence": confidence,
                    "read_type": best[2].read_type,
                    "read_layout": best[2].read_layout,
                    "read_id": best[2].read_id,
                    "read_path": str(best[2].source_path),
                    "candidate_count": str(candidate_count),
                    "score": str(best[0]),
                    "match_reason": ";".join(best[1]),
                }
            )
        else:
            # Long reads: include all equal-top-score candidates when requested.
            top_candidates = [item for item in scored if item[0] == top_score]
            for best in top_candidates:
                links.append(
                    {
                        "assembly_id": row["assembly_id"],
                        "sample_id": row["sample_id"],
                        "platform": row["platform"],
                        "assembler": row["assembler"],
                        "match_status": match_status,
                        "confidence": confidence,
                        "read_type": best[2].read_type,
                        "read_layout": best[2].read_layout,
                        "read_id": best[2].read_id,
                        "read_path": str(best[2].source_path),
                        "candidate_count": str(candidate_count),
                        "score": str(best[0]),
                        "match_reason": ";".join(best[1]),
                    }
                )
    return links


def choose_links_chunk(
    assembly_rows: list[dict[str, str]],
    fastq_records: list[FastqRecord],
    link_mode: str,
) -> list[dict[str, str]]:
    return choose_links(assembly_rows, fastq_records, link_mode=link_mode)


def chunk_rows(rows: list[dict[str, str]], chunk_count: int) -> list[list[dict[str, str]]]:
    if not rows:
        return []
    chunk_count = max(1, min(chunk_count, len(rows)))
    chunk_size = (len(rows) + chunk_count - 1) // chunk_count
    return [rows[index : index + chunk_size] for index in range(0, len(rows), chunk_size)]


def choose_links_parallel(
    assembly_rows: list[dict[str, str]],
    fastq_records: list[FastqRecord],
    *,
    link_mode: str,
    jobs: int,
) -> list[dict[str, str]]:
    if jobs <= 1 or len(assembly_rows) <= 1:
        return choose_links(assembly_rows, fastq_records, link_mode=link_mode)

    chunks = chunk_rows(assembly_rows, jobs)
    if len(chunks) == 1:
        return choose_links(assembly_rows, fastq_records, link_mode=link_mode)

    links: list[dict[str, str]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=min(jobs, len(chunks))) as executor:
        results = executor.map(
            choose_links_chunk,
            chunks,
            [fastq_records] * len(chunks),
            [link_mode] * len(chunks),
        )
        for chunk_links in results:
            links.extend(chunk_links)
    return links


def build_read_manifest_rows(fastq_root: Path, records: list[FastqRecord]) -> list[dict[str, str]]:
    return [
        {
            "read_id": record.read_id,
            "basename": record.basename,
            "source_path": str(record.source_path),
            "relative_path": record.relative_path,
            "top_group": record.top_group,
            "read_type": record.read_type,
            "read_layout": record.read_layout,
            "mate_label": record.mate_label,
            "token": record.token,
        }
        for record in records
    ]


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        print("--jobs must be at least 1", file=sys.stderr)
        return 1

    assembly_manifest = Path(args.assembly_manifest).resolve() if args.assembly_manifest else None
    assembly_root = Path(args.assembly_root).resolve() if args.assembly_root else None
    fastq_root = Path(args.fastq_root).resolve()
    output_root = Path(args.output_root).resolve()
    extra_fastq_lists = [Path(path).resolve() for path in args.extra_fastq_list]

    if not fastq_root.is_dir():
        print(f"Missing FASTQ root: {fastq_root}", file=sys.stderr)
        return 1

    default_manifest = (Path.cwd() / "work" / "assembly_manifest.tsv").resolve()
    if assembly_manifest is None and assembly_root is None:
        if default_manifest.is_file():
            assembly_manifest = default_manifest
        else:
            assembly_root = (Path.cwd() / "final_assemblies").resolve()

    if assembly_manifest is not None and not assembly_manifest.is_file():
        print(f"Missing assembly manifest: {assembly_manifest}", file=sys.stderr)
        return 1
    if assembly_root is not None and not assembly_root.is_dir():
        print(f"Missing assembly root: {assembly_root}", file=sys.stderr)
        return 1
    for manifest_path in extra_fastq_lists:
        if not manifest_path.is_file():
            print(f"Missing extra FASTQ list: {manifest_path}", file=sys.stderr)
            return 1

    default_extra_list = fastq_root / "fastq" / "all_kennel_fastq_files.txt"
    if default_extra_list.is_file() and default_extra_list.resolve() not in extra_fastq_lists:
        extra_fastq_lists.append(default_extra_list.resolve())

    metadata_dir = output_root / "metadata"
    logs_dir = output_root / "logs"
    ensure_dir(metadata_dir)
    ensure_dir(logs_dir)

    read_manifest_path = output_root / "read_manifest.tsv"
    link_manifest_path = output_root / "assembly_read_links.tsv"
    log_path = logs_dir / "read_linking.log"

    if not args.force and not args.refine:
        existing = [path for path in (read_manifest_path, link_manifest_path) if path.exists()]
        if existing and not args.dry_run:
            print(
                "Output file(s) already exist. Use --force to overwrite: "
                + ", ".join(str(path) for path in existing),
                file=sys.stderr,
            )
            return 1

    if assembly_manifest is not None:
        assembly_rows = read_tsv(assembly_manifest)
        assembly_source_desc = str(assembly_manifest)
    else:
        assembly_rows = discover_assembly_rows(assembly_root)
        assembly_source_desc = str(assembly_root)
        if not assembly_rows:
            print(f"No assembly files found under: {assembly_root}", file=sys.stderr)
            return 1
    fastq_records = discover_fastqs(fastq_root)
    seen_ids = {record.read_id for record in fastq_records}
    seen_sources = {str(record.source_path).lower() for record in fastq_records}
    for manifest_path in extra_fastq_lists:
        fastq_records.extend(
            discover_fastqs_from_list(
                manifest_path,
                seen_ids=seen_ids,
                seen_sources=seen_sources,
            )
        )
    read_manifest_rows = build_read_manifest_rows(fastq_root, fastq_records)
    link_rows = choose_links_parallel(
        assembly_rows,
        fastq_records,
        link_mode=args.link_mode,
        jobs=args.jobs,
    )

    if args.refine and link_manifest_path.exists():
        existing_link_rows = read_tsv(link_manifest_path)
        link_rows = merge_refined_links(existing_link_rows, link_rows)

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_message(
            log_handle,
            (
                f"[START] link_reads_to_assemblies.py assembly_source={assembly_source_desc} "
                f"fastq_root={fastq_root} output_root={output_root} dry_run={args.dry_run} "
                f"refine={args.refine} link_mode={args.link_mode} jobs={args.jobs} "
                f"extra_fastq_lists={','.join(str(path) for path in extra_fastq_lists) or 'none'}"
            ),
        )
        log_message(log_handle, f"[INFO] FASTQ records discovered: {len(fastq_records)}")
        log_message(log_handle, f"[INFO] Assemblies loaded: {len(assembly_rows)}")

        matched = 0
        review = 0
        unmatched = 0
        skipped = 0
        for row in link_rows:
            status = row["match_status"]
            if status == "matched":
                matched += 1
            elif status == "review":
                review += 1
            elif status == "unmatched":
                unmatched += 1
            else:
                skipped += 1

        log_message(log_handle, f"[INFO] Linked rows: {len(link_rows)}")
        log_message(log_handle, f"[INFO] Matched rows: {matched}")
        log_message(log_handle, f"[INFO] Review rows: {review}")
        log_message(log_handle, f"[INFO] Unmatched rows: {unmatched}")
        log_message(log_handle, f"[INFO] Skipped rows: {skipped}")

        preview = link_rows[:20]
        for row in preview:
            if row["read_path"]:
                log_message(
                    log_handle,
                    (
                        f"[LINK] {row['assembly_id']} -> {row['read_path']} "
                        f"(status={row['match_status']}, confidence={row['confidence']}, score={row['score']})"
                    ),
                )
            else:
                log_message(
                    log_handle,
                    f"[LINK] {row['assembly_id']} -> <none> (status={row['match_status']}, reason={row['match_reason']})",
                )
        if len(link_rows) > len(preview):
            log_message(log_handle, f"[INFO] Preview truncated at {len(preview)} rows")

        if args.dry_run:
            log_message(log_handle, f"[DRY-RUN] Would write {read_manifest_path}")
            log_message(log_handle, f"[DRY-RUN] Would write {link_manifest_path}")
        else:
            write_tsv(
                read_manifest_path,
                [
                    "read_id",
                    "basename",
                    "source_path",
                    "relative_path",
                    "top_group",
                    "read_type",
                    "read_layout",
                    "mate_label",
                    "token",
                ],
                read_manifest_rows,
            )
            write_tsv(
                link_manifest_path,
                [
                    "assembly_id",
                    "sample_id",
                    "platform",
                    "assembler",
                    "match_status",
                    "confidence",
                    "read_type",
                    "read_layout",
                    "read_id",
                    "read_path",
                    "candidate_count",
                    "score",
                    "match_reason",
                ],
                link_rows,
            )
            log_message(log_handle, f"[WRITE] {read_manifest_path}")
            log_message(log_handle, f"[WRITE] {link_manifest_path}")

    print(f"Log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
