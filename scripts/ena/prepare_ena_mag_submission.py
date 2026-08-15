#!/usr/bin/env python3
"""
Prepare DogMAG ANI95 MAG files for ENA submission.

This helper does NOT submit anything. It prepares reproducible staging files for:

1. Derived MAG sample registration planning.
2. Sanitised MAG FASTA copies with ENA-safe sequence names.
3. One Webin-CLI genome-context assembly manifest per MAG.

Default public MAG aliases are DOGMAG_ANI95_MAG001, DOGMAG_ANI95_MAG002, ... .
The corresponding sanitised contig names are DOGMAG_ANI95_MAG001_c000001, etc.

ENA treats MAGs as derived assemblies. Each MAG should have its own derived MAG
sample, normally using the GSC MIMAGS checklist, with a "sample derived from"
attribute pointing back to the source faecal/environmental sample(s) or run(s).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

IUPAC_DNA = set("ACGTUBDHKMNRSVWY")
FASTA_SUFFIXES = (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa")


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def first_existing(row: dict, names: Sequence[str], default: str = "") -> str:
    for name in names:
        if name in row and str(row[name]).strip() != "":
            return str(row[name]).strip()
    return default


def basename_any_path(value: str) -> str:
    return value.replace("\\", "/").rstrip("/").split("/")[-1]


def strip_fasta_suffixes(value: str) -> str:
    changed = True
    out = value
    while changed:
        changed = False
        for suffix in FASTA_SUFFIXES + (".gz",):
            if out.endswith(suffix):
                out = out[: -len(suffix)]
                changed = True
                break
    return out


def metadata_key_variants(value: str) -> List[str]:
    value = str(value).strip()
    if not value:
        return []
    value_slash = value.replace("\\", "/")
    base = basename_any_path(value)
    variants = [
        value,
        value_slash,
        base,
        strip_fasta_suffixes(value),
        strip_fasta_suffixes(value_slash),
        strip_fasta_suffixes(base),
    ]
    out: List[str] = []
    seen = set()
    for v in variants:
        if v and v not in seen:
            out.append(v)
            seen.add(v)
    return out


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    name = None
    seq_chunks: List[str] = []
    with open_text_maybe_gzip(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if name is not None:
        yield name, "".join(seq_chunks)


def write_fasta_gz(path: Path, records: Sequence[Tuple[str, str]], width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024 * 32), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_accessions(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for row in read_tsv(path):
        if row.get("TYPE") == "SAMPLE":
            alias = row.get("ALIAS", "").strip()
            acc = row.get("ACCESSION", "").strip()
            if alias and acc:
                mapping[alias] = acc
    return mapping


def build_assembly_sources(reads_to_assemblies: Path, sample_accessions: Path) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, List[str]]]:
    alias_to_ers = parse_accessions(sample_accessions)
    asm_to_ers: Dict[str, set] = defaultdict(set)
    asm_to_alias: Dict[str, set] = defaultdict(set)
    asm_to_platform: Dict[str, set] = defaultdict(set)

    for row in read_tsv(reads_to_assemblies):
        asm = first_existing(row, ["public_assembly_id", "assembly_id", "assembly", "assembly_file_public"])
        sample_alias = first_existing(row, ["public_sample_id", "sample", "sample_alias"])
        platform = first_existing(row, ["instrument_model", "read_platform", "platform"])
        if not asm or not sample_alias:
            continue
        asm_to_alias[asm].add(sample_alias)
        if sample_alias in alias_to_ers:
            asm_to_ers[asm].add(alias_to_ers[sample_alias])
        if platform:
            asm_to_platform[asm].add(platform)

    return (
        {k: sorted(v) for k, v in asm_to_ers.items()},
        {k: sorted(v) for k, v in asm_to_alias.items()},
        {k: sorted(v) for k, v in asm_to_platform.items()},
    )


def collect_mag_fastas(mag_fasta_dir: Path) -> List[Path]:
    return sorted(p for p in mag_fasta_dir.iterdir() if p.is_file() and any(str(p).endswith(s) for s in FASTA_SUFFIXES))


def read_mag_metadata(path: Optional[Path]) -> Dict[str, dict]:
    if path is None:
        return {}
    rows = read_tsv(path)
    out: Dict[str, dict] = {}
    key_columns = [
        "mag_alias",
        "mag_id",
        "genome_id",
        "genome",
        "genome_path",
        "user_genome",
        "bin_id",
        "Name",
        "name",
        "fasta",
        "filename",
        "file_name",
        "source_fasta",
        "source_basalt_quality_key",
    ]
    for row in rows:
        for col in key_columns:
            for key in metadata_key_variants(row.get(col, "")):
                out.setdefault(key, row)
    return out


def find_metadata_for_fasta(metadata: Dict[str, dict], fasta: Path) -> Tuple[dict, str]:
    candidates = []
    for value in [str(fasta), fasta.name, fasta.stem]:
        candidates.extend(metadata_key_variants(value))
    for candidate in candidates:
        if candidate in metadata:
            return metadata[candidate], candidate
    return {}, ""


def infer_source_assembly(row: dict) -> str:
    return first_existing(
        row,
        [
            "public_assembly_id",
            "assembly_id",
            "source_assembly_id",
            "parent_assembly_id",
            "assembly_file_public",
            "basalt_assembly",
            "assembly_file",
        ],
    )


def infer_source_basalt_quality_key(row: dict) -> str:
    return first_existing(row, ["source_basalt_quality_key", "quality_key", "bin_quality_key"])


def infer_taxonomy_text(row: dict) -> str:
    return first_existing(row, ["classification", "gtdb_taxonomy", "taxonomy", "gtdbtk_classification", "GTDB_taxonomy"])


def infer_gtdb_species(taxonomy: str) -> str:
    if not taxonomy:
        return ""
    for part in reversed(taxonomy.split(";")):
        part = part.strip()
        if part.startswith("s__"):
            return part[3:]
    return ""


def infer_assembly_program(row: dict, default_program: str) -> str:
    return first_existing(row, ["assembler", "assembly_program", "program", "PROGRAM"], default_program)


def infer_coverage(row: dict, default_coverage: str) -> str:
    return first_existing(row, ["coverage", "mean_coverage", "cov", "COVERAGE"], default_coverage)


def infer_completeness(row: dict) -> str:
    return first_existing(row, ["completeness", "Completeness", "checkm2_completeness", "completeness_percent"])


def infer_contamination(row: dict) -> str:
    return first_existing(row, ["contamination", "Contamination", "checkm2_contamination", "contamination_percent"])


def infer_n50(row: dict) -> str:
    return first_existing(row, ["n50", "N50"])


def infer_genome_size(row: dict) -> str:
    return first_existing(row, ["genome_size", "Genome_Size", "size"])


@dataclass
class SanitiseResult:
    mag_alias: str
    source_fasta: str
    sanitized_fasta: str
    contig_count: int
    total_bp: int
    md5: str
    warnings: List[str]


def sanitise_fasta(source: Path, dest: Path, mag_alias: str, trim_terminal_ns: bool) -> SanitiseResult:
    records: List[Tuple[str, str]] = []
    warnings: List[str] = []
    total_bp = 0
    seen = set()

    for idx, (_old_name, seq) in enumerate(iter_fasta(source), start=1):
        seq = seq.upper().replace(" ", "").replace("\t", "")
        if trim_terminal_ns:
            trimmed = seq.strip("N")
            if len(trimmed) != len(seq):
                warnings.append("terminal_Ns_trimmed")
            seq = trimmed
        elif seq.startswith("N") or seq.endswith("N"):
            warnings.append("terminal_Ns_present")

        invalid = sorted(set(seq) - IUPAC_DNA)
        if invalid:
            raise ValueError(f"{source}: invalid bases {invalid}")
        if len(seq) < 20:
            warnings.append("short_sequence_lt20bp_removed")
            continue

        name = f"{mag_alias}_c{idx:06d}"
        if len(name) >= 33:
            raise ValueError(f"Sequence name too long for ENA-style validation: {name}")
        if name in seen:
            raise ValueError(f"Duplicate generated sequence name: {name}")
        seen.add(name)
        records.append((name, seq))
        total_bp += len(seq)

    if len(records) == 0:
        raise ValueError(f"{source}: no valid FASTA records")
    if len(records) == 1:
        warnings.append("single_contig_MAG_check_ENA_chromosome_or_helpdesk")

    write_fasta_gz(dest, records)
    return SanitiseResult(
        mag_alias=mag_alias,
        source_fasta=str(source),
        sanitized_fasta=str(dest),
        contig_count=len(records),
        total_bp=total_bp,
        md5=md5_file(dest),
        warnings=sorted(set(warnings)),
    )


def write_manifest(path: Path, values: Dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    order = [
        "STUDY",
        "SAMPLE",
        "RUN_REF",
        "ASSEMBLYNAME",
        "ASSEMBLY_TYPE",
        "COVERAGE",
        "PROGRAM",
        "PLATFORM",
        "MINGAPLENGTH",
        "MOLECULETYPE",
        "DESCRIPTION",
        "FASTA",
    ]
    with path.open("w", encoding="utf-8") as out:
        for key in order:
            value = values.get(key, "")
            if value != "":
                out.write(f"{key}\t{value}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mag-fasta-dir", required=True, type=Path, help="Directory containing dereplicated MAG FASTAs.")
    parser.add_argument("--output-dir", type=Path, default=Path("ENA/PRJEB115259/mag_submission/prepared_ani95"))
    parser.add_argument("--study", default="PRJEB115259")
    parser.add_argument("--reads-to-assemblies", type=Path, default=Path("documents/supplementary_tables/supplementary_table_4_reads_to_assemblies.tsv"))
    parser.add_argument("--sample-accessions", type=Path, default=Path("ENA/PRJEB115259/sample_registration_ERC000013/Webin-accessions-2026-06-26T16_02_54.611+01_00.txt"))
    parser.add_argument("--mag-metadata", type=Path, default=None, help="Optional GTDB/CheckM/BASALT/dRep MAG metadata TSV.")
    parser.add_argument("--prefix", default="DOGMAG_ANI95_MAG")
    parser.add_argument("--default-program", default="BASALT, dRep")
    parser.add_argument("--default-coverage", default="1")
    parser.add_argument("--default-platform", default="Oxford Nanopore, Illumina")
    parser.add_argument("--mingaplength", default="10")
    parser.add_argument("--trim-terminal-ns", action="store_true")
    args = parser.parse_args()

    if not args.mag_fasta_dir.exists():
        raise SystemExit(f"MAG FASTA directory does not exist: {args.mag_fasta_dir}")
    if not args.reads_to_assemblies.exists():
        raise SystemExit(f"Missing reads-to-assemblies table: {args.reads_to_assemblies}")
    if not args.sample_accessions.exists():
        raise SystemExit(f"Missing source sample accession table: {args.sample_accessions}")

    fastas = collect_mag_fastas(args.mag_fasta_dir)
    if not fastas:
        raise SystemExit(f"No FASTA files found in {args.mag_fasta_dir}")

    metadata = read_mag_metadata(args.mag_metadata)
    asm_to_ers, asm_to_alias, asm_to_platform = build_assembly_sources(args.reads_to_assemblies, args.sample_accessions)

    fasta_out = args.output_dir / "fasta"
    manifests_out = args.output_dir / "assembly_manifests"
    review_rows: List[dict] = []
    sample_rows: List[dict] = []
    manifest_rows: List[dict] = []
    unmatched_metadata = 0

    for i, fasta in enumerate(fastas, start=1):
        mag_alias = f"{args.prefix}{i:03d}"
        row, metadata_match_key = find_metadata_for_fasta(metadata, fasta)
        if not row:
            unmatched_metadata += 1
        taxonomy = infer_taxonomy_text(row)
        gtdb_species = infer_gtdb_species(taxonomy)
        source_assembly = infer_source_assembly(row)
        source_basalt_quality_key = infer_source_basalt_quality_key(row)
        derived_from_ers = asm_to_ers.get(source_assembly, []) if source_assembly else []
        derived_from_aliases = asm_to_alias.get(source_assembly, []) if source_assembly else []
        platforms = asm_to_platform.get(source_assembly, []) if source_assembly else []
        platform = ", ".join(platforms) if platforms else args.default_platform
        program = infer_assembly_program(row, args.default_program)
        coverage = infer_coverage(row, args.default_coverage)
        completeness = infer_completeness(row)
        contamination = infer_contamination(row)
        n50 = infer_n50(row)
        genome_size = infer_genome_size(row)

        sanitized = fasta_out / f"{mag_alias}.fasta.gz"
        result = sanitise_fasta(fasta, sanitized, mag_alias, args.trim_terminal_ns)

        sample_description = (
            f"This sample represents a metagenome-assembled genome derived from canine faecal metagenomic sample(s): "
            f"{','.join(derived_from_ers) if derived_from_ers else 'REVIEW_SOURCE_SAMPLE'}"
        )
        sample_rows.append(
            {
                "mag_sample_alias": mag_alias,
                "sample_title": f"DogMAG ANI95 metagenome-assembled genome {mag_alias}",
                "sample_description": sample_description,
                "sample_derived_from": ",".join(derived_from_ers),
                "source_sample_aliases": ",".join(derived_from_aliases),
                "source_assembly_id": source_assembly,
                "source_basalt_quality_key": source_basalt_quality_key,
                "ncbi_tax_id": "REVIEW_NCBI_TAX_ID",
                "scientific_name": "REVIEW_NCBI_SCIENTIFIC_NAME",
                "gtdb_species": gtdb_species,
                "gtdb_taxonomy": taxonomy,
                "genome_size": genome_size,
                "completeness": completeness,
                "contamination": contamination,
                "n50": n50,
                "closest_genome_reference": first_existing(row, ["closest_genome_reference"]),
                "closest_genome_ani": first_existing(row, ["closest_genome_ani"]),
                "closest_genome_af": first_existing(row, ["closest_genome_af"]),
                "classification_method": first_existing(row, ["classification_method"]),
                "metadata_match_key": metadata_match_key,
                "binning_software": "BASALT",
                "dereplication_software": "dRep ANI95",
                "assembly_software": program,
                "notes": "Fill into ENA GSC MIMAGS checklist template; verify NCBI taxonomy and sample_derived_from before registration.",
            }
        )

        manifest_path = manifests_out / f"{mag_alias}.manifest.txt"
        write_manifest(
            manifest_path,
            {
                "STUDY": args.study,
                "SAMPLE": mag_alias,
                "ASSEMBLYNAME": mag_alias,
                "ASSEMBLY_TYPE": "Metagenome-Assembled Genome (MAG)",
                "COVERAGE": coverage,
                "PROGRAM": program,
                "PLATFORM": platform,
                "MINGAPLENGTH": args.mingaplength,
                "MOLECULETYPE": "genomic DNA",
                "DESCRIPTION": f"DogMAG ANI95 dereplicated canine faecal metagenome-assembled genome {mag_alias}.",
                "FASTA": str(sanitized.resolve()),
            },
        )

        review_rows.append(
            {
                "mag_alias": mag_alias,
                "source_fasta": str(fasta),
                "sanitized_fasta": str(sanitized),
                "sanitized_md5": result.md5,
                "contig_count": str(result.contig_count),
                "total_bp": str(result.total_bp),
                "genome_size": genome_size,
                "source_assembly_id": source_assembly,
                "source_basalt_quality_key": source_basalt_quality_key,
                "source_sample_accessions": ",".join(derived_from_ers),
                "source_sample_aliases": ",".join(derived_from_aliases),
                "platform": platform,
                "program": program,
                "coverage": coverage,
                "completeness": completeness,
                "contamination": contamination,
                "n50": n50,
                "gtdb_species": gtdb_species,
                "gtdb_taxonomy": taxonomy,
                "closest_genome_reference": first_existing(row, ["closest_genome_reference"]),
                "closest_genome_ani": first_existing(row, ["closest_genome_ani"]),
                "closest_genome_af": first_existing(row, ["closest_genome_af"]),
                "classification_method": first_existing(row, ["classification_method"]),
                "metadata_match_key": metadata_match_key,
                "warnings": ",".join(result.warnings),
                "manifest": str(manifest_path),
            }
        )
        manifest_rows.append({"mag_alias": mag_alias, "manifest": str(manifest_path), "fasta": str(sanitized)})

    sample_fields = [
        "mag_sample_alias",
        "sample_title",
        "sample_description",
        "sample_derived_from",
        "source_sample_aliases",
        "source_assembly_id",
        "source_basalt_quality_key",
        "ncbi_tax_id",
        "scientific_name",
        "gtdb_species",
        "gtdb_taxonomy",
        "genome_size",
        "completeness",
        "contamination",
        "n50",
        "closest_genome_reference",
        "closest_genome_ani",
        "closest_genome_af",
        "classification_method",
        "metadata_match_key",
        "binning_software",
        "dereplication_software",
        "assembly_software",
        "notes",
    ]
    review_fields = [
        "mag_alias",
        "source_fasta",
        "sanitized_fasta",
        "sanitized_md5",
        "contig_count",
        "total_bp",
        "genome_size",
        "source_assembly_id",
        "source_basalt_quality_key",
        "source_sample_accessions",
        "source_sample_aliases",
        "platform",
        "program",
        "coverage",
        "completeness",
        "contamination",
        "n50",
        "gtdb_species",
        "gtdb_taxonomy",
        "closest_genome_reference",
        "closest_genome_ani",
        "closest_genome_af",
        "classification_method",
        "metadata_match_key",
        "warnings",
        "manifest",
    ]
    write_tsv(args.output_dir / "mag_derived_sample_staging.tsv", sample_rows, sample_fields)
    write_tsv(args.output_dir / "mag_submission_review.tsv", review_rows, review_fields)
    write_tsv(args.output_dir / "mag_manifest_index.tsv", manifest_rows, ["mag_alias", "manifest", "fasta"])

    print(f"Prepared {len(fastas)} MAGs in {args.output_dir}")
    print(f"Metadata rows unmatched: {unmatched_metadata}")
    print(f"Review: {args.output_dir / 'mag_submission_review.tsv'}")
    print(f"MAG sample staging: {args.output_dir / 'mag_derived_sample_staging.tsv'}")
    print(f"Assembly manifests: {manifests_out}")
    print("NEXT: register derived MAG samples with the ENA GSC MIMAGS checklist, then validate manifests with Webin-CLI.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
