#!/usr/bin/env python3
"""Build a clean, auditable DogMAG article repository from the working project.

The builder deliberately uses a whitelist. It does not copy raw reads, BAM/SAM
files, databases, the BASALT workspace, private crosswalks, or run logs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
from pathlib import Path


MANUSCRIPT = [
    "documents/dogmag_scientific_data_dogfirst_rewrite_20260815_reconciled.md",
    "documents/dogmag_scientific_data_dogfirst_rewrite_20260815_reconciled.docx",
]

FIGURES = [
    "figure1_dogfirst_workflow_20260727.svg",
    "figure2_assembly_metrics_status_20260727.svg",
    "figure3_basalt_reselection_drep_20260727.svg",
    "figure4_mag_attrition_quality_20260727.svg",
    "figure5_gtdb_taxonomy_20260727.svg",
    "figure5_gtdb_taxonomy_sankey_links_20260727.tsv",
    "figure5_gtdb_taxonomy_sankey_nodes_20260727.tsv",
    "figure6_refseq_vs_dogmag_mapping_20260729.svg",
    "figure6_refseq_vs_dogmag_source_data_20260729.tsv",
    "figure7_viral_candidates_20260727.svg",
    "README_figure_sources_20260727.md",
]

SUPPLEMENTARY = [
    "supplementary_table_1_dogfirst_bin_attrition_20260727.tsv",
    "supplementary_table_2_mag_candidate_representative_metadata_20260727.tsv",
    "supplementary_table_3_wetlab_provenance_20260727.tsv",
    "supplementary_table_4_reads_to_assemblies_20260727.tsv",
    "supplementary_table_5_viral_candidate_metadata_20260727.tsv",
    "supplementary_table_6_viral_candidate_quality_summary_20260727.tsv",
    "supplementary_table_7_refseq_vs_dogmag_read_recruitment_20260727.tsv",
    "supplementary_table_8_final_dogfirst_assembly_metrics_20260806.tsv",
    "supplementary_table_9_ena_mag_deposition_status_20260815.tsv",
]

WORKFLOW_SCRIPTS = {
    "assembly": [
        "coassembly_dogfirst_20260505/run_flye_lr_only_dogfirst.sh",
        "opera_ms_dogfirst_20260505/run_opera_ms_named_lrs_dogs.sh",
        "scripts/collect_flye_lr_only_dogfirst_assemblies.sh",
        "scripts/collect_hybrid_dogfirst_assemblies.sh",
        "scripts/prepare_dogfirst_final_assemblies_for_basalt.sh",
        "scripts/normalize_final_assemblies.py",
        "scripts/link_reads_to_assemblies.py",
        "scripts/assembly_metrics_fast.py",
    ],
    "basalt": [
        "scripts/patch_basalt_hybrid_reassembly_resume_guard.sh",
        "scripts/patch_basalt_hybrid_failed_bin_skip.sh",
        "scripts/build_basalt_hybrid_failed_bin_skip_list.sh",
        "scripts/monitor_basalt_hybrid_reassembly.sh",
        "scripts/refresh_basalt_checkm2.py",
        "scripts/resume_basalt_step7_partial_lr_mapping.sh",
        "scripts/run_basalt_manual_final_drep_dogfirst.sh",
        "scripts/make_basalt_rerun_self_contained.sh",
        "scripts/stop_basalt_run_processes.sh",
        "scripts/normalize_basalt_workspace.sh",
        "scripts/rebuild_bins_folder.py",
    ],
    "catalogue": [
        "scripts/summarize_basalt_binset.py",
        "scripts/run_catalog_unit_dereplication.py",
        "scripts/build_canmag_depletion_panels.sh",
        "scripts/package_dogfirst_article_update_files.sh",
    ],
    "taxonomy": [
        "scripts/run_mag_taxonomy_workflow.sh",
        "scripts/prepare_mag_taxonomy_inputs.py",
        "scripts/merge_gtdbtk_taxonomy.py",
        "scripts/make_drep95_gtdbtk_summary_from_full_run.py",
        "scripts/taxonomy_sankey.R",
    ],
    "viral": [
        "scripts/prepare_viral_contig_inputs.py",
        "scripts/run_viral_contig_workflow.sh",
        "scripts/extract_viral_candidate_set.py",
        "scripts/summarize_viral_candidates.py",
    ],
    "mapping": [
        "scripts/run_genome_16s_minitax_db_workflow.sh",
        "scripts/build_genome_16s_minitax_db.py",
        "waltham_mapping_pilot/minitax_DogMAG/scripts/plot_final_balanced46_svg.py",
    ],
    "ena": [
        "scripts/prepare_ena_mag_submission.py",
        "scripts/build_ena_fastq1_manifest_prjeb115259.py",
    ],
    "release": [
        "scripts/build_dogmag_public_repo.py",
    ],
}

ACCESSIONS = [
    "ENA/PRJEB115259/fastq1_PRJEB115259_long_reads_public_ids.tsv",
    "ENA/PRJEB115259/sample_registration_ERC000013/Webin-accessions-2026-06-26T16_02_54.611+01_00.txt",
    "ENA/PRJEB115259/mag_submission/sample_registration_ERC000047/mag_sample_accessions.PROD.tsv",
    "ENA/PRJEB115259/mag_submission/prepared_ani95/mag_assembly_accessions.PROD_PASS122.tsv",
]

FORBIDDEN_SUFFIXES = {
    ".fastq", ".fq", ".bam", ".sam", ".cram", ".bt2", ".bt2l", ".mmi",
    ".dmnd", ".tar", ".gz", ".zip", ".7z",
}
FORBIDDEN_TEXT = ("/home/mdbio/", "D:\\data\\CanMAG", "C:\\Users\\SZTE")
PUBLIC_PATH_REPLACEMENTS = (
    ("/path/to/DogMAG_workdir", "/path/to/DogMAG_workdir"),
    ("/path/to/BASALT", "/path/to/BASALT"),
    ("/path/to/basalt_env/bin/python", "/path/to/basalt_env/bin/python"),
    ("/path/to/flye_work", "/path/to/flye_work"),
    ("/path/to/minitax", "/path/to/minitax"),
    (r"/path/to/DogMAG_workdir", "/path/to/DogMAG_workdir"),
)


def copy_required(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def sanitize_tsv(source: Path, target: Path, kind: str) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r", encoding="utf-8", newline="") as src:
        reader = csv.DictReader(src, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No TSV header in {source}")
        with target.open("w", encoding="utf-8", newline="") as dst:
            writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if kind == "mag":
                    row["reselected_fasta_path"] = f"reselected_mq_candidates/{row['selected_fasta_name']}"
                    row["source_path"] = f"BASALT_final_binset/{row['source_fasta_name']}"
                elif kind == "viral":
                    assembly_id = row.get("context_assembly_id", "")
                    row["context_assembly_path"] = f"assemblies/{assembly_id}.fasta" if assembly_id else ""
                writer.writerow(row)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def make_script_portable(path: Path) -> None:
    """Remove workstation-specific defaults from a public script copy."""
    text = path.read_text(encoding="utf-8", errors="strict")
    for old, new in PUBLIC_PATH_REPLACEMENTS:
        text = text.replace(old, new)
    if path.name == "plot_final_balanced46_svg.py":
        text = text.replace(
            'ROOT = Path(r"/path/to/DogMAG_workdir")',
            "ROOT = Path(__file__).resolve().parents[2]",
        )
    path.write_text(text, encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_readme() -> str:
    return """# DogMAG

DogMAG is a reusable canine faecal metagenome resource comprising long-read
assemblies, quality-controlled and dereplicated metagenome-assembled genome
(MAG) catalogues, taxonomic metadata, viral/proviral candidate metadata, and
reproducible workflow scripts.

This repository accompanies the DogMAG Scientific Data manuscript. It contains
the manuscript, figures, supplementary tables, public accession mappings, and
the scripts required to reproduce the reported processing and summaries.

## Key catalogue outputs

- 30,556 polished or reassembled candidate versions evaluated by BASALT.
- 11,276 BASALT-selected bin/version records.
- 3,418 reselected medium-quality-or-better MAG candidates.
- 792 representatives at 99% ANI and 135 representatives at 95% ANI.
- GTDB-Tk taxonomy for all 792 strain-like representatives.
- 22,068 viral/proviral candidate rows before final quality filtering.

The sequence files themselves are not stored in Git. Raw reads, assemblies,
MAG FASTAs, and deposited genome records are available through the ENA projects
and article data record described in the manuscript and `accessions/`.

## Repository map

- `manuscript/`: reconciled Markdown and Word manuscript.
- `figures/`: publication figures and the source tables used for Figures 5 and 6.
- `supplementary_tables/`: the nine manuscript supplementary tables.
- `scripts/`: public workflow and reporting scripts grouped by analysis stage.
- `accessions/`: public ENA sample, read, and MAG accession mappings.
- `docs/DATA_SHARING_SCOPE.md`: what is and is not distributed here.
- `checksums/SHA256SUMS.tsv`: SHA-256 inventory of repository files.

## Reuse

The 95% ANI catalogue is intended for species/SGB-like catalogue analyses and
the 99% ANI catalogue for strain-like analyses. DogMAG can also be used with
miniTax to obtain rapid, detailed taxonomic profiles from canine metagenomes.

## Citation

Please cite the DogMAG data descriptor once published. Interim citation
metadata are provided in `CITATION.cff`.
"""


def build_scope() -> str:
    return """# Data sharing scope

## Included in this Git repository

- Final manuscript source and rendered DOCX.
- Publication figures and figure source tables.
- Public supplementary metadata tables.
- Public ENA accession mappings.
- Workflow, catalogue, taxonomy, viral, mapping, and deposition scripts.

## Available through ENA or the article data record

- Raw sequencing reads.
- The 81 long-read assemblies and linked assembly products.
- The 11,276-member final BASALT binset.
- The 3,418 reselected medium-quality-or-better candidate MAG FASTAs.
- The 792 dRep99 and 135 dRep95 representative MAG FASTAs.
- Viral/proviral candidate sequence FASTAs and depletion-panel FASTAs.

## Deliberately excluded from Git

- Raw FASTQ, BAM, SAM, and intermediate mapping files.
- The full BASALT working directory and all 30,556 candidate FASTA versions.
- Bowtie2, minimap2, DIAMOND, GTDB, CheckM2, and other databases or indexes.
- Temporary work directories, logs, checkpoints, and restart snapshots.
- Private dog-name crosswalks, internal barcode mappings, credentials, and test receipts.

Supplementary Table 2 preserves candidate-level provenance but uses portable,
repository-relative sequence locations. Supplementary Table 5 likewise uses
portable assembly paths rather than local server paths.
"""


def build_accessions_readme() -> str:
    return """# Public accession mappings

This directory contains public-safe mappings exported from ENA submission work.

- `long_read_manifest.tsv`: public coded sample and FASTQ manifest for PRJEB115259.
- `source_sample_accessions.txt`: production ENA source-sample accessions.
- `mag_sample_accessions.tsv`: production MAG BioSample accessions.
- `mag_assembly_accessions_pass122.tsv`: 122 accepted multi-contig MAG assemblies.
- Supplementary Table 9 records the complete 135-representative deposition status,
  including the 13 single-contig records requiring the separate ENA route.

Submission receipts, test-service records, credentials, and private crosswalks
are intentionally excluded.
"""


def write_checksums(target: Path) -> None:
    checksum_path = target / "checksums" / "SHA256SUMS.tsv"
    rows = []
    for path in sorted(p for p in target.rglob("*") if p.is_file() and p != checksum_path):
        rows.append((sha256(path), path.relative_to(target).as_posix(), path.stat().st_size))
    checksum_path.parent.mkdir(parents=True, exist_ok=True)
    with checksum_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sha256", "relative_path", "size_bytes"))
        writer.writerows(rows)


def audit(target: Path) -> list[str]:
    errors: list[str] = []
    files = [p for p in target.rglob("*") if p.is_file()]
    if len(list((target / "figures").glob("*.svg"))) != 7:
        errors.append("Expected exactly 7 SVG figures")
    if len(list((target / "supplementary_tables").glob("supplementary_table_*.tsv"))) != 9:
        errors.append("Expected exactly 9 supplementary TSV files")
    for path in files:
        rel = path.relative_to(target).as_posix()
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"Forbidden data/archive suffix: {rel}")
        if path.stat().st_size > 25 * 1024 * 1024:
            errors.append(f"File exceeds 25 MiB: {rel}")
        if "private" in path.name.lower() or "crosswalk" in path.name.lower():
            errors.append(f"Private/crosswalk filename: {rel}")
        if path.suffix.lower() in {".md", ".tsv", ".txt", ".py", ".sh", ".r", ".cff", ".svg"}:
            text = path.read_text(encoding="utf-8", errors="replace")
            for needle in FORBIDDEN_TEXT:
                # The manuscript records exact historical commands and execution
                # paths; workflow scripts themselves must remain portable.
                if (
                    needle in text
                    and "manuscript/" not in rel
                    and rel != "scripts/release/build_dogmag_public_repo.py"
                ):
                    errors.append(f"Local path '{needle}' remains in {rel}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--target", type=Path, default=Path(r"D:\data\DogMAG"))
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Finish or refresh only the builder-managed files in a non-empty target",
    )
    args = parser.parse_args()
    source = args.source.resolve()
    target = args.target.resolve()

    if target.exists() and any(target.iterdir()) and not args.resume:
        raise SystemExit(f"Refusing non-empty target: {target}")
    target.mkdir(parents=True, exist_ok=True)

    for relative in MANUSCRIPT:
        copy_required(source / relative, target / "manuscript" / Path(relative).name)
    for name in FIGURES:
        copy_required(source / "figures" / "dogmag_final_20260727" / name, target / "figures" / name)
    for name in SUPPLEMENTARY:
        src = source / "documents" / "supplementary_tables" / name
        dst = target / "supplementary_tables" / name
        if name.startswith("supplementary_table_2_"):
            sanitize_tsv(src, dst, "mag")
        elif name.startswith("supplementary_table_5_"):
            sanitize_tsv(src, dst, "viral")
        else:
            copy_required(src, dst)

    copied_scripts = []
    for group, paths in WORKFLOW_SCRIPTS.items():
        for relative in paths:
            src = source / relative
            if not src.is_file():
                print(f"WARNING: optional script missing: {relative}", file=sys.stderr)
                continue
            dst = target / "scripts" / group / src.name
            copy_required(src, dst)
            make_script_portable(dst)
            copied_scripts.append((group, src.name, relative))

    accession_targets = [
        "long_read_manifest.tsv",
        "source_sample_accessions.txt",
        "mag_sample_accessions.tsv",
        "mag_assembly_accessions_pass122.tsv",
    ]
    for relative, name in zip(ACCESSIONS, accession_targets):
        copy_required(source / relative, target / "accessions" / name)

    copy_required(source / "LICENSE", target / "LICENSE")
    write_text(target / "README.md", build_readme())
    write_text(target / "docs" / "DATA_SHARING_SCOPE.md", build_scope())
    write_text(target / "accessions" / "README.md", build_accessions_readme())
    write_text(target / ".gitignore", """# Generated/local files
__pycache__/
*.py[cod]
.Rhistory
.RData
.DS_Store
Thumbs.db
*.log
*.tmp
*.bak

# Sequence and alignment data belong in ENA/the article data record
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.bam
*.sam
*.cram
*.bt2
*.bt2l
*.mmi
*.dmnd
""")
    write_text(target / "CITATION.cff", """cff-version: 1.2.0
message: "If you use DogMAG, please cite the accompanying Scientific Data article."
title: "DogMAG: a reusable canine gut metagenome resource"
type: dataset
authors:
  - family-names: Kakuk
    given-names: Balazs
repository-code: "https://github.com/Balays/DogMAG"
license: MIT
""")
    script_rows = ["group\tpublic_filename\tsource_worktree_path"]
    script_rows.extend(f"{group}\t{name}\t{relative}" for group, name, relative in copied_scripts)
    write_text(target / "scripts" / "SCRIPT_INDEX.tsv", "\n".join(script_rows))

    old_audit = target / "checksums" / "BUILD_AUDIT.tsv"
    if old_audit.exists():
        old_audit.unlink()
    write_checksums(target)
    errors = audit(target)
    audit_lines = [
        "DogMAG public repository audit",
        f"files\t{sum(1 for p in target.rglob('*') if p.is_file())}",
        f"svg_figures\t{len(list((target / 'figures').glob('*.svg')))}",
        f"supplementary_tables\t{len(list((target / 'supplementary_tables').glob('supplementary_table_*.tsv')))}",
        f"workflow_scripts\t{len(copied_scripts)}",
        f"status\t{'FAIL' if errors else 'PASS'}",
    ]
    audit_lines.extend(f"error\t{error}" for error in errors)
    write_text(target / "checksums" / "BUILD_AUDIT.tsv", "\n".join(audit_lines))
    print("\n".join(audit_lines))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
