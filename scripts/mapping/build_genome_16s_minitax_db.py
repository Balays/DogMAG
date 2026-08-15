#!/usr/bin/env python3
"""Build a minitax-compatible 16S database from barrnap calls.

The companion shell script runs barrnap on each genome FASTA. This script reads
those GFF files, extracts 16S rRNA intervals from the original FASTAs, and
attaches taxonomy in the same table format used by minitax. It is intended for
any genome set with one FASTA per genome: raw BASALT MAGs, dRep95 panels,
dRep99 panels, or another curated genome collection.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
DB_FIELDS = ["seqnames", "taxid", *RANKS]
FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


CONFIG_TEMPLATE = [
    ("platform", "ONT", "both", "Either: 'Illumina', 'PacBio' or 'ONT'"),
    ("db", "{db_name}", "both", "16S rRNA reference database"),
    ("db.dir", "{db_dir}", "both", "minitax-compatible database directory"),
    ("project", "{db_name}", "optional", "project identifier"),
    ("Vregion", "16S", "optional", "16S rRNA reference sequences"),
    ("indir", "{indir}", "minimap2", "input FASTQ directory; edit before running minitax"),
    ("outdir", "{outdir}", "both", "minitax output directory"),
    ("debug", "TRUE", "both", "print planned samples while mapping"),
    ("mm2_path", "mm2-fast", "minimap2", "mm2-fast executable in the active minitax environment"),
    ("mm2_index", "{idx_name}", "minimap2", "minimap2 index filename relative to db.dir"),
    ("mm2_ref", "{fa_name}", "minimap2", "reference FASTA filename relative to db.dir"),
    ("mapper_backend", "mm2-fast", "minimap2", "force CPU mapper for reproducible first runs"),
    ("parabricks_image", "nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1", "minimap2", "Parabricks Docker image for NVIDIA GPU mapping"),
    ("parabricks_num_gpus", "1", "minimap2", "number of GPUs to request for Parabricks"),
    ("parabricks_extra_flags", "NA", "minimap2", "extra flags passed to pbrun minimap2"),
    ("mm2_index_batch", "16G", "minimap2", "minimap2 -I batch size"),
    ("fastq_pair_pattern", "_L001", "minimap2", "unused for ONT"),
    ("fastq_suffix", ".fastq.gz", "minimap2", "extension of ONT FASTQ files"),
    ("reads", "merged", "minimap2", "unused for ONT"),
    ("Nsec", "20", "minimap2", "number of secondary alignments to keep in minimap2"),
    ("nproc", "16", "both", "number of cores to use"),
    ("minitax.dir", "/path/to/minitax", "minitax", "minitax repository root"),
    ("misc.dir", "/path/to/minitax", "minitax", "helper-function directory"),
    ("metadata", "NA", "minitax", "metadata table; optional"),
    ("keep.highest.mapq.aln.only", "T", "minitax", "keep one highest-MAPQ alignment per read before CIGAR scoring"),
    ("crop.na.tax", "F", "minitax", "keep rank columns even when some ranks are NA"),
    ("multicore", "T", "minitax", "run minitax in multicore mode"),
    ("saveRAM", "F", "minitax", "keep intermediate objects in memory where possible"),
    ("mapq.filt", "NA", "minitax", "do not filter by MAPQ by default"),
    ("outputs", "bam.sum;best_alignments_w_taxa", "minitax", "per-sample mapping summary and cached best-alignment taxonomy"),
    ("methods", "BestAln;SpeciesEstimate", "minitax", "taxonomic annotation methods"),
    ("pardir", "NA", "minitax", "legacy option; minitax uses outdir/bam"),
    ("pattern", ".bam", "minitax", "legacy BAM suffix"),
    ("keep.max.cigar", "T", "minitax", "keep only alignments with the highest CIGAR score for each read"),
    (
        "CIGAR_points",
        "match_score = 1; mismatch_score = -3; insertion_score = -2; deletion_score = -2; gap_opening_penalty = -2; gap_extension_penalty = -1",
        "minitax",
        "rules for scoring alignments",
    ),
    ("BestAln_thresholds", "0.6", "minitax", "BestAln support threshold"),
    ("best.mapq", "T", "minitax", "use highest MAPQ alignment before CIGAR scoring"),
]


def clean(value: object) -> str:
    if value is None:
        return "NA"
    text = str(value).strip()
    return text if text and text.upper() not in {"N/A", "NA"} else "NA"


def safe_id(value: str) -> str:
    value = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", value, flags=re.IGNORECASE)
    value = re.sub(r"[^A-Za-z0-9_.:-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "unnamed"


def strip_fasta_suffix(name: str) -> str:
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", name, flags=re.IGNORECASE)


def genome_id_from_path(path: Path) -> str:
    stem = strip_fasta_suffix(path.name)
    if "__" in stem:
        left, right = stem.split("__", 1)
        right = strip_fasta_suffix(right)
        return safe_id(f"{left}_{right}")
    return safe_id(stem)


def candidate_taxonomy_keys(path_or_id: str) -> list[str]:
    base = Path(path_or_id).name
    no_ext = strip_fasta_suffix(base)
    candidates = [path_or_id, base, no_ext, safe_id(no_ext)]
    if "__" in no_ext:
        left, right = no_ext.split("__", 1)
        right = strip_fasta_suffix(right)
        candidates.extend([f"{left}_{right}", safe_id(f"{left}_{right}"), right, safe_id(right)])
    return list(dict.fromkeys(candidates))


def parse_gtdb_classification(classification: str) -> dict[str, str]:
    values = {rank: "NA" for rank in RANKS}
    prefix_to_rank = {
        "d__": "superkingdom",
        "p__": "phylum",
        "c__": "class",
        "o__": "order",
        "f__": "family",
        "g__": "genus",
        "s__": "species",
    }
    for part in clean(classification).split(";"):
        for prefix, rank in prefix_to_rank.items():
            if part.startswith(prefix):
                values[rank] = clean(part[len(prefix) :])
                break
    return values


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_taxonomy(path: Path) -> dict[str, dict[str, str]]:
    taxonomy: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        classification = row.get("classification", "")
        ranks = parse_gtdb_classification(classification)
        taxid = clean(row.get("genome_id") or row.get("user_genome") or row.get("file_name"))
        tax_row = {"taxid": safe_id(taxid), **ranks}
        keys = []
        for field in ("genome_id", "file_name", "genome_path", "user_genome"):
            value = row.get(field)
            if value:
                keys.extend(candidate_taxonomy_keys(value))
        for key in keys:
            taxonomy.setdefault(key, tax_row)
    return taxonomy


def fasta_paths(genome_dir: Path) -> list[Path]:
    files = [p for p in genome_dir.iterdir() if p.is_file() and p.name.lower().endswith(FASTA_SUFFIXES)]
    return sorted(files, key=lambda p: p.name)


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line)
    return {name: "".join(parts).upper() for name, parts in records.items()}


def parse_attributes(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in text.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def find_gff(path: Path, genome_id: str) -> Path | None:
    candidates = [
        path / f"{genome_id}.gff",
        path / f"{genome_id}.bac.gff",
        path / f"{genome_id}.arc.gff",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = sorted(path.glob(f"{genome_id}*.gff"))
    return matches[0] if matches else None


def parse_16s_hits(gff_path: Path, min_length: int) -> list[dict[str, str | int]]:
    hits = []
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs = parts
            attr = parse_attributes(attrs)
            feature_text = " ".join([feature_type, attrs]).lower()
            if "16s" not in feature_text:
                continue
            start_i = int(start)
            end_i = int(end)
            length = abs(end_i - start_i) + 1
            if length < min_length:
                continue
            hits.append(
                {
                    "contig": seqid,
                    "source": source,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "length": length,
                    "product": clean(attr.get("product") or attr.get("Name") or attr.get("note")),
                }
            )
    return hits


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
    return seq.translate(table)[::-1].upper()


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_config(path: Path, *, db_name: str, db_dir: Path, outdir: Path, fa_name: str, idx_name: str, indir: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variable", "value", "used_by", "description"])
        for variable, value, used_by, description in CONFIG_TEMPLATE:
            writer.writerow(
                [
                    variable,
                    value.format(
                        db_name=db_name,
                        db_dir=db_dir.as_posix(),
                        outdir=outdir.as_posix(),
                        fa_name=fa_name,
                        idx_name=idx_name,
                        indir=indir,
                    ),
                    used_by,
                    description,
                ]
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-dir", required=True, type=Path, help="Directory containing one genome FASTA per MAG/genome.")
    parser.add_argument(
        "--taxonomy",
        required=False,
        type=Path,
        help=(
            "Taxonomy TSV for the same genome set. GTDB-Tk summary-style tables are supported. "
            "If omitted, taxonomy columns are filled with NA."
        ),
    )
    parser.add_argument("--barrnap-dir", required=True, type=Path, help="Directory containing per-genome barrnap GFF files.")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output root for minitax database files.")
    parser.add_argument("--db-name", default="GenomeSet_16S", help="Database basename.")
    parser.add_argument("--min-length", type=int, default=1000, help="Minimum extracted 16S length to retain.")
    parser.add_argument("--indir", default="/path/to/fastq", help="Placeholder input FASTQ path for the minitax config.")
    args = parser.parse_args()

    genome_dir = args.genome_dir
    taxonomy_path = args.taxonomy
    barrnap_dir = args.barrnap_dir
    out_dir = args.out_dir
    db_dir = out_dir / "db"
    config_dir = out_dir / "configs"
    db_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)

    taxonomy = load_taxonomy(taxonomy_path) if taxonomy_path else {}
    fasta_files = fasta_paths(genome_dir)
    if not fasta_files:
        raise SystemExit(f"No FASTA files found in {genome_dir}")

    fasta_name = f"{args.db_name}.fa"
    idx_name = f"{args.db_name}.idx"
    fasta_out = db_dir / fasta_name
    db_rows: list[dict[str, str]] = []
    hit_rows: list[dict[str, object]] = []
    genome_summary: list[dict[str, object]] = []
    missing_taxonomy: list[dict[str, object]] = []
    missing_gff: list[dict[str, object]] = []
    skipped_hits: list[dict[str, object]] = []
    total_hits = 0

    with fasta_out.open("w", encoding="utf-8", newline="\n") as fasta_handle:
        for genome_path in fasta_files:
            genome_id = genome_id_from_path(genome_path)
            gff_path = find_gff(barrnap_dir, genome_id)
            if gff_path is None:
                missing_gff.append({"genome_id": genome_id, "file_name": genome_path.name})
                genome_summary.append(
                    {"genome_id": genome_id, "file_name": genome_path.name, "status": "missing_barrnap_gff", "retained_16s": 0}
                )
                continue

            tax_row = None
            for key in candidate_taxonomy_keys(genome_path.name) + candidate_taxonomy_keys(genome_id):
                tax_row = taxonomy.get(key)
                if tax_row:
                    break
            if tax_row is None:
                missing_taxonomy.append({"genome_id": genome_id, "file_name": genome_path.name})
                tax_row = {"taxid": genome_id, **{rank: "NA" for rank in RANKS}}

            contigs = read_fasta(genome_path)
            hits = parse_16s_hits(gff_path, args.min_length)
            retained = 0
            for hit in hits:
                contig = str(hit["contig"])
                seq = contigs.get(contig)
                if seq is None:
                    skipped_hits.append(
                        {
                            "genome_id": genome_id,
                            "file_name": genome_path.name,
                            "gff_file": gff_path.name,
                            "contig": contig,
                            "reason": "contig_not_found_in_fasta",
                        }
                    )
                    continue
                start = int(hit["start"])
                end = int(hit["end"])
                piece = seq[start - 1 : end]
                if str(hit["strand"]) == "-":
                    piece = revcomp(piece)
                retained += 1
                total_hits += 1
                seq_id = safe_id(f"{genome_id}__16S_{retained:02d}")
                fasta_handle.write(
                    f">{seq_id} genome={genome_id} contig={contig} start={start} end={end} strand={hit['strand']} length={len(piece)}\n"
                )
                fasta_handle.write(wrap(piece) + "\n")
                db_rows.append({"seqnames": seq_id, **tax_row})
                hit_rows.append(
                    {
                        "seqnames": seq_id,
                        "genome_id": genome_id,
                        "file_name": genome_path.name,
                        "gff_file": gff_path.name,
                        "contig": contig,
                        "start": start,
                        "end": end,
                        "strand": hit["strand"],
                        "length": len(piece),
                        "product": hit["product"],
                    }
                )
            genome_summary.append(
                {
                    "genome_id": genome_id,
                    "file_name": genome_path.name,
                    "status": "ok" if retained else "no_16s_retained",
                    "retained_16s": retained,
                    "raw_16s_hits_min_length": len(hits),
                }
            )

    write_tsv(db_dir / "db_data.tsv", db_rows, DB_FIELDS)
    write_tsv(db_dir / "MAG.db.tsv", db_rows, DB_FIELDS)
    write_tsv(db_dir / "MAG.db.uni.tsv", db_rows, DB_FIELDS)
    write_tsv(
        db_dir / "rrna_hits.tsv",
        hit_rows,
        ["seqnames", "genome_id", "file_name", "gff_file", "contig", "start", "end", "strand", "length", "product"],
    )
    write_tsv(db_dir / "genome_16s_summary.tsv", genome_summary, ["genome_id", "file_name", "status", "retained_16s", "raw_16s_hits_min_length"])
    write_tsv(db_dir / "missing_taxonomy.tsv", missing_taxonomy, ["genome_id", "file_name"])
    write_tsv(db_dir / "missing_barrnap_gff.tsv", missing_gff, ["genome_id", "file_name"])
    write_tsv(db_dir / "skipped_16s_hits.tsv", skipped_hits, ["genome_id", "file_name", "gff_file", "contig", "reason"])

    status_counts = Counter(str(row["status"]) for row in genome_summary)
    manifest = [
        {"key": "db_name", "value": args.db_name},
        {"key": "genome_dir", "value": genome_dir.as_posix()},
        {"key": "taxonomy", "value": taxonomy_path.as_posix() if taxonomy_path else "NA"},
        {"key": "barrnap_dir", "value": barrnap_dir.as_posix()},
        {"key": "fasta", "value": fasta_out.as_posix()},
        {"key": "db_data", "value": (db_dir / "db_data.tsv").as_posix()},
        {"key": "genomes_seen", "value": str(len(fasta_files))},
        {"key": "genomes_with_16s", "value": str(sum(1 for row in genome_summary if int(row["retained_16s"]) > 0))},
        {"key": "total_16s_sequences", "value": str(total_hits)},
        {"key": "min_length", "value": str(args.min_length)},
        {"key": "missing_taxonomy", "value": str(len(missing_taxonomy))},
        {"key": "missing_barrnap_gff", "value": str(len(missing_gff))},
    ]
    for status, count in sorted(status_counts.items()):
        manifest.append({"key": f"genome_status_{status}", "value": str(count)})
    write_tsv(db_dir / "build_manifest.tsv", manifest, ["key", "value"])

    config_path = config_dir / f"minitax_config_{args.db_name}.tsv"
    write_config(
        config_path,
        db_name=args.db_name,
        db_dir=db_dir,
        outdir=out_dir / "out",
        fa_name=fasta_name,
        idx_name=idx_name,
        indir=args.indir,
    )

    readme = out_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                f"# {args.db_name}",
                "",
                "Minitax-compatible 16S rRNA database built from genome FASTAs.",
                "",
                "Inputs:",
                f"- Genome FASTAs: `{genome_dir.as_posix()}`",
                f"- Taxonomy table: `{taxonomy_path.as_posix() if taxonomy_path else 'not supplied; taxonomy columns are NA'}`",
                f"- barrnap GFFs: `{barrnap_dir.as_posix()}`",
                "",
                "Main outputs:",
                f"- `{(db_dir / fasta_name).as_posix()}`",
                f"- `{(db_dir / 'MAG.db.tsv').as_posix()}`",
                f"- `{config_path.as_posix()}`",
                "",
                "Build the minimap2/mm2-fast index from the database directory:",
                "",
                "```bash",
                f"cd {db_dir.as_posix()}",
                f"mm2-fast -I 16G -d {idx_name} {fasta_name}",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Genomes scanned: {len(fasta_files)}")
    print(f"Genomes with retained 16S: {sum(1 for row in genome_summary if int(row['retained_16s']) > 0)}")
    print(f"16S sequences written: {total_hits}")
    print(f"Missing barrnap GFFs: {len(missing_gff)}")
    print(f"Missing taxonomy rows: {len(missing_taxonomy)}")
    print(f"Wrote {fasta_out}")
    print(f"Wrote {config_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
