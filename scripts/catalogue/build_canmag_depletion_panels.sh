#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# build_canmag_depletion_panels.sh
#
# Purpose
#   1. Re-select the best assembly version for every BASALT bin group from the
#      FULL candidate table (e.g. ~30,556 polished/reassembled candidates).
#   2. Keep only medium-quality-or-better genomes:
#         completeness >= 50
#         contamination <= 10
#   3. Prefer high-quality candidates:
#         completeness >= 90
#         contamination <= 5
#   4. Rank candidates within each original bin group by:
#         quality tier -> completeness - 5*contamination -> N50 -> size
#   5. Build independent 95% and 99% ANI dRep catalogues.
#   6. Produce ready-to-use FASTA panels for ONT adaptive depletion.
#
# Recommended use
#   - Use the 95% ANI panel as the primary discovery-oriented depletion set.
#   - Keep the 99% ANI panel for offline sensitivity analysis or a more
#     aggressive depletion experiment.
#
# Dependencies
#   python3, dRep, seqkit, minimap2
#
# Candidate table requirements
#   One row per candidate assembly version, with columns corresponding to:
#     candidate ID/name, completeness, contamination, N50, genome size,
#     and FASTA path or FASTA filename.
#
#   Common column names are auto-detected. Override them with the environment
#   variables shown below when needed.
#
# Example
#   CANDIDATE_TABLE=all_candidate_versions.tsv \
#   FASTA_ROOT=/path/to/all_candidate_fastas \
#   OUTDIR=CanMAG_depletion_build \
#   THREADS=52 \
#   bash build_canmag_depletion_panels.sh
###############################################################################

CANDIDATE_TABLE="${CANDIDATE_TABLE:-}"
FASTA_ROOT="${FASTA_ROOT:-}"
OUTDIR="${OUTDIR:-CanMAG_depletion_build}"
THREADS="${THREADS:-32}"
LINK_MODE="${LINK_MODE:-symlink}"       # symlink or copy
REUSE_EXISTING="${REUSE_EXISTING:-0}"   # 1 skips completed expensive stages
PANELS_ONLY="${PANELS_ONLY:-0}"         # 1 builds panels from existing dRep outputs only
FORCE_PANELS="${FORCE_PANELS:-0}"       # 1 rebuilds panel FASTA/MMI/BED even if present
MINIMAP2_BIN="${MINIMAP2_BIN:-${MM2FAST_BIN:-}}"
MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-1000}"
MIN_COMPLETENESS="${MIN_COMPLETENESS:-50}"
MAX_CONTAMINATION="${MAX_CONTAMINATION:-10}"
HQ_COMPLETENESS="${HQ_COMPLETENESS:-90}"
HQ_CONTAMINATION="${HQ_CONTAMINATION:-5}"

# Optional explicit column names. Leave empty to auto-detect.
ID_COLUMN="${ID_COLUMN:-}"
COMPLETENESS_COLUMN="${COMPLETENESS_COLUMN:-}"
CONTAMINATION_COLUMN="${CONTAMINATION_COLUMN:-}"
N50_COLUMN="${N50_COLUMN:-}"
SIZE_COLUMN="${SIZE_COLUMN:-}"
FASTA_COLUMN="${FASTA_COLUMN:-}"

usage() {
  cat <<EOF
Usage:
  CANDIDATE_TABLE=... FASTA_ROOT=... OUTDIR=... bash scripts/build_canmag_depletion_panels.sh [options]

Options:
  --reuse-existing        Reuse completed re-selection, dRep, and panel outputs when present.
  --panels-only           Do not reselect or run dRep; build panels from existing drep95/drep99.
  --force-panels          Rebuild panels even if panel FASTA/MMI/BED already exist.
  --minimap2-bin PATH     minimap2-compatible binary for panel indexing; use mm2-fast here.
  --threads INT           Threads for dRep/minimap2 index.
  --copy                  Copy selected FASTAs instead of symlinking.
  --help                  Show this message.

Environment equivalents:
  REUSE_EXISTING=1, PANELS_ONLY=1, FORCE_PANELS=1,
  MINIMAP2_BIN=/path/to/minimap2, MM2FAST_BIN=/path/to/mm2-fast/minimap2,
  THREADS, LINK_MODE, MIN_CONTIG_LEN, MIN_COMPLETENESS, MAX_CONTAMINATION.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reuse-existing)
      REUSE_EXISTING=1
      shift
      ;;
    --panels-only)
      PANELS_ONLY=1
      REUSE_EXISTING=1
      shift
      ;;
    --force-panels)
      FORCE_PANELS=1
      shift
      ;;
    --minimap2-bin)
      MINIMAP2_BIN="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --copy)
      LINK_MODE="copy"
      shift
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

find_minimap2_bin() {
  if [[ -n "$MINIMAP2_BIN" ]]; then
    printf '%s\n' "$MINIMAP2_BIN"
    return
  fi
  for candidate in \
    mm2-fast \
    minimap2-fast \
    /mnt/c/ubuntu/programs/mm2-fast/minimap2 \
    /mnt/c/programs/mm2-fast/minimap2 \
    /mnt/d/programs/mm2-fast/minimap2 \
    minimap2
  do
    if command -v "$candidate" >/dev/null 2>&1; then
      command -v "$candidate"
      return
    fi
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return
    fi
  done
  return 1
}

MINIMAP2_BIN="$(find_minimap2_bin)" || {
  echo "ERROR: Could not find minimap2/mm2-fast. Set MINIMAP2_BIN or MM2FAST_BIN." >&2
  exit 1
}

if [[ "$PANELS_ONLY" != "1" ]]; then
  : "${CANDIDATE_TABLE:?Set CANDIDATE_TABLE to the full candidate-version TSV}"
  : "${FASTA_ROOT:?Set FASTA_ROOT to the directory containing candidate FASTAs}"
fi

mkdir -p \
  "$OUTDIR"/{tables,reselected_fastas,drep95,drep99,panels,logs,tmp}

echo "OUTDIR: $OUTDIR"
echo "THREADS: $THREADS"
echo "MINIMAP2_BIN: $MINIMAP2_BIN"
echo "REUSE_EXISTING: $REUSE_EXISTING"
echo "PANELS_ONLY: $PANELS_ONLY"

if [[ "$PANELS_ONLY" != "1" ]]; then
  if [[ "$REUSE_EXISTING" == "1" \
      && -s "$OUTDIR/tables/reselected_fasta_manifest.tsv" ]] \
      && find "$OUTDIR/reselected_fastas" -maxdepth 1 \( -type f -o -type l \) \
        \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) -print -quit | grep -q .; then
    echo "Reusing existing reselected FASTAs: $OUTDIR/reselected_fastas"
  else
    python3 - "$CANDIDATE_TABLE" "$FASTA_ROOT" "$OUTDIR" \
      "$LINK_MODE" "$MIN_COMPLETENESS" "$MAX_CONTAMINATION" \
      "$HQ_COMPLETENESS" "$HQ_CONTAMINATION" \
      "$ID_COLUMN" "$COMPLETENESS_COLUMN" "$CONTAMINATION_COLUMN" \
      "$N50_COLUMN" "$SIZE_COLUMN" "$FASTA_COLUMN" <<'PY'
import csv
import ast
import os
import re
import shutil
import sys
from pathlib import Path

FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")

(
    table_path,
    fasta_root,
    outdir,
    link_mode,
    min_comp,
    max_cont,
    hq_comp,
    hq_cont,
    explicit_id,
    explicit_comp,
    explicit_cont,
    explicit_n50,
    explicit_size,
    explicit_fasta,
) = sys.argv[1:]

table_path = Path(table_path).resolve()
fasta_root = Path(fasta_root).resolve()
outdir = Path(outdir).resolve()
min_comp = float(min_comp)
max_cont = float(max_cont)
hq_comp = float(hq_comp)
hq_cont = float(hq_cont)

table_out = outdir / "tables"
fasta_out = outdir / "reselected_fastas"
table_out.mkdir(parents=True, exist_ok=True)
fasta_out.mkdir(parents=True, exist_ok=True)
FASTA_INDEX = None

ALIASES = {
    "id": [
        "candidate_id", "candidate", "genome", "genome_id", "bin_name",
        "bin_id", "name", "assembly", "assembly_id", "file", "filename"
    ],
    "completeness": [
        "completeness", "completeness_percent", "checkm_completeness",
        "checkm2_completeness", "complete"
    ],
    "contamination": [
        "contamination", "contamination_percent", "checkm_contamination",
        "checkm2_contamination", "contam"
    ],
    "n50": ["n50", "N50", "contig_n50", "genome_n50"],
    "size": [
        "genome_size", "size", "length", "total_length", "assembly_size",
        "genome_length"
    ],
    "fasta": [
        "fasta", "fasta_path", "path", "file_path", "genome_path",
        "assembly_path", "filename", "file"
    ],
}

def norm(s):
    return re.sub(r"[^a-z0-9]+", "_", s.strip().lower()).strip("_")

def is_fasta_name(name):
    lower = name.lower()
    return any(lower.endswith(suffix) for suffix in FASTA_SUFFIXES)

def build_fasta_index():
    index = {}
    for path in fasta_root.rglob("*"):
        if path.is_file() and is_fasta_name(path.name):
            resolved = path.resolve()
            keys = {path.name, path.stem}
            if path.name.lower().endswith(".gz"):
                without_gz = Path(path.stem)
                keys.add(without_gz.name)
                keys.add(without_gz.stem)
            for key in keys:
                index.setdefault(key, []).append(resolved)
    for key in list(index):
        index[key] = sorted(set(index[key]))
    return index

def lookup_fasta_index(names):
    global FASTA_INDEX
    if FASTA_INDEX is None:
        print(f"Indexing FASTA files under {fasta_root}", flush=True)
        FASTA_INDEX = build_fasta_index()
        print(f"Indexed {sum(len(v) for v in FASTA_INDEX.values()):,} FASTA name entries", flush=True)
    hits = []
    for name in names:
        hits.extend(FASTA_INDEX.get(name, []))
    hits = sorted(set(hits))
    if hits:
        return hits[0]
    return None

def detect_column(fieldnames, explicit, aliases, required=True):
    if explicit:
        if explicit not in fieldnames:
            raise SystemExit(
                f"Explicit column {explicit!r} not found. Available columns: "
                + ", ".join(fieldnames)
            )
        return explicit

    normalized = {norm(x): x for x in fieldnames}
    for alias in aliases:
        key = norm(alias)
        if key in normalized:
            return normalized[key]

    if required:
        raise SystemExit(
            "Could not auto-detect a required column. Available columns: "
            + ", ".join(fieldnames)
        )
    return None

def sniff_delimiter(path):
    sample = path.read_text(errors="replace")[:10000]
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return "\t"

def first_nonempty_line(path):
    with path.open(errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                return stripped
    return ""

def parse_basalt_dict_candidate_table(path):
    """
    Parse BASALT Reassembled_bins_comparison.txt style lines, for example:
      bin642_mag_polished.fa {'N50': 645326, 'Completeness': 22.09, ...}
    """
    parsed_rows = []
    parse_errors = []
    pattern = re.compile(r"^(\S+)\s+(\{.*\})\s*$")
    with path.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            match = pattern.match(stripped)
            if not match:
                parse_errors.append(f"line {line_no}: not BASALT dict format: {stripped[:200]}")
                continue
            candidate_id, payload = match.groups()
            try:
                values = ast.literal_eval(payload)
            except (SyntaxError, ValueError) as exc:
                parse_errors.append(f"line {line_no}: could not parse quality dict: {exc}")
                continue
            values = {norm(str(key)): value for key, value in values.items()}
            parsed_rows.append({
                "__candidate_id": candidate_id,
                "__completeness": values.get("completeness"),
                "__contamination": values.get("contamination"),
                "__n50": values.get("n50"),
                "__size": values.get("genome_size") or values.get("genome_length") or values.get("size"),
                "__fasta": candidate_id,
                "__source_line": line_no,
            })
    return parsed_rows, parse_errors

def to_float(value, default=0.0):
    if value is None:
        return default
    value = str(value).strip().replace("%", "").replace(",", "")
    if value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default

def original_bin_group(candidate_id):
    """
    Extract the original BASALT bin group.

    Expected examples:
      bin10203_mag_polished.fa      -> bin10203
      bin6068_spades_hybrid.fasta   -> bin6068
      projectX_bin4398_IDBA.fa      -> bin4398

    Modify this regex if your original bin identifiers use another pattern.
    """
    match = re.search(r"(?i)(?:^|[^A-Za-z0-9])(bin[._-]?\d+)(?:[^0-9]|$)",
                      candidate_id)
    if not match:
        match = re.search(r"(?i)(bin\d+)", candidate_id)
    if not match:
        raise ValueError(
            f"Could not extract original bin group from candidate ID: "
            f"{candidate_id!r}"
        )
    return re.sub(r"[._-]", "", match.group(1).lower())

def resolve_fasta(row, candidate_id, fasta_column):
    candidates = []

    if fasta_column and row.get(fasta_column):
        raw = Path(str(row[fasta_column]).strip())
        candidates.append(raw)
        if not raw.is_absolute():
            candidates.append(fasta_root / raw)

    raw_id = Path(candidate_id)
    candidates.extend([
        raw_id,
        fasta_root / raw_id,
        fasta_root / f"{candidate_id}.fa",
        fasta_root / f"{candidate_id}.fna",
        fasta_root / f"{candidate_id}.fasta",
        fasta_root / f"{candidate_id}.fa.gz",
        fasta_root / f"{candidate_id}.fna.gz",
        fasta_root / f"{candidate_id}.fasta.gz",
    ])

    for path in candidates:
        path = path.expanduser()
        if path.exists() and path.is_file():
            return path.resolve()

    # Last-resort lookup by basename/stem. Build the recursive index once
    # instead of repeatedly walking a huge BASALT run directory.
    names = {
        raw_id.name,
        raw_id.stem,
        f"{candidate_id}.fa",
        f"{candidate_id}.fna",
        f"{candidate_id}.fasta",
        f"{candidate_id}.fa.gz",
        f"{candidate_id}.fna.gz",
        f"{candidate_id}.fasta.gz",
    }
    indexed_hit = lookup_fasta_index(names)
    if indexed_hit:
        return indexed_hit
    raise FileNotFoundError(f"No FASTA found for candidate {candidate_id!r}")

rows = []
errors = []

first_line = first_nonempty_line(table_path)
if re.match(r"^\S+\s+\{.*\}\s*$", first_line):
    raw_rows, errors = parse_basalt_dict_candidate_table(table_path)
    id_col = "__candidate_id"
    comp_col = "__completeness"
    cont_col = "__contamination"
    n50_col = "__n50"
    size_col = "__size"
    fasta_col = "__fasta"
    raw_iter = ((int(row.get("__source_line", index)), row) for index, row in enumerate(raw_rows, start=1))
else:
    delimiter = sniff_delimiter(table_path)
    with table_path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise SystemExit("Candidate table has no header.")
        fields = reader.fieldnames

        id_col = detect_column(fields, explicit_id, ALIASES["id"])
        comp_col = detect_column(fields, explicit_comp, ALIASES["completeness"])
        cont_col = detect_column(fields, explicit_cont, ALIASES["contamination"])
        n50_col = detect_column(fields, explicit_n50, ALIASES["n50"], required=False)
        size_col = detect_column(fields, explicit_size, ALIASES["size"], required=False)
        fasta_col = detect_column(fields, explicit_fasta, ALIASES["fasta"],
                                  required=False)
        raw_iter = list(enumerate(reader, start=2))

for line_no, row in raw_iter:
    candidate_id = str(row.get(id_col, "")).strip()
    if not candidate_id:
        errors.append(f"line {line_no}: missing candidate ID")
        continue

    try:
        group = original_bin_group(candidate_id)
    except ValueError as exc:
        errors.append(f"line {line_no}: {exc}")
        continue

    completeness = to_float(row.get(comp_col))
    contamination = to_float(row.get(cont_col))
    n50 = to_float(row.get(n50_col)) if n50_col else 0.0
    size = to_float(row.get(size_col)) if size_col else 0.0

    eligible = (
        completeness >= min_comp and contamination <= max_cont
    )
    high_quality = (
        completeness >= hq_comp and contamination <= hq_cont
    )
    quality_score = completeness - 5.0 * contamination

    try:
        fasta = resolve_fasta(row, candidate_id, fasta_col)
    except (ValueError, FileNotFoundError) as exc:
        errors.append(f"line {line_no}: {exc}")
        continue

    rows.append({
        "group": group,
        "candidate_id": candidate_id,
        "completeness": completeness,
        "contamination": contamination,
        "n50": n50,
        "genome_size": size,
        "eligible": eligible,
        "high_quality": high_quality,
        "quality_score": quality_score,
        "fasta": str(fasta),
        "source_line": line_no,
    })

if errors:
    error_file = table_out / "candidate_parsing_errors.txt"
    error_file.write_text("\n".join(errors) + "\n")
    print(
        f"WARNING: {len(errors)} candidate rows could not be used. "
        f"See {error_file}",
        file=sys.stderr,
    )

groups = {}
for row in rows:
    groups.setdefault(row["group"], []).append(row)

selected = []
no_eligible = []

for group, candidates in sorted(groups.items()):
    eligible = [x for x in candidates if x["eligible"]]
    if not eligible:
        best_any = max(
            candidates,
            key=lambda x: (
                x["quality_score"],
                x["completeness"],
                -x["contamination"],
                x["n50"],
                x["genome_size"],
            ),
        )
        no_eligible.append({
            "group": group,
            "candidate_count": len(candidates),
            "best_candidate_id": best_any["candidate_id"],
            "best_completeness": best_any["completeness"],
            "best_contamination": best_any["contamination"],
            "best_quality_score": best_any["quality_score"],
            "reason": "No candidate passed completeness/contamination thresholds",
        })
        continue

    # Ranking:
    #   1. HQ candidates before MQ-only candidates.
    #   2. completeness - 5*contamination.
    #   3. lower contamination.
    #   4. higher completeness.
    #   5. higher N50.
    #   6. larger assembly.
    #   7. deterministic candidate ID.
    winner = max(
        eligible,
        key=lambda x: (
            int(x["high_quality"]),
            x["quality_score"],
            -x["contamination"],
            x["completeness"],
            x["n50"],
            x["genome_size"],
            x["candidate_id"],
        ),
    )
    winner = dict(winner)
    winner["candidate_count_in_group"] = len(candidates)
    winner["eligible_candidate_count"] = len(eligible)
    selected.append(winner)

def write_tsv(path, data, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(data)

selected_columns = [
    "group", "candidate_id", "candidate_count_in_group",
    "eligible_candidate_count", "high_quality", "completeness",
    "contamination", "quality_score", "n50", "genome_size",
    "fasta", "source_line"
]
write_tsv(table_out / "reselected_best_candidates.tsv",
          selected, selected_columns)

write_tsv(
    table_out / "groups_without_medium_quality_candidate.tsv",
    no_eligible,
    [
        "group", "candidate_count", "best_candidate_id",
        "best_completeness", "best_contamination",
        "best_quality_score", "reason"
    ],
)

manifest = []
used_names = set()

for row in selected:
    src = Path(row["fasta"])
    suffix = "".join(src.suffixes)
    if suffix.lower().endswith(".gz"):
        # dRep generally expects uncompressed FASTA. Keep source linkage here;
        # decompression is handled by the shell stage.
        suffix = suffix[:-3]
    if suffix.lower() not in {".fa", ".fna", ".fasta"}:
        suffix = ".fa"

    dest_name = f"{row['group']}__{row['candidate_id']}{suffix}"
    dest_name = re.sub(r"[^A-Za-z0-9._-]+", "_", dest_name)
    if dest_name in used_names:
        dest_name = f"{row['group']}__line{row['source_line']}{suffix}"
    used_names.add(dest_name)

    dest = fasta_out / dest_name
    src_is_gz = str(src).lower().endswith(".gz")

    if src_is_gz:
        import gzip
        with gzip.open(src, "rt") as inp, dest.open("w") as out:
            shutil.copyfileobj(inp, out)
    elif link_mode == "copy":
        shutil.copy2(src, dest)
    else:
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        dest.symlink_to(src)

    manifest.append({
        "genome": dest.name,
        "bin_id": dest.name,
        "group": row["group"],
        "candidate_id": row["candidate_id"],
        "source_fasta": str(src),
        "selected_fasta": str(dest),
        "completeness": row["completeness"],
        "contamination": row["contamination"],
        "quality_score": row["quality_score"],
        "high_quality": row["high_quality"],
        "n50": row["n50"],
        "genome_size": row["genome_size"],
    })

write_tsv(
    table_out / "reselected_fasta_manifest.tsv",
    manifest,
    [
        "genome", "bin_id", "group", "candidate_id", "source_fasta",
        "selected_fasta", "completeness", "contamination", "quality_score",
        "high_quality", "n50", "genome_size"
    ],
)

print(f"Candidate rows parsed: {len(rows):,}")
print(f"Original bin groups: {len(groups):,}")
print(f"MQ-or-better groups selected: {len(selected):,}")
print(f"Groups without an eligible candidate: {len(no_eligible):,}")
print(f"HQ selected candidates: {sum(x['high_quality'] for x in selected):,}")
PY
  fi
else
  echo "PANELS_ONLY=1; skipping candidate re-selection and dRep."
fi

RESELECTED_GLOB="$OUTDIR/reselected_fastas/*"

# Confirm FASTAs exist before starting expensive dRep jobs.
shopt -s nullglob
RESELECTED_FASTAS=( $RESELECTED_GLOB )
shopt -u nullglob

if (( ${#RESELECTED_FASTAS[@]} == 0 )); then
  if [[ "$PANELS_ONLY" == "1" ]]; then
    echo "No reselected FASTAs found, but PANELS_ONLY=1; continuing to existing dRep outputs."
  else
    echo "ERROR: No reselected FASTA files were generated." >&2
    exit 1
  fi
else
  echo "Reselected FASTAs: ${#RESELECTED_FASTAS[@]}"
fi

###############################################################################
# Independent 95% and 99% dRep catalogues
#
# -pa: primary ANI threshold
# -sa: secondary ANI threshold
#
# The same quality-reselected input set is used for both catalogues.
###############################################################################

drep_done() {
  local dir="$1"
  [[ -d "$dir/dereplicated_genomes" ]] || return 1
  [[ -s "$dir/data_tables/Widb.csv" ]] || return 1
  find "$dir/dereplicated_genomes" -maxdepth 1 -type f \
    \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) -print -quit | grep -q .
}

if [[ "$PANELS_ONLY" != "1" ]]; then
  if [[ "$REUSE_EXISTING" == "1" ]] && drep_done "$OUTDIR/drep95"; then
    echo "Reusing completed dRep95: $OUTDIR/drep95"
  else
    dRep dereplicate "$OUTDIR/drep95" \
      -g "${RESELECTED_FASTAS[@]}" \
      -p "$THREADS" \
      -comp "$MIN_COMPLETENESS" \
      -con "$MAX_CONTAMINATION" \
      -pa 0.90 \
      -sa 0.95 \
      --S_algorithm fastANI \
      --cov_thresh 0.10 \
      --clusterAlg average \
      2>&1 | tee "$OUTDIR/logs/drep95.log"
  fi

  if [[ "$REUSE_EXISTING" == "1" ]] && drep_done "$OUTDIR/drep99"; then
    echo "Reusing completed dRep99: $OUTDIR/drep99"
  else
    dRep dereplicate "$OUTDIR/drep99" \
      -g "${RESELECTED_FASTAS[@]}" \
      -p "$THREADS" \
      -comp "$MIN_COMPLETENESS" \
      -con "$MAX_CONTAMINATION" \
      -pa 0.95 \
      -sa 0.99 \
      --S_algorithm fastANI \
      --cov_thresh 0.10 \
      --clusterAlg average \
      2>&1 | tee "$OUTDIR/logs/drep99.log"
  fi
else
  drep_done "$OUTDIR/drep95" || {
    echo "ERROR: PANELS_ONLY requested but completed dRep95 output is missing: $OUTDIR/drep95" >&2
    exit 1
  }
  drep_done "$OUTDIR/drep99" || {
    echo "ERROR: PANELS_ONLY requested but completed dRep99 output is missing: $OUTDIR/drep99" >&2
    exit 1
  }
fi

###############################################################################
# Build clean adaptive-sampling FASTAs
#
# We remove short contigs and exact duplicate sequences. FASTA headers are
# prefixed by representative genome name so all contig IDs are globally unique.
###############################################################################

build_panel () {
  local derep_dir="$1"
  local panel_name="$2"
  local winners_dir="$derep_dir/dereplicated_genomes"
  local raw="$OUTDIR/tmp/${panel_name}.raw.fa"
  local clean="$OUTDIR/panels/${panel_name}.fa"
  local index="$OUTDIR/panels/${panel_name}.mmi"
  local bed="$OUTDIR/panels/${panel_name}.full.bed"
  local manifest="$OUTDIR/tables/${panel_name}.genomes.tsv"

  if [[ "$REUSE_EXISTING" == "1" && "$FORCE_PANELS" != "1" \
      && -s "$clean" && -s "$index" && -s "$bed" && -s "$manifest" ]]; then
    echo "Reusing existing panel: $panel_name"
    return
  fi

  if [[ ! -d "$winners_dir" ]]; then
    echo "ERROR: Missing dRep winners directory: $winners_dir" >&2
    exit 1
  fi

  : > "$raw"
  printf "representative_genome\tfasta\n" > "$manifest"

  shopt -s nullglob
  local genomes=( "$winners_dir"/*.fa "$winners_dir"/*.fna "$winners_dir"/*.fasta )
  shopt -u nullglob

  if (( ${#genomes[@]} == 0 )); then
    echo "ERROR: No representative genomes found in $winners_dir" >&2
    exit 1
  fi

  for genome in "${genomes[@]}"; do
    local genome_id
    genome_id="$(basename "$genome")"
    genome_id="${genome_id%.fasta}"
    genome_id="${genome_id%.fna}"
    genome_id="${genome_id%.fa}"

    printf "%s\t%s\n" "$genome_id" "$genome" >> "$manifest"

    # Prefix contig headers with genome ID and remove descriptions.
    seqkit replace \
      -p '^(\S+).*$' \
      -r "${genome_id}__\${1}" \
      "$genome" >> "$raw"
  done

  seqkit seq -m "$MIN_CONTIG_LEN" "$raw" \
    | seqkit rmdup -s \
    > "$clean"

  # FASTA index and minimap2 index for offline testing.
  seqkit faidx "$clean"
  "$MINIMAP2_BIN" -t "$THREADS" -d "$index" "$clean"

  # Whole-reference BED accepted by adaptive-sampling workflows that request
  # target intervals. Coordinates are zero-based, half-open.
  seqkit fx2tab -n -l "$clean" \
    | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' \
    > "$bed"

  seqkit stats -a "$clean" \
    > "$OUTDIR/tables/${panel_name}.stats.txt"

  echo "Built $panel_name with ${#genomes[@]} representative genomes."
}

build_panel "$OUTDIR/drep95" "CanMAG_depletion_95ANI"
build_panel "$OUTDIR/drep99" "CanMAG_depletion_99ANI"

###############################################################################
# Checksums and summary
###############################################################################

(
  cd "$OUTDIR/panels"
  sha256sum \
    CanMAG_depletion_95ANI.fa \
    CanMAG_depletion_95ANI.mmi \
    CanMAG_depletion_95ANI.full.bed \
    CanMAG_depletion_99ANI.fa \
    CanMAG_depletion_99ANI.mmi \
    CanMAG_depletion_99ANI.full.bed \
    > SHA256SUMS
)

cat <<EOF

Finished.

Primary discovery-oriented panel:
  $OUTDIR/panels/CanMAG_depletion_95ANI.fa

Aggressive strain-aware comparison panel:
  $OUTDIR/panels/CanMAG_depletion_99ANI.fa

Selection audit:
  $OUTDIR/tables/reselected_best_candidates.tsv
  $OUTDIR/tables/groups_without_medium_quality_candidate.tsv
  $OUTDIR/tables/reselected_fasta_manifest.tsv

Before sequencing, map several independent canine metagenomes against both
panels and compare primary mapped reads, mapped bases, MAPQ, and the taxonomy
of reads captured only by the 99% panel.
EOF
