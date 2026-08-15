#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Build a 16S-only minitax database from any genome FASTA set.

This workflow:
  1. runs barrnap on each genome FASTA,
  2. extracts retained 16S rRNA intervals from the MAG FASTAs,
  3. joins each 16S sequence to taxonomy when a taxonomy table is supplied,
  4. writes minitax-compatible FASTA and MAG.db.tsv/db_data.tsv files.

Use it for raw BASALT MAG folders, dRep95, dRep99, or another one-FASTA-per-genome set.

Usage:
  bash scripts/run_genome_16s_minitax_db_workflow.sh --genome-dir PATH --out-dir PATH --db-name NAME [options]

Options:
  --project-dir PATH       CanMAG project directory
                          [default: /path/to/DogMAG_workdir if present, otherwise /mnt/d/data/CanMAG]
  --genome-dir PATH        Directory with one genome FASTA per MAG/genome [required unless GENOME_DIR is set]
  --taxonomy PATH          Optional taxonomy TSV, preferably GTDB-Tk summary style
  --out-dir PATH           Output root [default: PROJECT_DIR/minitax_16s_db]
  --db-name NAME           Database basename [default: GenomeSet_16S]
  --threads N              Threads for barrnap/index commands [default: 16]
  --min-length N           Minimum 16S length retained [default: 1000]
  --kingdom NAME           barrnap kingdom: bac, arc, mito, euk [default: bac]
  --force-barrnap          Re-run barrnap even when a GFF already exists
  --build-index            Also build DB_NAME.idx with mm2-fast/minimap2
  --mapper-bin PATH        Index builder executable [default: mm2-fast if found, else minimap2]
  --indir PATH             Placeholder FASTQ directory written into minitax config
  -h, --help               Show this help

Environment overrides:
  PROJECT_DIR GENOME_DIR TAXONOMY OUT_DIR DB_NAME THREADS MIN_16S_LENGTH
  BARRNAP_KINGDOM BUILD_INDEX MAPPER_BIN INDIR

Example:
  mamba activate minitax-mm2fast
  bash scripts/run_genome_16s_minitax_db_workflow.sh \
    --genome-dir /path/to/DogMAG_workdir/final_drep_all_final_20260505 \
    --taxonomy /path/to/DogMAG_workdir/work/mag_taxonomy/mag_taxonomy_gtdbtk_summary.tsv \
    --out-dir /path/to/DogMAG_workdir/minitax_final_drep_16s \
    --db-name FinalDrep_16S \
    --build-index
USAGE
}

if [[ -d /path/to/DogMAG_workdir ]]; then
  DEFAULT_PROJECT_DIR=/path/to/DogMAG_workdir
else
  DEFAULT_PROJECT_DIR=/mnt/d/data/CanMAG
fi

PROJECT_DIR=${PROJECT_DIR:-$DEFAULT_PROJECT_DIR}
GENOME_DIR=${GENOME_DIR:-}
TAXONOMY=${TAXONOMY:-}
OUT_DIR=${OUT_DIR:-"$PROJECT_DIR/minitax_16s_db"}
DB_NAME=${DB_NAME:-GenomeSet_16S}
THREADS=${THREADS:-16}
MIN_16S_LENGTH=${MIN_16S_LENGTH:-1000}
BARRNAP_KINGDOM=${BARRNAP_KINGDOM:-bac}
BUILD_INDEX=${BUILD_INDEX:-0}
MAPPER_BIN=${MAPPER_BIN:-}
INDIR=${INDIR:-/path/to/fastq}
FORCE_BARRNAP=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-dir)
      PROJECT_DIR=$2
      shift 2
      ;;
    --genome-dir)
      GENOME_DIR=$2
      shift 2
      ;;
    --taxonomy)
      TAXONOMY=$2
      shift 2
      ;;
    --out-dir)
      OUT_DIR=$2
      shift 2
      ;;
    --db-name)
      DB_NAME=$2
      shift 2
      ;;
    --threads)
      THREADS=$2
      shift 2
      ;;
    --min-length)
      MIN_16S_LENGTH=$2
      shift 2
      ;;
    --kingdom)
      BARRNAP_KINGDOM=$2
      shift 2
      ;;
    --force-barrnap)
      FORCE_BARRNAP=1
      shift
      ;;
    --build-index)
      BUILD_INDEX=1
      shift
      ;;
    --mapper-bin)
      MAPPER_BIN=$2
      shift 2
      ;;
    --indir)
      INDIR=$2
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BUILD_SCRIPT="$SCRIPT_DIR/build_genome_16s_minitax_db.py"
BARRNAP_DIR="$OUT_DIR/barrnap_gff"
DB_DIR="$OUT_DIR/db"

require_file() {
  local path=$1
  local label=$2
  if [[ ! -f "$path" ]]; then
    echo "ERROR: missing $label: $path" >&2
    exit 1
  fi
}

require_dir() {
  local path=$1
  local label=$2
  if [[ ! -d "$path" ]]; then
    echo "ERROR: missing $label: $path" >&2
    exit 1
  fi
}

if [[ -z "$GENOME_DIR" ]]; then
  echo "ERROR: --genome-dir is required unless GENOME_DIR is set." >&2
  usage >&2
  exit 2
fi

require_dir "$GENOME_DIR" "genome directory"
if [[ -n "$TAXONOMY" ]]; then
  require_file "$TAXONOMY" "taxonomy table"
fi
require_file "$BUILD_SCRIPT" "Python build script"

if ! command -v barrnap >/dev/null 2>&1; then
  echo "ERROR: barrnap is not on PATH. Install it in the active environment, e.g.:" >&2
  echo "  mamba install -c conda-forge -c bioconda barrnap" >&2
  exit 1
fi

mkdir -p "$BARRNAP_DIR" "$DB_DIR" "$OUT_DIR/logs"

echo "Project dir: $PROJECT_DIR"
echo "Genome dir:  $GENOME_DIR"
echo "Taxonomy:    ${TAXONOMY:-not supplied; ranks will be NA}"
echo "Output dir:  $OUT_DIR"
echo "Database:    $DB_NAME"
echo "barrnap:     kingdom=$BARRNAP_KINGDOM threads=$THREADS"

find "$GENOME_DIR" -maxdepth 1 -type f \
  \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' -o -name '*.fa.gz' -o -name '*.fna.gz' -o -name '*.fasta.gz' \) \
  | sort > "$OUT_DIR/genome_fastas.list"

if [[ ! -s "$OUT_DIR/genome_fastas.list" ]]; then
  echo "ERROR: no genome FASTAs found in $GENOME_DIR" >&2
  exit 1
fi

while IFS= read -r fasta; do
  base=$(basename "$fasta")
  stem=$base
  stem=${stem%.gz}
  stem=${stem%.fa}
  stem=${stem%.fna}
  stem=${stem%.fasta}
  if [[ "$stem" == *"__"* ]]; then
    left=${stem%%__*}
    right=${stem#*__}
    right=${right%.fa}
    right=${right%.fna}
    right=${right%.fasta}
    genome_id="${left}_${right}"
  else
    genome_id=$stem
  fi
  genome_id=$(printf '%s' "$genome_id" | sed -E 's/[^A-Za-z0-9_.:-]+/_/g; s/_+/_/g; s/^_//; s/_$//')
  gff="$BARRNAP_DIR/${genome_id}.gff"
  log="$OUT_DIR/logs/${genome_id}.barrnap.log"
  if [[ -s "$gff" && "$FORCE_BARRNAP" -eq 0 ]]; then
    continue
  fi
  echo "Running barrnap: $base"
  barrnap --kingdom "$BARRNAP_KINGDOM" --threads "$THREADS" "$fasta" > "$gff" 2> "$log"
done < "$OUT_DIR/genome_fastas.list"

BUILD_ARGS=(
  --genome-dir "$GENOME_DIR"
  --barrnap-dir "$BARRNAP_DIR"
  --out-dir "$OUT_DIR"
  --db-name "$DB_NAME"
  --min-length "$MIN_16S_LENGTH"
  --indir "$INDIR"
)

if [[ -n "$TAXONOMY" ]]; then
  BUILD_ARGS+=(--taxonomy "$TAXONOMY")
fi

python3 "$BUILD_SCRIPT" "${BUILD_ARGS[@]}"

if [[ "$BUILD_INDEX" -eq 1 ]]; then
  if [[ -z "$MAPPER_BIN" ]]; then
    if command -v mm2-fast >/dev/null 2>&1; then
      MAPPER_BIN=mm2-fast
    elif command -v minimap2 >/dev/null 2>&1; then
      MAPPER_BIN=minimap2
    else
      echo "ERROR: neither mm2-fast nor minimap2 is on PATH; cannot build index." >&2
      exit 1
    fi
  fi
  echo "Building minimap2 index with $MAPPER_BIN"
  (
    cd "$DB_DIR"
    "$MAPPER_BIN" -I 16G -d "${DB_NAME}.idx" "${DB_NAME}.fa"
  )
fi

echo "Done."
echo "Minitax DB:     $DB_DIR"
echo "Minitax config: $OUT_DIR/configs/minitax_config_${DB_NAME}.tsv"
