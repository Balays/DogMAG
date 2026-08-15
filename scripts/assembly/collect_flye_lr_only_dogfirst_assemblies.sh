#!/usr/bin/env bash
set -euo pipefail

CANMAG="${CANMAG:-$HOME/CanMAG}"
FLYE_BASE="${FLYE_BASE:-$HOME/flye_work/flye_lr_only_dogfirst_20260513}"
MANIFEST_TSV="${MANIFEST_TSV:-$CANMAG/metadata/flye_lr_only_dog_manifest_20260513.tsv}"
OUTDIR="${OUTDIR:-$CANMAG/final_flye_lr_only_dogfirst_assemblies_20260518}"
DRY_RUN="${DRY_RUN:-0}"
FORCE="${FORCE:-0}"
EXCLUDE_DOG_IDS="${EXCLUDE_DOG_IDS:-}"

usage() {
  cat <<'EOF'
Usage:
  collect_flye_lr_only_dogfirst_assemblies.sh [check|collect]

Env:
  CANMAG=/path/to/DogMAG_workdir
  FLYE_BASE=/path/to/flye_work/flye_lr_only_dogfirst_20260513
  MANIFEST_TSV=/path/to/DogMAG_workdir/metadata/flye_lr_only_dog_manifest_20260513.tsv
  OUTDIR=/path/to/DogMAG_workdir/final_flye_lr_only_dogfirst_assemblies_20260518
  DRY_RUN=1
  FORCE=1
  EXCLUDE_DOG_IDS=DOG_FANCSI

Behavior:
  Uses the LR-only dog manifest, not a blind directory scan, so backup/failed
  Flye directories are ignored. Copies each final assembly.fasta to:
    lr_flye_meta__dogfirst_<dog_slug>.fa
EOF
}

mode="${1:-check}"

is_excluded_dog_id() {
  local dog_id="$1" token
  [[ -n "$EXCLUDE_DOG_IDS" ]] || return 1
  for token in ${EXCLUDE_DOG_IDS//,/ }; do
    [[ "$token" == "$dog_id" ]] && return 0
  done
  return 1
}

sanitize_name() {
  echo "$1" | tr '[:upper:]' '[:lower:]' | tr ' ' '_' | tr '.-' '__' | tr -cd '[:alnum:]_'
}

read_manifest_rows() {
  python3 - "$MANIFEST_TSV" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        dog_id = (row.get("dog_id") or "").strip()
        dog = (row.get("canonical_dog_name") or "").strip()
        if dog_id and dog:
            print(f"{dog_id}\t{dog}")
PY
}

copy_or_plan() {
  local src="$1" dest="$2"
  if [[ "$DRY_RUN" == "1" || "$mode" == "check" ]]; then
    echo "[PLAN] $src -> $dest"
    return 0
  fi
  mkdir -p "$(dirname "$dest")"
  if [[ -e "$dest" && "$FORCE" != "1" ]]; then
    if cmp -s "$src" "$dest" 2>/dev/null; then
      echo "[SKIP] identical existing: $dest"
      return 0
    fi
    echo "[SKIP] existing destination, use FORCE=1 to overwrite: $dest" >&2
    return 1
  fi
  cp -a "$src" "$dest"
}

case "$mode" in
  check|collect)
    ;;
  -h|--help|help)
    usage
    exit 0
    ;;
  *)
    usage >&2
    exit 1
    ;;
esac

[[ -s "$MANIFEST_TSV" ]] || { echo "[ERROR] Missing manifest: $MANIFEST_TSV" >&2; exit 1; }
[[ -d "$FLYE_BASE/results" ]] || { echo "[ERROR] Missing Flye results dir: $FLYE_BASE/results" >&2; exit 1; }
echo "[INFO] EXCLUDE_DOG_IDS=${EXCLUDE_DOG_IDS:-<none>}"

if [[ "$mode" == "collect" && "$DRY_RUN" != "1" ]]; then
  mkdir -p "$OUTDIR"
  printf 'dog_id\tdog\toutput_name\tsource_path\n' > "$OUTDIR/assembly_sources.tsv"
  printf 'dog_id\tdog\treason\n' > "$OUTDIR/assembly_exclusions.tsv"
fi

found=0
missing=0
excluded=0

while IFS=$'\t' read -r dog_id dog; do
  [[ -n "$dog_id" ]] || continue
  slug="$(sanitize_name "$dog")"
  src="$FLYE_BASE/results/lrs_flye_meta__${slug}/assembly.fasta"
  out_name="lr_flye_meta__dogfirst_${slug}.fa"
  dest="$OUTDIR/$out_name"

  if is_excluded_dog_id "$dog_id"; then
    echo "[EXCLUDE] $dog_id $dog"
    if [[ "$mode" == "collect" && "$DRY_RUN" != "1" ]]; then
      printf '%s\t%s\t%s\n' "$dog_id" "$dog" "explicitly excluded via EXCLUDE_DOG_IDS" >> "$OUTDIR/assembly_exclusions.tsv"
    fi
    excluded=$((excluded + 1))
    continue
  fi

  if [[ ! -s "$src" ]]; then
    echo "[MISSING] $dog_id $dog: $src" >&2
    missing=$((missing + 1))
    continue
  fi

  copy_or_plan "$src" "$dest"
  if [[ "$mode" == "collect" && "$DRY_RUN" != "1" ]]; then
    printf '%s\t%s\t%s\t%s\n' "$dog_id" "$dog" "$out_name" "$src" >> "$OUTDIR/assembly_sources.tsv"
  fi
  found=$((found + 1))
done < <(read_manifest_rows)

echo "[INFO] Flye assemblies found/planned: $found"
echo "[INFO] Missing: $missing"
echo "[INFO] Excluded: $excluded"
echo "[INFO] OUTDIR=$OUTDIR"

if [[ "$missing" -gt 0 ]]; then
  exit 1
fi
