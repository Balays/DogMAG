#!/usr/bin/env bash
set -euo pipefail

CANMAG="${CANMAG:-$HOME/CanMAG}"
HYBRID_DIR="${HYBRID_DIR:-$CANMAG/final_hybrid_dogfirst_assemblies_20260515}"
FLYE_DIR="${FLYE_DIR:-$CANMAG/final_flye_lr_only_dogfirst_assemblies_20260518}"
MERGED_DIR="${MERGED_DIR:-$CANMAG/final_dogfirst_assemblies_for_basalt_20260518}"
FILTERED_DIR="${FILTERED_DIR:-$CANMAG/final_dogfirst_assemblies_for_basalt_20260518_min1500}"
MIN_LENGTH="${MIN_LENGTH:-1500}"
FORCE="${FORCE:-0}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  prepare_dogfirst_final_assemblies_for_basalt.sh [check|prepare]

Env:
  CANMAG=/path/to/DogMAG_workdir
  HYBRID_DIR=/path/to/DogMAG_workdir/final_hybrid_dogfirst_assemblies_20260515
  FLYE_DIR=/path/to/DogMAG_workdir/final_flye_lr_only_dogfirst_assemblies_20260518
  MERGED_DIR=/path/to/DogMAG_workdir/final_dogfirst_assemblies_for_basalt_20260518
  FILTERED_DIR=/path/to/DogMAG_workdir/final_dogfirst_assemblies_for_basalt_20260518_min1500
  MIN_LENGTH=1500
  FORCE=1
  EXCLUDE_DOG_IDS=DOG_FANCSI

Behavior:
  1. Collects manifest-approved Flye LR-only assemblies into FLYE_DIR.
  2. Symlinks HYBRID_DIR + FLYE_DIR FASTAs into MERGED_DIR.
  3. Filters merged assemblies to >= MIN_LENGTH into FILTERED_DIR.
EOF
}

mode="${1:-check}"
case "$mode" in
  check|prepare)
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

count_fastas() {
  local dir="$1"
  if [[ ! -d "$dir" ]]; then
    echo 0
    return
  fi
  find "$dir" -maxdepth 1 \( -type f -o -type l \) \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | wc -l
}

echo "[INFO] CANMAG=$CANMAG"
echo "[INFO] HYBRID_DIR=$HYBRID_DIR"
echo "[INFO] FLYE_DIR=$FLYE_DIR"
echo "[INFO] MERGED_DIR=$MERGED_DIR"
echo "[INFO] FILTERED_DIR=$FILTERED_DIR"
echo "[INFO] MIN_LENGTH=$MIN_LENGTH"

[[ -d "$HYBRID_DIR" ]] || { echo "[ERROR] Missing HYBRID_DIR: $HYBRID_DIR" >&2; exit 1; }

echo "[INFO] Hybrid FASTAs: $(count_fastas "$HYBRID_DIR")"

if [[ "$mode" == "check" ]]; then
  DRY_RUN=1 bash "$script_dir/collect_flye_lr_only_dogfirst_assemblies.sh" check
  exit 0
fi

bash "$script_dir/collect_flye_lr_only_dogfirst_assemblies.sh" collect
echo "[INFO] Flye FASTAs: $(count_fastas "$FLYE_DIR")"

if [[ -d "$MERGED_DIR" && "$FORCE" == "1" ]]; then
  find "$MERGED_DIR" -maxdepth 1 -type l -delete
fi
mkdir -p "$MERGED_DIR"

summary="$MERGED_DIR/assembly_sources.tsv"
printf 'source_group\tfile_name\tsource_path\n' > "$summary"

link_fastas() {
  local source_group="$1" dir="$2" src base dest
  shopt -s nullglob
  for src in "$dir"/*.fa "$dir"/*.fasta "$dir"/*.fna; do
    base="$(basename "$src")"
    dest="$MERGED_DIR/$base"
    if [[ -e "$dest" && ! -L "$dest" ]]; then
      echo "[ERROR] Refusing to overwrite non-symlink in MERGED_DIR: $dest" >&2
      exit 1
    fi
    if [[ -e "$dest" && "$FORCE" != "1" ]]; then
      echo "[SKIP] existing merged symlink: $dest"
    else
      ln -sfn "$src" "$dest"
    fi
    printf '%s\t%s\t%s\n' "$source_group" "$base" "$src" >> "$summary"
  done
}

link_fastas "hybrid" "$HYBRID_DIR"
link_fastas "flye_lr_only" "$FLYE_DIR"

echo "[INFO] Merged FASTAs: $(count_fastas "$MERGED_DIR")"

if [[ "$MIN_LENGTH" =~ ^[0-9]+$ && "$MIN_LENGTH" -gt 0 ]]; then
  if [[ -d "$FILTERED_DIR" && "$(count_fastas "$FILTERED_DIR")" -gt 0 && "$FORCE" != "1" ]]; then
    echo "[ERROR] FILTERED_DIR already has FASTAs; use FORCE=1 to rebuild: $FILTERED_DIR" >&2
    exit 1
  fi
  if [[ -d "$FILTERED_DIR" && "$FORCE" == "1" ]]; then
    find "$FILTERED_DIR" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) -delete
  fi
  bash "$script_dir/filter_assemblies_by_length.sh" \
    --input-dir "$MERGED_DIR" \
    --output-dir "$FILTERED_DIR" \
    --min-length "$MIN_LENGTH"
  cp -f "$summary" "$FILTERED_DIR/assembly_sources.unfiltered.tsv"
  echo "[DONE] BASALT assembly input: $FILTERED_DIR"
else
  echo "[DONE] BASALT assembly input: $MERGED_DIR"
fi
