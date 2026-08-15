#!/usr/bin/env bash
set -euo pipefail

opera_root="${OPERA_ROOT:-/path/to/DogMAG_workdir/opera_ms_dogfirst_20260505}"
metaspades_root="${METASPADES_ROOT:-/path/to/DogMAG_workdir/metaspades_dogfirst_20260511}"
output_dir="${OUTDIR:-/path/to/DogMAG_workdir/final_hybrid_dogfirst_assemblies_20260515}"
dry_run="false"
force_overwrite="false"
gzip_output="false"

dogs=(
  boszi
  brios
  cefre
  csuzli
  kek
  loki
  pink
  piros
  sugo
  toti
  zold
)

usage() {
  cat <<'EOF'
Usage:
  collect_hybrid_dogfirst_assemblies.sh [options]

Options:
  --opera-root PATH       OPERA-MS project root.
                          Default: /path/to/DogMAG_workdir/opera_ms_dogfirst_20260505
  --metaspades-root PATH  metaSPAdes project root.
                          Default: /path/to/DogMAG_workdir/metaspades_dogfirst_20260511
  -o, --output-dir PATH   Distinct flattened output directory.
                          Default: /path/to/DogMAG_workdir/final_hybrid_dogfirst_assemblies_20260515
  -gzip                   Write .fa.gz files instead of .fa files.
  -force                  Overwrite existing outputs.
  --dry-run               Print planned copies without writing files.
  -h, --help              Show this help.

Environment equivalents:
  OPERA_ROOT, METASPADES_ROOT, OUTDIR

Behavior:
  Collects one OPERA-MS assembly and one metaSPAdes hybrid assembly per dog when present.
  Output names are distinct and BASALT-friendly:
    lrs_srs_opera_ms_hybrid__dmd_dog_<dog>.fa
    lrs_srs_metaspades_hybrid__dmd_dog_<dog>.fa

Preferred sources:
  OPERA-MS:
    OPERA_ROOT/work/opera_ms_dog_manifest_20260509/results/<dog>/contigs.fasta
    fallback: .../intermediate_files/opera_long_read/scaffoldSeq.fasta.filled
    fallback: .../intermediate_files/opera_long_read/scaffoldSeq.fasta

  metaSPAdes:
    METASPADES_ROOT/results/<dog>_metaspades_hybrid/scaffolds.fasta
    fallback: .../contigs.fasta

Examples:
  bash scripts/collect_hybrid_dogfirst_assemblies.sh --dry-run

  bash scripts/collect_hybrid_dogfirst_assemblies.sh \
    --opera-root /path/to/DogMAG_workdir/opera_ms_dogfirst_20260505 \
    --metaspades-root /path/to/DogMAG_workdir/metaspades_dogfirst_20260511 \
    --output-dir /path/to/DogMAG_workdir/final_hybrid_dogfirst_assemblies_20260515
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --opera-root)
      opera_root="$2"
      shift 2
      ;;
    --metaspades-root)
      metaspades_root="$2"
      shift 2
      ;;
    -o|--output-dir)
      output_dir="$2"
      shift 2
      ;;
    -gzip)
      gzip_output="true"
      shift
      ;;
    -force)
      force_overwrite="true"
      shift
      ;;
    --dry-run)
      dry_run="true"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

gzip_cmd() {
  if command -v pigz >/dev/null 2>&1; then
    pigz -c
  else
    gzip -c
  fi
}

first_existing_file() {
  local candidate
  for candidate in "$@"; do
    if [[ -s "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

copy_one() {
  local source_file="$1"
  local dest_file="$2"

  if [[ "$gzip_output" == "true" ]]; then
    dest_file="${dest_file}.gz"
  fi

  if [[ "$dry_run" == "true" ]]; then
    if [[ -e "$dest_file" ]]; then
      printf '[dry-run] existing destination: %s\n' "$dest_file"
    fi
    printf '[dry-run] %s -> %s\n' "$source_file" "$dest_file"
    return 0
  fi

  mkdir -p "$(dirname "$dest_file")"

  if [[ -e "$dest_file" && "$force_overwrite" != "true" ]]; then
    if [[ "$gzip_output" == "false" ]] && cmp -s "$source_file" "$dest_file"; then
      printf 'Skipping existing identical file: %s\n' "$dest_file"
      return 0
    fi
    printf 'Skipping existing destination (use -force to overwrite): %s\n' "$dest_file"
    return 0
  fi

  if [[ "$gzip_output" == "true" ]]; then
    printf 'Gzipping %s -> %s\n' "$source_file" "$dest_file"
    gzip_cmd < "$source_file" > "$dest_file"
  else
    printf 'Copying %s -> %s\n' "$source_file" "$dest_file"
    cp -a "$source_file" "$dest_file"
  fi
}

if [[ "$dry_run" != "true" ]]; then
  mkdir -p "$output_dir"
fi

summary_file="${output_dir%/}/assembly_sources.tsv"
if [[ "$dry_run" != "true" ]]; then
  printf 'assembler\tdog\toutput_name\tsource_path\n' > "$summary_file"
fi

copied=0
missing=0

for dog in "${dogs[@]}"; do
  opera_source="$(
    first_existing_file \
      "${opera_root%/}/work/opera_ms_dog_manifest_20260509/results/${dog}/contigs.fasta" \
      "${opera_root%/}/work/opera_ms_dog_manifest_20260509/results/${dog}/intermediate_files/opera_long_read/scaffoldSeq.fasta.filled" \
      "${opera_root%/}/work/opera_ms_dog_manifest_20260509/results/${dog}/intermediate_files/opera_long_read/scaffoldSeq.fasta" \
      2>/dev/null || true
  )"

  if [[ -n "$opera_source" ]]; then
    opera_name="lrs_srs_opera_ms_hybrid__dmd_dog_${dog}.fa"
    copy_one "$opera_source" "${output_dir%/}/${opera_name}"
    if [[ "$dry_run" != "true" ]]; then
      printf 'opera_ms\t%s\t%s\t%s\n' "$dog" "$opera_name" "$opera_source" >> "$summary_file"
    fi
    ((copied+=1))
  else
    printf 'Missing OPERA-MS assembly for dog: %s\n' "$dog" >&2
    ((missing+=1))
  fi

  metaspades_source="$(
    first_existing_file \
      "${metaspades_root%/}/results/${dog}_metaspades_hybrid/scaffolds.fasta" \
      "${metaspades_root%/}/results/${dog}_metaspades_hybrid/contigs.fasta" \
      2>/dev/null || true
  )"

  if [[ -n "$metaspades_source" ]]; then
    metaspades_name="lrs_srs_metaspades_hybrid__dmd_dog_${dog}.fa"
    copy_one "$metaspades_source" "${output_dir%/}/${metaspades_name}"
    if [[ "$dry_run" != "true" ]]; then
      printf 'metaspades\t%s\t%s\t%s\n' "$dog" "$metaspades_name" "$metaspades_source" >> "$summary_file"
    fi
    ((copied+=1))
  else
    printf 'Missing metaSPAdes assembly for dog: %s\n' "$dog" >&2
    ((missing+=1))
  fi
done

echo "Collected/planned $copied assemblies"
echo "Missing $missing expected assemblies"
echo "Output: $output_dir"
if [[ "$dry_run" != "true" ]]; then
  echo "Source manifest: $summary_file"
fi
