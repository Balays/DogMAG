#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-58}"
CANMAG="${CANMAG:-$(pwd)}"
BASE="${BASE:-$CANMAG/opera_ms_dogfirst}"
SRS_MANIFEST="${SRS_MANIFEST:-$CANMAG/coassembly_dogfirst_20260505/metadata/dmd_dogfirst_srs_samples.tsv}"
LONG_DIR="${LONG_DIR:-$CANMAG/BASALT_v2}"
OPERA_MS="${OPERA_MS:-$(command -v OPERA-MS.pl || true)}"
OPERA_MINIMAP2="${OPERA_MINIMAP2:-$(command -v mm2-fast || command -v minimap2 || true)}"
OPERA_SAMTOOLS="${OPERA_SAMTOOLS:-$(command -v samtools || true)}"
LONG_READ_MAPPER="${LONG_READ_MAPPER:-minimap2}"
REF_CLUSTERING="${REF_CLUSTERING:-NO}"
STRAIN_CLUSTERING="${STRAIN_CLUSTERING:-YES}"
DECOMPRESS_INPUTS="${DECOMPRESS_INPUTS:-1}"
NOPOLISHING="${NOPOLISHING:-YES}"

mkdir -p "$BASE"/{configs,inputs,logs,results}

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERROR] Missing command in PATH: $1" >&2
    exit 1
  fi
}

if [[ ! -f "$OPERA_MS" ]]; then
  echo "[ERROR] OPERA-MS.pl not found: $OPERA_MS" >&2
  echo "Set OPERA_MS=/path/to/OPERA-MS.pl" >&2
  exit 1
fi

need_cmd perl
need_cmd zcat
if [[ "$LONG_READ_MAPPER" == "minimap2" ]]; then
  if [[ ! -x "$OPERA_MINIMAP2" ]]; then
    echo "[ERROR] OPERA_MINIMAP2 is not executable: $OPERA_MINIMAP2" >&2
    exit 1
  fi
fi

ensure_opera_minimap2_wrapper() {
  local opera_dir
  local tools_dir
  local bundled
  opera_dir="$(cd "$(dirname "$OPERA_MS")" && pwd)"
  tools_dir="$opera_dir/tools_opera_ms"
  bundled="$tools_dir/minimap2"

  if [[ ! -d "$tools_dir" ]]; then
    echo "[ERROR] OPERA tools directory not found: $tools_dir" >&2
    exit 1
  fi

  if [[ -e "$bundled" ]] && ! grep -qF "$OPERA_MINIMAP2" "$bundled" 2>/dev/null; then
    mv "$bundled" "$bundled.bundled_$(date +%Y%m%d_%H%M%S)"
  fi

  cat > "$bundled" <<EOF
#!/usr/bin/env bash
exec "$OPERA_MINIMAP2" "\$@"
EOF
  chmod +x "$bundled"
  echo "[INFO] OPERA minimap2 wrapper: $bundled -> $OPERA_MINIMAP2" >&2
}

ensure_opera_perl_wrapper() {
  local opera_dir
  local tools_dir
  local bundled
  local real_perl
  opera_dir="$(cd "$(dirname "$OPERA_MS")" && pwd)"
  tools_dir="$opera_dir/tools_opera_ms"
  bundled="$tools_dir/perl"
  real_perl="$(command -v perl)"

  if [[ -z "$real_perl" || ! -x "$real_perl" ]]; then
    echo "[ERROR] perl is not executable in PATH" >&2
    exit 1
  fi

  cat > "$bundled" <<EOF
#!/usr/bin/env bash
exec "$real_perl" "\$@"
EOF
  chmod +x "$bundled"
  echo "[INFO] OPERA perl wrapper: $bundled -> $real_perl" >&2
}

ensure_opera_perl_wrapper

if [[ "$LONG_READ_MAPPER" == "minimap2" ]]; then
  ensure_opera_minimap2_wrapper
fi

ensure_opera_samtools_wrapper() {
  local opera_dir
  local tools_dir
  local bundled
  opera_dir="$(cd "$(dirname "$OPERA_MS")" && pwd)"
  tools_dir="$opera_dir/tools_opera_ms"
  bundled="$tools_dir/samtools"

  if [[ -z "$OPERA_SAMTOOLS" || ! -x "$OPERA_SAMTOOLS" ]]; then
    echo "[ERROR] OPERA_SAMTOOLS is not executable. Install samtools or set OPERA_SAMTOOLS=/path/to/samtools" >&2
    exit 1
  fi

  if [[ -e "$bundled" ]] && ! grep -qF "$OPERA_SAMTOOLS" "$bundled" 2>/dev/null; then
    mv "$bundled" "$bundled.bundled_$(date +%Y%m%d_%H%M%S)"
  fi

  cat > "$bundled" <<EOF
#!/usr/bin/env bash
if [[ "\${1:-}" == "sort" && "\${2:-}" == "-" && -n "\${3:-}" ]]; then
  prefix="\$3"
  shift 3
  exec "$OPERA_SAMTOOLS" sort -o "\${prefix}.bam" -T "\${prefix}.tmp" - "\$@"
fi
exec "$OPERA_SAMTOOLS" "\$@"
EOF
  chmod +x "$bundled"
  echo "[INFO] OPERA samtools wrapper: $bundled -> $OPERA_SAMTOOLS" >&2
}

ensure_opera_samtools_wrapper

dog_long_reads() {
  case "$1" in
    DMD_B_PF1)
      printf '%s\n' \
        "$LONG_DIR/DMD_B_PF1_210504_MN.fastq.gz" \
        "$LONG_DIR/DMD_B_PF1_210504_HMW.fastq.gz"
      ;;
    DMD_B_PM2)
      printf '%s\n' \
        "$LONG_DIR/DMD_B_PM2_210505_MN.fastq.gz" \
        "$LONG_DIR/DMD_B_PM2_210505_HMW.fastq.gz"
      ;;
    DMD_0_AM1)
      printf '%s\n' \
        "$LONG_DIR/DMD_0_AM1_210513_MN.fastq.gz" \
        "$LONG_DIR/DMD_0_AM1_210513_HMW.fastq.gz"
      ;;
    DMD_0_AF1)
      printf '%s\n' \
        "$LONG_DIR/DMD_0_AF1_210508_MN.fastq.gz" \
        "$LONG_DIR/DMD_0_AF1_210508_HMW.fastq.gz"
      ;;
    *)
      echo "[ERROR] No explicit LRS mapping for dog: $1" >&2
      return 1
      ;;
  esac
}

check_files() {
  local missing=0
  for f in "$@"; do
    if [[ ! -s "$f" ]]; then
      echo "[MISSING] $f" >&2
      missing=1
    fi
  done
  return "$missing"
}

collect_srs_paths() {
  local dog="$1"
  local col="$2"
  awk -F '\t' -v dog="$dog" -v col="$col" 'NR>1 && $1 == dog {print $col}' "$SRS_MANIFEST" |
    sed "s#^/path/to/DogMAG_workdir#$CANMAG#" |
    sed 's#_trimmed_host_removed_R1\.fq\.gz#_R1_filt.fq.gz#' |
    sed 's#_trimmed_host_removed_R2\.fq\.gz#_R2_filt.fq.gz#'
}

concat_fastq_inputs() {
  local out_plain="$1"
  shift
  if [[ -s "$out_plain" ]]; then
    echo "[SKIP] exists: $out_plain" >&2
    echo "$out_plain"
    return 0
  fi
  if [[ "$DECOMPRESS_INPUTS" == "1" ]]; then
    echo "[ZCAT] $out_plain" >&2
    zcat "$@" > "$out_plain"
    echo "$out_plain"
  else
    local out_gz="${out_plain}.gz"
    echo "[CAT] $out_gz" >&2
    cat "$@" > "$out_gz"
    echo "$out_gz"
  fi
}

run_dog() {
  local dog="$1"
  local dog_id
  dog_id="$(echo "$dog" | tr '[:upper:]' '[:lower:]')"
  local input_dir="$BASE/inputs/$dog_id"
  local out_dir="$BASE/results/$dog_id"
  mkdir -p "$input_dir" "$out_dir"

  mapfile -t r1s < <(collect_srs_paths "$dog" 9)
  mapfile -t r2s < <(collect_srs_paths "$dog" 10)
  mapfile -t lrs < <(dog_long_reads "$dog")

  if [[ "${#r1s[@]}" -eq 0 || "${#r2s[@]}" -eq 0 ]]; then
    echo "[ERROR] No SRS rows found for dog $dog in $SRS_MANIFEST" >&2
    return 1
  fi
  check_files "${r1s[@]}" "${r2s[@]}" "${lrs[@]}"

  local r1_input
  local r2_input
  local lr_input
  r1_input="$(concat_fastq_inputs "$input_dir/${dog_id}_R1.fq" "${r1s[@]}")"
  r2_input="$(concat_fastq_inputs "$input_dir/${dog_id}_R2.fq" "${r2s[@]}")"
  lr_input="$(concat_fastq_inputs "$input_dir/${dog_id}_long_reads.fastq" "${lrs[@]}")"

  local config="$BASE/configs/${dog_id}.config"
  cat > "$config" <<EOF
ILLUMINA_READ_1 $r1_input
ILLUMINA_READ_2 $r2_input
LONG_READ $lr_input
OUTPUT_DIR $out_dir
NUM_PROCESSOR $THREADS
LONG_READ_MAPPER $LONG_READ_MAPPER
REF_CLUSTERING $REF_CLUSTERING
STRAIN_CLUSTERING $STRAIN_CLUSTERING
NOPOLISHING $NOPOLISHING
EOF

  echo "[RUN] OPERA-MS $dog"
  perl "$OPERA_MS" "$config" \
    > "$BASE/logs/${dog_id}.stdout.log" \
    2> "$BASE/logs/${dog_id}.stderr.log"
}

if [[ "${1:-all}" == "check" ]]; then
  echo "[INFO] OPERA-MS: $OPERA_MS"
  echo "[INFO] SRS manifest: $SRS_MANIFEST"
  echo "[INFO] LONG_DIR: $LONG_DIR"
  echo "[INFO] DECOMPRESS_INPUTS: $DECOMPRESS_INPUTS"
  for dog in DMD_B_PF1 DMD_B_PM2 DMD_0_AM1 DMD_0_AF1; do
    echo "[CHECK] $dog"
    mapfile -t r1s < <(collect_srs_paths "$dog" 9)
    mapfile -t r2s < <(collect_srs_paths "$dog" 10)
    mapfile -t lrs < <(dog_long_reads "$dog")
    check_files "${r1s[@]}" "${r2s[@]}" "${lrs[@]}" || true
    echo "  SRS pairs: ${#r1s[@]}"
    echo "  LRS files: ${#lrs[@]}"
  done
  exit 0
fi

case "${1:-all}" in
  DMD_B_PF1|DMD_B_PM2|DMD_0_AM1|DMD_0_AF1)
    run_dog "$1"
    ;;
  pilot|all)
    for dog in DMD_B_PF1 DMD_B_PM2 DMD_0_AM1 DMD_0_AF1; do
      run_dog "$dog"
    done
    ;;
  *)
    echo "Usage: bash run_opera_ms_named_lrs_dogs.sh [check|pilot|all|DMD_B_PF1|DMD_B_PM2|DMD_0_AM1|DMD_0_AF1]" >&2
    exit 1
    ;;
esac

echo "[DONE] OPERA-MS workflow finished"
