#!/usr/bin/env bash
set -euo pipefail

# Dog-first Flye metagenome assemblies for dogs that currently have long reads
# only. The default manifest is curated from the final coassembly metadata plus
# the corrected Serteperti six-barcode mapping.
#
# Modes:
#   list     Print dogs in the queue.
#   check    Resolve inputs and write resolved manifest; do not stage/run.
#   prepare  Resolve inputs and stage per-dog symlinks; do not run Flye.
#   run      Stage symlinks and run one Flye --meta assembly per dog.

THREADS="${THREADS:-20}"
BASE="${BASE:-$HOME/flye_work/flye_lr_only_dogfirst_$(date +%Y%m%d)}"
CANMAG="${CANMAG:-$HOME/CanMAG}"
MANIFEST_TSV="${MANIFEST_TSV:-$CANMAG/metadata/flye_lr_only_dog_manifest_20260513.tsv}"
LONG_DIR="${LONG_DIR:-$CANMAG/fastq/All_Kennel_ONT_WGS}"
EXTRA_LONG_DIR="${EXTRA_LONG_DIR:-$CANMAG/BASALT_v2}"
COASSEMBLY_SCRIPT="${COASSEMBLY_SCRIPT:-$CANMAG/scripts/run_flye_coassembly.sh}"

READ_TYPE="${READ_TYPE:---nano-hq}"
META="${META:-1}"
FORCE="${FORCE:-0}"
DOG_FILTER="${DOG_FILTER:-}"
CONTINUE_ON_DOG_ERROR="${CONTINUE_ON_DOG_ERROR:-0}"
VALIDATE_GZIP="${VALIDATE_GZIP:-0}"
RESUME_PARTIAL="${RESUME_PARTIAL:-1}"
RESTART_FAILED_RESUME="${RESTART_FAILED_RESUME:-1}"
RETRY_NO_POLISH_ON_FAILURE="${RETRY_NO_POLISH_ON_FAILURE:-1}"
FLYE_EXTRA_ARGS="${FLYE_EXTRA_ARGS:-}"

mkdir -p "$BASE"/{inputs,logs,results,status,metadata,tmp}
export TMPDIR="${TMPDIR:-$BASE/tmp}"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERROR] Missing command in PATH: $1" >&2
    exit 1
  fi
}

sanitize_name() {
  echo "$1" | tr '[:upper:]' '[:lower:]' | tr ' ' '_' | tr '.-' '__' | tr -cd '[:alnum:]_'
}

split_fastq_cell() {
  local cell="$1"
  echo "$cell" \
    | tr ';' '\n' \
    | tr ',' '\n' \
    | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' \
    | awk 'NF && $0 != "nan" && $0 != "NA"'
}

barcode_key() {
  local path="$1" base prefix num
  base="$(basename "$path")"
  if [[ "$base" =~ ^([AB])([45])?_barcode0*([0-9]+)\.f(ast)?q\.gz$ ]]; then
    prefix="${BASH_REMATCH[1]}"
    num="${BASH_REMATCH[3]}"
    printf '%s:%02d' "$prefix" "$num"
    return 0
  fi
  if [[ "$base" =~ barcode0*([0-9]+)\.f(ast)?q\.gz$ ]]; then
    num="${BASH_REMATCH[1]}"
    printf 'serteperti_six:%02d' "$num"
    return 0
  fi
  printf '%s' "$base"
}

candidate_paths() {
  local raw="$1" p base prefix kit num other_kit
  p="${raw# }"
  p="${p% }"
  p="${p/#\/home\/mdbio\/CanMAG/$CANMAG}"
  p="${p/#\/home\/kakuk\/CanMAG/$CANMAG}"
  [[ "$p" != /* ]] && p="$CANMAG/$p"
  base="$(basename "$p")"

  printf '%s\n' "$p"
  printf '%s\n' "$CANMAG/fastq/long_reads/$base"
  printf '%s\n' "$CANMAG/fastq/All_Kennel_ONT_WGS/$base"
  printf '%s\n' "$LONG_DIR/$base"
  printf '%s\n' "$EXTRA_LONG_DIR/$base"

  if [[ "$base" =~ ^([AB])([45])_barcode0*([0-9]+)\.f(ast)?q\.gz$ ]]; then
    prefix="${BASH_REMATCH[1]}"
    kit="${BASH_REMATCH[2]}"
    num="${BASH_REMATCH[3]}"
    [[ "$kit" == "4" ]] && other_kit="5" || other_kit="4"
    printf '%s\n' "$CANMAG/fastq/All_Kennel_ONT_WGS/${prefix}_barcode$(printf '%02d' "$num").fastq.gz"
    printf '%s\n' "$CANMAG/fastq/All_Kennel_ONT_WGS/${prefix}${other_kit}_barcode$(printf '%02d' "$num").fastq.gz"
    printf '%s\n' "$CANMAG/fastq/long_reads/${prefix}_barcode$(printf '%02d' "$num").fastq.gz"
    printf '%s\n' "$CANMAG/fastq/long_reads/${prefix}${other_kit}_barcode$(printf '%02d' "$num").fastq.gz"
  fi

  if [[ "$base" =~ ^([AB])_barcode0*([0-9]+)\.f(ast)?q\.gz$ ]]; then
    prefix="${BASH_REMATCH[1]}"
    num="${BASH_REMATCH[2]}"
    printf '%s\n' "$CANMAG/fastq/All_Kennel_ONT_WGS/${prefix}4_barcode$(printf '%02d' "$num").fastq.gz"
    printf '%s\n' "$CANMAG/fastq/All_Kennel_ONT_WGS/${prefix}5_barcode$(printf '%02d' "$num").fastq.gz"
    printf '%s\n' "$CANMAG/fastq/long_reads/${prefix}4_barcode$(printf '%02d' "$num").fastq.gz"
    printf '%s\n' "$CANMAG/fastq/long_reads/${prefix}5_barcode$(printf '%02d' "$num").fastq.gz"
  fi

  if [[ "$base" =~ barcode0*([0-9]+)\.f(ast)?q\.gz$ ]]; then
    num="${BASH_REMATCH[1]}"
    if [[ "$base" =~ Serte[Pp]erti ]]; then
      printf '%s\n' "$EXTRA_LONG_DIR/concatenated_barcode$(printf '%02d' "$num").fastq.gz"
      printf '%s\n' "$CANMAG/BASALT_v2/concatenated_barcode$(printf '%02d' "$num").fastq.gz"
      printf '%s\n' "$CANMAG/fastq/long_reads/SertePerti_ONT_WGS_concatenated_barcode$(printf '%02d' "$num").fastq.gz"
    fi
  fi
}

resolve_one() {
  local raw="$1" candidate
  while IFS= read -r candidate; do
    [[ -n "$candidate" ]] || continue
    if [[ -s "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  done < <(candidate_paths "$raw" | awk '!seen[$0]++')

  echo "[MISSING] $raw" >&2
  candidate_paths "$raw" | awk '!seen[$0]++' | sed 's/^/          tried: /' >&2
  return 1
}

validate_fastq_gz() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    echo "[ERROR] Missing or empty FASTQ: $f" >&2
    return 1
  fi
  if [[ "$VALIDATE_GZIP" == "1" && "$f" == *.gz ]]; then
    gzip -t "$f" || { echo "[ERROR] gzip validation failed: $f" >&2; return 1; }
  fi
}

read_manifest_rows() {
  python3 - "$MANIFEST_TSV" "$DOG_FILTER" <<'PY'
import csv
import sys

path, dog_filter = sys.argv[1], sys.argv[2]
wanted = {x.strip() for x in dog_filter.replace(";", ",").split(",") if x.strip()}

with open(path, newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"dog_id", "canonical_dog_name", "long_read_fastqs"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"Missing required manifest columns: {sorted(missing)}")
    for row in reader:
        dog_id = (row.get("dog_id") or "").strip()
        dog = (row.get("canonical_dog_name") or "").strip()
        long_reads = (row.get("long_read_fastqs") or "").strip()
        if not dog_id or not dog or not long_reads:
            continue
        if dog_id.startswith("DOG_SERTEPERTI_LR_BC"):
            continue
        if (row.get("short_read_R1_fastqs") or "").strip() or (row.get("short_read_R2_fastqs") or "").strip():
            continue
        if wanted and dog_id not in wanted and dog not in wanted:
            continue
        set_name = (row.get("set") or "lr_only_manifest").strip()
        notes = (row.get("notes") or "").strip()
        values = [dog_id, dog, set_name, long_reads, notes]
        print("\t".join(v.replace("\n", " ").replace("\r", " ") for v in values))
PY
}

resolve_dog_reads() {
  local dog_id="$1" dog="$2" lrs_cell="$3" out_list="$4"
  local raw resolved key
  declare -A seen=()
  : > "$out_list"

  while IFS= read -r raw; do
    [[ -n "$raw" ]] || continue
    key="$(barcode_key "$raw")"
    if [[ -v seen["$key"] ]]; then
      echo "[SKIP] duplicate barcode alias for $dog: $raw (key=$key)" >&2
      continue
    fi
    if resolved="$(resolve_one "$raw")"; then
      validate_fastq_gz "$resolved" || return 1
      seen["$key"]=1
      printf '%s\n' "$resolved" >> "$out_list"
    else
      return 1
    fi
  done < <(split_fastq_cell "$lrs_cell")

  [[ -s "$out_list" ]] || { echo "[ERROR] No resolved long reads for $dog_id $dog" >&2; return 1; }
}

stage_links() {
  local dog="$1" resolved_list="$2" input_dir="$3"
  local i=0 src dest
  mkdir -p "$input_dir"
  while IFS= read -r src; do
    i=$((i + 1))
    dest="$input_dir/$(printf '%03d' "$i")__$(basename "$src")"
    if [[ -L "$dest" || -e "$dest" ]]; then
      if [[ "$FORCE" == "1" ]]; then
        rm -f "$dest"
      else
        continue
      fi
    fi
    ln -s "$src" "$dest"
  done < "$resolved_list"
  echo "[INFO] staged $(find "$input_dir" -maxdepth 1 -type l -name '*.f*q*' | wc -l) links for $dog"
}

run_coassembly_once() {
  local input_dir="$1" out_dir="$2" stdout_log="$3" stderr_log="$4" extra_args="$5"
  local args=(
    --input-dir "$input_dir"
    --output-dir "$out_dir"
    --threads "$THREADS"
    --flye-type "$READ_TYPE"
  )

  if [[ "$META" != "1" ]]; then
    args+=(--no-meta)
  fi
  if [[ -n "$extra_args" ]]; then
    args+=(--flye-extra-args "$extra_args")
  fi

  bash "$COASSEMBLY_SCRIPT" "${args[@]}" > "$stdout_log" 2> "$stderr_log"
}

retry_without_polishing_if_needed() {
  local dog="$1" dog_slug="$2" input_dir="$3" out_dir="$4"
  local backup_dir no_polish_extra

  [[ "$RETRY_NO_POLISH_ON_FAILURE" == "1" ]] || return 1
  [[ -s "$out_dir/flye_coassembly.log" ]] || return 1
  grep -q 'Bubble format error' "$out_dir/flye_coassembly.log" || return 1

  backup_dir="${out_dir}.failed_polishing_$(date +%Y%m%d_%H%M%S)"
  echo "[WARN] Flye polishing failed for $dog with Bubble format error; moving failed polished attempt to $backup_dir"
  mv "$out_dir" "$backup_dir"

  no_polish_extra="$FLYE_EXTRA_ARGS"
  if [[ " $no_polish_extra " != *" --iterations "* ]]; then
    no_polish_extra="${no_polish_extra:+$no_polish_extra }--iterations 0"
  fi

  echo "[RETRY] Starting fresh no-polish Flye run for $dog"
  run_coassembly_once \
    "$input_dir" \
    "$out_dir" \
    "$BASE/logs/${dog_slug}.no_polish.stdout.log" \
    "$BASE/logs/${dog_slug}.no_polish.stderr.log" \
    "$no_polish_extra"
}

run_dog() {
  local dog_id="$1" dog="$2" set_name="$3" lrs="$4" notes="$5" mode="$6"
  local dog_slug input_dir out_dir status_done status_failed status_prepared resolved_list count
  dog_slug="$(sanitize_name "$dog")"
  input_dir="$BASE/inputs/$dog_slug"
  out_dir="$BASE/results/lrs_flye_meta__${dog_slug}"
  status_done="$BASE/status/$dog_slug.done"
  status_failed="$BASE/status/$dog_slug.failed"
  status_prepared="$BASE/status/$dog_slug.prepared"
  resolved_list="$BASE/metadata/${dog_slug}.long_reads.resolved.list"

  echo "[DOG] $dog ($dog_id, $set_name)"
  rm -f "$status_failed"

  if [[ "$FORCE" != "1" && "$mode" == "run" && -s "$out_dir/assembly.fasta" ]]; then
    date '+%Y-%m-%d %H:%M:%S' > "$status_done"
    echo "[SKIP] complete assembly exists: $out_dir/assembly.fasta"
    return 0
  fi

  if ! resolve_dog_reads "$dog_id" "$dog" "$lrs" "$resolved_list"; then
    echo "[FAILED] input resolution failed for $dog" | tee "$status_failed" >&2
    return 1
  fi

  count="$(wc -l < "$resolved_list" | tr -d ' ')"
  printf '%s\t%s\t%s\t%s\t%s\n' "$dog_id" "$dog" "$set_name" "$count" "$(paste -sd ';' "$resolved_list")" >> "$BASE/metadata/flye_lr_only_resolved_queue.tsv"

  if [[ "$mode" == "check" ]]; then
    echo "[OK] resolved $count long-read FASTQ(s) for $dog"
    return 0
  fi

  stage_links "$dog" "$resolved_list" "$input_dir"
  date '+%Y-%m-%d %H:%M:%S' > "$status_prepared"

  if [[ "$mode" == "prepare" ]]; then
    echo "[PREPARED] $dog"
    return 0
  fi

  log "[RUN] Flye --meta $dog"
  date '+%Y-%m-%d %H:%M:%S' > "$BASE/status/$dog_slug.started"

  local run_ok=0
  if [[ "$RESUME_PARTIAL" == "1" && -d "$out_dir" && -n "$(find "$out_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" && ! -s "$out_dir/assembly.fasta" ]]; then
    if [[ "$RETRY_NO_POLISH_ON_FAILURE" == "1" && -s "$out_dir/flye_coassembly.log" ]] && grep -q 'Bubble format error' "$out_dir/flye_coassembly.log"; then
      echo "[RESCUE] Existing partial Flye directory already has Bubble format error for $dog"
      if retry_without_polishing_if_needed "$dog" "$dog_slug" "$input_dir" "$out_dir"; then
        run_ok=1
      fi
    else
      echo "[RESUME] Existing partial Flye directory found for $dog: $out_dir"
      if flye --resume \
          -o "$out_dir" \
          -t "$THREADS" \
          > "$BASE/logs/${dog_slug}.resume.stdout.log" \
          2> "$BASE/logs/${dog_slug}.resume.stderr.log"; then
        run_ok=1
      elif [[ "$RESTART_FAILED_RESUME" == "1" ]]; then
        local backup_dir
        backup_dir="${out_dir}.failed_resume_$(date +%Y%m%d_%H%M%S)"
        echo "[WARN] Flye resume failed for $dog; moving partial directory to $backup_dir"
        mv "$out_dir" "$backup_dir"
        echo "[RESTART] Starting fresh Flye run for $dog"
        if run_coassembly_once \
            "$input_dir" \
            "$out_dir" \
            "$BASE/logs/${dog_slug}.restart.stdout.log" \
            "$BASE/logs/${dog_slug}.restart.stderr.log" \
            "$FLYE_EXTRA_ARGS"; then
          run_ok=1
        fi
      fi
    fi
  else
    if run_coassembly_once \
        "$input_dir" \
        "$out_dir" \
        "$BASE/logs/${dog_slug}.stdout.log" \
        "$BASE/logs/${dog_slug}.stderr.log" \
        "$FLYE_EXTRA_ARGS"; then
      run_ok=1
    fi
  fi

  if [[ "$run_ok" != "1" ]]; then
    if retry_without_polishing_if_needed "$dog" "$dog_slug" "$input_dir" "$out_dir"; then
      run_ok=1
    fi
  fi

  if [[ "$run_ok" == "1" ]]; then
    if [[ -s "$out_dir/assembly.fasta" ]]; then
      date '+%Y-%m-%d %H:%M:%S' > "$status_done"
      log "[DONE] $dog -> $out_dir/assembly.fasta"
    else
      echo "[FAILED] Flye exited but assembly.fasta is missing for $dog" | tee "$status_failed" >&2
      return 1
    fi
  else
    echo "[FAILED] $dog; see $BASE/logs/${dog_slug}.stderr.log, $BASE/logs/${dog_slug}.resume.stderr.log, $BASE/logs/${dog_slug}.restart.stderr.log, or $BASE/logs/${dog_slug}.no_polish.stderr.log" | tee "$status_failed" >&2
    return 1
  fi
}

run_queue() {
  local mode="$1"
  local dog_id dog set_name lrs notes

  : > "$BASE/metadata/flye_lr_only_resolved_queue.tsv"
  printf 'dog_id\tdog\tset\tresolved_long_read_count\tresolved_long_reads\n' >> "$BASE/metadata/flye_lr_only_resolved_queue.tsv"
  cp "$MANIFEST_TSV" "$BASE/metadata/input_flye_lr_only_manifest.tsv"

  read_manifest_rows | while IFS=$'\t' read -r dog_id dog set_name lrs notes; do
    [[ -n "$dog_id" ]] || continue
    if ! run_dog "$dog_id" "$dog" "$set_name" "$lrs" "$notes" "$mode"; then
      [[ "$CONTINUE_ON_DOG_ERROR" == "1" ]] && continue || exit 1
    fi
  done
}

case "${1:-run}" in
  list)
    read_manifest_rows | cut -f1-3
    ;;
  check)
    need_cmd python3
    [[ -f "$MANIFEST_TSV" ]] || { echo "[ERROR] Missing MANIFEST_TSV: $MANIFEST_TSV" >&2; exit 1; }
    echo "[INFO] BASE=$BASE"
    echo "[INFO] CANMAG=$CANMAG"
    echo "[INFO] MANIFEST_TSV=$MANIFEST_TSV"
    echo "[INFO] LONG_DIR=$LONG_DIR"
    echo "[INFO] EXTRA_LONG_DIR=$EXTRA_LONG_DIR"
    echo "[INFO] DOG_FILTER=${DOG_FILTER:-<none>}"
    echo "[INFO] VALIDATE_GZIP=$VALIDATE_GZIP"
    echo "[INFO] RESUME_PARTIAL=$RESUME_PARTIAL"
    echo "[INFO] RESTART_FAILED_RESUME=$RESTART_FAILED_RESUME"
    echo "[INFO] RETRY_NO_POLISH_ON_FAILURE=$RETRY_NO_POLISH_ON_FAILURE"
    echo "[INFO] FLYE_EXTRA_ARGS=${FLYE_EXTRA_ARGS:-<none>}"
    run_queue check
    ;;
  prepare)
    need_cmd python3
    need_cmd ln
    [[ -f "$MANIFEST_TSV" ]] || { echo "[ERROR] Missing MANIFEST_TSV: $MANIFEST_TSV" >&2; exit 1; }
    run_queue prepare
    ;;
  run)
    need_cmd python3
    need_cmd ln
    need_cmd flye
    [[ -f "$COASSEMBLY_SCRIPT" ]] || { echo "[ERROR] Missing COASSEMBLY_SCRIPT: $COASSEMBLY_SCRIPT" >&2; exit 1; }
    [[ -f "$MANIFEST_TSV" ]] || { echo "[ERROR] Missing MANIFEST_TSV: $MANIFEST_TSV" >&2; exit 1; }
    run_queue run
    ;;
  *)
    echo "Usage: bash run_flye_lr_only_dogfirst.sh [list|check|prepare|run]" >&2
    echo "Optional env: BASE=... CANMAG=... MANIFEST_TSV=... LONG_DIR=... EXTRA_LONG_DIR=... THREADS=20 READ_TYPE=--nano-hq DOG_FILTER=DOG_AJSA,Ajsa CONTINUE_ON_DOG_ERROR=1 VALIDATE_GZIP=1 RESUME_PARTIAL=1 RESTART_FAILED_RESUME=1 RETRY_NO_POLISH_ON_FAILURE=1 FLYE_EXTRA_ARGS='--iterations 0' FORCE=1" >&2
    exit 1
    ;;
esac

echo "[DONE] Flye LR-only dog-first workflow mode=${1:-run}"
