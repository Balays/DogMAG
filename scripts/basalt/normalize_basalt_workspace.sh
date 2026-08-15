#!/usr/bin/env bash
set -euo pipefail

workspace=""
nested_name="basalt_run"
apply=0

usage() {
  cat <<'EOF'
Usage:
  bash scripts/normalize_basalt_workspace.sh --workspace /path/to/workspace [--nested-name basalt_run] [--apply]

Default behavior is dry-run.
With --apply, nested workspace contents are copied into the workspace root with rsync --ignore-existing.
Trailing carriage-return characters in root filenames and assembly filenames are also normalized away.
EOF
}

strip_cr_filenames() {
  local target_dir="$1"
  local renamed=0

  [[ -d "$target_dir" ]] || return 0

  while IFS= read -r -d '' path; do
    local base clean dir
    base="$(basename "$path")"
    clean="${base%$'\r'}"
    dir="$(dirname "$path")"
    if [[ "$base" != "$clean" ]]; then
      if [[ "$apply" -eq 1 ]]; then
        mv -- "$path" "$dir/$clean"
      else
        echo "[PLAN] Would rename: $path -> $dir/$clean"
      fi
      renamed=$((renamed + 1))
    fi
  done < <(find "$target_dir" -maxdepth 1 -name "*"$'\r' -print0)

  if [[ "$renamed" -gt 0 ]]; then
    echo "[INFO] $([[ "$apply" -eq 1 ]] && echo Renamed || echo Found) $renamed carriage-return-suffixed entr$( [[ "$renamed" -eq 1 ]] && echo y || echo ies ) in $target_dir"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workspace)
      workspace="${2:-}"
      shift 2
      ;;
    --nested-name)
      nested_name="${2:-}"
      shift 2
      ;;
    --apply)
      apply=1
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

if [[ -z "$workspace" ]]; then
  echo "Missing required --workspace" >&2
  exit 1
fi

if [[ ! -d "$workspace" ]]; then
  echo "Workspace does not exist: $workspace" >&2
  exit 1
fi

nested_dir="${workspace}/${nested_name}"
manifest="${workspace}/metadata/basalt_input_manifest.tsv"

expected_assemblies=0
expected_srs=0
expected_lrs=0

if [[ -f "$manifest" ]]; then
  expected_assemblies=$(awk -F'\t' 'NR>1 {count++} END {print count+0}' "$manifest")
  expected_srs=$(awk -F'\t' 'NR>1 && $3=="SRS" {count++} END {print count+0}' "$manifest")
  expected_lrs=$(awk -F'\t' 'NR>1 && $3=="LRS" {count++} END {print count+0}' "$manifest")
fi

echo "[INFO] Workspace: $workspace"
echo "[INFO] Nested candidate: $nested_dir"
echo "[INFO] Mode: $([[ "$apply" -eq 1 ]] && echo apply || echo dry-run)"

if [[ -d "$nested_dir" ]]; then
  echo "[INFO] Nested workspace exists."
  echo "[INFO] Entries inside nested workspace:"
  find "$nested_dir" -mindepth 1 -maxdepth 1 -printf '  %f\n' | sort

  if [[ "$apply" -eq 1 ]]; then
    if ! command -v rsync >/dev/null 2>&1; then
      echo "rsync is required for --apply" >&2
      exit 1
    fi
    echo "[INFO] Copying nested contents into workspace root"
    rsync -a --ignore-existing "$nested_dir"/ "$workspace"/
  else
    echo "[PLAN] Would run:"
    echo "       rsync -a --ignore-existing \"$nested_dir/\" \"$workspace/\""
  fi
else
  echo "[INFO] No nested workspace named '$nested_name' found."
fi

strip_cr_filenames "$workspace"
strip_cr_filenames "$workspace/assemblies"

top_assemblies=0
top_sr_r1=0
top_lr=0

if [[ -d "$workspace/assemblies" ]]; then
  top_assemblies=$(find "$workspace/assemblies" -maxdepth 1 -type f | wc -l | tr -d ' ')
fi

top_sr_r1=$(find "$workspace" -maxdepth 1 -type f \( -name 'sr__*R1*.fq' -o -name 'sr__*R1*.fastq' \) | wc -l | tr -d ' ')
top_lr=$(find "$workspace" -maxdepth 1 -type f \( -name 'lr__*.fq' -o -name 'lr__*.fastq' \) | wc -l | tr -d ' ')

echo "[INFO] Expected from manifest: assemblies=$expected_assemblies srs=$expected_srs lrs=$expected_lrs"
echo "[INFO] Found at workspace root: assemblies=$top_assemblies sr_r1=$top_sr_r1 lr=$top_lr"

if [[ "$expected_assemblies" -gt 0 && "$top_assemblies" -eq "$expected_assemblies" && "$top_sr_r1" -eq "$expected_srs" && "$top_lr" -eq "$expected_lrs" ]]; then
  echo "[OK] Workspace layout looks runnable."
else
  echo "[WARN] Workspace root counts do not yet match the manifest."
  echo "[WARN] Do not remove the nested workspace until the counts match and run_basalt.sh is verified."
fi
