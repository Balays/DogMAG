#!/usr/bin/env bash
set -euo pipefail

basalt_run=""
assembly_dir=""
output_dir="work/viral_contigs"
min_length="1500"
threads="16"
prefix_headers="false"
run_callers="false"

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_viral_contig_workflow.sh --basalt-run DIR --assembly-dir DIR [options]

Options:
  --basalt-run DIR       BASALT run directory. Repeat by passing a comma-separated list.
  --assembly-dir DIR     Directory with assembly FASTAs to screen.
  --output-dir DIR       Output directory. Default: work/viral_contigs
  --min-length INT       Minimum contig length. Default: 1500
  --threads INT          Threads for external viral callers. Default: 16
  --prefix-headers       Prefix FASTA headers with assembly_id|.
  --run-callers          Also run geNomad and CheckV if available in the active environment.
  --help                 Show this help.

Outputs:
  basalt_contig_assignments.tsv
  basalt_binned_contigs.txt
  viral_screen_all_contigs.fna
  viral_screen_unbinned_contigs.fna
  viral_screen_binned_contigs_for_prophages.fna
  viral_contig_context.tsv
  commands_viral_callers.sh
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --basalt-run)
      basalt_run="$2"
      shift 2
      ;;
    --assembly-dir)
      assembly_dir="$2"
      shift 2
      ;;
    --output-dir)
      output_dir="$2"
      shift 2
      ;;
    --min-length)
      min_length="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --prefix-headers)
      prefix_headers="true"
      shift
      ;;
    --run-callers)
      run_callers="true"
      shift
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$basalt_run" || -z "$assembly_dir" ]]; then
  echo "Missing required --basalt-run or --assembly-dir." >&2
  usage >&2
  exit 2
fi

mkdir -p "$output_dir"

mask_args=()
IFS=',' read -r -a basalt_runs <<< "$basalt_run"
for run_dir in "${basalt_runs[@]}"; do
  mask_args+=(--basalt-run "$run_dir")
done

python scripts/make_basalt_contig_mask.py \
  "${mask_args[@]}" \
  --output "$output_dir/basalt_contig_assignments.tsv" \
  --contig-list "$output_dir/basalt_binned_contigs.txt"

prepare_args=()
if [[ "$prefix_headers" == "true" ]]; then
  prepare_args+=(--prefix-headers)
fi

python scripts/prepare_viral_contig_inputs.py \
  --assembly-dir "$assembly_dir" \
  --recursive \
  --basalt-assignments "$output_dir/basalt_contig_assignments.tsv" \
  --output-dir "$output_dir" \
  --min-length "$min_length" \
  "${prepare_args[@]}"

cat > "$output_dir/commands_viral_callers.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail

# Run viral discovery on all contigs so binned prophages are not lost.
genomad end-to-end --cleanup --threads "$threads" "$output_dir/viral_screen_all_contigs.fna" "$output_dir/genomad_all" "\${GENOMAD_DB:?set GENOMAD_DB}"

# Optional: run unbinned contigs separately for a cleaner free-virus subset.
genomad end-to-end --cleanup --threads "$threads" "$output_dir/viral_screen_unbinned_contigs.fna" "$output_dir/genomad_unbinned" "\${GENOMAD_DB:?set GENOMAD_DB}"

# CheckV expects a viral candidate FASTA. Replace the input below with geNomad/VirSorter2
# candidate sequences after reviewing the caller output.
checkv end_to_end "$output_dir/viral_screen_all_contigs.fna" "$output_dir/checkv_all" -t "$threads" -d "\${CHECKVDB:?set CHECKVDB}"
EOF

chmod +x "$output_dir/commands_viral_callers.sh" 2>/dev/null || true

if [[ "$run_callers" == "true" ]]; then
  bash "$output_dir/commands_viral_callers.sh"
fi

echo "Prepared viral contig workflow inputs in $output_dir"
