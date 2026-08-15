#!/usr/bin/env bash
set -euo pipefail

genome_dir="final_drep_all_current_20260504/dereplicated_genomes"
quality_report="final_drep_input_all_current_bestbins_20260504/combined_quality_report.tsv"
output_dir="work/mag_taxonomy"
min_completeness="50"
max_contamination="10"
threads="16"
run_gtdbtk="false"
include_low_quality="false"

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_mag_taxonomy_workflow.sh [options]

Options:
  --genome-dir DIR              Dereplicated MAG FASTA directory.
  --quality-report TSV          Combined BASALT quality report.
  --output-dir DIR              Output directory. Default: work/mag_taxonomy
  --min-completeness FLOAT      Default: 50
  --max-contamination FLOAT     Default: 10
  --include-low-quality         Include genomes below quality cutoffs.
  --threads INT                 GTDB-Tk threads. Default: 16
  --run-gtdbtk                  Run GTDB-Tk after preparing inputs.
  --help                        Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome-dir)
      genome_dir="$2"
      shift 2
      ;;
    --quality-report)
      quality_report="$2"
      shift 2
      ;;
    --output-dir)
      output_dir="$2"
      shift 2
      ;;
    --min-completeness)
      min_completeness="$2"
      shift 2
      ;;
    --max-contamination)
      max_contamination="$2"
      shift 2
      ;;
    --include-low-quality)
      include_low_quality="true"
      shift
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --run-gtdbtk)
      run_gtdbtk="true"
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

prep_args=()
if [[ "$include_low_quality" == "true" ]]; then
  prep_args+=(--include-low-quality)
fi

python scripts/prepare_mag_taxonomy_inputs.py \
  --genome-dir "$genome_dir" \
  --quality-report "$quality_report" \
  --output-dir "$output_dir" \
  --min-completeness "$min_completeness" \
  --max-contamination "$max_contamination" \
  --threads "$threads" \
  "${prep_args[@]}"

if [[ "$run_gtdbtk" == "true" ]]; then
  bash "$output_dir/run_gtdbtk_classify_wf.sh"
  python scripts/merge_gtdbtk_taxonomy.py \
    --metadata "$output_dir/mag_taxonomy_input_metadata.tsv" \
    --gtdbtk-dir "$output_dir/gtdbtk_out" \
    --prefix canmag \
    --output "$output_dir/mag_taxonomy_gtdbtk_summary.tsv"
fi

echo "Prepared MAG taxonomy workflow in $output_dir"
