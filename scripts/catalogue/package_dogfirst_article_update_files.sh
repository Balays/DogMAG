#!/usr/bin/env bash
set -euo pipefail

# Package the small metadata/log/summary files needed to update the article.
# This deliberately avoids large FASTA/MMI/BAM/SAM files.

PROJECT_DIR="${PROJECT_DIR:-/path/to/DogMAG_workdir}"
DEPLETION_OUT="${DEPLETION_OUT:-${PROJECT_DIR}/CanMAG_depletion_panels_20260721_indexed}"
GTDB_OUT="${GTDB_OUT:-${PROJECT_DIR}/work/dogfirst_gtdbtk_reselected_drep99_20260722}"
VIRAL_OUT="${VIRAL_OUT:-${PROJECT_DIR}/work/dogfirst_viral_contigs_final_assemblies_20260722}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/article_update_package_$(date +%Y%m%d_%H%M%S)}"
TARBALL="${TARBALL:-${OUTDIR}.tar.gz}"
INCLUDE_VIRAL="${INCLUDE_VIRAL:-auto}" # auto, 0, or 1

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

relpath() {
  python3 - "$PROJECT_DIR" "$1" <<'PY'
from pathlib import Path
import os
import sys

root = Path(sys.argv[1]).resolve()
path = Path(sys.argv[2]).resolve()
try:
    print(path.relative_to(root))
except ValueError:
    print(Path("external") / path.name)
PY
}

add_file() {
  local src="$1"
  if [[ -s "$src" ]]; then
    local rel
    rel="$(relpath "$src")"
    mkdir -p "$OUTDIR/$(dirname "$rel")"
    cp -a "$src" "$OUTDIR/$rel"
    printf "present\t%s\t%s\n" "$src" "$rel" >> "$OUTDIR/package_manifest.tsv"
  else
    printf "missing\t%s\t\n" "$src" >> "$OUTDIR/package_manifest.tsv"
  fi
}

add_dir_files() {
  local src_dir="$1"
  if [[ ! -d "$src_dir" ]]; then
    printf "missing_dir\t%s\t\n" "$src_dir" >> "$OUTDIR/package_manifest.tsv"
    return
  fi

  while IFS= read -r -d '' src; do
    add_file "$src"
  done < <(
    find "$src_dir" -maxdepth 1 -type f \
      \( -name '*.tsv' -o -name '*.txt' -o -name '*.csv' -o -name '*.log' -o -name '*.json' \) \
      -print0
  )
}

mkdir -p "$OUTDIR"
printf "status\tsource_path\tpackaged_relative_path\n" > "$OUTDIR/package_manifest.tsv"

log "PROJECT_DIR=$PROJECT_DIR"
log "DEPLETION_OUT=$DEPLETION_OUT"
log "GTDB_OUT=$GTDB_OUT"
log "VIRAL_OUT=$VIRAL_OUT"
log "OUTDIR=$OUTDIR"

# Depletion/reselection summaries.
add_file "$DEPLETION_OUT/tables/reselected_best_candidates.tsv"
add_file "$DEPLETION_OUT/tables/groups_without_medium_quality_candidate.tsv"
add_file "$DEPLETION_OUT/tables/reselected_fasta_manifest.tsv"
add_file "$DEPLETION_OUT/tables/CanMAG_depletion_95ANI.genomes.tsv"
add_file "$DEPLETION_OUT/tables/CanMAG_depletion_99ANI.genomes.tsv"
add_file "$DEPLETION_OUT/tables/CanMAG_depletion_95ANI.stats.txt"
add_file "$DEPLETION_OUT/tables/CanMAG_depletion_99ANI.stats.txt"
add_file "$DEPLETION_OUT/panels/SHA256SUMS"
add_file "$DEPLETION_OUT/run.indexed_fasta_lookup.log"
add_file "$DEPLETION_OUT/run.panels_only_mm2fast.log"

# Catalogue-unit summaries from the comparison workflow.
add_file "$DEPLETION_OUT/catalog_unit_summary/catalog_unit_summary.tsv"
add_file "$DEPLETION_OUT/catalog_unit_summary/quality_summary.tsv"
add_file "$DEPLETION_OUT/catalog_unit_summary/quality_counts_by_group.tsv"
add_file "$DEPLETION_OUT/catalog_unit_summary/comparison_context.tsv"
add_file "$DEPLETION_OUT/catalog_unit_summary/manifest_quality_pass_bins.tsv"
add_file "$DEPLETION_OUT/catalog_unit_summary/manifest_all_bins.tsv"

# dRep winner/member tables and warnings.
add_file "$DEPLETION_OUT/drep95/data_tables/Widb.csv"
add_file "$DEPLETION_OUT/drep95/data_tables/Cdb.csv"
add_file "$DEPLETION_OUT/drep95/data_tables/Ndb.csv"
add_file "$DEPLETION_OUT/drep95/data_tables/Wdb.csv"
add_file "$DEPLETION_OUT/drep95/log/warnings.txt"
add_file "$DEPLETION_OUT/drep99/data_tables/Widb.csv"
add_file "$DEPLETION_OUT/drep99/data_tables/Cdb.csv"
add_file "$DEPLETION_OUT/drep99/data_tables/Ndb.csv"
add_file "$DEPLETION_OUT/drep99/data_tables/Wdb.csv"
add_file "$DEPLETION_OUT/drep99/log/warnings.txt"

# GTDB-Tk taxonomy deliverables and metadata.
add_file "$GTDB_OUT/mag_taxonomy_gtdbtk_summary.tsv"
add_file "$GTDB_OUT/mag_taxonomy_input_metadata.tsv"
add_file "$GTDB_OUT/gtdbtk_batchfile.tsv"
add_file "$GTDB_OUT/excluded_genomes.tsv"
add_file "$GTDB_OUT/deliverables/gtdbtk_taxonomy.tsv"
add_file "$GTDB_OUT/deliverables/README.txt"
add_file "$GTDB_OUT/logs/gtdbtk_classify_wf.log"

# Viral workflow deliverables, if present.
if [[ "$INCLUDE_VIRAL" == "1" || ( "$INCLUDE_VIRAL" == "auto" && -d "$VIRAL_OUT/deliverables" ) ]]; then
  add_dir_files "$VIRAL_OUT/deliverables"
  add_file "$VIRAL_OUT/logs/prepare_viral_contigs.log"
  add_file "$VIRAL_OUT/logs/genomad_unbinned.log"
  add_file "$VIRAL_OUT/logs/checkv_unbinned.log"
  add_file "$VIRAL_OUT/logs/genomad_binned_prophages.log"
  add_file "$VIRAL_OUT/logs/checkv_binned_prophages.log"
fi

# Lightweight inventory for confidence.
{
  echo "package_created_at	$(date '+%Y-%m-%d %H:%M:%S')"
  echo "project_dir	$PROJECT_DIR"
  echo "depletion_out	$DEPLETION_OUT"
  echo "gtdb_out	$GTDB_OUT"
  echo "viral_out	$VIRAL_OUT"
  echo "include_viral	$INCLUDE_VIRAL"
  echo "present_files	$(awk -F'\t' 'NR>1 && $1=="present"{n++} END{print n+0}' "$OUTDIR/package_manifest.tsv")"
  echo "missing_files	$(awk -F'\t' 'NR>1 && $1 ~ /^missing/{n++} END{print n+0}' "$OUTDIR/package_manifest.tsv")"
} > "$OUTDIR/package_summary.tsv"

tmp_tar="${TARBALL}.tmp"
rm -f "$tmp_tar" "$TARBALL"
tar -C "$(dirname "$OUTDIR")" -czf "$tmp_tar" "$(basename "$OUTDIR")"
mv "$tmp_tar" "$TARBALL"

log "Package directory: $OUTDIR"
log "Tarball: $TARBALL"
log "Summary:"
cat "$OUTDIR/package_summary.tsv" >&2
