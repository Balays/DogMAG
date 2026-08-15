#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Make a continued BASALT run independent of the earlier run directory.

Required environment variables:
  OLD       Earlier BASALT run currently providing symlink targets
  FINAL     Continued/completed BASALT run that must be retained

Optional environment variables:
  STATE_DIR  Audit directory outside OLD and FINAL
  RAW_MODE   direct-link (default) or copy

Modes:
  audit        Inventory and validate FINAL symlinks into OLD; no changes
  materialize Replace those symlinks without changing or deleting OLD
  verify       Verify every recorded dependency and core FINAL outputs

For RAW_MODE=direct-link, raw-read aliases are rewired directly to their
ultimate canonical files outside OLD. All other OLD-backed dependencies are
copied into their existing paths under FINAL. This script never deletes OLD.

Example:
  OLD=/path/to/DogMAG_workdir/basalt_run_dogfirst_full_20260518 \
  FINAL=/path/to/DogMAG_workdir/basalt_run_dogfirst_full_20260518_wanted_only_copybins_checkm2_padidx_20260608 \
  bash scripts/make_basalt_rerun_self_contained.sh audit
USAGE
}

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"
}

die() {
  log "ERROR: $*" >&2
  exit 1
}

line_count() {
  awk 'END { print NR + 0 }' "$1"
}

fasta_count() {
  find "$1" -maxdepth 1 -type f \
    \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) \
    -size +0c -printf '.' | wc -c
}

classify_target() {
  local target="$1"
  local name="${target##*/}"

  case "$name" in
    lr__*.fastq.gz|lr__*.fq.gz|lr__*.fastq|lr__*.fq) printf 'long_read_fastq' ;;
    sr__*.fastq.gz|sr__*.fq.gz|sr__*.fastq|sr__*.fq) printf 'short_read_fastq' ;;
    *_assembly.depth.txt) printf 'assembly_depth' ;;
    Bins_total_connections_*.txt) printf 'bin_connections' ;;
    condense_connections_*.txt) printf 'condensed_connections' ;;
    *.fa|*.fna|*.fasta) printf 'source_assembly_fasta' ;;
    *)
      if [[ "$target" == *checkm* ]]; then
        printf 'checkm_output'
      else
        printf 'other'
      fi
      ;;
  esac
}

check_core_outputs() {
  local final_bin="$FINAL/BestBinset_outlier_refined_MAGs_polished_re-assembly"
  local quality_report="$final_bin/Best_binset_quality_report.tsv"
  local comparison="$FINAL/Reassembled_bins_comparison.txt"
  local bin_map="$FINAL/BestBinset_outlier_refined_mod/Bin_name_mod.txt"
  local bin_map_gf_lr="$FINAL/BestBinset_outlier_refined_gf_lr_mod/Bin_name_mod.txt"

  [[ -d "$final_bin" ]] || die "Missing final bin directory: $final_bin"
  [[ -f "$quality_report" ]] || die "Missing quality report: $quality_report"
  [[ -f "$comparison" ]] || die "Missing comparison table: $comparison"
  [[ -f "$bin_map" ]] || die "Missing bin map: $bin_map"
  [[ -f "$bin_map_gf_lr" ]] || die "Missing GF/LR bin map: $bin_map_gf_lr"

  local fastas links quality_rows comparison_rows map_rows gf_lr_rows
  fastas="$(fasta_count "$final_bin")"
  links="$(find "$final_bin" -maxdepth 1 -type l -printf '.' | wc -c)"
  quality_rows="$(line_count "$quality_report")"
  comparison_rows="$(line_count "$comparison")"
  map_rows="$(line_count "$bin_map")"
  gf_lr_rows="$(line_count "$bin_map_gf_lr")"

  printf 'final_fastas\t%s\n' "$fastas"
  printf 'final_fasta_symlinks\t%s\n' "$links"
  printf 'quality_report_rows\t%s\n' "$quality_rows"
  printf 'comparison_rows\t%s\n' "$comparison_rows"
  printf 'bin_map_rows\t%s\n' "$map_rows"
  printf 'bin_map_gf_lr_rows\t%s\n' "$gf_lr_rows"

  [[ "$fastas" == 11276 ]] || die "Expected 11276 final FASTAs, found $fastas"
  [[ "$links" == 0 ]] || die "Final FASTA directory contains $links symlinks"
  [[ "$quality_rows" == 11277 ]] || die "Expected 11277 quality-report rows, found $quality_rows"
  [[ "$comparison_rows" == 30556 ]] || die "Expected 30556 comparison rows, found $comparison_rows"
  [[ "$map_rows" == 11276 ]] || die "Expected 11276 bin-map rows, found $map_rows"
  [[ "$gf_lr_rows" == 11276 ]] || die "Expected 11276 GF/LR map rows, found $gf_lr_rows"
}

check_no_active_processes() {
  local active
  active="$(pgrep -af 'BASALT.py|spades.py|spades-core|checkm2|diamond|bowtie2|bwa|samtools' || true)"
  if [[ -n "$active" ]] && grep -F -e "$OLD" -e "$FINAL" <<<"$active" >/dev/null; then
    printf '%s\n' "$active" >&2
    die "A BASALT-related process still references OLD or FINAL"
  fi
}

audit_links() {
  local manifest="$STATE_DIR/symlinks_to_old.current.tsv"
  local original_manifest="$STATE_DIR/symlinks_to_old.tsv"
  local tmp="$manifest.tmp.$$"
  local link target category target_kind ultimate bytes link_rel target_rel same_relative
  local links=0 missing=0 mismatched=0 uncovered=0 raw_inside_old=0

  printf 'link_path\ttarget_path\tcategory\ttarget_kind\tbytes\tultimate_target\tsame_relative_path\n' > "$tmp"

  while IFS= read -r -d '' link; do
    target="$(readlink -- "$link")"
    [[ "$target" == "$OLD/"* ]] || continue
    links=$((links + 1))
    category="$(classify_target "$target")"
    ultimate="$(readlink -f -- "$target" 2>/dev/null || true)"

    if [[ -L "$target" ]]; then
      target_kind=symlink
    elif [[ -f "$target" ]]; then
      target_kind=file
    elif [[ -d "$target" ]]; then
      target_kind=directory
    else
      target_kind=missing
      missing=$((missing + 1))
    fi

    if [[ -n "$ultimate" && -e "$ultimate" ]]; then
      bytes="$(stat -Lc '%s' -- "$ultimate" 2>/dev/null || du -sb -- "$ultimate" | awk '{print $1}')"
    else
      bytes=0
    fi

    link_rel="${link#"$FINAL"/}"
    target_rel="${target#"$OLD"/}"
    if [[ "$link_rel" == "$target_rel" ]]; then
      same_relative=yes
    else
      same_relative=no
      mismatched=$((mismatched + 1))
    fi

    if [[ "$category" == long_read_fastq || "$category" == short_read_fastq ]]; then
      if [[ -z "$ultimate" || "$ultimate" == "$OLD" || "$ultimate" == "$OLD/"* ]]; then
        raw_inside_old=$((raw_inside_old + 1))
      fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$link" "$target" "$category" "$target_kind" "$bytes" "$ultimate" "$same_relative" >> "$tmp"
  done < <(find "$FINAL" -type l -print0)

  {
    head -n 1 "$tmp"
    tail -n +2 "$tmp" | LC_ALL=C sort -t $'\t' -k1,1
  } > "$manifest"
  rm -f -- "$tmp"

  # Freeze the first complete dependency inventory. Later audits describe only
  # links that still point into OLD, while materialize/verify retain all rows.
  if [[ ! -f "$original_manifest" ]]; then
    cp -a -- "$manifest" "$original_manifest"
  fi

  # Some BASALT connection files are intentionally linked at two FINAL paths.
  # Alternate-path links are safe only when the OLD-relative canonical path is
  # also represented in FINAL, so historical absolute paths remain recoverable.
  while IFS= read -r target; do
    [[ -n "$target" ]] || continue
    target_rel="${target#"$OLD"/}"
    link="$FINAL/$target_rel"
    if [[ -L "$link" ]]; then
      [[ "$(readlink -- "$link")" == "$target" ]] || uncovered=$((uncovered + 1))
    elif [[ ! -e "$link" ]]; then
      uncovered=$((uncovered + 1))
    fi
  done < <(tail -n +2 "$manifest" | cut -f2 | LC_ALL=C sort -u)

  {
    printf 'metric\tvalue\n'
    printf 'links_to_old\t%s\n' "$links"
    printf 'missing_targets\t%s\n' "$missing"
    printf 'relative_path_mismatches\t%s\n' "$mismatched"
    printf 'target_paths_without_final_counterpart\t%s\n' "$uncovered"
    printf 'raw_targets_still_inside_old\t%s\n' "$raw_inside_old"
    awk -F '\t' 'NR > 1 { n[$3]++ } END { for (k in n) print "category_" k "\t" n[k] }' "$manifest" | LC_ALL=C sort
  } | tee "$STATE_DIR/audit_summary.tsv"

  [[ "$links" -gt 0 ]] || die "No FINAL symlinks into OLD were found"
  [[ "$missing" == 0 ]] || die "$missing OLD-backed targets are missing"
  [[ "$uncovered" == 0 ]] || die "$uncovered OLD-relative target paths have no counterpart in FINAL"
  if [[ "$RAW_MODE" == direct-link ]]; then
    [[ "$raw_inside_old" == 0 ]] || die "$raw_inside_old raw-read targets cannot be rewired outside OLD"
  fi

  log "Current-link audit: $manifest"
  log "Original dependency manifest: $original_manifest"
}

materialize_links() {
  local manifest="$STATE_DIR/symlinks_to_old.tsv"
  local actions="$STATE_DIR/materialization_actions.tsv"
  local link target category target_kind bytes ultimate same_relative tmp action
  local copied=0 rewired=0 skipped=0

  audit_links
  printf 'link_path\toriginal_target\tcategory\taction\tfinal_target\n' > "$actions"

  while IFS=$'\t' read -r link target category target_kind bytes ultimate same_relative; do
    [[ "$link" != link_path ]] || continue
    if [[ ! -L "$link" ]]; then
      [[ -e "$link" ]] || die "Recorded dependency is now missing: $link"
      if [[ -f "$target" && -f "$link" ]]; then
        cmp -s -- "$target" "$link" || die "Existing materialized file differs from OLD: $link"
      elif [[ -d "$target" && -d "$link" ]]; then
        diff -qr -- "$target" "$link" >/dev/null || die "Existing materialized directory differs from OLD: $link"
      else
        die "Existing materialized path has the wrong type: $link"
      fi
      skipped=$((skipped + 1))
      printf '%s\t%s\t%s\talready_materialized\t%s\n' "$link" "$target" "$category" "$link" >> "$actions"
      continue
    fi

    if [[ "$(readlink -- "$link")" != "$target" ]]; then
      if [[ "$category" == long_read_fastq || "$category" == short_read_fastq ]] && \
         [[ "$(readlink -f -- "$link" 2>/dev/null || true)" == "$ultimate" ]]; then
        skipped=$((skipped + 1))
        printf '%s\t%s\t%s\trewired_to_canonical_raw\t%s\n' "$link" "$target" "$category" "$ultimate" >> "$actions"
        continue
      fi
      die "Recorded symlink changed unexpectedly: $link"
    fi

    tmp="${link}.basalt-materialize.$$"
    [[ ! -e "$tmp" && ! -L "$tmp" ]] || die "Temporary path already exists: $tmp"

    if [[ "$category" == long_read_fastq || "$category" == short_read_fastq ]] && [[ "$RAW_MODE" == direct-link ]]; then
      [[ -n "$ultimate" && -e "$ultimate" ]] || die "Missing canonical raw input for $link"
      [[ "$ultimate" != "$OLD" && "$ultimate" != "$OLD/"* ]] || die "Raw input still resolves inside OLD: $link"
      ln -s -- "$ultimate" "$tmp"
      mv -Tf -- "$tmp" "$link"
      action=rewired_to_canonical_raw
      rewired=$((rewired + 1))
      printf '%s\t%s\t%s\t%s\t%s\n' "$link" "$target" "$category" "$action" "$ultimate" >> "$actions"
    else
      cp -aL --reflink=auto -- "$target" "$tmp"
      if [[ -d "$tmp" && ! -L "$tmp" ]]; then
        # GNU mv will not replace a symlink-to-directory with a directory.
        # Remove only the verified symlink, then install the prepared copy.
        [[ -L "$link" && "$(readlink -- "$link")" == "$target" ]] || \
          die "Refusing to replace unexpected directory path: $link"
        rm -- "$link"
        mv -T -- "$tmp" "$link"
      else
        mv -Tf -- "$tmp" "$link"
      fi
      action=copied_into_final
      copied=$((copied + 1))
      printf '%s\t%s\t%s\t%s\t%s\n' "$link" "$target" "$category" "$action" "$link" >> "$actions"
    fi
  done < "$manifest"

  printf 'copied_into_final\t%s\nrewired_to_canonical_raw\t%s\nalready_materialized\t%s\n' \
    "$copied" "$rewired" "$skipped" | tee "$STATE_DIR/materialization_summary.tsv"
  log "Materialization completed; OLD was not modified"
}

verify_materialization() {
  local manifest="$STATE_DIR/symlinks_to_old.tsv"
  local actions="$STATE_DIR/materialization_actions.tsv"
  local link target category action final_target resolved
  local checked=0 failed=0 old_links broken_links

  [[ -f "$manifest" ]] || die "Missing audit manifest: $manifest"
  [[ -f "$actions" ]] || die "Missing action manifest: $actions"

  while IFS=$'\t' read -r link target category action final_target; do
    [[ "$link" != link_path ]] || continue
    checked=$((checked + 1))

    case "$action" in
      rewired_to_canonical_raw)
        resolved="$(readlink -f -- "$link" 2>/dev/null || true)"
        if [[ ! -L "$link" || -z "$resolved" || "$resolved" != "$final_target" || ! -e "$resolved" ]]; then
          printf 'FAIL raw link: %s\n' "$link" >&2
          failed=$((failed + 1))
        fi
        ;;
      copied_into_final|already_materialized)
        if [[ -L "$link" || ! -e "$link" ]]; then
          printf 'FAIL materialized path: %s\n' "$link" >&2
          failed=$((failed + 1))
        elif [[ -f "$target" && -f "$link" ]]; then
          if ! cmp -s -- "$target" "$link"; then
            printf 'FAIL file comparison: %s\n' "$link" >&2
            failed=$((failed + 1))
          fi
        elif [[ -d "$target" && -d "$link" ]]; then
          if ! diff -qr -- "$target" "$link" >/dev/null; then
            printf 'FAIL directory comparison: %s\n' "$link" >&2
            failed=$((failed + 1))
          fi
        fi
        ;;
      *)
        printf 'FAIL unknown action %s for %s\n' "$action" "$link" >&2
        failed=$((failed + 1))
        ;;
    esac
  done < "$actions"

  old_links="$(find "$FINAL" -type l -lname "$OLD/*" -printf '.' | wc -c)"
  broken_links="$(find -L "$FINAL" -type l -printf '.' | wc -c)"

  {
    printf 'metric\tvalue\n'
    printf 'recorded_dependencies_checked\t%s\n' "$checked"
    printf 'failed_dependency_checks\t%s\n' "$failed"
    printf 'remaining_links_to_old\t%s\n' "$old_links"
    printf 'broken_links_in_final\t%s\n' "$broken_links"
  } | tee "$STATE_DIR/verification_summary.tsv"

  check_core_outputs | tee "$STATE_DIR/core_output_verification.tsv"

  [[ "$failed" == 0 ]] || die "$failed materialized dependencies failed validation"
  [[ "$old_links" == 0 ]] || die "$old_links symlinks from FINAL still point into OLD"
  [[ "$broken_links" == 0 ]] || die "$broken_links broken symlinks remain in FINAL"
  log "SELF_CONTAINED_RERUN=PASS"
}

[[ $# -le 1 ]] || {
  usage >&2
  exit 2
}

mode="${1:-audit}"
[[ "$mode" != -h && "$mode" != --help ]] || {
  usage
  exit 0
}

[[ -n "${OLD:-}" ]] || die "OLD is required"
[[ -n "${FINAL:-}" ]] || die "FINAL is required"
[[ -d "$OLD" ]] || die "OLD does not exist as a directory: $OLD"
[[ -d "$FINAL" ]] || die "FINAL does not exist: $FINAL"

OLD="$(realpath "$OLD")"
FINAL="$(realpath "$FINAL")"
PROJECT_DIR="$(dirname "$OLD")"
STATE_DIR="${STATE_DIR:-$PROJECT_DIR/basalt_rerun_self_contained_state_20260810}"
RAW_MODE="${RAW_MODE:-direct-link}"
STATE_DIR="$(realpath -m "$STATE_DIR")"

[[ "$RAW_MODE" == direct-link || "$RAW_MODE" == copy ]] || die "RAW_MODE must be direct-link or copy"
case "$STATE_DIR/" in
  "$OLD/"*|"$FINAL/"*) die "STATE_DIR must be outside OLD and FINAL" ;;
esac

mkdir -p "$STATE_DIR"
check_no_active_processes
log "MODE=$mode"
log "OLD=$OLD"
log "FINAL=$FINAL"
log "STATE_DIR=$STATE_DIR"
log "RAW_MODE=$RAW_MODE"
check_core_outputs | tee "$STATE_DIR/core_output_precheck.tsv"

case "$mode" in
  audit) audit_links ;;
  materialize) materialize_links ;;
  verify) verify_materialization ;;
  *)
    usage >&2
    die "Unknown mode: $mode"
    ;;
esac
