#!/usr/bin/env bash
set -Eeuo pipefail

mode="${1:-preflight}"

RUN="${RUN:-/path/to/DogMAG_workdir/basalt_run_dogfirst_full_20260518_wanted_only_copybins_checkm2_padidx_20260608}"
FINAL="${FINAL:-BestBinset_outlier_refined_MAGs_polished_re-assembly}"
BASALT_DIR="${BASALT_DIR:-/path/to/BASALT}"
PYTHON_BIN="${PYTHON_BIN:-/path/to/basalt_env/bin/python}"
THREADS="${THREADS:-32}"
MIN_MOD_MATCH_FRAC="${MIN_MOD_MATCH_FRAC:-0.90}"
SNAPSHOT="${SNAPSHOT:-1}"

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

die() {
  log "ERROR: $*"
  exit 1
}

usage() {
  cat <<EOF
Usage:
  RUN=/path/to/basalt_run bash scripts/run_basalt_manual_final_drep_dogfirst.sh preflight
  RUN=/path/to/basalt_run bash scripts/run_basalt_manual_final_drep_dogfirst.sh run

Modes:
  inventory   Print relevant files/folders.
  preflight   Validate final folder, CheckM2 quality report, and BASALT mod mapping.
  run         Snapshot final FASTAs, then run CheckM2-aware BASALT final_drep directly.

Important env vars:
  RUN         BASALT workspace. Default: ${RUN}
  FINAL       Post-reassembly selected bin folder. Default: ${FINAL}
  BASALT_DIR  BASALT Python module folder. Default: ${BASALT_DIR}
  PYTHON_BIN  Python in basalt_env. Default avoids stale mamba activation.
  THREADS     Threads for final_drep/checkm2 fallback. Default: ${THREADS}
EOF
}

if [[ "${mode}" == "-h" || "${mode}" == "--help" ]]; then
  usage
  exit 0
fi

cd "${RUN}" || die "Cannot cd to RUN=${RUN}"

[[ -d "${FINAL}" ]] || die "Missing final folder: ${RUN}/${FINAL}"
[[ -d "${BASALT_DIR}" ]] || die "Missing BASALT_DIR=${BASALT_DIR}"
[[ -x "${PYTHON_BIN}" ]] || die "Missing executable PYTHON_BIN=${PYTHON_BIN}; set PYTHON_BIN to the basalt_env python"

expected_mod="${FINAL%%_re-assembly*}_mod"
[[ "${expected_mod}" != "${FINAL}_mod" ]] || die "Could not derive expected mod folder from FINAL=${FINAL}"

count_final_fastas() {
  find "${FINAL}" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) -size +0c | wc -l | tr -d '[:space:]'
}

quality_files() {
  find "${FINAL}" -maxdepth 1 -type f -name '*quality_report.tsv' -size +0c -printf '%f\t%s bytes\n' | sort
}

print_inventory() {
  log "RUN=${RUN}"
  log "FINAL=${FINAL}"
  log "BASALT_DIR=${BASALT_DIR}"
  log "PYTHON_BIN=${PYTHON_BIN}"
  log "THREADS=${THREADS}"
  log "Expected mod folder: ${expected_mod}"
  printf 'final_fastas\t%s\n' "$(count_final_fastas)"
  printf '\n== quality reports in FINAL ==\n'
  quality_files || true
  printf '\n== Bin_name_mod candidates ==\n'
  find . -maxdepth 3 -type f -name 'Bin_name_mod.txt' -printf '%p\t%s bytes\n' | sort || true
  printf '\n== checkpoints ==\n'
  tail -30 Basalt_checkpoint.txt 2>/dev/null || true
  printf '\n== active BASALT-like processes ==\n'
  pgrep -af 'BASALT.py|checkm2|diamond|spades.py|spades-core|bowtie2|samtools|minimap2' || true
}

validate_mod_candidate() {
  local candidate="$1"
  "${PYTHON_BIN}" - "$RUN" "$FINAL" "$candidate" "$MIN_MOD_MATCH_FRAC" <<'PY'
from pathlib import Path
import re
import sys

run = Path(sys.argv[1])
final = sys.argv[2]
candidate = Path(sys.argv[3])
min_frac = float(sys.argv[4])

final_dir = run / final
mod_file = candidate / "Bin_name_mod.txt"
if not mod_file.exists() or mod_file.stat().st_size == 0:
    print(f"status\tbad\treason\tmissing_or_empty\tcandidate\t{candidate}")
    raise SystemExit(2)

final_prefixes = set()
for p in final_dir.iterdir():
    if not p.is_file() or p.stat().st_size == 0:
        continue
    if p.suffix.lower() not in {".fa", ".fasta", ".fna"}:
        continue
    m = re.match(r"^(bin[0-9]+)", p.name)
    if m:
        final_prefixes.add(m.group(1))

mapped = set()
original_with_genomes = 0
line_count = 0
with mod_file.open() as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        line_count += 1
        original, mod = parts[0], parts[1]
        if "_genomes." in original:
            original_with_genomes += 1
        mapped.add(mod.rsplit(".fa", 1)[0].rsplit(".fasta", 1)[0].rsplit(".fna", 1)[0])

matches = len(final_prefixes & mapped)
denom = max(1, len(final_prefixes))
frac = matches / denom
print(f"candidate\t{candidate}")
print(f"line_count\t{line_count}")
print(f"original_with__genomes\t{original_with_genomes}")
print(f"final_bin_prefixes\t{len(final_prefixes)}")
print(f"mapped_bin_ids\t{len(mapped)}")
print(f"matches\t{matches}")
print(f"match_fraction\t{frac:.6f}")
if line_count == 0:
    print("status\tbad\treason\tno_mapping_rows")
    raise SystemExit(2)
if original_with_genomes == 0:
    print("status\tbad\treason\tno_original_genomes_names")
    raise SystemExit(2)
if final_prefixes and frac < min_frac:
    print("status\tbad\treason\tlow_match_fraction")
    raise SystemExit(2)
print("status\tok")
PY
}

select_or_link_mod_folder() {
  if [[ -s "${expected_mod}/Bin_name_mod.txt" ]]; then
    log "Expected mod folder exists: ${expected_mod}"
    validate_mod_candidate "${expected_mod}"
    return 0
  fi

  local candidates=()
  local unpolished="${expected_mod/_MAGs_polished/}"
  candidates+=("${unpolished}")
  candidates+=("BestBinset_outlier_refined_mod")
  candidates+=("BestBinset_outlier_refined_filtrated_retrieved_mod")

  local candidate
  for candidate in "${candidates[@]}"; do
    [[ -s "${candidate}/Bin_name_mod.txt" ]] || continue
    log "Testing mod candidate: ${candidate}"
    if validate_mod_candidate "${candidate}"; then
      if [[ "${mode}" == "run" ]]; then
        if [[ -e "${expected_mod}" ]]; then
          die "Expected mod path exists but lacks Bin_name_mod.txt: ${expected_mod}; move/inspect it manually"
        fi
        log "Linking expected mod folder ${expected_mod} -> ${candidate}"
        ln -s "${candidate}" "${expected_mod}"
      else
        log "Preflight would link ${expected_mod} -> ${candidate} during run mode"
      fi
      return 0
    fi
  done

  die "No safe Bin_name_mod candidate found for expected folder ${expected_mod}"
}

preflight() {
  local n_fastas
  n_fastas="$(count_final_fastas)"
  [[ "${n_fastas}" -gt 0 ]] || die "No non-empty FASTA bins in ${FINAL}"
  log "Final FASTAs: ${n_fastas}"

  if ! quality_files | grep -q 'quality_report.tsv'; then
    die "No non-empty *quality_report.tsv in ${FINAL}; do not run final_drep until CheckM2 quality is present"
  fi
  log "Quality report(s) found in ${FINAL}:"
  quality_files >&2

  select_or_link_mod_folder
  log "Preflight OK"
}

run_final_drep() {
  preflight

  local stamp snap stdout stderr
  stamp="$(date +%Y%m%d_%H%M%S)"
  snap="${RUN}/final_drep_checkm2_snapshot_${stamp}"
  stdout="${RUN}/stdout.manual_final_drep_checkm2_${stamp}.log"
  stderr="${RUN}/stderr.manual_final_drep_checkm2_${stamp}.log"

  mkdir -p "${snap}"
  cp -a Basalt_checkpoint.txt Hybrid_re-assembly_status.txt Reassembled_bins_comparison.txt BASALT_command.txt run_wanted_only_initial_drep.command.sh "${snap}/" 2>/dev/null || true
  if [[ "${SNAPSHOT}" == "1" ]]; then
    log "Creating hardlink snapshot: ${snap}/${FINAL}.hardlink_snapshot"
    cp -al "${FINAL}" "${snap}/${FINAL}.hardlink_snapshot" || die "Hardlink snapshot failed; free inode/permission issue? Set SNAPSHOT=0 only if you accept the risk."
  fi

  log "Running CheckM2-aware BASALT final_drep directly"
  PYTHONPATH="${BASALT_DIR}" FINAL="${FINAL}" THREADS="${THREADS}" "${PYTHON_BIN}" - <<'PY' >"${stdout}" 2>"${stderr}"
import os
from pathlib import Path
from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_10242023 import final_binset_comparitor

pwd = os.getcwd()
final_folder = os.environ["FINAL"]
threads = int(os.environ.get("THREADS", "32"))

manifest = Path("metadata/dogfirst_full_srs_read_pairs.tsv")
if not manifest.exists():
    manifest = Path("/path/to/DogMAG_workdir/basalt_run_dogfirst_full_20260518/metadata/dogfirst_full_srs_read_pairs.tsv")

datasets = {}
if manifest.exists():
    with manifest.open() as fh:
        next(fh, None)
        for i, line in enumerate(fh, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                datasets[str(i)] = [Path(parts[3]).name, Path(parts[4]).name]

print(f"[INFO] final_folder={final_folder}", flush=True)
print(f"[INFO] datasets={len(datasets)}", flush=True)
print(f"[INFO] threads={threads}", flush=True)

final_binset_comparitor(final_folder, [], datasets, threads, pwd, "final_drep")
print("[DONE] direct CheckM2-aware BASALT final_drep finished", flush=True)
PY

  log "stdout=${stdout}"
  log "stderr=${stderr}"

  if grep -q '\[DONE\] direct CheckM2-aware BASALT final_drep finished' "${stdout}"; then
    grep -q '12th final de-replication done!' Basalt_checkpoint.txt || echo '12th final de-replication done!' >> Basalt_checkpoint.txt
    log "Final dRep completed"
  else
    die "Final dRep did not print DONE; inspect ${stdout} and ${stderr}"
  fi

  printf 'post_final_fastas\t%s\n' "$(count_final_fastas)"
  tail -40 "${stdout}" || true
  tail -80 "${stderr}" || true
}

case "${mode}" in
  inventory)
    print_inventory
    ;;
  preflight)
    print_inventory
    preflight
    ;;
  run)
    run_final_drep
    ;;
  *)
    usage
    die "Unknown mode: ${mode}"
    ;;
esac
