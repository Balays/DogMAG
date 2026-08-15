#!/usr/bin/env bash
set -euo pipefail

# Resume BASALT step 7 after the long-read mapping/splitting loop stopped
# before writing "Long-read mapping done!". This keeps completed split outputs
# from being duplicated and reuses an existing lrN.sam if one was produced.

RUN="${RUN:-/path/to/DogMAG_workdir/basalt_run_dogfirst_full_20260518_wanted_only_copybins_checkm2_padidx_20260608}"
BASALT_S7="${BASALT_S7:-/path/to/BASALT/S7lr_finding_sr_contigs_basing_lr_and_polishing_11022023.py}"
MIN_FREE_GB="${MIN_FREE_GB:-150}"
RESUME_STDOUT="${RESUME_STDOUT:-stdout.resume24_step7_partial_continue.log}"
RESUME_STDERR="${RESUME_STDERR:-stderr.resume24_step7_partial_continue.log}"
MODE="${1:-audit}"

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

require_run() {
  [[ -d "${RUN}" ]] || die "RUN does not exist: ${RUN}"
  [[ -s "${RUN}/run_wanted_only_initial_drep.command.sh" ]] || die "Missing run command in ${RUN}"
}

free_gb_for_run() {
  df -BG "${RUN}" | awk 'NR==2 {gsub(/G/,"",$4); print $4}'
}

audit() {
  require_run
  cd "${RUN}"

  log "RUN=${RUN}"
  log "BASALT_S7=${BASALT_S7}"
  log "Free GB on RUN filesystem: $(free_gb_for_run)"

  log "Checkpoint:"
  tail -20 Basalt_checkpoint.txt 2>/dev/null || true

  log "Polishing status:"
  ls -lh 1_polishing_status.txt 2>/dev/null || true
  cat 1_polishing_status.txt 2>/dev/null || true

  local mapped split_a split_b lrfq sams
  mapped="$(grep -c 'Real time:' stderr.resume23_gzipfastqfix.log 2>/dev/null || true)"
  split_a="$(find . -maxdepth 1 -type f -name 'Bin_long_read*.txt' | wc -l)"
  split_b="$(find . -maxdepth 1 -type f -name 'Long_read_bin*.txt' | wc -l)"
  lrfq="$(find . -maxdepth 1 -type f -name '*_lr.fq' | wc -l)"
  sams="$(find . -maxdepth 1 -type f -name 'lr*.sam' | wc -l)"

  log "Completed minimap2 mappings in resume23 stderr: ${mapped}"
  log "Bin_long_read*.txt: ${split_a}"
  log "Long_read_bin*.txt: ${split_b}"
  log "*_lr.fq files: ${lrfq}"
  log "lr*.sam files: ${sams}"

  log "Last split outputs:"
  find . -maxdepth 1 -type f -name 'Bin_long_read*.txt' -printf '%f\t%s\t%TY-%Tm-%Td %TH:%TM\n' | sort -V | tail -10 || true
  find . -maxdepth 1 -type f -name 'Long_read_bin*.txt' -printf '%f\t%s\t%TY-%Tm-%Td %TH:%TM\n' | sort -V | tail -10 || true

  log "SAM files:"
  find . -maxdepth 1 -type f -name 'lr*.sam' -printf '%f\t%s\t%TY-%Tm-%Td %TH:%TM\n' | sort -V || true

  if [[ "${split_a}" -gt 0 && "${split_a}" -eq "${split_b}" ]]; then
    log "Audit OK: completed split outputs are paired."
  else
    die "Split output counts differ; inspect before resume."
  fi
}

patch_s7() {
  [[ -s "${BASALT_S7}" ]] || die "Missing BASALT S7 script: ${BASALT_S7}"

  if grep -q 'Skipping long reads .* split outputs already exist' "${BASALT_S7}"; then
    log "BASALT S7 already has partial long-read resume patch."
    return 0
  fi

  local backup="${BASALT_S7}.before_partial_lr_resume.$(date +%Y%m%d_%H%M%S).bak"
  cp -p "${BASALT_S7}" "${backup}"
  log "Backed up BASALT S7 script: ${backup}"

  python3 - "$BASALT_S7" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text()

old = """                if x == 0:
                    n, bin_lr=0, {}
                    for lrs in long_read:
                        n+=1
                        os.system('minimap2 -t '+str(num_threads)+' -ax map-ont Total_bins.fa '+str(lrs)+' > lr'+str(n)+'.sam')
                        print('Splitting long reads '+str(n))
                        bin_lr.update(parse_lr_sam('lr'+str(n)+'.sam', lrs, n))
                        os.system('rm lr'+str(n)+'.sam')
"""

new = """                if x == 0:
                    n, bin_lr=0, {}
                    # Resume-safe path: keep already split long reads and avoid
                    # appending them into *_lr.fq a second time after an interrupted run.
                    for existing_lr_fq in os.listdir('.'):
                        if existing_lr_fq.endswith('_lr.fq'):
                            bin_lr[existing_lr_fq[:-6]]={}
                    for lrs in long_read:
                        n+=1
                        bin_long_file='Bin_long_read'+str(n)+'.txt'
                        long_bin_file='Long_read_bin'+str(n)+'.txt'
                        sam_file='lr'+str(n)+'.sam'
                        if os.path.exists(bin_long_file) and os.path.getsize(bin_long_file) > 0 and os.path.exists(long_bin_file) and os.path.getsize(long_bin_file) > 0:
                            print('Skipping long reads '+str(n)+' because split outputs already exist')
                            continue
                        if os.path.exists(sam_file) and os.path.getsize(sam_file) > 0:
                            print('Reusing existing '+str(sam_file)+' for long reads '+str(n))
                        else:
                            os.system('minimap2 -t '+str(num_threads)+' -ax map-ont Total_bins.fa '+str(lrs)+' > '+str(sam_file))
                        print('Splitting long reads '+str(n))
                        bin_lr.update(parse_lr_sam(sam_file, lrs, n))
                        os.system('rm '+str(sam_file))
"""

if old not in text:
    raise SystemExit("Could not find the expected long-read mapping loop. Refusing to patch.")

path.write_text(text.replace(old, new, 1))
PY

  log "Patched BASALT S7 for partial long-read mapping resume."
}

pre_resume_checks() {
  require_run
  cd "${RUN}"

  local free_gb
  free_gb="$(free_gb_for_run)"
  if [[ "${free_gb}" -lt "${MIN_FREE_GB}" ]]; then
    die "Only ${free_gb} GB free on RUN filesystem; require at least ${MIN_FREE_GB} GB."
  fi

  local active_processes
  active_processes="$(pgrep -af 'BASALT.py|minimap2|samtools|checkm2|pilon' | grep -v -E 'pgrep|resume_basalt_step7_partial_lr_mapping' || true)"
  if [[ -n "${active_processes}" ]]; then
    log "Active BASALT-related processes:"
    printf '%s\n' "${active_processes}"
    die "Stop/verify active processes before restarting."
  fi

  [[ -s "${BASALT_S7}" ]] || die "Missing BASALT S7 script: ${BASALT_S7}"
  grep -q 'Skipping long reads .* split outputs already exist' "${BASALT_S7}" || die "BASALT S7 is not patched yet; run: $0 patch"
}

resume_run() {
  pre_resume_checks
  cd "${RUN}"
  log "Starting BASALT resume. Logs: ${RUN}/${RESUME_STDOUT}, ${RUN}/${RESUME_STDERR}"
  bash run_wanted_only_initial_drep.command.sh >> "${RESUME_STDOUT}" 2>> "${RESUME_STDERR}"
}

case "${MODE}" in
  audit)
    audit
    ;;
  patch)
    patch_s7
    ;;
  preflight)
    audit
    patch_s7
    pre_resume_checks
    log "Preflight OK. Resume with: RUN='${RUN}' bash $0 resume"
    ;;
  resume)
    resume_run
    ;;
  all)
    audit
    patch_s7
    pre_resume_checks
    resume_run
    ;;
  *)
    die "Unknown mode: ${MODE}. Use audit, patch, preflight, resume, or all."
    ;;
esac
