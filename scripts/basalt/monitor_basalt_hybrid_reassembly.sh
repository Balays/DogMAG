#!/usr/bin/env bash
set -euo pipefail

# Monitor BASALT step-7/8 hybrid SPAdes reassembly progress.
#
# Usage:
#   RUN=/path/to/basalt_run bash scripts/monitor_basalt_hybrid_reassembly.sh
#
# Optional:
#   RECENT_MIN=120 RUN=/path/to/basalt_run bash scripts/monitor_basalt_hybrid_reassembly.sh

RUN="${RUN:-}"
if [[ -z "$RUN" ]]; then
  echo "[ERROR] Set RUN=/path/to/basalt_run" >&2
  exit 1
fi

RUN="$(readlink -f "$RUN")"
RECENT_MIN="${RECENT_MIN:-120}"

if [[ ! -d "$RUN" ]]; then
  echo "[ERROR] RUN directory not found: $RUN" >&2
  exit 1
fi

python - "$RUN" "$RECENT_MIN" <<'PY'
from pathlib import Path
from datetime import datetime, timezone
import os
import re
import sys
import time

run = Path(sys.argv[1])
recent_min = float(sys.argv[2])
now = time.time()

long_dir = run / "BestBinset_long_read"
final_dir = run / "BestBinset_outlier_refined_MAGs_polished_re-assembly_binset"
skip_file = run / "dogfirst_hybrid_reassembly_skip_bins.txt"
skip_report = run / "dogfirst_hybrid_reassembly_skip_bins.report.tsv"

def bin_sort_key(bin_id):
    try:
        return int(bin_id[3:])
    except Exception:
        return bin_id

def count_lines(logs, needle):
    total = 0
    latest_name = ""
    latest_total = 0
    latest_mtime = -1
    for log in logs:
        try:
            mtime = log.stat().st_mtime
            with log.open("r", errors="replace") as handle:
                n = sum(1 for line in handle if needle in line)
        except OSError:
            continue
        total += n
        if mtime > latest_mtime:
            latest_mtime = mtime
            latest_name = log.name
            latest_total = n
    return total, latest_name, latest_total

candidate_bins = set()
if long_dir.is_dir():
    pat = re.compile(r"^(bin\d+)_lr\.fq(?:\.gz)?$")
    for item in long_dir.iterdir():
        if not item.is_file():
            continue
        match = pat.match(item.name)
        if match:
            candidate_bins.add(match.group(1))

success_bins = {}
if final_dir.is_dir():
    pat = re.compile(r"^(bin\d+)_SPAdes_hybrid_re-assembly_contigs\.fa$")
    for item in final_dir.iterdir():
        if not item.is_file():
            continue
        match = pat.match(item.name)
        if match and item.stat().st_size > 0:
            success_bins[match.group(1)] = item.stat()

stderr_logs = sorted(run.glob("stderr*.log")) + sorted(run.glob("stderr.resume*.log"))
stdout_logs = sorted(run.glob("stdout*.log")) + sorted(run.glob("stdout.resume*.log"))
stderr_logs = sorted(set(stderr_logs), key=lambda p: p.name)
stdout_logs = sorted(set(stdout_logs), key=lambda p: p.name)

explicit_failed = {}
fail_pat = re.compile(r"(bin\d+)_spades_hybrid_reassembly/contigs\.fasta")
for log in stderr_logs:
    try:
        with log.open("r", errors="replace") as handle:
            for lineno, line in enumerate(handle, 1):
                if "cannot stat" not in line or "spades_hybrid_reassembly/contigs.fasta" not in line:
                    continue
                match = fail_pat.search(line)
                if not match:
                    continue
                explicit_failed.setdefault(match.group(1), []).append((log.name, lineno))
    except OSError:
        continue

explicit_failed_no_success = set(explicit_failed) - set(success_bins)

skip_bins = set()
if skip_file.is_file():
    try:
        with skip_file.open("r", errors="replace") as handle:
            skip_bins = {line.strip().split()[0] for line in handle if line.strip() and not line.lstrip().startswith("#")}
    except OSError:
        pass

spades_dirs = []
recent_dirs = []
dir_pat = re.compile(r"^(bin\d+)_spades_hybrid_reassembly$")
for item in run.iterdir():
    if not item.is_dir():
        continue
    match = dir_pat.match(item.name)
    if not match:
        continue
    stat = item.stat()
    spades_dirs.append((match.group(1), stat))
    if now - stat.st_mtime <= recent_min * 60:
        recent_dirs.append((match.group(1), stat))

completed_skip_total, completed_skip_latest_log, completed_skip_latest = count_lines(
    stdout_logs, "BASALT_DOGFIRST_RESUME_GUARD skip existing hybrid reassembly"
)
failed_skip_total, failed_skip_latest_log, failed_skip_latest = count_lines(
    stdout_logs, "BASALT_DOGFIRST_FAILED_GUARD skip known failed hybrid reassembly"
)
spades_start_total, spades_start_latest_log, spades_start_latest = count_lines(
    stdout_logs, "======= SPAdes pipeline started"
)
spades_finish_total, spades_finish_latest_log, spades_finish_latest = count_lines(
    stdout_logs, "SPAdes pipeline finished"
)

unresolved = candidate_bins - set(success_bins) - explicit_failed_no_success
if skip_bins:
    unresolved_without_skip = unresolved - skip_bins
else:
    unresolved_without_skip = unresolved

print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] RUN={run}")
print(f"Candidates with LR reads:              {len(candidate_bins):>8}")
print(f"Successful final hybrid FASTAs:        {len(success_bins):>8}")
print(f"Explicit failed bins, no final FASTA:  {len(explicit_failed_no_success):>8}")
print(f"Failed-bin skip-list entries:          {len(skip_bins):>8}")
print(f"Unresolved candidates:                 {len(unresolved):>8}")
print(f"Unresolved not in skip-list:           {len(unresolved_without_skip):>8}")
print(f"Hybrid SPAdes work dirs:               {len(spades_dirs):>8}")
print(f"Hybrid SPAdes dirs modified <= {recent_min:g} min: {len(recent_dirs):>8}")
print()
print("Cumulative log counters:")
print(f"  completed-output skips: {completed_skip_total}")
print(f"  failed-bin skips:       {failed_skip_total}")
print(f"  SPAdes starts:          {spades_start_total}")
print(f"  SPAdes finished lines:  {spades_finish_total}")
print()
print("Latest stdout-log counters:")
if stdout_logs:
    print(f"  latest stdout log:      {max(stdout_logs, key=lambda p: p.stat().st_mtime).name}")
print(f"  completed-output skips: {completed_skip_latest} ({completed_skip_latest_log})")
print(f"  failed-bin skips:       {failed_skip_latest} ({failed_skip_latest_log})")
print(f"  SPAdes starts:          {spades_start_latest} ({spades_start_latest_log})")
print(f"  SPAdes finished lines:  {spades_finish_latest} ({spades_finish_latest_log})")
print()

latest_success = sorted(success_bins.items(), key=lambda kv: kv[1].st_mtime, reverse=True)[:8]
if latest_success:
    print("Latest successful final FASTAs:")
    for bin_id, stat in latest_success:
        ts = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M")
        print(f"  {ts}\t{stat.st_size}\t{bin_id}")
    print()

latest_fail = sorted(
    ((bin_id, evidence[0]) for bin_id, evidence in explicit_failed.items() if bin_id in explicit_failed_no_success),
    key=lambda x: (x[1][0], x[1][1]),
    reverse=True,
)[:8]
if latest_fail:
    print("Example explicit failed bins without final FASTA:")
    for bin_id, (log_name, line_no) in latest_fail:
        print(f"  {bin_id}\t{log_name}:{line_no}")
    print()

if skip_file.exists():
    print(f"Skip list:  {skip_file}")
if skip_report.exists():
    print(f"Report:     {skip_report}")
PY
