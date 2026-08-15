#!/usr/bin/env bash
set -euo pipefail

# Stop only BASALT/SPAdes processes whose cwd or command line points at RUN.
# This avoids printing the enormous BASALT command line to the terminal.
#
# Usage:
#   RUN=/path/to/basalt_run bash scripts/stop_basalt_run_processes.sh

RUN="${RUN:-}"
if [[ -z "$RUN" ]]; then
  echo "[ERROR] Set RUN=/path/to/basalt_run" >&2
  exit 1
fi

RUN="$(readlink -f "$RUN")"
if [[ ! -d "$RUN" ]]; then
  echo "[ERROR] RUN directory not found: $RUN" >&2
  exit 1
fi

SNAP="${SNAP:-$RUN/restart_guard_stop_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$SNAP"

pid_file="$SNAP/run_processes.pid"
details_file="$SNAP/run_processes.details.tsv"
: > "$pid_file"
: > "$details_file"

for pid in $(pgrep -f 'BASALT.py|spades.py|spades-core' || true); do
  cwd="$(readlink -f "/proc/$pid/cwd" 2>/dev/null || true)"
  cmd="$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null || true)"
  if [[ "$cwd" == "$RUN"* || "$cmd" == *"$RUN"* ]]; then
    printf '%s\n' "$pid" >> "$pid_file"
    printf '%s\t%s\t%s\n' "$pid" "$cwd" "$cmd" >> "$details_file"
  fi
done

count="$(wc -l < "$pid_file" | tr -d ' ')"
echo "[INFO] Matched run-related BASALT/SPAdes processes: $count"
echo "[INFO] Details: $details_file"

if [[ "$count" -eq 0 ]]; then
  echo "[OK] Nothing to stop."
  exit 0
fi

echo "[INFO] Sending TERM..."
xargs -r kill < "$pid_file"
sleep "${TERM_WAIT_SECONDS:-20}"

still_file="$SNAP/run_processes_after_term.pid"
: > "$still_file"
while read -r pid; do
  if [[ -n "$pid" && -d "/proc/$pid" ]]; then
    printf '%s\n' "$pid" >> "$still_file"
  fi
done < "$pid_file"

still="$(wc -l < "$still_file" | tr -d ' ')"
if [[ "$still" -gt 0 ]]; then
  echo "[WARN] Still alive after TERM: $still; sending KILL..."
  xargs -r kill -9 < "$still_file"
  sleep 5
fi

final_file="$SNAP/run_processes_after_kill.pid"
: > "$final_file"
while read -r pid; do
  if [[ -n "$pid" && -d "/proc/$pid" ]]; then
    printf '%s\n' "$pid" >> "$final_file"
  fi
done < "$pid_file"

final="$(wc -l < "$final_file" | tr -d ' ')"
if [[ "$final" -gt 0 ]]; then
  echo "[WARN] Some processes are still present, likely uninterruptible disk sleep:"
  while read -r pid; do
    ps -p "$pid" -o pid,ppid,stat,etime,wchan:32,comm || true
  done < "$final_file"
  echo "[WARN] If STAT contains D, only WSL restart usually clears them."
  exit 2
fi

echo "[OK] Run-related BASALT/SPAdes processes stopped."
