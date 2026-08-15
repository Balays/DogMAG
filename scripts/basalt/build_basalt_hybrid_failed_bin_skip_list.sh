#!/usr/bin/env bash
set -euo pipefail

# Build a skip list for BASALT hybrid SPAdes bins that already failed to
# produce contigs.fasta. This is intentionally conservative: it only skips
# bins with explicit "mv: cannot stat .../bin*_spades_hybrid_reassembly/contigs.fasta"
# evidence, and removes bins that already have a final non-empty hybrid FASTA.
#
# Usage:
#   RUN=/path/to/basalt_run bash scripts/build_basalt_hybrid_failed_bin_skip_list.sh

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

OUT="${OUT:-$RUN/dogfirst_hybrid_reassembly_skip_bins.txt}"
REPORT="${REPORT:-$RUN/dogfirst_hybrid_reassembly_skip_bins.report.tsv}"

python - "$RUN" "$OUT" "$REPORT" <<'PY'
from pathlib import Path
import re
import sys

run = Path(sys.argv[1])
out = Path(sys.argv[2])
report = Path(sys.argv[3])

failed = {}
pattern = re.compile(r"(bin\d+)_spades_hybrid_reassembly/contigs\.fasta")

for log in sorted(run.glob("stderr.resume*.log")) + sorted(run.glob("stderr*.log")):
    if not log.is_file():
        continue
    try:
        with log.open("r", errors="replace") as handle:
            for lineno, line in enumerate(handle, 1):
                if "cannot stat" not in line or "spades_hybrid_reassembly/contigs.fasta" not in line:
                    continue
                match = pattern.search(line)
                if not match:
                    continue
                failed.setdefault(match.group(1), []).append((log.name, lineno))
    except OSError:
        continue

final_dir = run / "BestBinset_outlier_refined_MAGs_polished_re-assembly_binset"
skippable = []
kept_completed = []

for bin_id in sorted(failed, key=lambda x: int(x[3:]) if x[3:].isdigit() else x):
    final_fa = final_dir / f"{bin_id}_SPAdes_hybrid_re-assembly_contigs.fa"
    if final_fa.is_file() and final_fa.stat().st_size > 0:
        kept_completed.append(bin_id)
    else:
        skippable.append(bin_id)

out.write_text("\n".join(skippable) + ("\n" if skippable else ""))

with report.open("w") as rep:
    rep.write("bin_id\tstatus\tevidence_count\tfirst_evidence\n")
    for bin_id in sorted(failed, key=lambda x: int(x[3:]) if x[3:].isdigit() else x):
        status = "already_completed_not_skipped" if bin_id in kept_completed else "skip_known_failed"
        first = failed[bin_id][0]
        rep.write(f"{bin_id}\t{status}\t{len(failed[bin_id])}\t{first[0]}:{first[1]}\n")

print(f"[INFO] explicit failed bins found: {len(failed)}")
print(f"[INFO] bins already completed and not added: {len(kept_completed)}")
print(f"[INFO] bins written to skip list: {len(skippable)}")
print(f"[INFO] skip list: {out}")
print(f"[INFO] report: {report}")
PY

