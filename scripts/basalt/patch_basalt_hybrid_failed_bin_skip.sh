#!/usr/bin/env bash
set -euo pipefail

# Patch BASALT hybrid reassembly so known failed bins can be skipped using
# a run-local skip list. This complements the completed-output resume guard.
#
# Usage:
#   bash scripts/patch_basalt_hybrid_failed_bin_skip.sh
#
# Optional:
#   BASALT_S9P=/path/to/S9p_Hybrid_Reassembly_10262023.py bash scripts/patch_basalt_hybrid_failed_bin_skip.sh
#
# Runtime skip list path:
#   BASALT_HYBRID_FAILED_SKIP_BINS=/path/to/list.txt bash run_wanted_only_initial_drep.command.sh
# If the env var is omitted, BASALT looks for:
#   <current BASALT run>/dogfirst_hybrid_reassembly_skip_bins.txt

TARGET="${BASALT_S9P:-/path/to/BASALT/S9p_Hybrid_Reassembly_10262023.py}"
MARKER="BASALT_DOGFIRST_FAILED_GUARD"

if [[ ! -f "$TARGET" ]]; then
  echo "[ERROR] Target BASALT script not found: $TARGET" >&2
  exit 1
fi

if grep -q "$MARKER" "$TARGET"; then
  echo "[OK] Failed-bin skip guard already present in: $TARGET"
  exit 0
fi

backup="${TARGET}.before_dogfirst_failed_skip_guard.$(date +%Y%m%d_%H%M%S).bak"
cp -a "$TARGET" "$backup"
echo "[INFO] Backup written: $backup"

python - "$TARGET" <<'PY'
from pathlib import Path
import sys

target = Path(sys.argv[1])
text = target.read_text()
marker = "BASALT_DOGFIRST_FAILED_GUARD"

if marker in text:
    print(f"[OK] Failed-bin skip guard already present in: {target}")
    raise SystemExit(0)

func = "def hybrid_assembly_mul(sr_folder, bin_seq, item, bin_lr_reads, lr_folder,"
gzip_block = "    try:\n        os.system('gzip -d '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')"
completed_guard_tail = (
    "    if os.path.isfile(existing_hybrid_out) and os.path.getsize(existing_hybrid_out) > 0:\n"
    "        print('BASALT_DOGFIRST_RESUME_GUARD skip existing hybrid reassembly '+str(item)+': '+existing_hybrid_out)\n"
    "        return\n\n"
)

try:
    func_pos = text.index(func)
except ValueError:
    raise SystemExit(f"[ERROR] Could not find hybrid_assembly_mul() in {target}")

try:
    insert_pos = text.index(completed_guard_tail, func_pos) + len(completed_guard_tail)
except ValueError:
    try:
        insert_pos = text.index(gzip_block, func_pos)
    except ValueError:
        raise SystemExit(f"[ERROR] Could not find expected insertion point inside hybrid_assembly_mul() in {target}")

guard = (
    "    failed_skip_file = os.environ.get('BASALT_HYBRID_FAILED_SKIP_BINS', pwd+'/dogfirst_hybrid_reassembly_skip_bins.txt')\n"
    "    try:\n"
    "        with open(failed_skip_file, 'r') as _skip_handle:\n"
    "            failed_skip_bins = set([line.strip().split()[0] for line in _skip_handle if line.strip() and not line.lstrip().startswith('#')])\n"
    "    except:\n"
    "        failed_skip_bins = set()\n"
    "    if str(item) in failed_skip_bins:\n"
    "        print('BASALT_DOGFIRST_FAILED_GUARD skip known failed hybrid reassembly '+str(item)+': '+failed_skip_file)\n"
    "        return\n"
    "\n"
)

target.write_text(text[:insert_pos] + guard + text[insert_pos:])
print(f"[OK] Patched failed-bin skip guard into: {target}")
PY

python -m py_compile "$TARGET"
echo "[OK] Python syntax check passed."
echo "[DONE] BASALT will now skip bins listed in dogfirst_hybrid_reassembly_skip_bins.txt."
