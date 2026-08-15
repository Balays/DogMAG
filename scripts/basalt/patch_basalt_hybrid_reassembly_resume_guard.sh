#!/usr/bin/env bash
set -euo pipefail

# Patch BASALT hybrid reassembly so interrupted restarts do not rerun bins
# whose final hybrid SPAdes FASTA already exists.
#
# Usage:
#   bash scripts/patch_basalt_hybrid_reassembly_resume_guard.sh
#
# Optional:
#   BASALT_S9P=/path/to/S9p_Hybrid_Reassembly_10262023.py bash scripts/patch_basalt_hybrid_reassembly_resume_guard.sh

TARGET="${BASALT_S9P:-/path/to/BASALT/S9p_Hybrid_Reassembly_10262023.py}"
MARKER="BASALT_DOGFIRST_RESUME_GUARD"

if [[ ! -f "$TARGET" ]]; then
  echo "[ERROR] Target BASALT script not found: $TARGET" >&2
  exit 1
fi

if grep -q "$MARKER" "$TARGET"; then
  echo "[OK] Resume guard already present in: $TARGET"
  exit 0
fi

backup="${TARGET}.before_dogfirst_hybrid_resume_guard.$(date +%Y%m%d_%H%M%S).bak"
cp -a "$TARGET" "$backup"
echo "[INFO] Backup written: $backup"

python - "$TARGET" <<'PY'
from pathlib import Path
import sys

target = Path(sys.argv[1])
text = target.read_text()
marker = "BASALT_DOGFIRST_RESUME_GUARD"

if marker in text:
    print(f"[OK] Resume guard already present in: {target}")
    raise SystemExit(0)

func = "def hybrid_assembly_mul(sr_folder, bin_seq, item, bin_lr_reads, lr_folder,"
try_line = "    try:\n        os.system('gzip -d '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')"

try:
    func_pos = text.index(func)
except ValueError:
    raise SystemExit(f"[ERROR] Could not find hybrid_assembly_mul() in {target}")

try:
    insert_pos = text.index(try_line, func_pos)
except ValueError:
    raise SystemExit(f"[ERROR] Could not find expected gzip block inside hybrid_assembly_mul() in {target}")

guard = (
    "    existing_hybrid_out = pwd+'/'+reassembly_bin_folder+'/'+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa'\n"
    "    if os.path.isfile(existing_hybrid_out) and os.path.getsize(existing_hybrid_out) > 0:\n"
    "        print('BASALT_DOGFIRST_RESUME_GUARD skip existing hybrid reassembly '+str(item)+': '+existing_hybrid_out)\n"
    "        return\n"
    "\n"
)

target.write_text(text[:insert_pos] + guard + text[insert_pos:])
print(f"[OK] Patched hybrid reassembly resume guard into: {target}")
PY

python -m py_compile "$TARGET"
echo "[OK] Python syntax check passed."
echo "[DONE] BASALT will now skip bins with existing non-empty *_SPAdes_hybrid_re-assembly_contigs.fa outputs."
