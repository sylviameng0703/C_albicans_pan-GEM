#!/usr/bin/env bash
set -euo pipefail

BASE=${BASE:-$PWD}
RUN=${RUN:-$BASE/results/pangenome_relaxed_80}
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3)}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

MODEL=${MODEL:-$RUN/342_final_panGEM_detail_closed_release/Calb_panGEM_final_release.xml}
PROJ_OUT=${PROJ_OUT:-$RUN/344_final_panGEM_detail_closed_reaction_projection_for_ssGEM}
OUT=${OUT:-$RUN/345_regenerate_ssGEMs_from_final_detail_closed_panGEM}
SCRIPT=${SCRIPT:-$SCRIPT_DIR/build_ssgems.py}
REACTION_PA=${REACTION_PA:-$PROJ_OUT/locked_sample_reaction_presence_absence.tsv}
REACTION_CATALOG=${REACTION_CATALOG:-$PROJ_OUT/locked_reaction_projection_catalog.tsv}

for f in "$MODEL" "$REACTION_PA" "$REACTION_CATALOG" "$SCRIPT"; do
  [[ -s "$f" ]] || { echo "[ERROR] Missing required input: $f" >&2; exit 1; }
done
mkdir -p "$OUT" "$RUN/logs"

"$PYTHON_BIN" "$SCRIPT" \
  --model "$MODEL" \
  --reaction-presence-absence "$REACTION_PA" \
  --reaction-catalog "$REACTION_CATALOG" \
  --outdir "$OUT" \
  --samples all \
  --removable-categories "shell_gpr;absent_gpr_projection" \
  --rescue-categories "shell_gpr" \
  --atpm-id calb_ATPM__c \
  --model-label "Calb"

"$PYTHON_BIN" - "$OUT" "$MODEL" <<'PY'
import csv
import hashlib
import sys
from pathlib import Path

outdir, model = map(Path, sys.argv[1:])
summary_path = outdir / "qc" / "v18_afumigatus_style_generation_summary.tsv"
with summary_path.open(newline="") as handle:
    summary = {row["metric"]: row.get("value", "") for row in csv.DictReader(handle, delimiter="\t")}

required = {
    "samples_written": 80,
    "final_growth_samples_after_rescue": 80,
    "final_zero_or_failed_growth_samples": 0,
    "egc_pass_samples": 80,
}
failures = []
for metric, expected in required.items():
    try:
        observed = int(float(summary.get(metric, "")))
    except ValueError:
        observed = None
    if observed != expected:
        failures.append(f"{metric}={observed}, expected {expected}")

if len(list((outdir / "models").glob("*.xml"))) != 80:
    failures.append("model_count is not 80")
if summary.get("input_model_sha256") != hashlib.sha256(model.read_bytes()).hexdigest():
    failures.append("input_model_sha256 does not match the frozen 342 pan-GEM")
if failures:
    raise SystemExit("[ERROR] final ssGEM release gate failed: " + "; ".join(failures))
print("[PASS] final 342-derived ssGEM gate: 80 models, 80 growth, 80 ATPM-cycle QC passes")
PY

echo "[DONE] Final 342-derived ssGEM generation"
echo "  models:  $OUT/models"
echo "  summary: $OUT/qc/v18_afumigatus_style_generation_summary.tsv"
