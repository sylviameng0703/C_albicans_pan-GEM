#!/usr/bin/env bash
set -euo pipefail

BASE=${BASE:-$PWD}
RUN=${RUN:-$BASE/results/pangenome_relaxed_80}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

MODEL=${MODEL:-$RUN/342_final_panGEM_detail_closed_release/Calb_panGEM_final_release.xml}
RELEASE_OVERVIEW=${RELEASE_OVERVIEW:-$RUN/342_final_panGEM_detail_closed_release/final_panGEM_detail_closed_release_overview.tsv}
OUT=${OUT:-$RUN/344_final_panGEM_detail_closed_reaction_projection_for_ssGEM}

MODEL="$MODEL" \
RELEASE_OVERVIEW="$RELEASE_OVERVIEW" \
OUT="$OUT" \
EXPECTED_RELEASE_GATE=pass_final_detail_closed_panGEM_release \
bash "$SCRIPT_DIR/project_reactions.sh"
