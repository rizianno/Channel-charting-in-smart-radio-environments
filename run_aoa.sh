#!/usr/bin/env bash
# Run the AoA experiment headless (e.g. on the lab server).
# Usage examples:
#   PHI_OPTIONS_DEG="0:20:40" RUN_TAG=coarse_3x3 ./run_aoa.sh
#   PHI_OPTIONS_DEG="0:4:40"  RUN_TAG=full_11x11 ./run_aoa.sh
# Outputs: results_aoa/<RUN_TAG>/ (csv, png, notes.txt; results.mat stays local)
# Log:     logs/<RUN_TAG>.log
set -euo pipefail
cd "$(dirname "$0")"
: "${RUN_TAG:=aoa_$(date +%Y%m%d_%H%M)}"
export RUN_TAG
mkdir -p logs results_aoa
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Estimation" 2>&1 | tee "logs/${RUN_TAG}.log"
