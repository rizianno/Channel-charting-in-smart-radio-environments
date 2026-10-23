#!/usr/bin/env bash
# Run the single-EMS deployment experiment headless (see experiments/Main_AoA_Serving_SingleEMS.m).
#   SGL_WORKERS=48 ./run_single_ems.sh
# Log: logs/single_ems.log
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Serving_SingleEMS" > logs/single_ems.log 2>&1
