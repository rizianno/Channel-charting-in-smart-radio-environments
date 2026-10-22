#!/usr/bin/env bash
# Run the EMS beam selection vs sweeps experiment headless (see experiments/Main_AoA_Serving_Sweeps.m).
#   SRV_WORKERS=48 ./run_serving.sh
# Log: logs/serving_sweeps.log
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Serving_Sweeps" > logs/serving_sweeps.log 2>&1
