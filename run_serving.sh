#!/usr/bin/env bash
# Run the AoA serving experiment headless (see experiments/Main_AoA_Serving.m).
#   SRV_WORKERS=48 ./run_serving.sh
# Log: logs/serving.log
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Serving" > logs/serving.log 2>&1
