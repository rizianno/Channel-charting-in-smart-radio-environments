#!/usr/bin/env bash
# Run the single-EMS deployment experiment headless (see experiments/Main_AoA_Serving_Deployments.m).
#   SDP_WORKERS=48 ./run_deployments.sh
# Log: logs/deployments.log
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Serving_Deployments" > logs/deployments.log 2>&1
