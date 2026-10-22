#!/usr/bin/env bash
# Run the full-NLOS surface-size experiment headless (see experiments/Main_AoA_Serving_NLoS_Size.m).
#   SRV_WORKERS=48 ./run_serving.sh
# Log: logs/serving_nlos_size.log
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Serving_NLoS_Size" > logs/serving_nlos_size.log 2>&1
