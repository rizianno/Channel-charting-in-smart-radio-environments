#!/usr/bin/env bash
# Run the AoA ablation study headless (see experiments/Main_AoA_Ablations.m).
#   ABL_STAGE=calibration ABL_WORKERS=32 ./run_ablations.sh
#   ABL_STAGE=matrix ABL_NOISE=est_G30dB ABL_WORKERS=48 ./run_ablations.sh
# Log: logs/abl_<stage>[_<noise>].log
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
TAG="${ABL_STAGE:-calibration}${ABL_NOISE:+_$ABL_NOISE}"
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Ablations" > "logs/abl_${TAG}.log" 2>&1
