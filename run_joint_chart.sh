#!/usr/bin/env bash
# Run the joint 2-D chart (K x K training) experiment headless (see experiments/Main_AoA_Serving_JointChart.m).
#   SJC_WORKERS=40 ./run_joint_chart.sh
# Log: logs/joint_chart.log
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
if [ -z "${MATLAB_BIN:-}" ]; then
    MATLAB_BIN=$(command -v matlab || echo "$HOME/matlab/bin/matlab")
fi
"$MATLAB_BIN" -batch "addpath('experiments'); Main_AoA_Serving_JointChart" > logs/joint_chart.log 2>&1
