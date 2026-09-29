#!/usr/bin/env bash
# setup_env.sh
# Creates the Python environment for the AE experiment.
#
# Usage:
#   chmod +x setup_env.sh
#   ./setup_env.sh venv  env_tensorflow_ae    # creates venv folder env_tensorflow_ae/
#   ./setup_env.sh conda env_tensorflow_ae    # creates conda env named env_tensorflow_ae
#
# Arguments:
#   $1  MODE      : 'venv' or 'conda'  (default: venv)
#   $2  ENV_NAME  : name of the environment  (default: env_tensorflow_ae)

set -e   # stop immediately if any command fails

MODE=${1:-venv}
ENV_NAME=${2:-env_tensorflow_ae}

# Always work relative to the folder where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -------------------------------------------------------------------------
#  Helper: install all required Python packages into the active environment
# -------------------------------------------------------------------------
install_packages() {
    pip install --upgrade pip

    ARCH=$(uname -m)
    if [ "$ARCH" = "arm64" ]; then
        echo "Apple Silicon (arm64) detected — installing tensorflow-macos + tensorflow-metal"
        pip install tensorflow-macos tensorflow-metal
    else
        echo "Intel/x86_64 detected — installing tensorflow"
        pip install tensorflow
    fi

    pip install numpy scipy scikit-learn
}

# -------------------------------------------------------------------------
#  venv mode
# -------------------------------------------------------------------------
if [ "$MODE" = "venv" ]; then
    echo "=== Creating venv: ${SCRIPT_DIR}/${ENV_NAME} ==="
    cd "$SCRIPT_DIR"
    python3 -m venv "$ENV_NAME"
    source "${ENV_NAME}/bin/activate"
    install_packages
    deactivate
    echo "=== Done. Python at: ${SCRIPT_DIR}/${ENV_NAME}/bin/python ==="

# -------------------------------------------------------------------------
#  conda mode
# -------------------------------------------------------------------------
elif [ "$MODE" = "conda" ]; then
    echo "=== Creating conda env: ${ENV_NAME} ==="
    conda create -y -n "$ENV_NAME" python=3.10
    conda run -n "$ENV_NAME" pip install --upgrade pip

    ARCH=$(uname -m)
    if [ "$ARCH" = "arm64" ]; then
        echo "Apple Silicon (arm64) detected — installing tensorflow-macos + tensorflow-metal"
        conda run -n "$ENV_NAME" pip install tensorflow-macos tensorflow-metal
    else
        echo "Intel/x86_64 detected — installing tensorflow"
        conda run -n "$ENV_NAME" pip install tensorflow
    fi

    conda run -n "$ENV_NAME" pip install numpy scipy scikit-learn
    echo "=== Done. Conda env ready: ${ENV_NAME} ==="

else
    echo "ERROR: Unknown mode '${MODE}'. Use 'venv' or 'conda'."
    exit 1
fi
