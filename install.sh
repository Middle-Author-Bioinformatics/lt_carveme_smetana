#!/usr/bin/env bash
set -euo pipefail

# Run once: bash install.sh [environment-name]
# Requires an existing conda/mamba installation; does not use sudo.
ENV_NAME="${1:-lt-metabolism}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if command -v mamba >/dev/null 2>&1; then
    MANAGER=mamba
elif command -v conda >/dev/null 2>&1; then
    MANAGER=conda
else
    echo 'Install Miniforge first, then rerun this installer: https://github.com/conda-forge/miniforge' >&2
    exit 1
fi
# create deliberately refuses to replace an existing named environment.
"$MANAGER" create -y -n "$ENV_NAME" --override-channels -c conda-forge -c bioconda \
    python=3.12 pip diamond=2.1.11
"$MANAGER" run -n "$ENV_NAME" python -m pip install \
    'carveme==1.6.6' 'smetana==1.2.1' 'reframed==1.6.0' 'pyscipopt==6.3.0'
"$MANAGER" run -n "$ENV_NAME" diamond version
"$MANAGER" run -n "$ENV_NAME" carve --help >/dev/null
"$MANAGER" run -n "$ENV_NAME" python "$SCRIPT_DIR/smoke_test.py" scip
printf '\nInstallation checked. Activate with: conda activate %s\n' "$ENV_NAME"
