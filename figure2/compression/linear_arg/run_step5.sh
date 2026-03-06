#!/bin/bash
job_metadata=$1

set -euo pipefail

curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv tool install git+https://github.com/quattro/linear-dag.git

kodama multi-step-compress step5 \
    --job-metadata $job_metadata \