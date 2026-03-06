#!/bin/bash
job_metadata=$1
large_job_id=$2

set -euo pipefail

curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv tool install git+https://github.com/quattro/linear-dag.git

kodama multi-step-compress step3 \
    --job-metadata $job_metadata \
    --large-job-id $large_job_id