#!/bin/bash
vcf_metadata=$1
partition_size=$2
n_small_blocks=$3
out=$4

set -euo pipefail

curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv tool install git+https://github.com/quattro/linear-dag.git

kodama multi-step-compress step0 \
    --vcf-metadata $vcf_metadata \
    --partition-size $partition_size \
    --n-small-blocks $n_small_blocks \
    --flip-minor-alleles \
    --mount-point "/mnt/project/linear_args/" \
    --out $out \