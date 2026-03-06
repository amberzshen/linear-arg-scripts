#!/bin/bash
set -euo pipefail
beta_prefix=$1
dataset=$2
node=$3
shift 3   # drop first 3 args, now $1 is the first chromosome
chromosomes=("$@")   # everything left goes into an array

sudo apt update
sudo apt install libffi-dev
curl https://pyenv.run | bash
export PYENV_ROOT="$HOME/.pyenv"
[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init -)"
pyenv install 3.10.13
pyenv global 3.10.13

pip install h5py
pip install polars
pip install memory_profiler
pip install git+https://github.com/quattro/linear-dag.git

if [[ "$dataset" == "all" ]]; then
    dataset_path="ukb20279_chr1-22.h5"
else
    dataset_path="ukb20279_maf_0.01_chr1-22.h5"
fi

score_cols=$(cat "${beta_prefix}_cols.txt" | tr '\n' ' ')

python /mnt/project/amber/scripts/make_parquet.py "${beta_prefix}_weights.tsv" "${beta_prefix}_weights"

bash /mnt/project/amber/scripts/profile.sh linarg_prs_${beta_prefix}_${node}.csv \
    kodama score \
        --linarg-path $dataset_path \
        --beta-path "${beta_prefix}_weights.parquet" \
        --score-cols $score_cols  \
        --chrom ${chromosomes[@]} \
        --out linarg_prs_${beta_prefix}_${node}