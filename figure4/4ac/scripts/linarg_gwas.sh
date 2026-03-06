#!/bin/bash
set -euo pipefail
dataset=$1
node=$2
shift 2   # drop first 2 args, now $1 is the first chromosome
chromosomes=("$@")   # everything left goes into an array

curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv tool install git+https://github.com/quattro/linear-dag.git

if [[ "$dataset" == "all" ]]; then
    linarg="ukb20279_chr1-22.h5"
else
    linarg="ukb20279_maf_0.01_chr1-22.h5"
fi

pheno_cols=$(cat phenotypes.txt | tr '\n' ' ')
covar_cols=$(cat covariates.txt | tr '\n' ' ')

bash /mnt/project/amber/scripts/profile.sh linarg_gwas_${dataset}_chr1-22_${node}.csv \
    kodama assoc \
        $linarg \
        phenotypes.tsv \
        --covar covariates.tsv \
        --pheno-name $pheno_cols \
        --covar-name $covar_cols  \
        --chrom ${chromosomes[@]}  \
        --out linarg_gwas_${dataset}_chr1-22_${node}

