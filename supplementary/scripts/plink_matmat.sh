#!/bin/bash
set -euo pipefail
beta_prefix=$1
dataset=$2
node=$3
chrom=$4

if [[ "$dataset" == "all" ]]; then
    pgen="ukb20279_c${chrom}_b0_v1_whitelist"
else
    pgen="ukb20279_c${chrom}_b0_v1_whitelist_maf_0.01"
fi

start_col=3
end_col=$(awk -F'\t' '{print NF; exit}' "${beta_prefix}_weights.tsv")
if [[ "$start_col" -eq "$end_col" ]]; then
    cols="$start_col"
else
    cols="${start_col}-${end_col}"
fi

bash /mnt/project/amber/scripts/profile.sh plink2_prs_${chrom}_${beta_prefix}_${node}.csv \
    plink2 \
        --pfile $pgen \
        --score "${beta_prefix}_weights.tsv" 1 2 header cols=+scoresums \
        --score-col-nums $cols \
        --out plink2_prs_${chrom}_${beta_prefix}_${node}