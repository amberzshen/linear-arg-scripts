#!/bin/bash
set -euo pipefail
dataset=$1
node=$2
chrom=$3

OUT_DIR=`pwd`

if [[ "$dataset" == "all" ]]; then
    pgen="ukb20279_c${chrom}_b0_v1_whitelist"
else
    pgen="ukb20279_c${chrom}_b0_v1_whitelist_maf_0.01"
fi

pheno_cols=$(grep -vwE "IID|FID" phenotypes.txt | tr '\n' ' ')
covar_cols=$(grep -vwE "IID|FID" covariates.txt | tr '\n' ' ')

bash /mnt/project/amber/scripts/profile.sh $OUT_DIR/plink_gwas_${dataset}_chr${chrom}_${node}.csv \
    plink2 \
        --pfile $pgen \
        --double-id \
        --no-input-missing-phenotype \
        --pheno phenotypes.tsv --pheno-name $pheno_cols \
        --covar covariates.tsv --covar-name $covar_cols \
        --glm hide-covar omit-ref \
        --out plink_gwas_${dataset}_chr${chrom}_${node}