import numpy as np
import polars as pl
from scipy.stats import norm
    
def get_covariates():
    covariates_path = '/mnt/project/phenotypes/age_sex_height_pcs.csv'
    whitelist_path = '/mnt/project/sample_metadata/ukb20279/250129_whitelist.txt'
    with open(whitelist_path, 'r') as f:
        whitelist = np.array([int(line.strip()) for line in f], dtype=np.int64)
    covariates = pl.read_csv(covariates_path)
    covariates = covariates.filter(pl.col('eid').is_in(whitelist))
    covariates = covariates.rename({"eid": "FID"})
    covariates = covariates.with_columns(pl.col("FID").cast(pl.Utf8))
    covariates = covariates.with_columns(pl.col("FID").alias("IID"))
    covariates = covariates.with_columns(
        pl.when(pl.col("p31") == "Male").then(1).otherwise(0).alias("p31")
    )
    covariate_ids = ['p21022', 'p31'] + [f'p22009_a{i}' for i in range(1, 41)]
    for col in covariate_ids:
        if (col == 'IID') or (col == 'FID'):
            continue
        mean = covariates[col].mean()
        std_dev = covariates[col].std()   
        covariates = covariates.with_columns(
            ((pl.col(col) - mean) / std_dev).alias(col)
        )
    covariates = covariates.select(["FID", "IID"] + covariate_ids)
    covariate_ids = ["IID"] + covariate_ids

    return covariates, covariate_ids


def get_phenotypes(quant_norm=True):
    phenotypes = pl.read_csv('/mnt/project/phenotypes/phenotypes_median_20250617.csv')
    whitelist_path = '/mnt/project/sample_metadata/ukb20279/250129_whitelist.txt'
    with open(whitelist_path, 'r') as f:
        whitelist = np.array([int(line.strip()) for line in f], dtype=np.int64)

    # filter and rename
    phenotypes = phenotypes.filter(pl.col("eid").is_in(whitelist))
    phenotypes = phenotypes.rename({"eid": "FID"})
    phenotypes = phenotypes.with_columns(pl.col("FID").cast(pl.Utf8))
    phenotypes = phenotypes.with_columns(pl.col("FID").alias("IID"))
    
    # select numeric columns, excluding iid
    numeric_types = {pl.Int8, pl.Int16, pl.Int32, pl.Int64,
                     pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64,
                     pl.Float32, pl.Float64}

    numeric_cols = [col for col, dtype in zip(phenotypes.columns, phenotypes.dtypes)
                    if (col != "IID") and (col != "FID") and (dtype in numeric_types)]

    phenotypes = phenotypes.select(["FID", "IID"] + numeric_cols)

    # handle missingness
    n_rows = phenotypes.height
    fraction_missing = phenotypes.select([
        (pl.col(col).is_null().sum() / n_rows).alias(col)
        for col in phenotypes.columns
    ])
    high_missing_cols = [
        col for col in fraction_missing.columns
        if fraction_missing.select(pl.col(col)).item() > 0.8
    ]
    phenotypes = phenotypes.drop(high_missing_cols)
    phenotype_ids = phenotypes.columns

    # quantile normalize
    for col in phenotype_ids:
        if (col == 'IID') or (col == 'FID'):
            continue
        if quant_norm:
            col_norm = quantile_norm(np.array(phenotypes[col]))
            phenotypes = phenotypes.with_columns(
                pl.Series(name=col, values=col_norm)
            )
        else:
            mean = phenotypes[col].mean()
            std_dev = phenotypes[col].std()   
            phenotypes = phenotypes.with_columns(
                ((pl.col(col) - mean) / std_dev).alias(col)
            )

    return phenotypes, phenotype_ids



def quantile_norm(x):
    n = len(x)
    x_norm = np.full(n, np.nan)
    valid_mask = ~np.isnan(x)
    valid_x = x[valid_mask]
    if len(valid_x) == 0:
        return x_norm
    ranks = (np.arange(1, len(valid_x) + 1) - 0.5) / len(valid_x)  # better empirical CDF
    normalized_values = norm.ppf(ranks)
    sort_indices = np.argsort(valid_x)
    x_norm_indices = np.where(valid_mask)[0]
    x_norm[x_norm_indices[sort_indices]] = normalized_values
    return x_norm


if __name__ == "__main__":
    
    phenotypes, phenotype_ids = get_phenotypes()
    covariates, covariate_ids = get_covariates()
        
    phenotypes.write_csv(file="phenotypes.tsv", separator='\t', null_value='NA')
    covariates.write_csv(file="covariates.tsv", separator='\t', null_value='NA')
    
    with open("phenotypes.txt", "w") as f:
        for x in phenotype_ids:
            f.write(f"{x}\n")
    
    with open("covariates.txt", "w") as f:
        for x in covariate_ids:
            f.write(f"{x}\n")