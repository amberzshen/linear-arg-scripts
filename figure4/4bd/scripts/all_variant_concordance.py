import pandas as pd
import numpy as np
import pyarrow.parquet as pq
import pyarrow as pa
import polars as pl


def read_parquet_row_range(file_path, columns, start, end):
    pf = pq.ParquetFile(file_path)
    
    rows_to_read = end - start
    accumulated_rows = 0
    row_groups_to_read = []
    row_group_start_offsets = []

    # Identify which row groups to read
    for rg_index in range(pf.num_row_groups):
        rg_num_rows = pf.metadata.row_group(rg_index).num_rows
        if accumulated_rows + rg_num_rows > start:
            row_groups_to_read.append(rg_index)
            row_group_start_offsets.append(accumulated_rows)
        if accumulated_rows >= end:
            break
        accumulated_rows += rg_num_rows

    # Read and slice
    tables = []
    for rg_index, rg_start in zip(row_groups_to_read, row_group_start_offsets):
        rg_table = pf.read_row_group(rg_index, columns=columns)
        
        # Determine slice within this row group
        slice_start = max(start - rg_start, 0)
        slice_end = min(end - rg_start, rg_table.num_rows)
        rg_table = rg_table.slice(slice_start, slice_end - slice_start)

        tables.append(rg_table)

    full_table = pa.concat_tables(tables)
    return full_table.to_pandas()


def compute_trait_r(trait, linarg_sumstats, phenotypes, result):
    
    c = 21
    plink_path = f'/mnt/project/amber/final_gwas_benchmark/plink/plink_gwas_all_chr{c}_mem3_ssd1_v2_x32.{trait}.glm.linear'
    plink_sumstats = pd.read_csv(plink_path, sep='\t')
    
    sumstats = pd.merge(linarg_sumstats, plink_sumstats, on='ID', how='inner')
    sumstats = sumstats.dropna()
    
    missingness = phenotypes[trait].isna().mean()
    r_beta = np.corrcoef(sumstats['BETA'], sumstats[f'{trait}_BETA'])[0, 1]
    r_se = np.corrcoef(sumstats['SE'], sumstats[f'{trait}_SE'])[0, 1]
    r_z = np.corrcoef(sumstats['BETA']/sumstats['SE'], sumstats[f'{trait}_BETA']/sumstats[f'{trait}_SE'])[0, 1]
    
    result.loc[len(result)] = [trait, r_z, r_beta, r_se, missingness, len(sumstats)]
    
    return result


def compute_all_trait_r(linarg_sumstats_path, out):
    
    c = 21
    
    phenotypes = pd.read_csv('/mnt/project/amber/final_gwas_benchmark/phenotypes/phenotypes.tsv', sep='\t')
    traits = [x for x in phenotypes.columns if (x!='FID') and (x!='IID')]
    
    metadata = pd.read_csv('/mnt/project/amber/ukb20279_chr1-22_variant_info.tsv', sep='\t')
    metadata = metadata[metadata.chrom.isin([1, 11, 21])]    
    start = np.sum(metadata[metadata.chrom < c].n_variants)
    end = np.sum(metadata[metadata.chrom < c+1].n_variants)
    cols = ['ID'] + [f'{trait}_BETA' for trait in traits] + [f'{trait}_SE' for trait in traits]
    
    linarg_sumstats = read_parquet_row_range(linarg_sumstats_path, cols, start, end)    
    linarg_sumstats.ID = [x.decode('utf-8')[3:] for x in linarg_sumstats.ID]

    result = pd.DataFrame(columns=['trait', 'r_z', 'r_beta', 'r_se', 'missingness', 'n_variants'])
    
    for trait in traits:
       result = compute_trait_r(trait, linarg_sumstats, phenotypes, result)
       
    result.to_csv(out, sep='\t', index=False)


if __name__ == "__main__":
    
    experiments = [
        'linarg_gwas_all_chr1-11-21_mem3_ssd1_v2_x64',
        'linarg_gwas_repeat_covar_all_chr1-11-21_mem3_ssd1_v2_x64',
        'linarg_gwas_no_hwe_all_chr1-11-21_mem3_ssd1_v2_x64',
        'linarg_gwas_no_hwe_repeat_covar_all_chr1-11-21_mem3_ssd1_v2_x64'        
    ]
    
    for exp in experiments:
        linarg_sumstats_path = f'/mnt/project/amber/final_gwas_benchmark/linarg/{exp}.parquet'
        out = f'{exp}_plink_concordance.tsv'
        compute_all_trait_r(linarg_sumstats_path, out)
