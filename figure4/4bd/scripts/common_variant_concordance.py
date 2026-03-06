import pandas as pd
import numpy as np

def compute_trait_r(trait, linarg_sumstats, phenotypes, result):
    
    plink_sumstats = []
    for c in range(1,23):
        plink_path = f'/mnt/project/amber/final_gwas_benchmark/plink/plink_gwas_common_chr{c}_mem3_ssd1_v2_x32.{trait}.glm.linear'
        df = pd.read_csv(plink_path, sep='\t')
        plink_sumstats.append(df)        
    plink_sumstats = pd.concat(plink_sumstats)
    
    sumstats = pd.merge(linarg_sumstats, plink_sumstats, on='ID', how='inner')
    sumstats = sumstats.dropna()
    
    missingness = phenotypes[trait].isna().mean()
    r_beta = np.corrcoef(sumstats['BETA'], sumstats[f'{trait}_BETA'])[0, 1]
    r_se = np.corrcoef(sumstats['SE'], sumstats[f'{trait}_SE'])[0, 1]
    r_z = np.corrcoef(sumstats['BETA']/sumstats['SE'], sumstats[f'{trait}_BETA']/sumstats[f'{trait}_SE'])[0, 1]
    
    result.loc[len(result)] = [trait, r_z, r_beta, r_se, missingness, len(sumstats)]
    
    return result


def compute_all_trait_r(linarg_sumstats_path, out):
    
    linarg_sumstats = pd.read_parquet(linarg_sumstats_path)
    linarg_sumstats.ID = [x.decode('utf-8')[3:] for x in linarg_sumstats.ID]
    phenotypes = pd.read_csv('/mnt/project/amber/final_gwas_benchmark/phenotypes/phenotypes.tsv', sep='\t')
    traits = [x for x in phenotypes.columns if (x!='FID') and (x!='IID')]
    result = pd.DataFrame(columns=['trait', 'r_z', 'r_beta', 'r_se', 'missingness', 'n_variants'])
    
    for trait in traits:
       result = compute_trait_r(trait, linarg_sumstats, phenotypes, result)
       
    result.to_csv(out, sep='\t', index=False)


if __name__ == "__main__":
    
    experiments = [
        'linarg_gwas_common_chr1-22_mem3_ssd1_v2_x32',
        'linarg_gwas_repeat_covar_common_chr1-22_mem3_ssd1_v2_x32',
        'linarg_gwas_no_hwe_common_chr1-22_mem3_ssd1_v2_x64',
        'linarg_gwas_no_hwe_repeat_covar_common_chr1-22_mem3_ssd1_v2_x64'        
    ]
    
    for exp in experiments:
        linarg_sumstats_path = f'/mnt/project/amber/final_gwas_benchmark/linarg/{exp}.parquet'
        out = f'{exp}_plink_concordance.tsv'
        compute_all_trait_r(linarg_sumstats_path, out)
