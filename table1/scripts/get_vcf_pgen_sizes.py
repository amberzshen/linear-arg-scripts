import os
import numpy as np
import polars as pl

def get_vcf_disk_size(chrom, maf=False):
    if maf:
        disk_size = os.path.getsize(f"/mnt/project/amber/vcfs/ukb20279_c{chrom}_b0_v1_whitelist_maf_0.01.vcf.gz")
    else:
        disk_size = os.path.getsize(f"/mnt/project/Bulk/Previous WGS releases/GATK and GraphTyper WGS/SHAPEIT Phased VCFs/ukb20279_c{chrom}_b0_v1.vcf.gz")
    return disk_size / 10**9


def get_pgen_disk_size(chrom, maf=False):
    disk_sizes = []
    for suff in ["pgen", "psam", "pvar"]:
        if maf:
            disk_sizes.append(os.path.getsize(f"/mnt/project/amber/pgens/ukb20279_c{chrom}_b0_v1_whitelist_maf_0.01.{suff}"))
        else:
            disk_sizes.append(os.path.getsize(f"/mnt/project/amber/pgens/ukb20279_c{chrom}_b0_v1_whitelist.{suff}"))
    return np.array(disk_sizes) / 10**9 
    

if __name__ == "__main__":
    rows = []
    for chrom in np.arange(1, 23):
        for dataset in ['common', 'all']:
            vcf_size = get_vcf_disk_size(chrom, maf=(dataset == 'common'))
            pgen_size, psam_size, pvar_size = get_pgen_disk_size(chrom, maf=(dataset == 'common'))
            rows.append((chrom, dataset, vcf_size, pgen_size, psam_size, pvar_size))
    df = pl.DataFrame(rows, schema=["chrom", "dataset", "vcf_size", "pgen_size", "psam_size", "pvar_size"])
    df.write_csv("ukb20279_vcf_pgen_sizes.csv")
