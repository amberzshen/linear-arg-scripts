import h5py
import os
import polars as pl

def get_linarg_disk_size(chrom):
    disk_size = 0
    linarg_dir = f'/mnt/project/linear_args/ukb20279/chr{chrom}/'
    for partition in os.listdir(linarg_dir):
        with h5py.File(f'{linarg_dir}{partition}/linear_arg.h5', 'r') as f:
            dataset_names = ['indptr', 'indices', 'data', 'variant_indices', 'flip'] 
            for name in dataset_names:
                dset = f[name]
                disk_size += dset.id.get_storage_size()
    return disk_size        


def get_grg_disk_size(chrom):
    return os.path.getsize(f'/mnt/project/methods_comparisons/grg/ukb20279_c{chrom}_b0_v1_250129_whitelist.grg')


def get_xsi_disk_size(chrom):
    return os.path.getsize(f'/mnt/project/methods_comparisons/xsi/ukb20279_c{chrom}_b0_v1_250129_whitelist.xsi')


def get_vcf_disk_size(chrom):
    return os.path.getsize(f'/mnt/project/amber/filtered_vcfs/ukb20279_c{chrom}_b0_v1_250129_whitelist.vcf.gz')


if __name__ == "__main__":
    rows = []
    for chrom in [1, 11, 21]:
        linarg_bytes = get_linarg_disk_size(chrom)
        grg_bytes = get_grg_disk_size(chrom)
        xsi_bytes = get_xsi_disk_size(chrom)
        vcf_bytes = get_vcf_disk_size(chrom)

        rows.append((
            chrom,
            linarg_bytes / (1024**3),
            grg_bytes / (1024**3),
            xsi_bytes / (1024**3),
            vcf_bytes / (1024**3),
        ))

    df = pl.DataFrame(
        rows,
        schema=["chr", "linear_arg_disk_size", "grg_disk_size", "xsi_disk_size", "vcf_disk_size"]
    )

    print(df)
    os.makedirs("amber/methods_comparisons/results/", exist_ok=True)
    df.write_csv("amber/methods_comparisons/results/disk_size_benchmark.csv")