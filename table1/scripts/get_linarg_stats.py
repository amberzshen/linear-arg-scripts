import numpy as np
import linear_dag as ld
from linear_dag.core.lineararg import list_blocks
import pandas as pd
import polars as pl
import h5py


def get_size_of_mtx(A):    
    if (np.min(A.data) > -127) and (np.max(A.data) < 128):
        data_size = A.data.astype(np.int8).nbytes
    else:
        data_size = A.data.nbytes
    size_in_memory = data_size + A.indices.astype(np.uint32).nbytes + A.indptr.astype(np.uint32).nbytes
    return size_in_memory


def get_size_of_linarg(linarg):
    mtx_size = get_size_of_mtx(linarg.A)
    indices_size = linarg.variant_indices.astype(np.uint32).nbytes + linarg.flip.nbytes
    return mtx_size + indices_size


def get_linarg_size(linarg_path, chrom):
    size_in_memory = 0
    block_metadata = list_blocks(linarg_path)
    block_metadata = block_metadata.filter(pl.col("chrom") == chrom)
    
    for block in block_metadata.get_column("block_name").to_list():
        linarg = ld.LinearARG.read(linarg_path, block=block)
        size_in_memory += get_size_of_linarg(linarg)
    return size_in_memory / 10**9


def get_linarg_disk_size(linarg_path, chrom):
    disk_size = 0
    block_metadata = list_blocks(linarg_path)
    block_metadata = block_metadata.filter(pl.col("chrom") == chrom)
    
    dataset_names = ['indptr', 'indices', 'data', 'variant_indices', 'flip'] 
    for block in block_metadata.get_column("block_name").to_list():
        with h5py.File(linarg_path, 'r') as f:
            for name in dataset_names:
                dset = f[block][name]
                disk_size += dset.id.get_storage_size()
    return disk_size / 10**9


def get_variant_metadata_disk_size(linarg_path, chrom):
    disk_size = 0
    block_metadata = list_blocks(linarg_path)
    block_metadata = block_metadata.filter(pl.col("chrom") == chrom)
    
    dataset_names = ["CHROM", "POS", "ID", "REF", "ALT"]
    for block in block_metadata.get_column("block_name").to_list():
        with h5py.File(linarg_path, 'r') as f:
            for name in dataset_names:
                dset = f[block][name]
                disk_size += dset.id.get_storage_size()
    return disk_size / 10**9


def get_nnz_ratio(linarg_path, chrom):
    linarg_nnz = 0
    genotypes_nnz = 0
    n_variants = 0
    block_metadata = list_blocks(linarg_path)
    block_metadata = block_metadata.filter(pl.col("chrom") == chrom)
    
    for block in block_metadata.get_column("block_name").to_list():
        with h5py.File(linarg_path, 'r') as f:
            linarg_nnz += f[block].attrs["n_entries"]
            n_variants += f[block].attrs["n_variants"]
            linarg = ld.LinearARG.read(linarg_path, block=block)
            v = np.ones(linarg.shape[0])
            allele_counts = v @ linarg
            genotypes_nnz += np.sum(allele_counts)
    return genotypes_nnz, linarg_nnz, n_variants


if __name__ == "__main__":
    
    linarg_dir = '/mnt/project/final_linear_args/'
    linargs = [
        'ukb20279_chr1-22',
        'ukb20279_chr1-22_individual',
        'ukb20279_maf_0.01_chr1-22',
        'ukb20279_maf_0.01_chr1-22_individual'
    ]
    chroms = [str(i) for i in range(1, 23)]
    
    result = pd.DataFrame(columns=["linarg", "chrom", "size_in_memory", "linarg_disk_size", "variant_disk_size", "genotypes_nnz", "linarg_nnz", "n_variants"])
    
    for linarg in linargs:
        
        print(linarg, flush=True)
        
        linarg_path = f'{linarg_dir}/{linarg}.h5'

        for chrom in chroms:
            print(chrom, flush=True)
            
            print('getting size in memory', flush=True)
            size_in_memory = get_linarg_size(linarg_path, chrom)
            print('getting linarg disk size', flush=True)
            linarg_disk_size = get_linarg_disk_size(linarg_path, chrom)
            print('getting variant disk size', flush=True)
            variant_disk_size = get_variant_metadata_disk_size(linarg_path, chrom)
            print('getting nnz ratio', flush=True)
            genotypes_nnz, linarg_nnz, n_variants = get_nnz_ratio(linarg_path, chrom)
            
            result.loc[len(result)] = [linarg, chrom, size_in_memory, linarg_disk_size, variant_disk_size, genotypes_nnz, linarg_nnz, n_variants]

    result.to_csv("ukb20279_linarg_stats.csv", index=False)