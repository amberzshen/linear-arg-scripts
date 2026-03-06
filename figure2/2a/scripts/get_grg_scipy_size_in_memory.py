import polars as pl
import numpy as np
from scipy.sparse import load_npz
import pygrgl
import psutil
import os

def get_size_of_mtx(A):    
    if (np.min(A.data) > -127) and (np.max(A.data) < 128):
        data_size = A.data.astype(np.int8).nbytes
    else:
        data_size = A.data.nbytes
    size_in_memory = data_size + A.indices.astype(np.uint32).nbytes + A.indptr.astype(np.uint32).nbytes
    return size_in_memory


def get_ukb_scipy_size(chroms):
    size_in_memory = []
    for chrom in chroms:        
        linarg_dir = f'/mnt/project/linear_args/ukb20279/chr{chrom}/'
        for large_partition in os.listdir(linarg_dir):
            for file in os.listdir(f'{linarg_dir}/{large_partition}/genotype_matrices/'):
                genotypes_path = f'{linarg_dir}/{large_partition}/genotype_matrices/{file}'
                A = load_npz(genotypes_path)
                size_in_memory.append(get_size_of_mtx(A))
    return np.sum(size_in_memory) / 10**9


def get_ukb_grg_size(chroms):
    size_in_memory = []
    for chrom in chroms:
        grg_path = f'/mnt/project/methods_comparisons/grg/ukb20279_c{chrom}_b0_v1_250129_whitelist.grg'
        size_in_memory.append(estimate_grg_memory_usage(grg_path))
    return np.sum(size_in_memory) / 10**9


def estimate_grg_memory_usage(grg_path):
    """
    Loads a GRG object from disk and estimates its memory usage in MB.
    """
    def get_process_memory():
        process = psutil.Process(os.getpid())
        return process.memory_info().rss  # in bytes

    mem_before = get_process_memory()
    grg = pygrgl.load_immutable_grg(grg_path)
    mem_after = get_process_memory()

    return mem_after - mem_before


if __name__ == "__main__":
    
    results = []
    
    for chrom in [1, 11, 21]:
        results.append((chrom, 'scipy', get_ukb_scipy_size([chrom])))
        results.append((chrom, 'grg', get_ukb_grg_size([chrom])))

    df = pl.DataFrame(
        results,
        schema=["chr", "method", "size_in_memory_gb"]
    )
    
    df.write_csv("grg_scipy_size_in_memory.tsv", separator='\t')