import polars as pl
import numpy as np
import linear_dag as ld
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


def get_size_of_linarg(linarg):
    mtx_size = get_size_of_mtx(linarg.A)
    indices_size = linarg.variant_indices.astype(np.uint32).nbytes + linarg.flip.nbytes
    return mtx_size + indices_size
    

def get_ukb_linarg_size(chroms):
    size_in_memory = []
    for chrom in chroms:        
        linarg_dir = f'/mnt/project/linear_args/ukb20279/chr{chrom}'
        files = os.listdir(f'{linarg_dir}')
        ind_arr = np.array([int(f.split('_')[0]) for f in files])
        order = ind_arr.argsort()
        files = np.array(files)[order].tolist() # sort files by index
        for f in files:
            linarg = ld.LinearARG.read(f'{linarg_dir}/{f}/linear_arg.h5')
            size_in_memory.append(get_size_of_linarg(linarg))
    return np.sum(size_in_memory) / 10**9


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
    
    results.append(('1_11_21', 'linarg', get_ukb_linarg_size([1, 11, 21])))
    results.append(('1_11_21', 'scipy', get_ukb_scipy_size([1, 11, 21])))
    results.append(('1_11_21', 'grg', get_ukb_grg_size([1, 11, 21])))

    df = pl.DataFrame(
        results,
        schema=["chr", "method", "size_in_memory_gb"]
    )
    
    os.makedirs("amber/methods_comparisons/results/", exist_ok=True)
    df.write_csv("amber/methods_comparisons/results/size_in_memory.csv")