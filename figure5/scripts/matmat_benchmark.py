import pygrgl
import numpy as np
from memory_profiler import memory_usage
import time
from scipy.sparse import load_npz
import polars as pl
from linear_dag import LinearARG
from linear_dag.core.lineararg import list_blocks
from linear_dag.cli import _filter_blocks
import os


def benchmark(func, *args, repeat=1, interval=0.001, **kwargs):
    total_time = 0.0
    peak_memories = []
    result = None
    for _ in range(repeat):
        start_time = time.perf_counter()
        mem_usage, result = memory_usage(
            (func, args, kwargs),
            retval=True,
            interval=interval,
            timeout=None,
            include_children=True,
            max_usage=True
        )
        end_time = time.perf_counter()
        total_time += (end_time - start_time)
        peak_memories.append(mem_usage)
    mean_time = total_time / repeat # seconds
    mean_peak_mem = sum(peak_memories) / repeat # MB
    return mean_time, mean_peak_mem


def grg_matmat(grg_path, n_vectors):
    grg = pygrgl.load_immutable_grg(grg_path)
    load_time, load_mem = benchmark(pygrgl.load_immutable_grg, grg_path)
    b = np.ones([grg.num_mutations, n_vectors]).T
    res = pygrgl.matmul(grg, b, pygrgl.TraversalDirection.DOWN)
    dp_time, dp_mem = benchmark(pygrgl.matmul, grg, b, pygrgl.TraversalDirection.DOWN)
    return load_time, load_mem, dp_time, dp_mem


def linarg_matmat(linarg_path, block_name, n_vectors):
    linarg = LinearARG.read(linarg_path, block_name)
    load_time, load_mem = benchmark(LinearARG.read, linarg_path, block_name)
    b = np.ones([linarg.shape[1], n_vectors])
    dp_time, dp_mem = benchmark(lambda: linarg @ b)
    return load_time, load_mem, dp_time, dp_mem


def scipy_matmat(genotypes_path, n_vectors):
    genotypes = load_npz(genotypes_path)
    load_time, load_mem = benchmark(load_npz, genotypes_path)
    b = np.ones([genotypes.shape[1], n_vectors])
    dp_time, dp_mem = benchmark(lambda: genotypes @ b)
    return load_time, load_mem, dp_time, dp_mem

    
if __name__ == "__main__":
    
    results = []
    
    for n_vectors in [1]:
        
        print(f'n_vectors: {n_vectors}')
        
        print('linarg')
        linarg_path = '/mnt/project/final_linear_args/ukb20279_chr1-22.h5'
        block_metadata = list_blocks(linarg_path)
        block_metadata = _filter_blocks(block_metadata, chromosomes=['1', '11', '21'])
        for block_name in block_metadata['block_name']:
            chrom = block_name.split('_')[0]
            load_time, load_mem, matmat_time, matmat_mem = linarg_matmat(linarg_path, block_name, n_vectors)
            results.append((chrom, block_name, n_vectors, 'linarg', load_time, load_mem, matmat_time, matmat_mem))
        
        print('grg')
        for chrom in [1, 11, 21]:
            grg_path = f'/mnt/project/methods_comparisons/grg/ukb20279_c{chrom}_b0_v1_250129_whitelist.grg'
            load_time, load_mem, matmat_time, matmat_mem = grg_matmat(grg_path, n_vectors)
            results.append((chrom, chrom, n_vectors, 'grg', load_time, load_mem, matmat_time, matmat_mem))
              
        print('scipy')
        for chrom in [1, 11, 21]:  
            genotypes_dir = f'/mnt/project/linear_args/ukb20279/chr{chrom}/'
            for large_partition in os.listdir(genotypes_dir):
                for file in os.listdir(f'{genotypes_dir}/{large_partition}/genotype_matrices/'):
                    genotypes_path = f'{genotypes_dir}/{large_partition}/genotype_matrices/{file}'
                    load_time, load_mem, matmat_time, matmat_mem = scipy_matmat(genotypes_path, n_vectors)
                    results.append((chrom, large_partition, n_vectors, 'scipy', load_time, load_mem, matmat_time, matmat_mem))
        
    df = pl.DataFrame(
        results,
        schema=["chr", "region", "n_vectors", "method", "load_time", "max_load_memory", "matmat_time", "max_matmat_memory"],
        orient="row",
    )
    
    df.write_csv("linarg_scipy_grg_matvec_serial_benchmark.csv")