import numpy as np
import linear_dag as ld
import os
import pandas as pd
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
    

def get_ukb_linarg_size(chroms):
    size_in_memory = []
    for chrom in chroms:        
        linarg_dir = f'/mnt/project/linear_args/ukb20279/chr{chrom}'
        files = os.listdir(f'{linarg_dir}')
        ind_arr = np.array([int(f.split('_')[0]) for f in files])
        order = ind_arr.argsort()
        files = np.array(files)[order].tolist() # sort files by index
        for f in files:
            if not os.path.exists(f'{linarg_dir}/{f}/linear_arg.h5'):
                continue
            linarg = ld.LinearARG.read(f'{linarg_dir}/{f}/linear_arg.h5')
            size_in_memory.append(get_size_of_linarg(linarg))
    return np.sum(size_in_memory) / 10**9


def get_linarg_disk_size(chroms):
    disk_size = 0
    for chrom in chroms:
        linarg_dir = f'/mnt/project/linear_args/ukb20279/chr{chrom}/'
        for partition in os.listdir(linarg_dir):
            if not os.path.exists(f'{linarg_dir}{partition}/linear_arg.h5'):
                continue
            with h5py.File(f'{linarg_dir}{partition}/linear_arg.h5', 'r') as f:
                dataset_names = ['indptr', 'indices', 'data', 'variant_indices', 'flip'] 
                for name in dataset_names:
                    dset = f[name]
                    disk_size += dset.id.get_storage_size()
    return disk_size / 10**9


def get_variant_metadata_disk_size(chroms):
    disk_size = 0
    for chrom in chroms:
        linarg_dir = f'/mnt/project/linear_args/ukb20279/chr{chrom}/'
        for partition in os.listdir(linarg_dir):
            if not os.path.exists(f'{linarg_dir}{partition}/linear_arg.h5'):
                continue
            with h5py.File(f'{linarg_dir}{partition}/linear_arg.h5', 'r') as f:
                dataset_names = ["CHROM", "POS", "ID", "REF", "ALT"]
                for name in dataset_names:
                    dset = f[name]
                    disk_size += dset.id.get_storage_size()
    return disk_size / 10**9


def get_vcf_disk_size(chroms):
    disk_size = 0
    for chrom in chroms:
        disk_size += os.path.getsize(f"/mnt/project/Bulk/Previous WGS releases/GATK and GraphTyper WGS/SHAPEIT Phased VCFs/ukb20279_c{chrom}_b0_v1.vcf.gz")
    return disk_size
        
        
def get_nnz_ratio(chroms):
    linarg_nnz = 0
    genotypes_nnz = 0
    n_variants = 0
    linarg_path = '/mnt/project/linear_args/ukb20279'
    for chrom in chroms:
        for p in os.listdir(f'{linarg_path}/chr{chrom}'):
            if not os.path.exists(f'{linarg_path}/chr{chrom}/{p}/linear_arg_stats.txt'):
                continue
            tmp = pd.read_csv(f'{linarg_path}/chr{chrom}/{p}/linear_arg_stats.txt', delim_whitespace=True, header=0)
            linarg_nnz += tmp.linarg_nnz
            genotypes_nnz += tmp.genotypes_nnz
            n_variants += tmp.m
    return genotypes_nnz / linarg_nnz, n_variants


if __name__ == "__main__":
    
    chroms = [str(i) for i in range(1, 23)] + ['X']
    
    size_in_memory = get_ukb_linarg_size(chroms)
    nnz_ratio, n_variants = get_nnz_ratio(chroms)
    linarg_disk_size = get_linarg_disk_size(chroms)
    variant_disk_size = get_variant_metadata_disk_size(chroms)
    vcf_disk_size = get_vcf_disk_size(chroms)
    
    results = {
        "size_in_memory": size_in_memory,
        "n_variants": n_variants,
        "nnz_ratio": nnz_ratio,
        "linarg_disk_size": linarg_disk_size,
        "variant_disk_size": variant_disk_size,
        "vcf_disk_size": vcf_disk_size,
    }

    with open("ukb20279_linarg_stats.txt", "w") as f:
        for name, value in results.items():
            f.write(f"{name}: {value}\n")
        