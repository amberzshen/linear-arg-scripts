# cython: language_level=3, boundscheck=False, wraparound=False

from scipy.sparse import csc_matrix, csr_matrix
import numpy as np
cimport numpy as np
import msprime
import tskit
import h5py
import time


def simulate_genotype_error(X, n_errors, seed):
    cdef int i, j, k
    cdef int n_rows = X.shape[0]
    cdef int n_cols = X.shape[1]
    np.random.seed(seed)
    cdef np.ndarray[np.int32_t, ndim=1] row_inds = np.random.randint(0, n_rows, size=n_errors, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] col_inds = np.random.randint(0, n_cols, size=n_errors, dtype=np.int32)
    X_lil = X.tolil()
    for k in range(n_errors):
        i = row_inds[k]
        j = col_inds[k]
        if X_lil[i, j] != 0:
            X_lil[i, j] = 0
        else:
            X_lil[i, j] = 1
    return X_lil.tocsc()


def simulate_position_switching(X, num_pairs, seed):
    np.random.seed(seed)
    n_cols = X.shape[1]
    variant_pairs = np.random.choice(n_cols, (num_pairs, 2), replace=False)
    X_lil = X.tolil()
    for i, j in variant_pairs:
        X_lil[:, [i, j]] = X_lil[:, [j, i]]
    return X_lil.tocsc()


def simulate_back_mutation(h5_path: str, n_mutations: int, seed: int) -> "csc_matrix":

    cdef int i, n_nodes, M
    cdef np.ndarray[np.float64_t, ndim=1] mutation_randoms = np.random.random(n_mutations) # pre-sample for speed
    cdef np.ndarray[np.float64_t, ndim=1] descendant_randoms = np.random.random(n_mutations)
    cdef double mutation_site
    cdef int mutation_node

    np.random.seed(seed)
    genotypes, mutation_ids, mutation_sites, mutation_nodes, left_intervals, right_intervals, parents, children, _ = read_h5(h5_path)
    
    mutation_nodes = mutation_nodes.astype(np.int32)
    cdef np.float64_t[:] mutation_sites_mv = mutation_sites
    cdef np.int32_t[:] mutation_nodes_mv = mutation_nodes
    M = len(mutation_ids)
    n_nodes = len(set(list(parents) + list(children)))

    for i in range(n_mutations):

        mut_idx = int(np.floor(mutation_randoms[i] * M))
        mutation_site = mutation_sites_mv[mut_idx]
        mutation_node = mutation_nodes_mv[mut_idx]

        marginal_tree = get_marginal_tree((mutation_site, mutation_site), left_intervals, right_intervals, parents, children, n_nodes)[0]
        mutation_descendants = get_descendants(mutation_node, marginal_tree)
        if len(mutation_descendants) == 0: # mutation has no descendants
            continue

        back_mutation_node = int(np.floor(descendant_randoms[i] * len(mutation_descendants)))
        samples_to_flip = get_descendants(back_mutation_node, marginal_tree)
        if len(samples_to_flip) == 0:
            continue
        for s in samples_to_flip:
            if s < genotypes.shape[0]:
                genotypes[s, mut_idx] = 0

    genotypes.eliminate_zeros()
    return genotypes


def simulate_recurrent_mutation(h5_path: str, n_mutations: int, seed: int) -> "csc_matrix":

    cdef int i, n_nodes, M
    cdef np.ndarray[np.float64_t, ndim=1] mutation_randoms = np.random.random(n_mutations) # pre-sample for speed
    cdef np.ndarray[np.float64_t, ndim=1] recurrent_randoms = np.random.random(n_mutations)
    cdef double mutation_site
    cdef int mutation_node

    np.random.seed(seed)
    genotypes, mutation_ids, mutation_sites, mutation_nodes, left_intervals, right_intervals, parents, children, _ = read_h5(h5_path)
    
    mutation_nodes = mutation_nodes.astype(np.int32)
    cdef np.float64_t[:] mutation_sites_mv = mutation_sites
    cdef np.int32_t[:] mutation_nodes_mv = mutation_nodes
    M = len(mutation_ids)
    n_nodes = len(set(list(parents) + list(children)))

    for i in range(n_mutations):

        mut_idx = int(np.floor(mutation_randoms[i] * M))
        mutation_site = mutation_sites_mv[mut_idx]
        mutation_node = mutation_nodes_mv[mut_idx]

        marginal_tree, nodes_in_interval = get_marginal_tree((mutation_site, mutation_site), left_intervals, right_intervals, parents, children, n_nodes)
        mutation_ancestor_descendants = set(get_descendants(mutation_node, marginal_tree)).union(get_descendants(mutation_node, marginal_tree.T)) # ancestor or descendant of variant
        recurrent_mutation_node_candidates = list(set(nodes_in_interval).difference(mutation_ancestor_descendants)) # nodes that are neither ancestors or descendants defined at the position
        if len(recurrent_mutation_node_candidates) == 0:
            continue

        recurrent_mutation_node = int(np.floor(recurrent_randoms[i] * len(recurrent_mutation_node_candidates)))
        samples_to_flip = get_descendants(recurrent_mutation_node, marginal_tree)
        if len(samples_to_flip) == 0:
            continue
        for s in samples_to_flip:
            if s < genotypes.shape[0]: # only flip samples
                genotypes[s, mut_idx] = 1

    return genotypes


def dfs(csr_indices: np.ndarray, csr_indptr: np.ndarray, node: int, visited: np.ndarray) -> None:
    cdef int idx, neighbor
    cdef list stack = [node]
    visited[node] = 1
    while stack:
        node = stack.pop()
        for idx in range(csr_indptr[node], csr_indptr[node + 1]):
            neighbor = csr_indices[idx]
            if not visited[neighbor]:
                visited[neighbor] = 1
                stack.append(neighbor)


def get_marginal_tree(interval: tuple, left_ints: np.ndarray, right_ints: np.ndarray, parents: np.ndarray, children: np.ndarray, N: int) -> tuple:
    cdef np.ndarray parents_in_int
    cdef np.ndarray children_in_int
    cdef set nodes_in_interval = set()

    mask = (left_ints <= interval[0]) & (right_ints >= interval[1])
    parents_in_int = parents[mask]
    children_in_int = children[mask]
    
    marginal_tree = csr_matrix((np.ones(len(parents_in_int)), (parents_in_int, children_in_int)), shape=(N, N))
    nodes_in_interval.update(parents_in_int.tolist())
    nodes_in_interval.update(children_in_int.tolist())
    
    return marginal_tree, nodes_in_interval


def get_descendants(node: int, G: "csr_matrix") -> np.ndarray:
    cdef np.ndarray csr_indices = G.indices
    cdef np.ndarray csr_indptr = G.indptr
    cdef np.ndarray visited = np.zeros(G.shape[0])
    dfs(csr_indices, csr_indptr, node, visited)
    decsendants = np.where(visited)[0]
    return decsendants


def extract_data_from_ts(ts):
    site_positions = {site.id: site.position for site in ts.sites()} # just in case mutation ids are not in genomic order
    mutation_ids = np.array([m.id for m in ts.mutations()])
    mutation_sites = np.array([site_positions[m.site] for m in ts.mutations()])
    mutation_nodes = np.array([m.node for m in ts.mutations()])
    left_intervals = ts.tables.edges.left
    right_intervals = ts.tables.edges.right
    parents = ts.tables.edges.parent
    children = ts.tables.edges.child
    observed_genotypes = ts.genotype_matrix().T
    population_labels = ts.tables.nodes.population[np.where(ts.tables.nodes.individual != -1)[0]] # subset to sample nodes
    return observed_genotypes, mutation_ids, mutation_sites, mutation_nodes, left_intervals, right_intervals, parents, children, population_labels


def ts_to_h5(ts, out):
    observed_genotypes, mutation_ids, mutation_sites, mutation_nodes, left_intervals, right_intervals, parents, children, population_labels = extract_data_from_ts(ts)
    observed_genotypes = csc_matrix(observed_genotypes)
    with h5py.File(out, "w") as f:
        f.create_dataset('geno_shape', data=observed_genotypes.shape, compression='gzip', shuffle=True)
        f.create_dataset('geno_data', data=observed_genotypes.data, compression='gzip', shuffle=True)
        f.create_dataset('geno_indptr', data=observed_genotypes.indptr, compression='gzip', shuffle=True)
        f.create_dataset('geno_indices', data=observed_genotypes.indices, compression='gzip', shuffle=True)
        f.create_dataset('mutation_ids', data=mutation_ids, compression='gzip', shuffle=True)
        f.create_dataset('mutation_sites', data=mutation_sites, compression='gzip', shuffle=True) 
        f.create_dataset('mutation_nodes', data=mutation_nodes, compression='gzip', shuffle=True) 
        f.create_dataset('left_intervals', data=left_intervals, compression='gzip', shuffle=True) 
        f.create_dataset('right_intervals', data=right_intervals, compression='gzip', shuffle=True) 
        f.create_dataset('parents', data=parents, compression='gzip', shuffle=True) 
        f.create_dataset('children', data=children, compression='gzip', shuffle=True) 
        f.create_dataset('population_labels', data=population_labels, compression='gzip', shuffle=True) 


def read_h5(h5_path):
    with h5py.File(h5_path, 'r') as f:
        observed_genotypes = csc_matrix((f['geno_data'][:], f['geno_indices'][:], f['geno_indptr'][:]), shape=f['geno_shape'][:])  
        mutation_ids = f['mutation_ids'][:]  
        mutation_sites = f['mutation_sites'][:]  
        mutation_nodes = f['mutation_nodes'][:]  
        left_intervals = f['left_intervals'][:]  
        right_intervals = f['right_intervals'][:]  
        parents = f['parents'][:]  
        children = f['children'][:]  
        population_labels = f['population_labels'][:]  
    return observed_genotypes, mutation_ids, mutation_sites, mutation_nodes, left_intervals, right_intervals, parents, children, population_labels
