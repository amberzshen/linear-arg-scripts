import argparse
import scipy.sparse as sp
import numpy as np
import pickle as pkl

from linear_dag.core.brick_graph import BrickGraph
from linear_dag.core.recombination import Recombination
from linear_dag.core.one_summed_cy import linearize_brick_graph


def infer_linarg(genotypes_path, simulation_name, mut_type, n_muts, linarg_out):
    print('loading genotypes...', flush=True)
    genotypes = sp.load_npz(f'{genotypes_path}/{mut_type}/{simulation_name}_{n_muts}.npz') # samples x variants matrix
    print('inferring linear ARG...', flush=True)
    brick_graph, sample_indices, variant_indices = BrickGraph.from_genotypes(genotypes)
    recom = Recombination.from_graph(brick_graph)
    recom.find_recombinations()
    linarg_adj_mat = linearize_brick_graph(recom)
    brick_graph_adj_mat = brick_graph.to_csr().T.tocsr()
    
    print('saving linear ARG...', flush=True)
    with open(f'{linarg_out}/{simulation_name}_{mut_type}_{n_muts}.pkl', 'wb') as file:
        pkl.dump([brick_graph_adj_mat, linarg_adj_mat, np.array(sample_indices), np.array(variant_indices)], file)
    
    nnzs = [genotypes.nnz, brick_graph_adj_mat.nnz, linarg_adj_mat.nnz]
    nnzs = [str(x) for x in nnzs]
    with open(f'{linarg_out}/{simulation_name}_{mut_type}_{n_muts}_stats.txt', 'w') as file:
        file.write(" ".join(['genotypes_nnz', 'brick_graph_nnz', 'linarg_nnz']) + "\n")
        file.write(" ".join(nnzs) + "\n")    
        
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('genotypes_path', type=str) # "/n/data1/hms/dbmi/oconnor/lab/amber/simulate_mutations/simulated_tree_sequences_w_simulated_mutations/"
    parser.add_argument('simulation_name', type=str) # f"100000_chr21-21990355-23324211_{seed}"
    parser.add_argument('mut_type', type=str)
    parser.add_argument('n_muts', type=str)
    parser.add_argument('linarg_out', type=str) # "/n/data1/hms/dbmi/oconnor/lab/amber/simulate_mutations/linear_args"
    args = parser.parse_args()
    
    infer_linarg(args.genotypes_path, args.simulation_name, args.mut_type, args.n_muts, args.linarg_out)

    
