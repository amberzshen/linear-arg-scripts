import h5py
import argparse
import simulate_mutations as sm

parser = argparse.ArgumentParser()
parser.add_argument('n_muts', type=int)
parser.add_argument('mut_type', type=str)
parser.add_argument('h5_dir', type=str)
parser.add_argument('simulation_name', type=str)
parser.add_argument('out_dir', type=str)
args = parser.parse_args()

seed = int(args.simulation_name.split('_')[-1])

if args.mut_type == 'back':
    X = sm.simulate_back_mutation(f'{args.h5_dir}/{args.simulation_name}.h5', args.n_muts, seed)
elif args.mut_type == 'recurrent':
    X = sm.simulate_recurrent_mutation(f'{args.h5_dir}/{args.simulation_name}.h5', args.n_muts, seed)
elif args.mut_type == 'error':
    X = sm.simulate_genotype_error(f'{args.h5_dir}/{args.simulation_name}.h5', args.n_muts, seed)
elif args.mut_type == 'position':
    X = sm.simulate_position_switching(f'{args.h5_dir}/{args.simulation_name}.h5', args.n_muts, seed)

out = f'{args.out_dir}/{args.mut_type}/{args.simulation_name}_{args.n_muts}.h5'
with h5py.File(out, "w") as f:
    f.create_dataset('shape', data=X.shape, compression='gzip', shuffle=True)
    f.create_dataset('indptr', data=X.indptr, compression='gzip', shuffle=True)
    f.create_dataset('indices', data=X.indices, compression='gzip', shuffle=True)
    f.create_dataset('data', data=X.data, compression='gzip', shuffle=True)