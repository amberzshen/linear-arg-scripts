import h5py
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("NUM_WEIGHTS", type=int)
parser.add_argument("dataset", type=str)
parser.add_argument("out", type=str)
parser.add_argument("chromosomes", type=str, nargs="+")
args = parser.parse_args()

if args.dataset == 'common':
    linarg_prefix = f'/mnt/project/linear_args/ukb20279_maf_0.01'
else:
    linarg_prefix = f'/mnt/project/linear_args/ukb20279'    

output_path = f"{args.out}_weights.tsv"
with open(output_path, "w") as out_f:
    
    header = "ID\tA1\t" + "\t".join([f"weight_{i+1}" for i in range(args.NUM_WEIGHTS)])
    out_f.write(header + "\n")
    
    for chrom in args.chromosomes:
        linarg_dir = f'{linarg_prefix}/chr{chrom}'
        for partition in os.listdir(linarg_dir):
            with h5py.File(f'{linarg_dir}/{partition}/linear_arg.h5', 'r') as f:
                POS = f['POS'][:]
                REF = [x.decode('utf-8') for x in f['REF'][:]]
                ALT = [x.decode('utf-8') for x in f['ALT'][:]]
                ID = [x.decode('utf-8') for x in f['ID'][:]]

                for pos, ref, alt, linarg_id in zip(POS, REF, ALT, ID):
                    alt1 = alt.split(",")[0]
                    variant_id = f'chr{chrom}:{pos}:{ref}:{alt1}'
                    weights = np.ones(args.NUM_WEIGHTS)
                    weights_str = "\t".join(f"{w:.6f}" for w in weights)
                    out_f.write(f"{variant_id}\t{alt1}\t{weights_str}\n")

with open(f"{args.out}_cols.txt", "w") as f:
    for i in range(args.NUM_WEIGHTS):
        f.write(f"weight_{i+1}\n")