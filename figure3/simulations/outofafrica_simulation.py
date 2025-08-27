import argparse
import stdpopsim

# region in the format of chrXX-XXX-XXX
def simulate_OutOfAfrica_3G09(sample_size, region, seed, out):
    n = int(sample_size/3)
    chrom, start, end = region.split('-')
    species = stdpopsim.get_species('HomSap')
    model = species.get_demographic_model('OutOfAfrica_3G09')
    contig = species.get_contig(chrom, mutation_rate=model.mutation_rate, inclusion_mask=[(start, end)])
    samples = {'YRI': n, 'CHB': n, 'CEU': n}
    engine = stdpopsim.get_engine('msprime')
    ts = engine.simulate(model, contig, samples, seed=seed)
    ts.dump(out)
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_size', type=int)
    parser.add_argument('region', type=str)
    parser.add_argument('seed', type=int)
    parser.add_argument('out', type=str)
    args = parser.parse_args()
    
    simulate_OutOfAfrica_3G09(args.sample_size, args.region, args.seed, f'{args.out}/{args.sample_size}_{args.region}_{args.seed}.trees')