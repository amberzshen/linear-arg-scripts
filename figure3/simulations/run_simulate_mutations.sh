#!/bin/bash
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH -p medium
#SBATCH --mem=32GB
#SBATCH --output=/n/data1/hms/dbmi/oconnor/lab/amber/slurm/250421_sim_muts_%a.out
#SBATCH --error=/n/data1/hms/dbmi/oconnor/lab/amber/slurm/250421_sim_muts_%a.err
#SBATCH --job-name="sim_muts"
#SBATCH --array=1-1200
#SBATCH --mail-user amber_shen@g.harvard.edu
#SBATCH --mail-type=ALL

params="sim_mut_params2.csv"
simulation_name=$(awk -F',' 'NR==i {print $j}' i=$SLURM_ARRAY_TASK_ID j=1 $params)
mut_type=$(awk -F',' 'NR==i {print $j}' i=$SLURM_ARRAY_TASK_ID j=2 $params)
n_muts=$(awk -F',' 'NR==i {print $j}' i=$SLURM_ARRAY_TASK_ID j=3 $params)


python run_simulate_mutations.py \
        $n_muts \
        $mut_type \
        "/n/data1/hms/dbmi/oconnor/lab/amber/msprime_simulations/small_ooa/" \
        $simulation_name \
        "/n/data1/hms/dbmi/oconnor/lab/amber/msprime_simulations/small_ooa_sim_mut/"
