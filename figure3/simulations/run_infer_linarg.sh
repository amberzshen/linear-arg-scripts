#!/bin/bash
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -p short
#SBATCH --mem=64GB
#SBATCH --output=/n/data1/hms/dbmi/oconnor/lab/amber/slurm/250204_sim_muts_infer_linarg_%a.out
#SBATCH --error=/n/data1/hms/dbmi/oconnor/lab/amber/slurm/250204_sim_muts_infer_linarg_%a.err
#SBATCH --job-name="sim_muts_infer_linarg"
#SBATCH --array=1-2000
#SBATCH --mail-user amber_shen@g.harvard.edu
#SBATCH --mail-type=ALL

params="../params/ooa_infer_linarg_params.csv"

simulation_name=$(awk -F',' 'NR==i {print $j}' i=$SLURM_ARRAY_TASK_ID j=1 $params)
mut_type=$(awk -F',' 'NR==i {print $j}' i=$SLURM_ARRAY_TASK_ID j=2 $params)
n_muts=$(awk -F',' 'NR==i {print $j}' i=$SLURM_ARRAY_TASK_ID j=3 $params)
genotypes_path="/n/data1/hms/dbmi/oconnor/lab/amber/simulate_mutations/simulated_tree_sequences_w_simulated_mutations/"
linarg_out="/n/data1/hms/dbmi/oconnor/lab/amber/simulate_mutations/linear_args/"

echo "${genotypes_path}, ${simulation_name}, ${mut_type}, ${n_muts}, ${linarg_out}"

python infer_linarg.py $genotypes_path $simulation_name $mut_type $n_muts $linarg_out
