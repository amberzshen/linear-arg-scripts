instance_type=$1
job_metadata="/mnt/project/linear_args/ukb30108/job_metadata.parquet"

dx download "/linear_args/ukb30108/job_metadata.parquet"
n_large_jobs=$(python -c "import polars as pl; print(pl.read_parquet('job_metadata.parquet')['large_job_id'].n_unique())")

echo $n_large_jobs

for ((i=0; i<n_large_jobs; i++)); do

    dx run app-swiss-army-knife \
        -icmd="bash /mnt/project/amber/scripts/run_step3.sh $job_metadata $i" \
        --destination "/linear_args/" \
        --instance-type $instance_type \
        --priority high \
        --name "ukb500k_step3_${i}" \
        --brief \
        -y

done