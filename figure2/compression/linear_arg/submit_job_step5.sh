instance_type=$1
job_metadata="/mnt/project/linear_args/ukb30108/job_metadata.parquet"

dx run app-swiss-army-knife \
    -icmd="bash /mnt/project/amber/scripts/run_step5.sh $job_metadata" \
    --destination "/linear_args/" \
    --instance-type $instance_type \
    --priority high \
    --name "ukb500k_step5" \
    --brief \
    --ignore-reuse \
    -y
