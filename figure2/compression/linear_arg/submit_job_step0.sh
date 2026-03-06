vcf_metadata="/mnt/project/amber/vcf_metadata/ukb30108_allchr_vcf_metadata.txt"
partition_size=20000000 # 20Mb large partitions
n_small_blocks=40 # 0.5Mb small partitions
out="ukb30108"

dx run app-swiss-army-knife \
    -icmd="bash /mnt/project/amber/scripts/run_step0.sh $vcf_metadata $partition_size $n_small_blocks $out" \
    --destination "/linear_args/" \
    --instance-type "mem1_ssd1_v2_x2" \
    --priority high \
    --name "ukb500k_step0" \
    --brief \
    --ignore-reuse \
    --extra-args '{"timeoutPolicyByExecutable": {"app-J2fv5P89f2ZFbj52533QZKPG": {"*": {"hours": 720}}}}' \
    -y