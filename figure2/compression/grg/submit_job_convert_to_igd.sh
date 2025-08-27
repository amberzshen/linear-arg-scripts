chrom=1

dx run swiss-army-knife \
    -iin="/amber/scripts/convert_to_igd.sh" \
    -icmd="bash convert_to_igd.sh $chrom" \
    --destination "/methods_comparisons/grg/" \
    --instance-type mem3_ssd1_v2_x2 \
    --priority high \
    --name "convert_to_igd_chr${chrom}" \
    --brief \
    --extra-args '{"timeoutPolicy": {"*": {"hours": 504}}}' \
    -y
