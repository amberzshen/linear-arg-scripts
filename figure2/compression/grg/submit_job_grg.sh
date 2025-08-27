# 72 cores, 144GB memory
dx run app-swiss-army-knife \
    -iin="/amber/scripts/run_grg.sh" \
    -icmd="bash run_grg.sh" \
    --destination "/methods_comparisons/grg/" \
    --instance-type mem1_ssd1_v2_x72 \
    --priority high \
    --name "grg_chr1" \
    --brief \
    --extra-args '{"timeoutPolicy": {"*": {"hours": 504}}}' \
    -y
