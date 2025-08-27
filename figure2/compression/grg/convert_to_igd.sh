#!/bin/bash
chrom=$1

vcf_path="/mnt/project/amber/filtered_vcfs/ukb20279_c${chrom}_b0_v1_250129_whitelist.vcf.gz"
igd_path="ukb20279_c${chrom}_b0_v1_250129_whitelist.igd"

# download python 3.10
sudo apt update
sudo apt install libffi-dev
curl https://pyenv.run | bash
export PYENV_ROOT="$HOME/.pyenv"
[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init -)"
pyenv install 3.10.13
pyenv global 3.10.13 

pip install pygrgl

grg convert "$vcf_path" $igd_path