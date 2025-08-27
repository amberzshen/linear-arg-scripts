#!/bin/bash
chrom=$1

set -euo pipefail

vcf_path="/mnt/project/amber/filtered_vcfs/ukb20279_c${chrom}_b0_v1_250129_whitelist.vcf.gz"
output_name="ukb20279_c${chrom}_b0_v1_250129_whitelist.xsi"

# Clone
OUT_DIR=`pwd`
git clone https://github.com/rwk-unil/xSqueezeIt.git $HOME/xSqueezeIt
cd $HOME/xSqueezeIt

# Clone and build htslib (if you already have htslib set Makefile accordingly and skip)
git submodule update --init htslib
cd htslib
git submodule update --init --recursive
autoheader
autoconf
./configure
make
sudo make install
sudo ldconfig
cd ..

# Clone, build, and install zstd (if you already have zstd set Makefile accordingly and skip)
git clone https://github.com/facebook/zstd.git
cd zstd
make
sudo make install
sudo ldconfig
cd ..

# Build application
make

/home/dnanexus/xSqueezeIt/xsqueezeit --compress -f "$vcf_path" -o $OUT_DIR/$output_name
