# linear-arg-scripts

Scripts and notebooks to reproduce the figures and tables.

## Directory structure

### Table 1: Compression performance
- `table1/`: Compression statistics for 1000 Genomes, UKB 200k, and UKB 500k
  - `notebooks/make_table1.ipynb`
  - `data/`: UKB linarg stats and vcf/pgen sizes
  - `scripts/`: Data collection scripts

### Figure 2: Comparison of compression methods
- `figure2/2a/`: **Fig 2a**: Size in memory on UKB 200k chr 1, 11, 21
  - `notebooks/make_figure2a.ipynb`
- `figure2/2b/`: **Fig 2b**: Size on disk on UKB 200k chr 1, 11, 21
  - `notebooks/make_figure2b.ipynb`
- `figure2/compression/`: Scripts to run compression methods (kodama, GRG, XSI) on UKB data
- `figure2/data/`: Size benchmarks and linarg stats for Figure 2

### Figure 3: Genealogical compression of simulated genotype data
- `figure3/3a/`: **Fig 3a**: NNZ comparison between simulated ARG, inferred linear ARG, and sparse genotype matrix
  - `make_figure3a.ipynb`
  - `data/`: Size comparison results on simulated data
- `figure3/3c/`: **Fig 3c**: Brick graph and linear ARG NNZ ratio vs. number of mutations/errors
  - `make_figure3c.ipynb`
  - `data/`: Robustness to error results
- `figure3/simulations/`: Simulation scripts (msprime Out of Africa model, mutation simulation, linear ARG inference)

### Figure 4: GWAS benchmark
- `figure4/4ac/`: **Fig 4a**: Common variant GWAS runtime; **Fig 4c**: All-variant GWAS runtime
  - `notebooks/make_figure4ac.ipynb`
  - `data/`: GWAS runtime benchmarks
  - `scripts/`: GWAS scripts to run kodama-gwas and PLINK 2.0
- `figure4/4bd/`: **Fig 4b**: Common variant z-score concordance; **Fig 4d**: All-variant z-score concordance
  - `notebooks/make_figure4bd.ipynb`
  - `data/`: GWAS concordance and correlation data (note: `kodama-gwas_plink_common_bilirubin_gwas_comparison_chr11.tsv` and `kodama-gwas-rac_plink_all_bilirubin_gwas_comparison_chr11.tsv` were excluded due to file size)
  - `scripts/`: GWAS concordance scripts

### Figure 5: Matrix-vector multiplication benchmark
- `figure5/`: **Fig 5**: Load time and matrix-vector CPU time on UKB 200k chr 1, 11, 21
  - `notebooks/make_figure5.ipynb`
  - `data/`: Matrix multiplication benchmarks
  - `scripts/`: Benchmark scripts

### Supplementary Figures
- `supplementary/`
  - `notebooks/make_supp_figure1.ipynb`: Common variant GWAS benchmark
  - `notebooks/make_supp_figure2.ipynb`: All variant GWAS benchmark
  - `notebooks/make_supp_figure3a.ipynb`: Serial right matrix multiplication
  - `notebooks/make_supp_figure3b.ipynb`: PRS runtime benchmark
  - `data/`: PRS runtime benchmarks
