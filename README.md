# Metagenomic Genome Quantification Pipeline

This pipeline implements and compares two approaches for quantifying closely-related genomes in metagenomic samples:

1. TwoPaCo-based colored de Bruijn graph (cDBG) with three quantification methods:
   - Unique-unitig RPKM
   - Poisson MLE
   - Non-negative Least Squares (NNLS)

2. A-Bruijn graph baseline with two quantification methods:
   - RPKM
   - Poisson MLE

## Requirements

- Python 3.9+
- 16 GB RAM minimum
- 4-8 CPU cores recommended
- Linux/Unix environment (for all tools)

## Installation

1. Create and activate the conda environment:
```bash
conda env create -f environment.yml
conda activate metagenome_quant
```

## Dataset

The pipeline uses synthetic data from 4 E. coli strains:
- E. coli K-12 MG1655 (NC_000913.3)
- E. coli BL21(DE3) (CP053602)
- E. coli REL606 (NC_012967.1)
- E. coli ATCC 25922 (CP009072)

Ground truth abundances: 10%, 20%, 30%, 40% respectively.

## Usage

1. Run the complete pipeline:
```bash
snakemake --cores 4 -p
```

2. Run specific steps:
```bash
# Download genomes only
snakemake download_genomes --cores 4

# Run TwoPaCo analysis only
snakemake build_twopaco --cores 4
```

3. Run tests:
```bash
# Unit tests
python tests/test_quantification.py

# Integration test
bash tests/integration_test.sh
```

## Pipeline Steps

1. **Data Preparation**
   - Download reference genomes
   - Simulate metagenomic reads using ART

2. **TwoPaCo Pipeline**
   - Build colored de Bruijn graph
   - Extract color information
   - Map reads to unitigs
   - Quantify abundances using three methods

3. **A-Bruijn Pipeline**
   - Analyze genome similarities
   - Build A-Bruijn graph
   - Map reads
   - Quantify abundances

4. **Evaluation**
   - Calculate accuracy metrics
   - Measure runtime and memory usage
   - Generate visualization plots

## Output

The pipeline generates:
- CSV files with abundance estimates
- Evaluation metrics and comparisons
- Performance benchmarks
- Visualization plots
  - Abundance comparisons
  - Performance metrics
  - Accuracy summaries

## Directory Structure

```
project/
├── config.yaml           # Pipeline configuration
├── environment.yml       # Conda environment specification
├── Snakefile            # Workflow definition
├── data/                # Input data
├── scripts/             # Analysis scripts
├── results/             # Output files
└── tests/               # Unit and integration tests
```


## Authors

This project was developed by
- [A ANANDAKRISHNA](https://github.com/azure-2k4)
- [RAMANAND BALAJI](https://github.com/Sky-Lance)
- [BALU M KRISHNA](https://github.com/yourgithubhandle)  

as part of the coursework for *22BIO211: Intelligence of Biological Systems II* within the BTech Computer Science Engineering (Artificial Intelligence) curriculum at *Amrita School of Computing, Amrita Vishwa Vidyapeetham, Amritapuri Campus*.  

## References

1. M. Wang, Y. Ye, and H. Tang, "**A de Bruijn graph approach to the quantification of closely-related genomes in a microbial community**", Journal of Computational Biology, vol. 19, no. 6, pp. 814–825, Jun. 2012, doi: 10.1089/cmb.2012.0058
2. I. Minkin, S. Pham, and P. Medvedev, "**TwoPaCo: an efficient algorithm to build the compacted de Bruijn graph from many complete genomes**", Bioinformatics, Volume 33, Issue 24, December 2017, Pages 4024–4032, doi: 10.1093/bioinformatics/btw609

