#!/bin/bash
# Integration test for metagenomic quantification pipeline

# Set up error handling
set -e
set -o pipefail

# Function to check if a file exists and is non-empty
check_file() {
    if [ ! -f "$1" ] || [ ! -s "$1" ]; then
        echo "Error: File $1 does not exist or is empty"
        exit 1
    fi
}

# Create test environment
echo "Creating test environment..."
conda env create -f environment.yml
conda activate metagenome_quant

# Run pipeline
echo "Running pipeline..."
snakemake --cores 4 -p

# Check output files
echo "Checking output files..."

# Check genome downloads
check_file "data/genomes/NC_000913.3.fna"
check_file "data/genomes/CP053602.fna"
check_file "data/genomes/NC_012967.1.fna"
check_file "data/genomes/CP009072.fna"

# Check simulated reads
check_file "data/reads/simulated_1.fq"
check_file "data/reads/simulated_2.fq"
check_file "data/ground_truth.json"

# Check TwoPaCo output
check_file "results/twopaco/ecoli_cdbg.gfa"
check_file "results/twopaco/unitigs.fasta"
check_file "results/twopaco/unitig_colors.json"
check_file "results/twopaco/alignments.bam"

# Check A-Bruijn output
check_file "results/abruijn/abruijn.gfa"
check_file "results/abruijn/consensus.fasta"
check_file "results/abruijn/alignments.bam"

# Check quantification results
check_file "results/evaluation/rpkm_abundances.csv"
check_file "results/evaluation/poisson_abundances.csv"
check_file "results/evaluation/nnls_abundances.csv"

# Check evaluation output
check_file "results/evaluation/summary_metrics.csv"
check_file "results/evaluation/runtime_memory.csv"

# Check plots
check_file "results/plots/abundance_comparison.pdf"
check_file "results/plots/performance_comparison.pdf"
check_file "results/plots/accuracy_metrics.pdf"

# Validate results
echo "Validating results..."
python tests/test_quantification.py

# Clean up
conda deactivate

echo "Integration test completed successfully!"
