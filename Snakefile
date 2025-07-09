import os
from snakemake.utils import min_version

# Set minimum Snakemake version
min_version("7.0")

# Load configuration
configfile: "config.yaml"

# Define final output files
rule all:
    input:
        "results/evaluation/summary_metrics.csv",
        "results/evaluation/runtime_memory.csv",
        "results/plots/abundance_comparison.pdf",
        "results/plots/performance_comparison.pdf",
        "results/plots/accuracy_metrics.pdf"

# Download reference genomes
rule download_genomes:
    output:
        directory("data/genomes")
    params:
        accessions=config["genomes"]
    script:
        "scripts/download_genomes.py"

# Simulate metagenomic reads
rule simulate_reads:
    input:
        genomes="data/genomes"
    output:
        reads=["data/reads/simulated_1.fq", "data/reads/simulated_2.fq"],
        truth="data/ground_truth.json"
    params:
        reads=config["reads"],
        accessions=config["genomes"]
    script:
        "scripts/simulate_reads.py"

# TwoPaCo Pipeline
rule build_twopaco:
    input:
        genomes="data/genomes"
    output:
        graph="results/twopaco/ecoli_cdbg.gfa",
        unitigs="results/twopaco/unitigs.fasta"
    params:
        twopaco=config["twopaco"]
    threads: config["threads"]
    script:
        "scripts/build_twopaco.py"

rule extract_colors:
    input:
        graph="results/twopaco/ecoli_cdbg.gfa"
    output:
        colors="results/twopaco/unitig_colors.json"
    script:
        "scripts/extract_colors.py"

rule map_reads_twopaco:
    input:
        unitigs="results/twopaco/unitigs.fasta",
        reads=["data/reads/simulated_1.fq", "data/reads/simulated_2.fq"]
    output:
        bam="results/twopaco/alignments.bam"
    params:
        config["bwa"]
    threads: config["threads"]
    script:
        "scripts/map_reads.py"

# A-Bruijn Pipeline
rule build_abruijn:
    input:
        genomes="data/genomes"
    output:
        graph="results/abruijn/abruijn.gfa",
        consensus="results/abruijn/consensus.fasta"
    params:
        config["blast"]
    threads: config["threads"]
    script:
        "scripts/build_abruijn.py"

rule map_reads_abruijn:
    input:
        consensus="results/abruijn/consensus.fasta",
        reads=["data/reads/simulated_1.fq", "data/reads/simulated_2.fq"]
    output:
        bam="results/abruijn/alignments.bam"
    params:
        config["bwa"]
    threads: config["threads"]
    script:
        "scripts/map_reads.py"

# Quantification Methods
rule quantify_rpkm:
    input:
        twopaco_bam="results/twopaco/alignments.bam",
        abruijn_bam="results/abruijn/alignments.bam",
        colors="results/twopaco/unitig_colors.json"
    output:
        "results/evaluation/rpkm_abundances.csv"
    script:
        "scripts/quantify_rpkm.py"

rule quantify_poisson:
    input:
        twopaco_bam="results/twopaco/alignments.bam",
        abruijn_bam="results/abruijn/alignments.bam",
        colors="results/twopaco/unitig_colors.json"
    output:
        "results/evaluation/poisson_abundances.csv"
    script:
        "scripts/quantify_poisson.py"

rule quantify_nnls:
    input:
        twopaco_bam="results/twopaco/alignments.bam",
        colors="results/twopaco/unitig_colors.json"
    output:
        "results/evaluation/nnls_abundances.csv"
    script:
        "scripts/quantify_nnls.py"

# Evaluation and Visualization
rule evaluate_methods:
    input:
        rpkm="results/evaluation/rpkm_abundances.csv",
        poisson="results/evaluation/poisson_abundances.csv",
        nnls="results/evaluation/nnls_abundances.csv",
        truth="data/ground_truth.json"
    output:
        metrics="results/evaluation/summary_metrics.csv",
        runtime="results/evaluation/runtime_memory.csv"
    script:
        "scripts/evaluate_methods.py"

rule visualize_results:
    input:
        metrics="results/evaluation/summary_metrics.csv",
        runtime="results/evaluation/runtime_memory.csv"
    output:
        abundance="results/plots/abundance_comparison.pdf",
        performance="results/plots/performance_comparison.pdf",
        accuracy="results/plots/accuracy_metrics.pdf"
    script:
        "scripts/visualize_results.py"
