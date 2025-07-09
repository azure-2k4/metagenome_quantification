#!/usr/bin/env python3
"""
Simulate metagenomic reads from reference genomes using ART.
"""
import os
import json
import logging
import subprocess
from pathlib import Path
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def calculate_read_counts(genomes, total_pairs):
    """
    Calculate number of read pairs per genome based on abundances.
    
    Args:
        genomes (list): List of genome dictionaries with abundances
        total_pairs (int): Total number of read pairs to simulate
    
    Returns:
        dict: Mapping of genome names to number of read pairs
    """
    read_counts = {}
    for genome in genomes:
        name = genome['name']
        abundance = genome['abundance']
        read_counts[name] = int(total_pairs * abundance)
    return read_counts

def run_art_simulation(
    genome_path,
    output_prefix,
    read_count,
    read_length,
    fragment_size,
    error_rate
):
    """
    Run ART read simulator for a single genome.
    
    Args:
        genome_path (Path): Path to reference genome
        output_prefix (str): Prefix for output files
        read_count (int): Number of read pairs to simulate
        read_length (int): Length of each read
        fragment_size (int): Mean fragment size
        error_rate (float): Sequencing error rate
    """
    try:
        cmd = [
            "art_illumina",
            "-ss", "HS25",  # Illumina HiSeq 2500
            "-i", str(genome_path),
            "-p",  # paired-end
            "-l", str(read_length),
            "-c", str(read_count),
            "-m", str(fragment_size),
            "-s", str(error_rate * 100),  # Convert to percentage
            "-o", output_prefix
        ]
        
        subprocess.run(cmd, check=True)
        logging.info(f"Successfully simulated reads for {genome_path}")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running ART for {genome_path}: {str(e)}")
        raise

def merge_fastq_files(individual_fastqs, output_fastq):
    """
    Merge multiple FASTQ files into a single file.
    
    Args:
        individual_fastqs (list): List of input FASTQ files
        output_fastq (Path): Output merged FASTQ file
    """
    with open(output_fastq, 'w') as outfile:
        for fastq in individual_fastqs:
            with open(fastq) as infile:
                outfile.write(infile.read())

def main(snakemake):
    """
    Main function to simulate metagenomic reads.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Create output directories
    reads_dir = Path(snakemake.output.reads[0]).parent
    reads_dir.mkdir(parents=True, exist_ok=True)
    
    # Get parameters
    params = snakemake.params.reads
    genomes = snakemake.params.accessions

    
    # Calculate read counts per genome
    read_counts = calculate_read_counts(genomes, params['total_pairs'])
    
    # Store ground truth
    ground_truth = {
        "abundances": {genome['name']: genome['abundance'] for genome in genomes},
        "read_counts": read_counts
    }
    
    with open(snakemake.output.truth, 'w') as f:
        json.dump(ground_truth, f, indent=2)
    
    # Simulate reads for each genome
    temp_files = []
    for genome in genomes:
        name = genome['name']
        genome_path = Path(snakemake.input.genomes) / f"{genome['accession']}.fna"

        # Create temporary output prefix
        temp_prefix = reads_dir / f"temp_{name}"

        # Run ART
        run_art_simulation(
            genome_path=genome_path,
            output_prefix=str(temp_prefix),
            read_count=read_counts[name],
            read_length=params['read_length'],
            fragment_size=params['fragment_size'],
            error_rate=params['error_rate']
        )

        # ART outputs files as temp_{name}1.fq and temp_{name}2.fq (no dot before 1/2)
        temp_files.extend([
            temp_prefix.parent / f"{temp_prefix.name}1.fq",
            temp_prefix.parent / f"{temp_prefix.name}2.fq"
        ])

    # Merge read files
    merge_fastq_files(
        [f for f in temp_files if f.name.endswith('1.fq')],
        snakemake.output.reads[0]
    )
    merge_fastq_files(
        [f for f in temp_files if f.name.endswith('2.fq')],
        snakemake.output.reads[1]
    )

    # Clean up temporary files
    for temp_file in temp_files:
        temp_file.unlink()

    logging.info("Read simulation completed successfully")

if __name__ == "__main__":
    main(snakemake)
