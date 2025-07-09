#!/usr/bin/env python3
"""
Quantify genome abundances using Non-Negative Least Squares method.
"""
import json
import logging
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def build_presence_matrix(unitigs, colors):
    """
    Build presence matrix A where A[i,j] = 1 if unitig i belongs to genome j.
    
    Args:
        unitigs (list): List of unitig IDs
        colors (dict): Mapping of unitigs to genome indices
    
    Returns:
        tuple: (presence matrix, list of genome IDs)
    """
    # Get unique genome IDs
    genomes = sorted(list(set(
        genome
        for color_list in colors.values()
        for genome in color_list
    )))
    
    # Build matrix
    n_unitigs = len(unitigs)
    n_genomes = len(genomes)
    A = np.zeros((n_unitigs, n_genomes))
    
    for i, unitig in enumerate(unitigs):
        if unitig in colors:
            for genome in colors[unitig]:
                j = genomes.index(genome)
                A[i, j] = 1
    
    return A, genomes

def get_read_counts(bam_file, unitigs):
    """
    Count reads mapped to each unitig.
    
    Args:
        bam_file (str): Path to BAM file
        unitigs (list): List of unitig IDs
    
    Returns:
        np.array: Read counts for each unitig
    """
    import pysam
    
    counts = np.zeros(len(unitigs))
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                unitig_idx = unitigs.index(read.reference_name)
                counts[unitig_idx] += 1
    
    return counts

def normalize_abundances(abundances):
    """
    Normalize abundance estimates to sum to 1.
    
    Args:
        abundances (np.array): Raw abundance estimates
    
    Returns:
        np.array: Normalized abundances
    """
    return abundances / np.sum(abundances)

def main(snakemake):
    """
    Main function to calculate abundances using NNLS.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Load color information
    with open(snakemake.input.colors) as f:
        colors = json.load(f)
    
    # Get unitig list from BAM file
    import pysam
    with pysam.AlignmentFile(snakemake.input.twopaco_bam, 'rb') as bam:
        unitigs = bam.references
    
    # Build presence matrix
    A, genomes = build_presence_matrix(unitigs, colors)
    
    # Get read counts
    b = get_read_counts(snakemake.input.twopaco_bam, unitigs)
    
    # Solve NNLS problem
    abundances, residuals = nnls(A, b)
    
    # Normalize abundances
    abundances = normalize_abundances(abundances)
    
    # Create results DataFrame
    results = pd.DataFrame({
        'genome': genomes,
        'nnls_abundance': abundances
    })
    
    # Save results
    results.to_csv(snakemake.output[0], index=False)
    
    # Log optimization quality
    relative_residual = np.sqrt(residuals) / np.sum(b)
    logging.info(f"NNLS optimization complete (relative residual: {relative_residual:.2e})")
    logging.info("NNLS quantification completed successfully")

if __name__ == "__main__":
    main(snakemake)
