#!/usr/bin/env python3
"""
Quantify genome abundances using Poisson MLE method.
"""
import json
import logging
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def poisson_loglik(abundances, read_counts, presence_matrix, unitig_lengths, total_reads):
    """
    Calculate Poisson log-likelihood for abundance estimation.
    
    Args:
        abundances (np.array): Current abundance estimates
        read_counts (np.array): Observed read counts per unitig
        presence_matrix (np.array): Matrix indicating which unitigs belong to which genomes
        unitig_lengths (np.array): Length of each unitig
        total_reads (int): Total number of mapped reads
    
    Returns:
        float: Negative log-likelihood (for minimization)
    """
    # Expected read counts
    expected = total_reads * np.dot(presence_matrix, abundances) * unitig_lengths
    
    # Calculate log-likelihood (ignoring constant terms)
    loglik = np.sum(read_counts * np.log(expected + 1e-10) - expected)
    
    # Return negative for minimization
    return -loglik

def optimize_abundances(read_counts, presence_matrix, unitig_lengths, total_reads):
    """
    Estimate abundances using Poisson MLE.
    
    Args:
        read_counts (np.array): Observed read counts per unitig
        presence_matrix (np.array): Matrix indicating which unitigs belong to which genomes
        unitig_lengths (np.array): Length of each unitig
        total_reads (int): Total number of mapped reads
    
    Returns:
        np.array: Estimated abundances
    """
    n_genomes = presence_matrix.shape[1]
    if n_genomes == 0:
        logging.warning("No genomes found for optimization. Returning empty array.")
        return np.array([])

    # Initial guess (equal abundances)
    x0 = np.ones(n_genomes) / n_genomes

    # Optimization bounds and constraints
    bounds = [(0, 1) for _ in range(n_genomes)] if n_genomes > 0 else None
    constraints = [{
        'type': 'eq',
        'fun': lambda x: np.sum(x) - 1
    }] if n_genomes > 0 else []

    # Run optimization
    result = minimize(
        poisson_loglik,
        x0,
        args=(read_counts, presence_matrix, unitig_lengths, total_reads),
        bounds=bounds,
        constraints=constraints,
        method='SLSQP'
    )

    if not result.success:
        logging.warning(f"Optimization did not converge: {result.message}")

    return result.x

def process_alignments(bam_file, colors=None):
    """
    Process BAM file to get read counts and unitig lengths.
    
    Args:
        bam_file (str): Path to BAM file
        colors (dict, optional): Color information for TwoPaCo unitigs
    
    Returns:
        tuple: (read counts array, unitig lengths array, total reads,
               list of unitig IDs, list of genome IDs)
    """
    import pysam
    
    # Get unitig information from BAM header
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        unitigs = bam.references
        lengths = bam.lengths
    
    # Count mapped reads
    read_counts = np.zeros(len(unitigs))
    total_reads = 0
    
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                unitig_idx = unitigs.index(read.reference_name)
                read_counts[unitig_idx] += 1
                total_reads += 1
    
    # Create presence matrix
    if colors is not None:
        # TwoPaCo: use color information
        genomes = sorted(list(set(
            genome
            for color_list in colors.values()
            for genome in color_list
        )))
        
        presence_matrix = np.zeros((len(unitigs), len(genomes)))
        for i, unitig in enumerate(unitigs):
            if unitig in colors:
                for genome in colors[unitig]:
                    j = genomes.index(genome)
                    presence_matrix[i, j] = 1
    else:
        # A-Bruijn: each unitig represents one genome
        genomes = unitigs
        presence_matrix = np.eye(len(unitigs))
    
    return (
        read_counts,
        np.array(lengths),
        total_reads,
        unitigs,
        genomes
    )

def main(snakemake):
    """
    Main function to calculate abundances using Poisson MLE.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Load color information for TwoPaCo
    with open(snakemake.input.colors) as f:
        colors = json.load(f)
    
    # Process TwoPaCo alignments
    (
        twopaco_counts,
        twopaco_lengths,
        twopaco_total,
        twopaco_unitigs,
        twopaco_genomes
    ) = process_alignments(snakemake.input.twopaco_bam, colors)

    # Build presence matrix for TwoPaCo
    n_unitigs = len(twopaco_unitigs)
    n_genomes = len(twopaco_genomes)
    presence_matrix = np.zeros((n_unitigs, n_genomes))
    for i, unitig in enumerate(twopaco_unitigs):
        if unitig in colors:
            for genome in colors[unitig]:
                j = twopaco_genomes.index(genome)
                presence_matrix[i, j] = 1

    # Estimate TwoPaCo abundances
    twopaco_abundances = optimize_abundances(
        twopaco_counts,
        presence_matrix,
        twopaco_lengths,
        twopaco_total
    )
    
    # Process A-Bruijn alignments
    (
        abruijn_counts,
        abruijn_lengths,
        abruijn_total,
        abruijn_unitigs,
        abruijn_genomes
    ) = process_alignments(snakemake.input.abruijn_bam)
    
    # Estimate A-Bruijn abundances
    abruijn_abundances = optimize_abundances(
        abruijn_counts,
        np.eye(len(abruijn_unitigs)),
        abruijn_lengths,
        abruijn_total
    )
    
    # Combine results
    results = pd.DataFrame({
        'genome': twopaco_genomes,
        'twopaco_poisson': twopaco_abundances,
        'abruijn_poisson': [
            abruijn_abundances[abruijn_genomes.index(g)]
            if g in abruijn_genomes else 0
            for g in twopaco_genomes
        ]
    })
    
    # Save results
    results.to_csv(snakemake.output[0], index=False)
    logging.info("Poisson MLE quantification completed successfully")

if __name__ == "__main__":
    main(snakemake)
