#!/usr/bin/env python3
"""
Quantify genome abundances using RPKM method.
"""
import json
import logging
import pandas as pd
import pysam
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def calculate_rpkm(reads_mapped, length, total_mapped_reads):
    """
    Calculate RPKM for a sequence.
    
    Args:
        reads_mapped (int): Number of reads mapped to sequence
        length (int): Length of sequence
        total_mapped_reads (int): Total number of mapped reads
    
    Returns:
        float: RPKM value
    """
    return (reads_mapped * 1e9) / (length * total_mapped_reads)

def get_unitig_lengths(bam_file):
    """
    Get lengths of all unitigs from BAM header.
    
    Args:
        bam_file (str): Path to BAM file
    
    Returns:
        dict: Mapping of unitig IDs to lengths
    """
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        return {ref: length for ref, length in zip(bam.references, bam.lengths)}

def count_mapped_reads(bam_file, unitigs):
    """
    Count reads mapped to each unitig.
    
    Args:
        bam_file (str): Path to BAM file
        unitigs (list): List of unitig IDs
    
    Returns:
        tuple: (dict of read counts per unitig, total mapped reads)
    """
    read_counts = {u: 0 for u in unitigs}
    total_mapped = 0

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                unitig = read.reference_name
                if unitig in read_counts:
                    read_counts[unitig] += 1
                else:
                    # Optionally log or skip unitigs not in the provided list
                    logging.warning(f"Read mapped to unknown unitig: {unitig}, skipping.")
                total_mapped += 1

    return read_counts, total_mapped

def calculate_genome_rpkm(
    read_counts,
    unitig_lengths,
    colors,
    total_mapped_reads
):
    """
    Calculate RPKM for each genome based on unique unitigs.
    
    Args:
        read_counts (dict): Read counts per unitig
        unitig_lengths (dict): Length of each unitig
        colors (dict): Genome membership of each unitig
        total_mapped_reads (int): Total mapped reads
    
    Returns:
        dict: RPKM values per genome
    """
    # Initialize counters
    genome_reads = {}
    genome_lengths = {}
    
    # Process each unitig
    for unitig, genomes in colors.items():
        # Only use unitigs unique to one genome
        if len(genomes) == 1:
            genome = genomes[0]
            if genome not in genome_reads:
                genome_reads[genome] = 0
                genome_lengths[genome] = 0
            
            genome_reads[genome] += read_counts[unitig]
            genome_lengths[genome] += unitig_lengths[unitig]
    
    # Calculate RPKM for each genome
    genome_rpkm = {
        genome: calculate_rpkm(reads, genome_lengths[genome], total_mapped_reads)
        for genome, reads in genome_reads.items()
    }
    
    return genome_rpkm

def normalize_abundances(rpkm_values):
    """
    Convert RPKM values to relative abundances.
    
    Args:
        rpkm_values (dict): RPKM values per genome
    
    Returns:
        dict: Relative abundances per genome
    """
    total_rpkm = sum(rpkm_values.values())
    return {
        genome: rpkm / total_rpkm
        for genome, rpkm in rpkm_values.items()
    }

def main(snakemake):
    """
    Main function to calculate abundances using RPKM.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Load color information
    with open(snakemake.input.colors) as f:
        colors = json.load(f)
    
    # Process TwoPaCo results
    twopaco_bam = snakemake.input.twopaco_bam
    
    # Get unitig lengths and read counts
    unitig_lengths = get_unitig_lengths(twopaco_bam)
    read_counts, total_mapped = count_mapped_reads(
        twopaco_bam,
        list(colors.keys())
    )
    
    # Calculate abundances
    twopaco_rpkm = calculate_genome_rpkm(
        read_counts,
        unitig_lengths,
        colors,
        total_mapped
    )
    twopaco_abundances = normalize_abundances(twopaco_rpkm)
    
    # Process A-Bruijn results
    abruijn_bam = snakemake.input.abruijn_bam
    
    # Get read counts for A-Bruijn graph
    abruijn_lengths = get_unitig_lengths(abruijn_bam)
    abruijn_counts, abruijn_total = count_mapped_reads(
        abruijn_bam,
        list(abruijn_lengths.keys())
    )
    
    # Calculate A-Bruijn abundances (simplified version)
    abruijn_rpkm = {
        ref: calculate_rpkm(count, abruijn_lengths[ref], abruijn_total)
        for ref, count in abruijn_counts.items()
    }
    abruijn_abundances = normalize_abundances(abruijn_rpkm)
    
    # Combine results
    results = pd.DataFrame({
        'genome': list(twopaco_abundances.keys()),
        'twopaco_rpkm': [twopaco_abundances[g] for g in twopaco_abundances],
        'abruijn_rpkm': [abruijn_abundances[g] for g in twopaco_abundances]
    })
    
    # Save results
    results.to_csv(snakemake.output[0], index=False)
    logging.info("RPKM quantification completed successfully")

if __name__ == "__main__":
    main(snakemake)
