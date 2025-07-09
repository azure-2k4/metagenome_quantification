#!/usr/bin/env python3
"""
Map reads to unitigs using BWA-MEM.
"""
import logging
import subprocess
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def run_bwa_index(reference):
    """
    Build BWA index for reference sequence.
    
    Args:
        reference (Path): Path to reference FASTA file
    """
    try:
        cmd = ["bwa", "index", str(reference)]
        subprocess.run(cmd, check=True)
        logging.info(f"BWA index built for {reference}")
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Error building BWA index: {str(e)}")
        raise

def run_bwa_mem(reference, reads, output_bam, threads, min_score):
    """
    Map reads using BWA-MEM.
    
    Args:
        reference (Path): Path to reference FASTA file
        reads (list): List of paths to FASTQ files
        output_bam (Path): Path for output BAM file
        threads (int): Number of threads to use
        min_score (int): Minimum alignment score
    """
    try:
        # Run BWA-MEM
        bwa_cmd = [
            "bwa", "mem",
            "-t", str(threads),
            "-T", str(min_score),
            str(reference)
        ]
        bwa_cmd.extend(str(r) for r in reads)
        
        # Pipe to samtools for BAM conversion and sorting
        samtools_cmd = [
            "samtools", "sort",
            "-@", str(threads),
            "-o", str(output_bam)
        ]
        
        # Run pipeline
        p1 = subprocess.Popen(
            bwa_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        p2 = subprocess.Popen(
            samtools_cmd,
            stdin=p1.stdout,
            stderr=subprocess.PIPE
        )
        
        # Close p1's stdout to allow p2 to receive SIGPIPE
        p1.stdout.close()
        
        # Wait for completion
        p2.communicate()
        
        if p2.returncode != 0:
            raise subprocess.CalledProcessError(p2.returncode, samtools_cmd)
        
        # Index BAM file
        subprocess.run(
            ["samtools", "index", str(output_bam)],
            check=True
        )
        
        logging.info(f"Reads mapped to {reference}")
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Error mapping reads: {str(e)}")
        raise

def main(snakemake):
    """
    Main function to map reads.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Get parameters
    params = snakemake.params[0]
    # Accept either 'unitigs' or 'reference' as input key, fallback to first input if needed
    if hasattr(snakemake.input, 'unitigs'):
        reference = Path(snakemake.input.unitigs)
    elif hasattr(snakemake.input, 'reference'):
        reference = Path(snakemake.input.reference)
    else:
        # fallback: use first input as reference
        reference = Path(snakemake.input[0])

    # Reads: use all other inputs as reads
    reads = [Path(r) for r in snakemake.input if str(r) != str(reference)]
    output_bam = Path(snakemake.output.bam)

    # Ensure output directory exists
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Build BWA index
    run_bwa_index(reference)

    # Map reads
    run_bwa_mem(
        reference=reference,
        reads=reads,
        output_bam=output_bam,
        threads=params['threads'],
        min_score=params['min_score']
    )

    logging.info("Read mapping completed successfully")

if __name__ == "__main__":
    main(snakemake)
