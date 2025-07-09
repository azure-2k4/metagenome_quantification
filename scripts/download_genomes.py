#!/usr/bin/env python3
"""
Download reference genomes from NCBI for metagenomic analysis.
"""
import os
import logging
from Bio import Entrez, SeqIO
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def setup_email():
    """Set up email for NCBI Entrez."""
    Entrez.email = "workspace.aanand@gmail.com"  # email

def download_genome(accession, output_dir):
    """
    Download a genome from NCBI using its accession number.
    
    Args:
        accession (str): NCBI accession number
        output_dir (Path): Directory to save the genome
    
    Returns:
        Path: Path to downloaded genome file
    """
    try:
        logging.info(f"Downloading genome: {accession}")
        
        # Download genome from NCBI
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta",
            retmode="text"
        )
        
        # Save genome to file
        output_file = output_dir / f"{accession}.fna"
        with open(output_file, 'w') as f:
            f.write(handle.read())
        
        return output_file
    
    except Exception as e:
        logging.error(f"Error downloading {accession}: {str(e)}")
        raise

def main(snakemake):
    """
    Main function to download all reference genomes.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Create output directory
    output_dir = Path(snakemake.output[0])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up NCBI access
    setup_email()
    
    # Download each genome
    for genome in snakemake.params.accessions:
        accession = genome['accession']
        download_genome(accession, output_dir)
    
    logging.info("All genomes downloaded successfully")

if __name__ == "__main__":
    main(snakemake)
