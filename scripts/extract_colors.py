#!/usr/bin/env python3
"""
Extract color information from TwoPaCo GFA output.
"""
import json
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def parse_gfa_colors(gfa_path):
    """
    Parse color (genome presence) information from GFA file.
    
    Args:
        gfa_path (Path): Path to GFA file
    
    Returns:
        dict: Mapping of unitig IDs to list of genome indices
    """
    colors = {}
    with open(gfa_path) as f:
        for line in f:
            if line.startswith('S'):  # Segment line
                fields = line.strip().split('\t')
                unitig_id = fields[1]
                
                # Find KC (color) tag
                for tag in fields[3:]:
                    if tag.startswith('KC:'):
                        # Parse color string
                        color_str = tag.split(':')[2]
                        # Convert to list of genome indices (1-based)
                        genome_indices = [
                            i + 1 for i, c in enumerate(color_str)
                            if c == '1'
                        ]
                        colors[unitig_id] = genome_indices
                        break
    
    return colors

def save_colors(colors, output_path):
    """
    Save color information to JSON file.
    
    Args:
        colors (dict): Mapping of unitig IDs to genome indices
        output_path (Path): Path to output JSON file
    """
    with open(output_path, 'w') as f:
        json.dump(colors, f, indent=2)

def main(snakemake):
    """
    Main function to extract color information.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Parse GFA file
    colors = parse_gfa_colors(Path(snakemake.input.graph))
    
    # Save results
    save_colors(colors, Path(snakemake.output.colors))
    
    logging.info(f"Extracted colors for {len(colors)} unitigs")

if __name__ == "__main__":
    main(snakemake)
