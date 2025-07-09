#!/usr/bin/env python3
"""
Build colored de Bruijn graph using TwoPaCo and print graphdump output to stdout for debugging.
"""
import logging
import subprocess
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def run_twopaco(input_dir: Path, output_prefix: Path, k: int, memory_gb: int, threads: int):
    tmp_dir = output_prefix.parent / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Build TwoPaCo command
    cmd_build = [
        "twopaco",
        "--filtermemory", str(memory_gb),
        "-k", str(k),
        "--tmpdir", str(tmp_dir),
        "-o", str(output_prefix) + ".bin",
        "--threads", str(threads)
    ]
    input_files = list(input_dir.glob("*.fna"))
    logging.info(f"Input directory: {input_dir.resolve()}")
    logging.info(f"Input genome files found: {[str(f) for f in input_files]}")
    cmd_build += [str(f) for f in input_files]

    logging.info(f"Running TwoPaCo: {' '.join(cmd_build)}")
    subprocess.run(cmd_build, check=True)


    # Generate GFA output and check correctness
    cmd_gfa = [
        "graphdump",
        str(output_prefix) + ".bin",
        "-f", "gfa1",
        "-k", str(k)
    ]
    for f in input_files:
        cmd_gfa += ["-s", str(f)]
    logging.info(f"Running graphdump for GFA: {' '.join(cmd_gfa)}")
    try:
        result = subprocess.run(cmd_gfa, check=True, capture_output=True, text=True)
        gfa_lines = result.stdout.splitlines()
        # Print first 40 lines for inspection
        print("\n--- graphdump GFA output (first 40 lines) ---")
        for line in gfa_lines[:40]:
            print(line)
        print("--- end of preview ---\n")
        # Write full GFA to file
        with open(str(output_prefix) + ".gfa", "w") as gf:
            gf.write(result.stdout)
        # Check GFA correctness: at least one S and one L line
        s_count = sum(1 for l in gfa_lines if l.startswith('S'))
        l_count = sum(1 for l in gfa_lines if l.startswith('L'))
        if s_count == 0 or l_count == 0:
            logging.error("GFA file may be incorrect: missing S or L lines.")
        else:
            logging.info(f"GFA file looks correct: {s_count} segments, {l_count} links.")
    except Exception as e:
        logging.error(f"graphdump failed for GFA: {e}")

    # Generate FASTA unitigs and check correctness
    cmd_fasta = [
        "graphdump",
        str(output_prefix) + ".bin",
        "-f", "fasta",
        "-k", str(k)
    ]
    for f in input_files:
        cmd_fasta += ["-s", str(f)]
    logging.info(f"Generating unitigs FASTA: {' '.join(cmd_fasta)}")
    try:
        result_fasta = subprocess.run(cmd_fasta, check=True, capture_output=True, text=True)
        fasta_lines = result_fasta.stdout.splitlines()
        with open(str(output_prefix) + "_unitigs.fasta", "w") as uf:
            uf.write(result_fasta.stdout)
        # Check FASTA correctness: at least one header and one sequence
        header_count = sum(1 for l in fasta_lines if l.startswith('>'))
        seq_count = sum(1 for l in fasta_lines if l and not l.startswith('>'))
        if header_count == 0 or seq_count == 0:
            logging.error("FASTA file may be incorrect: missing headers or sequences.")
        else:
            logging.info(f"FASTA file looks correct: {header_count} unitigs.")
    except Exception as e:
        logging.error(f"graphdump failed for FASTA: {e}")

    # Visualize the GFA with Bandage (if installed)
    try:
        import shutil
        bandage_path = shutil.which("Bandage")
        if bandage_path:
            bandage_img = str(output_prefix) + "_bandage.png"
            logging.info(f"Visualizing GFA with Bandage: Bandage image will be saved as {bandage_img}")
            subprocess.run([bandage_path, "image", str(output_prefix) + ".gfa", bandage_img], check=True)
        else:
            logging.warning("Bandage is not installed or not in PATH. Skipping visualization.")
    except Exception as e:
        logging.error(f"Bandage visualization failed: {e}")

    # Cleanup temporary files
    for temp_file in tmp_dir.iterdir():
        temp_file.unlink(missing_ok=True)
    tmp_dir.rmdir()


def main(snakemake):
    params = snakemake.params.twopaco
    input_dir = Path(snakemake.input.genomes)
    output_prefix = Path(snakemake.output.graph).with_suffix("")

    run_twopaco(
        input_dir=input_dir,
        output_prefix=output_prefix,
        k=params["k"],
        memory_gb=params["memory_gb"],
        threads=params["threads"]
    )

    logging.info("TwoPaCo graph construction completed successfully")


if __name__ == "__main__":
    main(snakemake)
