
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/anand/miniconda3/envs/metagenome_quant/lib/python3.9/site-packages', '/home/anand/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpk_2mes2_/file/mnt/c/Users/anand/Documents/Academics/S4/22BIO211- INTELLIGENCE OF BIOLOGICAL SYSTEMS II/Project/metagenome_quantification/scripts', '/mnt/c/Users/anand/Documents/Academics/S4/22BIO211- INTELLIGENCE OF BIOLOGICAL SYSTEMS II/Project/metagenome_quantification/scripts']); import pickle; snakemake = pickle.loads(b"\x80\x04\x95\xa9\x07\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c\x0cdata/genomes\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x07genomes\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x1eresults/twopaco/ecoli_cdbg.gfa\x94\x8c\x1dresults/twopaco/unitigs.fasta\x94e}\x94(h\x0c}\x94(\x8c\x05graph\x94K\x00N\x86\x94\x8c\x07unitigs\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh*h&h,h'ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(\x8c\x01k\x94Ke\x8c\tmemory_gb\x94K\x10\x8c\x07threads\x94K\x04ua}\x94(h\x0c}\x94\x8c\x07twopaco\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhAh;ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhcK\x01heK\x01hgh`ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x07genomes\x94]\x94(\x8c\x0bcollections\x94\x8c\x0bOrderedDict\x94\x93\x94)R\x94(\x8c\taccession\x94\x8c\x0bNC_000913.3\x94\x8c\x04name\x94\x8c\x10ecoli_k12_mg1655\x94\x8c\tabundance\x94G?\xb9\x99\x99\x99\x99\x99\x9auh\x87)R\x94(\x8c\taccession\x94\x8c\x08CP053602\x94\x8c\x04name\x94\x8c\x0eecoli_bl21_de3\x94\x8c\tabundance\x94G?\xc9\x99\x99\x99\x99\x99\x9auh\x87)R\x94(\x8c\taccession\x94\x8c\x0bNC_012967.1\x94\x8c\x04name\x94\x8c\x0cecoli_rel606\x94\x8c\tabundance\x94G?\xd3333333uh\x87)R\x94(\x8c\taccession\x94\x8c\x08CP009072\x94\x8c\x04name\x94\x8c\x0fecoli_atcc25922\x94\x8c\tabundance\x94G?\xd9\x99\x99\x99\x99\x99\x9aue\x8c\x05reads\x94}\x94(\x8c\x0btotal_pairs\x94M\xe8\x03\x8c\x0bread_length\x94K2\x8c\rfragment_size\x94K\xc8\x8c\nerror_rate\x94G?\x84z\xe1G\xae\x14{u\x8c\x07twopaco\x94h;\x8c\x03bwa\x94}\x94(\x8c\x07threads\x94K\x04\x8c\tmin_score\x94K\x1eu\x8c\x05blast\x94}\x94(\x8c\rperc_identity\x94Ka\x8c\x14min_raw_gapped_score\x94Kd\x8c\x07threads\x94K\x04u\x8c\x0boutput_dirs\x94}\x94(\x8c\x07genomes\x94\x8c\x0cdata/genomes\x94\x8c\x05reads\x94\x8c\ndata/reads\x94\x8c\x07twopaco\x94\x8c\x0fresults/twopaco\x94\x8c\x07abruijn\x94\x8c\x0fresults/abruijn\x94\x8c\nevaluation\x94\x8c\x12results/evaluation\x94\x8c\x05plots\x94\x8c\rresults/plots\x94u\x8c\x07threads\x94K\x04\x8c\tmemory_gb\x94K\x10u\x8c\x04rule\x94\x8c\rbuild_twopaco\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c\x83/mnt/c/Users/anand/Documents/Academics/S4/22BIO211- INTELLIGENCE OF BIOLOGICAL SYSTEMS II/Project/metagenome_quantification/scripts\x94ub."); from snakemake.logging import logger; logger.printshellcmds = True; __real_file__ = __file__; __file__ = '/mnt/c/Users/anand/Documents/Academics/S4/22BIO211- INTELLIGENCE OF BIOLOGICAL SYSTEMS II/Project/metagenome_quantification/scripts/build_twopaco.py';
######## snakemake preamble end #########
#!/usr/bin/env python3
"""
Build colored de Bruijn graph using TwoPaCo.
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
    """
    Run TwoPaCo to build a colored de Bruijn graph.

    Args:
        input_dir (Path): Directory containing input genome FASTA files.
        output_prefix (Path): Prefix for output files (without extension).
        k (int): k-mer size.
        memory_gb (int): Memory limit in GB for the Bloom filter.
        threads (int): Number of threads to use.
    """
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
    # Add all .fna input files
    input_files = list(input_dir.glob("*.fna"))
    cmd_build += [str(f) for f in input_files]

    logging.info(f"Running TwoPaCo: {' '.join(cmd_build)}")
    subprocess.run(cmd_build, check=True)

    # Generate GFA output
    cmd_gfa = [
        "graphdump",
        str(output_prefix) + ".bin",
        "-f", "gfa1",
        "-k", str(k)
    ]
    for f in input_files:
        cmd_gfa += ["-s", str(f)]
    logging.info(f"Generating GFA: {' '.join(cmd_gfa)}")
    with open(str(output_prefix) + ".gfa", "w") as gf:
        subprocess.run(cmd_gfa, check=True, stdout=gf)

    # Generate FASTA unitigs
    cmd_fasta = [
        "graphdump",
        str(output_prefix) + ".bin",
        "-f", "fasta",
        "-k", str(k)
    ]
    for f in input_files:
        cmd_fasta += ["-s", str(f)]
    logging.info(f"Generating unitigs FASTA: {' '.join(cmd_fasta)}")
    with open(str(output_prefix) + "_unitigs.fasta", "w") as uf:
        subprocess.run(cmd_fasta, check=True, stdout=uf)

    # Cleanup temporary files
    for temp_file in tmp_dir.iterdir():
        temp_file.unlink(missing_ok=True)
    tmp_dir.rmdir()


def main(snakemake):
    """
    Main entry point for Snakemake.

    Expects:
        - snakemake.params.twopaco to be a dict with keys 'k', 'memory_gb', 'threads'
        - snakemake.input.genomes to be the path to the genomes directory
        - snakemake.output.graph to be the desired .gfa file path
        - snakemake.output.unitigs to be the desired unitigs FASTA file path
    """
    params = snakemake.params.twopaco
    input_dir = Path(snakemake.input.genomes)
    # Remove suffix from graph output to get prefix
    output_prefix = Path(snakemake.output.graph).with_suffix("")

    run_twopaco(
        input_dir=input_dir,
        output_prefix=output_prefix,
        k=params["k"],
        memory_gb=params["memory_gb"],
        threads=params["threads"]
    )

    # Move/rename outputs (if necessary; Snakemake will already point to correct files)
    # But since we wrote directly to those paths, no further action is needed here.

    logging.info("TwoPaCo graph construction completed successfully")


if __name__ == "__main__":
    main(snakemake)
