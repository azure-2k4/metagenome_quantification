#!/usr/bin/env python3
"""
Build a de Bruijn graph from genome FASTA files and output GFA and consensus.
"""
import os
from pathlib import Path
from collections import defaultdict, Counter
def read_fasta(file_path):
    seq = ""
    with open(file_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq += line.strip()
    return seq

def build_debruijn_graph(sequences, k):
    edges = defaultdict(list)
    nodes = set()
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer1 = seq[i:i+k-1]
            kmer2 = seq[i+1:i+k]
            edges[kmer1].append(kmer2)
            nodes.add(kmer1)
            nodes.add(kmer2)
    return edges, nodes

def write_gfa(edges, nodes, gfa_path):
    with open(gfa_path, 'w') as f:
        for node in nodes:
            f.write(f"S\t{node}\t{node}\n")
        for src, dsts in edges.items():
            for dst in dsts:
                f.write(f"L\t{src}\t+\t{dst}\t+\t0M\n")

def get_consensus_path(edges, nodes):
    # Simple greedy walk for consensus (not optimal for complex graphs)
    in_deg = Counter()
    out_deg = Counter()
    for src, dsts in edges.items():
        out_deg[src] += len(dsts)
        for dst in dsts:
            in_deg[dst] += 1
    start_nodes = [n for n in nodes if in_deg[n] == 0]
    if not start_nodes:
        start_nodes = list(nodes)
    path = []
    stack = [start_nodes[0]]
    local_edges = {k: v.copy() for k, v in edges.items()}  # avoid mutating original
    while stack:
        node = stack[-1]
        if node in local_edges and local_edges[node]:
            nxt = local_edges[node].pop()
            stack.append(nxt)
        else:
            path.append(stack.pop())
    return path[::-1]

def write_consensus(path, fasta_path):
    if not path:
        return
    seq = path[0]
    for kmer in path[1:]:
        seq += kmer[-1]
    with open(fasta_path, 'w') as f:
        f.write(">consensus\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def main(snakemake):
    k = snakemake.params['perc_identity'] if 'perc_identity' in snakemake.params else 31  # fallback k
    input_dir = Path(snakemake.input.genomes)
    gfa_path = snakemake.output.graph
    fasta_path = snakemake.output.consensus
    # Use k=31 by default if not specified
    try:
        k = int(k)
    except Exception:
        k = 31
    # Read all genome sequences
    sequences = [read_fasta(f) for f in input_dir.glob("*.fna")]
    edges, nodes = build_debruijn_graph(sequences, k)
    write_gfa(edges, nodes, gfa_path)
    path = get_consensus_path(edges, nodes)
    write_consensus(path, fasta_path)
    print("A-Bruijn (de Bruijn) graph construction complete.")

if __name__ == "__main__":
    main(snakemake)
