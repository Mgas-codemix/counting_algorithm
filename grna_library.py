"""
CRISPR Guide RNA Reference Library Generator

This module generates a dummy gRNA reference library for testing purposes.
In real CRISPR screening experiments, this library would come from the
actual gRNA library design (e.g., Brunello, GeCKO, TKOv3).

Typical gRNA structure:
- 20 nucleotides targeting sequence
- Associated with target gene name
- May have additional metadata (PAM sequence, etc.)
"""

import random
from typing import List, Tuple, Optional
import pandas as pd


# Common gene names for CRISPR screening (oncogenes, tumor suppressors, etc.)
GENE_NAMES = [
    "TP53", "BRCA1", "BRCA2", "KRAS", "EGFR", "MYC", "PTEN", "RB1",
    "APC", "PIK3CA", "BRAF", "CDKN2A", "NRAS", "HRAS", "ATM", "VHL",
    "NF1", "NF2", "SMAD4", "STK11", "CDH1", "MLH1", "MSH2", "MSH6",
    "ARID1A", "CREBBP", "EP300", "KMT2A", "KMT2D", "SETD2", "EZH2",
    "BAP1", "PBRM1", "SMARCA4", "SMARCD1", "IDH1", "IDH2", "DNMT3A",
    "TET2", "ASXL1", "WT1", "NPM1", "FLT3", "KIT", "JAK2", "NOTCH1",
    "FBXW7", "RNF43", "ZNRF3", "CTNNB1", "AKT1", "MTOR", "TSC1", "TSC2",
    "KEAP1", "NFE2L2", "CUL3", "SOX9", "GATA3", "FOXA1", "ESR1", "AR",
    "SPOP", "NCOR1", "ATRX", "DAXX", "MEN1", "PDGFRA", "CDK4", "CDK6",
    "CCND1", "CCNE1", "RET", "ALK", "ROS1", "MET", "FGFR1", "FGFR2",
    "FGFR3", "ERBB2", "ERBB3", "MAP2K1", "MAP2K2", "RAF1", "ARAF",
    "POLE", "POLD1", "BRIP1", "PALB2", "RAD51C", "RAD51D", "CHEK2",
    "NBN", "MRE11", "RAD50", "BLM", "FANCA", "FANCC", "FANCD2", "FANCE"
]

NUCLEOTIDES = ['A', 'T', 'G', 'C']


def generate_random_grna_sequence(length: int = 20, seed: Optional[int] = None) -> str:
    """
    Generate a random gRNA sequence.

    Args:
        length: Length of the gRNA sequence (default 20 for SpCas9)
        seed: Random seed for reproducibility

    Returns:
        Random DNA sequence string
    """
    if seed is not None:
        random.seed(seed)
    return ''.join(random.choices(NUCLEOTIDES, k=length))


def generate_grna_library(
    num_guides: int = 1000,
    guides_per_gene: int = 4,
    grna_length: int = 20,
    seed: int = 42
) -> pd.DataFrame:
    """
    Generate a dummy gRNA reference library.

    In a real CRISPR screen, you would have:
    - Multiple guides per gene (typically 4-10)
    - Unique 20bp targeting sequences
    - Control guides (non-targeting, safe-targeting)

    Args:
        num_guides: Total number of guides to generate
        guides_per_gene: Number of guides per gene
        grna_length: Length of each gRNA sequence
        seed: Random seed for reproducibility

    Returns:
        DataFrame with columns: guide_id, gene_name, sequence
    """
    random.seed(seed)

    library_data = []
    sequences_used = set()  # Ensure unique sequences

    # Calculate how many genes we need
    num_genes = num_guides // guides_per_gene

    # Select genes (cycle through if we need more than available)
    genes = []
    while len(genes) < num_genes:
        remaining = num_genes - len(genes)
        genes.extend(GENE_NAMES[:min(remaining, len(GENE_NAMES))])

    guide_counter = 1

    for gene in genes:
        for guide_num in range(1, guides_per_gene + 1):
            # Generate unique sequence
            while True:
                sequence = generate_random_grna_sequence(grna_length)
                if sequence not in sequences_used:
                    sequences_used.add(sequence)
                    break

            guide_id = f"sg{gene}_{guide_num}"

            library_data.append({
                'guide_id': guide_id,
                'gene_name': gene,
                'sequence': sequence
            })

            guide_counter += 1

            if guide_counter > num_guides:
                break

        if guide_counter > num_guides:
            break

    # Add non-targeting control guides
    num_controls = max(10, num_guides // 100)
    for i in range(1, num_controls + 1):
        while True:
            sequence = generate_random_grna_sequence(grna_length)
            if sequence not in sequences_used:
                sequences_used.add(sequence)
                break

        library_data.append({
            'guide_id': f"sgNTC_{i}",
            'gene_name': 'Non-targeting',
            'sequence': sequence
        })

    return pd.DataFrame(library_data)


def save_library(library: pd.DataFrame, filepath: str) -> None:
    """
    Save the gRNA library to a CSV file.

    Args:
        library: DataFrame containing the gRNA library
        filepath: Output file path
    """
    library.to_csv(filepath, index=False)
    print(f"Library saved to {filepath}")
    print(f"Total guides: {len(library)}")
    print(f"Unique genes: {library['gene_name'].nunique()}")


def load_library(filepath: str) -> pd.DataFrame:
    """
    Load a gRNA library from a CSV file.

    Args:
        filepath: Input file path

    Returns:
        DataFrame containing the gRNA library
    """
    return pd.read_csv(filepath)


if __name__ == "__main__":
    # Generate a test library
    library = generate_grna_library(num_guides=400, guides_per_gene=4, seed=42)
    print("\nSample of generated library:")
    print(library.head(10))
    print(f"\nTotal guides: {len(library)}")
    print(f"Unique genes: {library['gene_name'].nunique()}")
