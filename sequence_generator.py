"""
DNA Sequence Generator for CRISPR Screening Simulation

This module generates simulated Illumina NovaSeq sequencing reads
for testing the gRNA counting algorithm.

Typical CRISPR screening read structure:
[5' adapter] - [stagger] - [gRNA 20bp] - [scaffold] - [3' adapter]

For simplicity, we simulate:
[flanking 5'] - [gRNA 20bp] - [flanking 3']
"""

import random
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np


NUCLEOTIDES = ['A', 'T', 'G', 'C']

# Typical flanking sequences in CRISPR vectors
# These would be the constant regions flanking the gRNA
FIVE_PRIME_FLANK = "TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG"  # U6 promoter end
THREE_PRIME_FLANK = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCC"  # Scaffold start


def introduce_sequencing_errors(
    sequence: str,
    error_rate: float = 0.001
) -> str:
    """
    Introduce random sequencing errors to simulate Illumina error profile.

    NovaSeq typically has error rates of ~0.1-0.5%

    Args:
        sequence: Original DNA sequence
        error_rate: Probability of error per base (default 0.1%)

    Returns:
        Sequence with introduced errors
    """
    if error_rate == 0:
        return sequence

    sequence_list = list(sequence)

    for i in range(len(sequence_list)):
        if random.random() < error_rate:
            # Substitution error (most common in Illumina)
            original_base = sequence_list[i]
            other_bases = [b for b in NUCLEOTIDES if b != original_base]
            sequence_list[i] = random.choice(other_bases)

    return ''.join(sequence_list)


def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence."""
    return ''.join(random.choices(NUCLEOTIDES, k=length))


def simulate_read(
    grna_sequence: str,
    read_length: int = 150,
    use_standard_flanks: bool = True,
    error_rate: float = 0.001
) -> str:
    """
    Simulate a single sequencing read containing a gRNA.

    Args:
        grna_sequence: The 20bp gRNA sequence
        read_length: Total read length (default 150 for NovaSeq)
        use_standard_flanks: Use realistic vector flanking sequences
        error_rate: Sequencing error rate

    Returns:
        Simulated sequencing read
    """
    grna_len = len(grna_sequence)

    if use_standard_flanks:
        # Use realistic flanking sequences
        five_prime = FIVE_PRIME_FLANK
        three_prime = THREE_PRIME_FLANK
    else:
        # Use random flanking sequences
        flank_len = (read_length - grna_len) // 2
        five_prime = generate_random_sequence(flank_len)
        three_prime = generate_random_sequence(flank_len)

    # Construct the full read
    full_construct = five_prime + grna_sequence + three_prime

    # If construct is longer than read length, take a random window
    if len(full_construct) > read_length:
        # Ensure the gRNA is included in the read
        grna_start = len(five_prime)
        grna_end = grna_start + grna_len

        # Calculate valid start positions that include the gRNA
        max_start = grna_start
        min_start = max(0, grna_end - read_length)

        start_pos = random.randint(min_start, max_start)
        read = full_construct[start_pos:start_pos + read_length]
    else:
        read = full_construct

    # Pad if necessary
    if len(read) < read_length:
        padding_needed = read_length - len(read)
        read = read + generate_random_sequence(padding_needed)

    # Introduce sequencing errors
    read = introduce_sequencing_errors(read, error_rate)

    return read


def generate_simulated_reads(
    grna_library: pd.DataFrame,
    num_reads: int = 100000,
    read_length: int = 150,
    abundance_distribution: str = 'lognormal',
    error_rate: float = 0.001,
    noise_fraction: float = 0.05,
    seed: int = 42
) -> Tuple[List[str], Dict[str, int]]:
    """
    Generate simulated sequencing reads from a gRNA library.

    This simulates a CRISPR screening experiment where:
    - Different gRNAs have different abundances (selection effects)
    - Some reads don't match any gRNA (noise)
    - Sequencing errors may occur

    Args:
        grna_library: DataFrame with 'sequence' column
        num_reads: Total number of reads to generate
        read_length: Length of each read
        abundance_distribution: Distribution of gRNA abundances ('uniform', 'lognormal')
        error_rate: Sequencing error rate
        noise_fraction: Fraction of reads that don't contain library gRNAs
        seed: Random seed for reproducibility

    Returns:
        Tuple of (list of read sequences, dict of true counts per gRNA)
    """
    random.seed(seed)
    np.random.seed(seed)

    sequences = grna_library['sequence'].tolist()
    num_grnas = len(sequences)

    # Calculate number of gRNA-containing vs noise reads
    num_grna_reads = int(num_reads * (1 - noise_fraction))
    num_noise_reads = num_reads - num_grna_reads

    # Generate abundance weights for each gRNA
    if abundance_distribution == 'uniform':
        weights = np.ones(num_grnas)
    elif abundance_distribution == 'lognormal':
        # Lognormal distribution mimics real screening data
        # Some guides are enriched, some depleted
        weights = np.random.lognormal(mean=0, sigma=1.5, size=num_grnas)
    else:
        weights = np.ones(num_grnas)

    # Normalize weights to probabilities
    probabilities = weights / weights.sum()

    # Sample gRNAs according to their abundance
    grna_indices = np.random.choice(
        num_grnas,
        size=num_grna_reads,
        p=probabilities
    )

    # Count true occurrences
    true_counts = {}
    for seq in sequences:
        true_counts[seq] = 0

    reads = []

    print(f"Generating {num_grna_reads} gRNA-containing reads...")
    for idx in grna_indices:
        grna_seq = sequences[idx]
        true_counts[grna_seq] += 1

        # Generate read containing this gRNA
        read = simulate_read(
            grna_seq,
            read_length=read_length,
            use_standard_flanks=True,
            error_rate=error_rate
        )
        reads.append(read)

    # Generate noise reads (random sequences not in library)
    print(f"Generating {num_noise_reads} noise reads...")
    for _ in range(num_noise_reads):
        noise_read = generate_random_sequence(read_length)
        reads.append(noise_read)

    # Shuffle reads
    random.shuffle(reads)

    print(f"Total reads generated: {len(reads)}")

    return reads, true_counts


def save_reads_fastq(reads: List[str], filepath: str, quality_score: int = 37) -> None:
    """
    Save reads in FASTQ format.

    Args:
        reads: List of read sequences
        filepath: Output file path
        quality_score: Phred quality score (default 37 = 'F')
    """
    quality_char = chr(quality_score + 33)

    with open(filepath, 'w') as f:
        for i, read in enumerate(reads):
            f.write(f"@READ_{i+1}\n")
            f.write(f"{read}\n")
            f.write("+\n")
            f.write(f"{quality_char * len(read)}\n")

    print(f"Saved {len(reads)} reads to {filepath}")


def save_reads_fasta(reads: List[str], filepath: str) -> None:
    """
    Save reads in FASTA format.

    Args:
        reads: List of read sequences
        filepath: Output file path
    """
    with open(filepath, 'w') as f:
        for i, read in enumerate(reads):
            f.write(f">READ_{i+1}\n")
            f.write(f"{read}\n")

    print(f"Saved {len(reads)} reads to {filepath}")


if __name__ == "__main__":
    from grna_library import generate_grna_library

    # Generate a test library
    library = generate_grna_library(num_guides=100, seed=42)

    # Generate simulated reads
    reads, true_counts = generate_simulated_reads(
        library,
        num_reads=10000,
        seed=42
    )

    print(f"\nGenerated {len(reads)} reads")
    print(f"Sample read: {reads[0][:50]}...")

    # Show some true counts
    print("\nTop 10 most abundant gRNAs (true counts):")
    sorted_counts = sorted(true_counts.items(), key=lambda x: x[1], reverse=True)
    for seq, count in sorted_counts[:10]:
        print(f"  {seq}: {count}")
