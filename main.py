#!/usr/bin/env python3
"""
CRISPR Guide RNA Counter - Main Entry Point

This script demonstrates the complete pipeline for counting guide RNAs
in Illumina NovaSeq sequencing data using the Aho-Corasick algorithm.

Pipeline:
1. Generate/Load gRNA reference library
2. Generate simulated reads (or load real FASTQ/FASTA)
3. Count gRNAs using Aho-Corasick algorithm
4. Output results table

Usage:
    python main.py                    # Run full demo
    python main.py --help            # Show help
    python main.py --library lib.csv --fastq reads.fq -o counts.csv
"""

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import numpy as np

from grna_library import generate_grna_library, save_library, load_library
from sequence_generator import (
    generate_simulated_reads,
    save_reads_fastq,
    save_reads_fasta
)
from grna_counter import GRNACounter, MismatchTolerantCounter, AhoCorasickExplainer
from grna_qc import (
    validate_library, generate_qc_report, detect_gc_bias,
    normalize_by_gc, GRNAException
)
from question_generator import ask_question, EXAMPLE_QUESTIONS


def run_demo(
    num_guides: int = 500,
    num_reads: int = 100000,
    guides_per_gene: int = 4,
    error_rate: float = 0.001,
    noise_fraction: float = 0.05,
    output_dir: str = "output",
    seed: int = 42,
    max_mismatches: int = 0,
    run_qc: bool = True
) -> pd.DataFrame:
    """
    Run a complete demonstration of the gRNA counting pipeline.

    Args:
        num_guides: Number of guides in the library
        num_reads: Number of reads to generate
        guides_per_gene: Number of guides per gene
        error_rate: Sequencing error rate
        noise_fraction: Fraction of noise reads
        output_dir: Directory for output files
        seed: Random seed
        max_mismatches: Maximum allowed mismatches for fuzzy matching
        run_qc: Whether to run quality control checks

    Returns:
        DataFrame with count results
    """
    print("="*70)
    print("CRISPR GUIDE RNA COUNTER - DEMONSTRATION")
    print("Using Aho-Corasick Algorithm for Efficient Pattern Matching")
    print("="*70)
    print()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    # Step 1: Generate gRNA library
    print("STEP 1: Generating gRNA Reference Library")
    print("-" * 50)
    library = generate_grna_library(
        num_guides=num_guides,
        guides_per_gene=guides_per_gene,
        seed=seed
    )
    library_file = output_path / "grna_library.csv"
    save_library(library, str(library_file))
    print()

    # Step 2: Generate simulated reads
    print("STEP 2: Generating Simulated NovaSeq Reads")
    print("-" * 50)
    reads, true_counts = generate_simulated_reads(
        library,
        num_reads=num_reads,
        error_rate=error_rate,
        noise_fraction=noise_fraction,
        seed=seed
    )

    # Save reads
    fastq_file = output_path / "simulated_reads.fastq"
    save_reads_fastq(reads, str(fastq_file))
    print()

    # Step 2.5: Run Library QC (if enabled)
    if run_qc:
        print("STEP 2.5: Library Quality Control")
        print("-" * 50)
        try:
            lib_qc = validate_library(library)
            print(f"  Valid guides: {lib_qc.valid_guides}/{lib_qc.total_guides}")
            print(f"  GC distribution: {lib_qc.gc_distribution}")
            if lib_qc.warnings:
                for w in lib_qc.warnings[:3]:
                    print(f"  Warning: {w}")
        except GRNAException as e:
            print(f"  QC Error: {e}")
        print()

    # Step 3: Count gRNAs using Aho-Corasick
    print("STEP 3: Counting gRNAs with Aho-Corasick Algorithm")
    print("-" * 50)
    if max_mismatches > 0:
        print(f"Using mismatch-tolerant counting (max mismatches: {max_mismatches})")
        counter = MismatchTolerantCounter(library, max_mismatches=max_mismatches)
    else:
        counter = GRNACounter(library)
    counted = counter.count_reads(reads)
    result_table = counter.get_count_table(counted)
    print()

    # Step 4: Save and display results
    print("STEP 4: Results Summary")
    print("-" * 50)

    # Save count table
    counts_file = output_path / "grna_counts.csv"
    result_table.to_csv(counts_file, index=False)
    print(f"Count table saved to: {counts_file}")

    # Add true counts for comparison (only in demo mode)
    result_table['true_count'] = result_table['sequence'].map(true_counts)
    result_table['count_error'] = result_table['count'] - result_table['true_count']

    # Calculate accuracy metrics
    total_true = sum(true_counts.values())
    total_counted = result_table['count'].sum()

    print(f"\nTotal true gRNA reads: {total_true}")
    print(f"Total counted reads: {total_counted}")
    print(f"Detection rate: {100*total_counted/total_true:.2f}%")

    # Check correlation between true and counted
    correlation = result_table[['count', 'true_count']].corr().iloc[0, 1]
    print(f"Correlation (counted vs true): {correlation:.4f}")

    # Show sample results
    print("\n" + "="*70)
    print("TOP 20 COUNTED gRNAs (OUTPUT TABLE)")
    print("="*70)
    display_cols = ['guide_id', 'gene_name', 'sequence', 'count']
    print(result_table[display_cols].head(20).to_string(index=False))

    # Gene-level summary
    print("\n" + "="*70)
    print("GENE-LEVEL SUMMARY (Top 20 by total counts)")
    print("="*70)
    gene_counts = result_table.groupby('gene_name')['count'].agg(['sum', 'mean', 'std'])
    gene_counts = gene_counts.sort_values('sum', ascending=False)
    gene_counts.columns = ['total_counts', 'mean_per_guide', 'std_per_guide']
    print(gene_counts.head(20).to_string())

    # Save gene-level summary
    gene_file = output_path / "gene_level_counts.csv"
    gene_counts.to_csv(gene_file)
    print(f"\nGene-level summary saved to: {gene_file}")

    # Step 5: Quality Control Report (if enabled)
    if run_qc:
        print("\n" + "="*70)
        print("STEP 5: Quality Control Report")
        print("="*70)
        qc_report = generate_qc_report(library, counted, len(reads))
        print(qc_report)

        # GC bias analysis
        gc_corr, gc_interp = detect_gc_bias(library, counted)
        print(f"\nGC Bias Analysis: {gc_interp}")

        # Save QC report
        qc_file = output_path / "qc_report.txt"
        with open(qc_file, 'w') as f:
            f.write(qc_report)
        print(f"\nQC report saved to: {qc_file}")

    print("\n" + "="*70)
    print("DEMO COMPLETE")
    print("="*70)

    return result_table


def count_from_files(
    library_path: str,
    reads_path: str,
    output_path: str,
    file_format: str = 'fastq'
) -> pd.DataFrame:
    """
    Count gRNAs from existing library and reads files.

    Args:
        library_path: Path to gRNA library CSV
        reads_path: Path to reads file (FASTQ or FASTA)
        output_path: Path for output counts CSV
        file_format: Format of reads file ('fastq' or 'fasta')

    Returns:
        DataFrame with count results
    """
    print("Loading gRNA library...")
    library = load_library(library_path)
    print(f"Loaded {len(library)} guides from {library_path}")

    print(f"\nLoading reads from {reads_path}...")

    # Load reads based on format
    reads = []
    if file_format.lower() == 'fastq':
        with open(reads_path, 'r') as f:
            line_num = 0
            for line in f:
                if line_num % 4 == 1:
                    reads.append(line.strip())
                line_num += 1
    else:  # FASTA
        with open(reads_path, 'r') as f:
            current_seq = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        reads.append(current_seq)
                    current_seq = ""
                else:
                    current_seq += line
            if current_seq:
                reads.append(current_seq)

    print(f"Loaded {len(reads)} reads")

    # Count
    print("\nCounting gRNAs...")
    counter = GRNACounter(library)
    counted = counter.count_reads(reads)
    result_table = counter.get_count_table(counted)

    # Save results
    result_table.to_csv(output_path, index=False)
    print(f"\nResults saved to {output_path}")

    return result_table


def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description='CRISPR Guide RNA Counter using Aho-Corasick Algorithm',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py                                    # Run demo with default settings
  python main.py --demo --num-guides 1000           # Demo with 1000 guides
  python main.py -l lib.csv -r reads.fq -o out.csv  # Count from files

For more information about the algorithm:
  python main.py --explain
        """
    )

    parser.add_argument('--demo', action='store_true', default=True,
                        help='Run demonstration mode (default)')
    parser.add_argument('--explain', action='store_true',
                        help='Explain the Aho-Corasick algorithm')
    parser.add_argument('--ask', type=str, metavar='QUESTION',
                        help='Ask a question about CRISPR screening (uses sub-question methodology)')
    parser.add_argument('--list-questions', action='store_true',
                        help='List example questions you can ask')

    # Demo mode options
    demo_group = parser.add_argument_group('Demo options')
    demo_group.add_argument('--num-guides', type=int, default=500,
                           help='Number of guides in library (default: 500)')
    demo_group.add_argument('--num-reads', type=int, default=100000,
                           help='Number of reads to generate (default: 100000)')
    demo_group.add_argument('--guides-per-gene', type=int, default=4,
                           help='Guides per gene (default: 4)')
    demo_group.add_argument('--error-rate', type=float, default=0.001,
                           help='Sequencing error rate (default: 0.001)')
    demo_group.add_argument('--noise-fraction', type=float, default=0.05,
                           help='Fraction of noise reads (default: 0.05)')
    demo_group.add_argument('--seed', type=int, default=42,
                           help='Random seed (default: 42)')
    demo_group.add_argument('--max-mismatches', type=int, default=0,
                           help='Maximum mismatches for fuzzy matching (default: 0=exact only)')
    demo_group.add_argument('--no-qc', action='store_true',
                           help='Skip quality control checks')
    demo_group.add_argument('--gc-normalize', action='store_true',
                           help='Apply GC content normalization to counts')

    # File mode options
    file_group = parser.add_argument_group('File options')
    file_group.add_argument('-l', '--library', type=str,
                           help='Path to gRNA library CSV file')
    file_group.add_argument('-r', '--reads', type=str,
                           help='Path to reads file (FASTQ/FASTA)')
    file_group.add_argument('-f', '--format', type=str, default='fastq',
                           choices=['fastq', 'fasta'],
                           help='Format of reads file (default: fastq)')
    file_group.add_argument('-o', '--output', type=str, default='grna_counts.csv',
                           help='Output file path (default: grna_counts.csv)')
    file_group.add_argument('--output-dir', type=str, default='output',
                           help='Output directory for demo mode (default: output)')

    args = parser.parse_args()

    # Handle explain mode
    if args.explain:
        AhoCorasickExplainer.explain_algorithm()
        return

    # Handle question mode
    if args.list_questions:
        print("Example questions you can ask with --ask:")
        print("-" * 50)
        for q in EXAMPLE_QUESTIONS:
            print(f"  {q}")
        return

    if args.ask:
        print("Analyzing your question using sub-question methodology...")
        print()
        ask_question(args.ask)
        return

    # Handle file mode
    if args.library and args.reads:
        result = count_from_files(
            library_path=args.library,
            reads_path=args.reads,
            output_path=args.output,
            file_format=args.format
        )
        print("\nTop 10 results:")
        print(result[['guide_id', 'gene_name', 'sequence', 'count']].head(10).to_string())
        return

    # Run demo mode
    run_demo(
        num_guides=args.num_guides,
        num_reads=args.num_reads,
        guides_per_gene=args.guides_per_gene,
        error_rate=args.error_rate,
        noise_fraction=args.noise_fraction,
        output_dir=args.output_dir,
        seed=args.seed,
        max_mismatches=args.max_mismatches,
        run_qc=not args.no_qc
    )


if __name__ == "__main__":
    main()
