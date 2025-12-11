"""
Quality Control and Bias Correction Module for CRISPR gRNA Counting

This module addresses common issues in CRISPR screening experiments:

1. EXCEPTIONS/EDGE CASES:
   - Empty or malformed input files
   - Invalid gRNA sequences (non-ATGC characters)
   - Duplicate guide sequences
   - Low-quality reads
   - Missing required columns

2. BIASES IN COUNTING:
   - GC content bias (guides with extreme GC are harder to amplify)
   - Position bias (guides at certain positions may have systematic effects)
   - Length bias (non-standard guide lengths)
   - PCR amplification bias

3. EXACT MATCHING IMPROVEMENTS:
   - Mismatch tolerance for sequencing errors
   - Hamming distance-based matching
   - Quality score-aware matching

Methodology:
-----------
For each issue type, we generate sub-questions to fully understand and address it:
- What causes this issue?
- How does it affect results?
- How can we detect it?
- How can we correct it?
"""

import re
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from collections import Counter
import pandas as pd
import numpy as np


# ============================================================================
# EXCEPTION CLASSES FOR CRISPR EXPERIMENTS
# ============================================================================

class GRNAException(Exception):
    """Base exception for gRNA-related errors."""
    pass


class InvalidSequenceError(GRNAException):
    """Raised when a sequence contains invalid characters."""
    def __init__(self, sequence: str, position: int = None, invalid_char: str = None):
        self.sequence = sequence
        self.position = position
        self.invalid_char = invalid_char
        msg = f"Invalid sequence: {sequence[:50]}..."
        if invalid_char and position is not None:
            msg += f" (character '{invalid_char}' at position {position})"
        super().__init__(msg)


class DuplicateGuideError(GRNAException):
    """Raised when duplicate guide sequences are detected."""
    def __init__(self, sequence: str, count: int):
        self.sequence = sequence
        self.count = count
        super().__init__(f"Duplicate guide sequence found {count} times: {sequence}")


class EmptyLibraryError(GRNAException):
    """Raised when the gRNA library is empty."""
    pass


class EmptyReadsError(GRNAException):
    """Raised when no reads are provided."""
    pass


class LowMappingRateError(GRNAException):
    """Raised when mapping rate is critically low."""
    def __init__(self, rate: float, threshold: float = 0.1):
        self.rate = rate
        self.threshold = threshold
        super().__init__(f"Mapping rate {rate:.2%} is below threshold {threshold:.2%}")


class MissingColumnError(GRNAException):
    """Raised when required column is missing from library."""
    def __init__(self, column: str, available: List[str]):
        self.column = column
        self.available = available
        super().__init__(f"Missing required column '{column}'. Available: {available}")


# ============================================================================
# DATA CLASSES FOR QC METRICS
# ============================================================================

@dataclass
class SequenceQC:
    """Quality control results for a single sequence."""
    sequence: str
    is_valid: bool
    gc_content: float
    length: int
    has_homopolymer: bool
    homopolymer_length: int
    issues: List[str]


@dataclass
class LibraryQC:
    """Quality control results for entire library."""
    total_guides: int
    valid_guides: int
    duplicate_count: int
    gc_distribution: Dict[str, int]  # bins
    length_distribution: Dict[int, int]
    homopolymer_guides: int
    issues: List[str]
    warnings: List[str]


@dataclass
class CountingQC:
    """Quality control results for counting operation."""
    total_reads: int
    mapped_reads: int
    unmapped_reads: int
    multi_mapped_reads: int
    mapping_rate: float
    guides_with_counts: int
    guides_without_counts: int
    zero_count_fraction: float
    gini_coefficient: float  # Measure of count inequality
    issues: List[str]
    warnings: List[str]


@dataclass
class BiasMetrics:
    """Metrics for detecting various biases."""
    gc_bias_correlation: float
    length_bias_present: bool
    position_bias_detected: bool
    extreme_outliers: List[str]  # guide IDs
    recommendations: List[str]


# ============================================================================
# SEQUENCE VALIDATION
# ============================================================================

def validate_sequence(sequence: str, allow_n: bool = False) -> SequenceQC:
    """
    Validate a DNA sequence and compute QC metrics.

    Sub-questions addressed:
    1. Is the sequence valid DNA? (only ATGC or ATGCN)
    2. What is the GC content?
    3. Are there problematic homopolymers?
    4. Is the length appropriate for gRNA?

    Args:
        sequence: DNA sequence string
        allow_n: Whether to allow N (ambiguous) bases

    Returns:
        SequenceQC object with validation results
    """
    issues = []
    is_valid = True

    # Check for valid characters
    valid_chars = set('ATGCatgc')
    if allow_n:
        valid_chars.update('Nn')

    sequence_upper = sequence.upper()

    for i, char in enumerate(sequence_upper):
        if char not in 'ATGCN':
            is_valid = False
            issues.append(f"Invalid character '{char}' at position {i}")
            break

    # Calculate GC content
    gc_count = sequence_upper.count('G') + sequence_upper.count('C')
    gc_content = gc_count / len(sequence) if sequence else 0

    # Check for extreme GC content
    if gc_content < 0.2:
        issues.append(f"Very low GC content ({gc_content:.1%})")
    elif gc_content > 0.8:
        issues.append(f"Very high GC content ({gc_content:.1%})")

    # Check length
    length = len(sequence)
    if length < 17:
        issues.append(f"Sequence too short ({length}bp, expected 17-23bp)")
    elif length > 25:
        issues.append(f"Sequence too long ({length}bp, expected 17-23bp)")

    # Check for homopolymers (4+ of same base)
    homopolymer_match = re.search(r'([ATGC])\1{3,}', sequence_upper)
    has_homopolymer = homopolymer_match is not None
    homopolymer_length = len(homopolymer_match.group()) if homopolymer_match else 0

    if homopolymer_length >= 5:
        issues.append(f"Long homopolymer run ({homopolymer_length}bp)")

    return SequenceQC(
        sequence=sequence,
        is_valid=is_valid,
        gc_content=gc_content,
        length=length,
        has_homopolymer=has_homopolymer,
        homopolymer_length=homopolymer_length,
        issues=issues
    )


def validate_library(library: pd.DataFrame, sequence_col: str = 'sequence') -> LibraryQC:
    """
    Validate an entire gRNA library.

    Sub-questions addressed:
    1. Are all sequences valid?
    2. Are there duplicates?
    3. What is the GC distribution?
    4. What is the length distribution?

    Args:
        library: DataFrame with gRNA library
        sequence_col: Column containing sequences

    Returns:
        LibraryQC object with validation results
    """
    issues = []
    warnings = []

    # Check for empty library
    if library.empty:
        raise EmptyLibraryError("Library DataFrame is empty")

    # Check for required column
    if sequence_col not in library.columns:
        raise MissingColumnError(sequence_col, library.columns.tolist())

    sequences = library[sequence_col].tolist()
    total_guides = len(sequences)

    # Check for duplicates
    sequence_counts = Counter(sequences)
    duplicates = {seq: count for seq, count in sequence_counts.items() if count > 1}
    duplicate_count = len(duplicates)

    if duplicate_count > 0:
        issues.append(f"Found {duplicate_count} duplicate sequences")
        for seq, count in list(duplicates.items())[:5]:
            warnings.append(f"  Duplicate: {seq[:20]}... ({count}x)")

    # Validate individual sequences
    valid_count = 0
    gc_values = []
    lengths = []
    homopolymer_count = 0

    for seq in sequences:
        qc = validate_sequence(seq)
        if qc.is_valid:
            valid_count += 1
        gc_values.append(qc.gc_content)
        lengths.append(qc.length)
        if qc.has_homopolymer and qc.homopolymer_length >= 4:
            homopolymer_count += 1

    # GC distribution
    gc_bins = {'low (<30%)': 0, 'normal (30-70%)': 0, 'high (>70%)': 0}
    for gc in gc_values:
        if gc < 0.3:
            gc_bins['low (<30%)'] += 1
        elif gc > 0.7:
            gc_bins['high (>70%)'] += 1
        else:
            gc_bins['normal (30-70%)'] += 1

    # Length distribution
    length_dist = Counter(lengths)

    # Add warnings
    invalid_count = total_guides - valid_count
    if invalid_count > 0:
        issues.append(f"{invalid_count} sequences have invalid characters")

    if gc_bins['low (<30%)'] > total_guides * 0.1:
        warnings.append(f"{gc_bins['low (<30%)']} guides have low GC content (<30%)")
    if gc_bins['high (>70%)'] > total_guides * 0.1:
        warnings.append(f"{gc_bins['high (>70%)']} guides have high GC content (>70%)")

    if homopolymer_count > total_guides * 0.05:
        warnings.append(f"{homopolymer_count} guides contain homopolymer runs (4+bp)")

    return LibraryQC(
        total_guides=total_guides,
        valid_guides=valid_count,
        duplicate_count=duplicate_count,
        gc_distribution=gc_bins,
        length_distribution=dict(length_dist),
        homopolymer_guides=homopolymer_count,
        issues=issues,
        warnings=warnings
    )


# ============================================================================
# COUNTING QC AND BIAS DETECTION
# ============================================================================

def calculate_gini_coefficient(values: List[int]) -> float:
    """
    Calculate Gini coefficient to measure inequality in counts.

    A high Gini coefficient (>0.8) indicates very unequal distribution,
    which could indicate strong selection or technical bias.

    Args:
        values: List of count values

    Returns:
        Gini coefficient (0 = equal, 1 = maximally unequal)
    """
    if not values or all(v == 0 for v in values):
        return 0.0

    values = sorted(values)
    n = len(values)
    cumulative = np.cumsum(values)

    gini = (2 * np.sum((np.arange(1, n + 1) * values))) / (n * cumulative[-1]) - (n + 1) / n
    return max(0, min(1, gini))


def evaluate_counting_results(
    counts: Dict[int, int],
    total_reads: int,
    library_size: int
) -> CountingQC:
    """
    Evaluate counting results for quality issues.

    Sub-questions addressed:
    1. What is the mapping rate?
    2. How many guides have zero counts?
    3. Is the count distribution reasonable?
    4. Are there concerning patterns?

    Args:
        counts: Dictionary of guide_idx -> count
        total_reads: Total number of reads processed
        library_size: Number of guides in library

    Returns:
        CountingQC object with evaluation results
    """
    issues = []
    warnings = []

    # Basic statistics
    count_values = list(counts.values())
    mapped_reads = sum(count_values)
    unmapped_reads = total_reads - mapped_reads
    mapping_rate = mapped_reads / total_reads if total_reads > 0 else 0

    guides_with_counts = sum(1 for c in count_values if c > 0)
    guides_without_counts = library_size - guides_with_counts
    zero_count_fraction = guides_without_counts / library_size if library_size > 0 else 0

    # Gini coefficient
    all_counts = [counts.get(i, 0) for i in range(library_size)]
    gini = calculate_gini_coefficient(all_counts)

    # Evaluate issues
    if mapping_rate < 0.1:
        issues.append(f"Critical: Very low mapping rate ({mapping_rate:.1%})")
    elif mapping_rate < 0.5:
        warnings.append(f"Low mapping rate ({mapping_rate:.1%})")

    if zero_count_fraction > 0.5:
        issues.append(f"High dropout: {zero_count_fraction:.1%} guides have zero counts")
    elif zero_count_fraction > 0.2:
        warnings.append(f"{zero_count_fraction:.1%} guides have zero counts")

    if gini > 0.9:
        warnings.append(f"Very unequal count distribution (Gini={gini:.3f})")

    return CountingQC(
        total_reads=total_reads,
        mapped_reads=mapped_reads,
        unmapped_reads=unmapped_reads,
        multi_mapped_reads=0,  # Would need additional tracking
        mapping_rate=mapping_rate,
        guides_with_counts=guides_with_counts,
        guides_without_counts=guides_without_counts,
        zero_count_fraction=zero_count_fraction,
        gini_coefficient=gini,
        issues=issues,
        warnings=warnings
    )


def detect_gc_bias(
    library: pd.DataFrame,
    counts: Dict[int, int],
    sequence_col: str = 'sequence'
) -> Tuple[float, str]:
    """
    Detect GC content bias in counting results.

    Sub-questions:
    1. Is there correlation between GC content and counts?
    2. Are high/low GC guides systematically under-counted?

    Args:
        library: gRNA library DataFrame
        counts: Counting results
        sequence_col: Column with sequences

    Returns:
        Tuple of (correlation coefficient, interpretation)
    """
    gc_contents = []
    count_values = []

    for idx, row in library.iterrows():
        seq = row[sequence_col].upper()
        gc = (seq.count('G') + seq.count('C')) / len(seq)
        gc_contents.append(gc)
        count_values.append(counts.get(idx, 0))

    # Calculate Spearman correlation
    if len(set(count_values)) < 2:
        return 0.0, "Cannot calculate correlation (no variation in counts)"

    # Use numpy for correlation
    gc_array = np.array(gc_contents)
    count_array = np.array(count_values)

    # Rank-based correlation (Spearman)
    gc_ranks = gc_array.argsort().argsort()
    count_ranks = count_array.argsort().argsort()

    correlation = np.corrcoef(gc_ranks, count_ranks)[0, 1]

    if abs(correlation) < 0.1:
        interpretation = "No significant GC bias detected"
    elif correlation > 0.3:
        interpretation = f"High GC guides may be over-represented (r={correlation:.3f})"
    elif correlation < -0.3:
        interpretation = f"Low GC guides may be over-represented (r={correlation:.3f})"
    else:
        interpretation = f"Weak GC bias detected (r={correlation:.3f})"

    return correlation, interpretation


# ============================================================================
# MISMATCH-TOLERANT MATCHING
# ============================================================================

def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between two sequences.

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Number of positions that differ
    """
    if len(seq1) != len(seq2):
        return max(len(seq1), len(seq2))
    return sum(c1 != c2 for c1, c2 in zip(seq1.upper(), seq2.upper()))


def find_best_match(
    query: str,
    library_sequences: List[str],
    max_mismatches: int = 1
) -> Tuple[Optional[int], int]:
    """
    Find the best matching guide for a query sequence with mismatch tolerance.

    Sub-questions:
    1. Is there an exact match?
    2. If not, is there a unique match within mismatch tolerance?
    3. How to handle ambiguous matches?

    Args:
        query: Query sequence
        library_sequences: List of library sequences
        max_mismatches: Maximum allowed mismatches

    Returns:
        Tuple of (best_match_index or None, number_of_mismatches)
    """
    best_idx = None
    best_distance = max_mismatches + 1
    match_count = 0

    query_upper = query.upper()

    for idx, lib_seq in enumerate(library_sequences):
        distance = hamming_distance(query_upper, lib_seq.upper())

        if distance == 0:
            return idx, 0  # Exact match

        if distance <= max_mismatches:
            match_count += 1
            if distance < best_distance:
                best_distance = distance
                best_idx = idx

    # Only return if unique match
    if match_count == 1 and best_idx is not None:
        return best_idx, best_distance
    elif match_count > 1:
        return None, best_distance  # Ambiguous
    else:
        return None, -1  # No match


def extract_grna_from_read(
    read: str,
    expected_length: int = 20,
    upstream_pattern: str = "ACCG",  # End of U6 promoter
    downstream_pattern: str = "GTTT"  # Start of scaffold
) -> Optional[str]:
    """
    Extract gRNA sequence from a read using flanking patterns.

    This improves exact matching by:
    1. Finding the expected position of gRNA
    2. Extracting only the gRNA region
    3. Avoiding partial matches with flanking sequences

    Args:
        read: Full sequencing read
        expected_length: Expected gRNA length
        upstream_pattern: Pattern immediately upstream of gRNA
        downstream_pattern: Pattern immediately downstream of gRNA

    Returns:
        Extracted gRNA sequence or None if not found
    """
    read_upper = read.upper()

    # Try to find upstream pattern
    upstream_pos = read_upper.find(upstream_pattern)
    if upstream_pos >= 0:
        start = upstream_pos + len(upstream_pattern)
        end = start + expected_length
        if end <= len(read):
            return read[start:end]

    # Try to find downstream pattern
    downstream_pos = read_upper.find(downstream_pattern)
    if downstream_pos >= expected_length:
        start = downstream_pos - expected_length
        return read[start:downstream_pos]

    return None


# ============================================================================
# BIAS CORRECTION
# ============================================================================

def normalize_by_gc(
    counts: Dict[int, int],
    library: pd.DataFrame,
    sequence_col: str = 'sequence'
) -> Dict[int, float]:
    """
    Normalize counts by GC content to correct for amplification bias.

    Sub-questions:
    1. What is the expected count per GC bin?
    2. How much does each guide deviate from expectation?
    3. What normalization factor to apply?

    Args:
        counts: Raw counts
        library: Library DataFrame
        sequence_col: Column with sequences

    Returns:
        GC-normalized counts
    """
    # Calculate GC for each guide
    gc_to_guides = {}
    for idx, row in library.iterrows():
        seq = row[sequence_col].upper()
        gc = round((seq.count('G') + seq.count('C')) / len(seq), 1)
        if gc not in gc_to_guides:
            gc_to_guides[gc] = []
        gc_to_guides[gc].append(idx)

    # Calculate median count per GC bin
    gc_medians = {}
    for gc, guides in gc_to_guides.items():
        bin_counts = [counts.get(g, 0) for g in guides]
        gc_medians[gc] = np.median(bin_counts) if bin_counts else 0

    # Global median
    global_median = np.median(list(counts.values())) if counts else 1

    # Normalize
    normalized = {}
    for idx, row in library.iterrows():
        seq = row[sequence_col].upper()
        gc = round((seq.count('G') + seq.count('C')) / len(seq), 1)

        gc_median = gc_medians.get(gc, global_median)
        scale_factor = global_median / gc_median if gc_median > 0 else 1

        raw_count = counts.get(idx, 0)
        normalized[idx] = raw_count * scale_factor

    return normalized


# ============================================================================
# COMPREHENSIVE QC REPORT
# ============================================================================

def generate_qc_report(
    library: pd.DataFrame,
    counts: Dict[int, int],
    total_reads: int,
    sequence_col: str = 'sequence'
) -> str:
    """
    Generate a comprehensive QC report for a counting run.

    Args:
        library: gRNA library
        counts: Counting results
        total_reads: Total reads processed
        sequence_col: Sequence column name

    Returns:
        Formatted QC report string
    """
    lines = []
    lines.append("=" * 70)
    lines.append("CRISPR gRNA COUNTING - QUALITY CONTROL REPORT")
    lines.append("=" * 70)

    # Library QC
    lines.append("\n1. LIBRARY QUALITY")
    lines.append("-" * 50)
    try:
        lib_qc = validate_library(library, sequence_col)
        lines.append(f"   Total guides: {lib_qc.total_guides}")
        lines.append(f"   Valid guides: {lib_qc.valid_guides}")
        lines.append(f"   Duplicates: {lib_qc.duplicate_count}")
        lines.append(f"   GC distribution: {lib_qc.gc_distribution}")
        if lib_qc.issues:
            lines.append("   Issues:")
            for issue in lib_qc.issues:
                lines.append(f"     - {issue}")
        if lib_qc.warnings:
            lines.append("   Warnings:")
            for warning in lib_qc.warnings:
                lines.append(f"     - {warning}")
    except GRNAException as e:
        lines.append(f"   ERROR: {e}")

    # Counting QC
    lines.append("\n2. COUNTING QUALITY")
    lines.append("-" * 50)
    counting_qc = evaluate_counting_results(counts, total_reads, len(library))
    lines.append(f"   Total reads: {counting_qc.total_reads:,}")
    lines.append(f"   Mapped reads: {counting_qc.mapped_reads:,} ({counting_qc.mapping_rate:.1%})")
    lines.append(f"   Guides with counts: {counting_qc.guides_with_counts}")
    lines.append(f"   Guides with zero counts: {counting_qc.guides_without_counts} ({counting_qc.zero_count_fraction:.1%})")
    lines.append(f"   Gini coefficient: {counting_qc.gini_coefficient:.3f}")
    if counting_qc.issues:
        lines.append("   Issues:")
        for issue in counting_qc.issues:
            lines.append(f"     - {issue}")
    if counting_qc.warnings:
        lines.append("   Warnings:")
        for warning in counting_qc.warnings:
            lines.append(f"     - {warning}")

    # GC Bias
    lines.append("\n3. BIAS ANALYSIS")
    lines.append("-" * 50)
    gc_corr, gc_interpretation = detect_gc_bias(library, counts, sequence_col)
    lines.append(f"   GC bias correlation: {gc_corr:.3f}")
    lines.append(f"   Interpretation: {gc_interpretation}")

    # Recommendations
    lines.append("\n4. RECOMMENDATIONS")
    lines.append("-" * 50)

    if counting_qc.mapping_rate < 0.5:
        lines.append("   - Consider checking read quality and trimming adapters")
        lines.append("   - Verify library and reads are from the same experiment")

    if counting_qc.zero_count_fraction > 0.3:
        lines.append("   - High dropout rate - check for selection effects or library issues")

    if abs(gc_corr) > 0.3:
        lines.append("   - Consider GC normalization for downstream analysis")

    if counting_qc.gini_coefficient > 0.8:
        lines.append("   - Very unequal distribution - may indicate strong selection")

    lines.append("\n" + "=" * 70)

    return "\n".join(lines)


# ============================================================================
# DEMONSTRATION
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("gRNA QC MODULE DEMONSTRATION")
    print("=" * 70)

    # Demo sequence validation
    print("\n1. SEQUENCE VALIDATION")
    print("-" * 40)

    test_sequences = [
        "ATGCATGCATGCATGCATGC",  # Valid
        "ATGCATGCATGCATGCNNNN",  # Contains N
        "AAAAAAATGCATGCATGCAT",  # Homopolymer
        "ATGCATGCATGCATGCATGCATGCATGCATGC",  # Too long
        "ATGX",  # Invalid character
    ]

    for seq in test_sequences:
        qc = validate_sequence(seq)
        status = "VALID" if qc.is_valid else "INVALID"
        print(f"\n{seq[:30]}...")
        print(f"  Status: {status}, GC: {qc.gc_content:.1%}, Length: {qc.length}")
        if qc.issues:
            print(f"  Issues: {qc.issues}")

    # Demo mismatch matching
    print("\n\n2. MISMATCH-TOLERANT MATCHING")
    print("-" * 40)

    library_seqs = [
        "ATGCATGCATGCATGCATGC",
        "GGCCGGCCGGCCGGCCGGCC",
        "AAAATTTTAAAATTTTAAAA",
    ]

    query = "ATGCATGCATGCATGCATGA"  # 1 mismatch from first

    idx, distance = find_best_match(query, library_seqs, max_mismatches=1)
    print(f"Query:   {query}")
    print(f"Library: {library_seqs}")
    print(f"Best match: index={idx}, mismatches={distance}")

    print("\n" + "=" * 70)
    print("Module ready for use!")
    print("=" * 70)
