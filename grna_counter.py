"""
CRISPR Guide RNA Counter using Aho-Corasick Algorithm

This module implements efficient gRNA counting in sequencing reads using
the Aho-Corasick algorithm, which is optimal for multiple pattern matching.

Why Aho-Corasick?
-----------------
The Aho-Corasick algorithm is ideal for CRISPR screening analysis because:

1. Multiple patterns: We need to search for thousands of gRNA sequences
   simultaneously (typical libraries have 50k-100k guides)

2. Efficiency: O(n + m + z) where:
   - n = total length of text (reads)
   - m = total length of patterns (gRNAs)
   - z = number of matches
   This is much faster than searching for each pattern separately O(n*k)

3. Trie-based: The algorithm builds a trie (prefix tree) from all gRNA
   sequences, then augments it with failure links for efficient traversal

How it works:
------------
1. BUILD PHASE: Construct an automaton from all gRNA sequences
   - Build a trie from all patterns
   - Add failure links (suffix links) for pattern matching
   - Add output links for patterns that are suffixes of other patterns

2. SEARCH PHASE: Process each read through the automaton
   - Traverse the automaton character by character
   - When reaching an accepting state, record the match
   - Use failure links to efficiently handle mismatches

The pyahocorasick library provides a highly optimized C implementation.
"""

import ahocorasick
from typing import Dict, List, Tuple, Optional
import pandas as pd
from collections import defaultdict
import time


class GRNACounter:
    """
    Efficient gRNA counter using Aho-Corasick algorithm.

    This class builds an Aho-Corasick automaton from a gRNA library
    and uses it to count occurrences in sequencing reads.
    """

    def __init__(self, library: pd.DataFrame, sequence_col: str = 'sequence'):
        """
        Initialize the counter with a gRNA library.

        Args:
            library: DataFrame containing gRNA sequences
            sequence_col: Name of the column containing sequences
        """
        self.library = library.copy()
        self.sequence_col = sequence_col
        self.automaton = None
        self.sequence_to_idx = {}

        # Build the automaton
        self._build_automaton()

    def _build_automaton(self) -> None:
        """
        Build the Aho-Corasick automaton from the gRNA library.

        The automaton is a finite state machine that can match
        all patterns simultaneously in a single pass through the text.
        """
        print("Building Aho-Corasick automaton...")
        start_time = time.time()

        self.automaton = ahocorasick.Automaton()

        # Add each gRNA sequence to the automaton
        for idx, row in self.library.iterrows():
            sequence = row[self.sequence_col]

            # Store sequence -> index mapping
            self.sequence_to_idx[sequence] = idx

            # Add to automaton: (pattern, value)
            # Value is the index in the library
            self.automaton.add_word(sequence, idx)

        # Finalize the automaton - this builds the failure links
        # Converting the trie to an Aho-Corasick automaton
        self.automaton.make_automaton()

        build_time = time.time() - start_time
        print(f"Automaton built in {build_time:.3f} seconds")
        print(f"  - Number of patterns: {len(self.sequence_to_idx)}")
        print(f"  - Automaton size: {self.automaton.get_stats()}")

    def count_reads(
        self,
        reads: List[str],
        count_multiple: bool = False,
        search_reverse_complement: bool = False
    ) -> Dict[int, int]:
        """
        Count gRNA occurrences in a list of reads.

        Args:
            reads: List of sequencing read strings
            count_multiple: If True, count multiple matches per read
                          If False, only count first match per read
            search_reverse_complement: Also search reverse complement

        Returns:
            Dictionary mapping library index to count
        """
        print(f"Counting gRNAs in {len(reads)} reads...")
        start_time = time.time()

        counts = defaultdict(int)
        matched_reads = 0
        total_matches = 0

        for read in reads:
            read_matched = False

            # Search forward strand
            for end_idx, library_idx in self.automaton.iter(read):
                if count_multiple or not read_matched:
                    counts[library_idx] += 1
                    total_matches += 1
                    read_matched = True

                if not count_multiple and read_matched:
                    break

            # Optionally search reverse complement
            if search_reverse_complement:
                rc_read = self._reverse_complement(read)
                for end_idx, library_idx in self.automaton.iter(rc_read):
                    if count_multiple or not read_matched:
                        counts[library_idx] += 1
                        total_matches += 1
                        read_matched = True

                    if not count_multiple and read_matched:
                        break

            if read_matched:
                matched_reads += 1

        count_time = time.time() - start_time
        print(f"Counting completed in {count_time:.3f} seconds")
        print(f"  - Reads processed: {len(reads)}")
        print(f"  - Reads with matches: {matched_reads} ({100*matched_reads/len(reads):.1f}%)")
        print(f"  - Total matches: {total_matches}")

        return dict(counts)

    def _reverse_complement(self, sequence: str) -> str:
        """
        Get the reverse complement of a DNA sequence.

        Args:
            sequence: DNA sequence string

        Returns:
            Reverse complement sequence
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                      'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence))

    def get_count_table(self, counts: Dict[int, int]) -> pd.DataFrame:
        """
        Generate a count table from counting results.

        Args:
            counts: Dictionary mapping library index to count

        Returns:
            DataFrame with sequence, gene_name, and count columns
        """
        result = self.library.copy()

        # Add counts column
        result['count'] = result.index.map(lambda x: counts.get(x, 0))

        # Sort by count (descending)
        result = result.sort_values('count', ascending=False)

        return result

    def get_stats(self) -> Dict:
        """Get statistics about the automaton."""
        return {
            'num_patterns': len(self.sequence_to_idx),
            'automaton_stats': self.automaton.get_stats() if self.automaton else None
        }


def count_grnas_in_fastq(
    fastq_path: str,
    library: pd.DataFrame,
    count_multiple: bool = False
) -> pd.DataFrame:
    """
    Count gRNAs in a FASTQ file.

    Args:
        fastq_path: Path to FASTQ file
        library: gRNA library DataFrame
        count_multiple: Count multiple matches per read

    Returns:
        Count table DataFrame
    """
    # Read FASTQ file
    reads = []
    with open(fastq_path, 'r') as f:
        line_num = 0
        for line in f:
            # FASTQ format: every 4 lines, line 2 (index 1) is the sequence
            if line_num % 4 == 1:
                reads.append(line.strip())
            line_num += 1

    print(f"Loaded {len(reads)} reads from {fastq_path}")

    # Create counter and count
    counter = GRNACounter(library)
    counts = counter.count_reads(reads, count_multiple=count_multiple)

    return counter.get_count_table(counts)


def count_grnas_in_fasta(
    fasta_path: str,
    library: pd.DataFrame,
    count_multiple: bool = False
) -> pd.DataFrame:
    """
    Count gRNAs in a FASTA file.

    Args:
        fasta_path: Path to FASTA file
        library: gRNA library DataFrame
        count_multiple: Count multiple matches per read

    Returns:
        Count table DataFrame
    """
    # Read FASTA file
    reads = []
    with open(fasta_path, 'r') as f:
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

    print(f"Loaded {len(reads)} reads from {fasta_path}")

    # Create counter and count
    counter = GRNACounter(library)
    counts = counter.count_reads(reads, count_multiple=count_multiple)

    return counter.get_count_table(counts)


class MismatchTolerantCounter:
    """
    Extended gRNA counter with mismatch tolerance for handling sequencing errors.

    Sub-questions addressed in design:
    1. How to handle reads with 1-2 sequencing errors? -> Hamming distance matching
    2. How to avoid ambiguous matches? -> Only count unique matches
    3. How to maintain efficiency? -> First try exact match, then mismatch
    4. How to extract gRNA from reads? -> Use flanking patterns

    This counter first attempts exact matching with Aho-Corasick,
    then falls back to mismatch-tolerant matching for unmatched reads.
    """

    def __init__(
        self,
        library: pd.DataFrame,
        sequence_col: str = 'sequence',
        max_mismatches: int = 1
    ):
        """
        Initialize counter with mismatch tolerance.

        Args:
            library: DataFrame containing gRNA sequences
            sequence_col: Name of the column containing sequences
            max_mismatches: Maximum mismatches to allow (0-2 recommended)
        """
        self.library = library.copy()
        self.sequence_col = sequence_col
        self.max_mismatches = max_mismatches
        self.sequences = library[sequence_col].tolist()
        self.sequence_length = len(self.sequences[0]) if self.sequences else 20

        # Build exact match automaton
        self.exact_counter = GRNACounter(library, sequence_col)

        # Pre-compute for mismatch matching
        self._build_mismatch_index()

    def _build_mismatch_index(self):
        """Build index structures for efficient mismatch matching."""
        print(f"Building mismatch index for {len(self.sequences)} sequences...")
        # Store sequences as uppercase for comparison
        self.sequences_upper = [s.upper() for s in self.sequences]

    def _hamming_distance(self, seq1: str, seq2: str) -> int:
        """Calculate Hamming distance between two sequences."""
        if len(seq1) != len(seq2):
            return max(len(seq1), len(seq2))
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

    def _find_mismatch_match(self, query: str) -> Tuple[Optional[int], int]:
        """
        Find best matching guide with mismatch tolerance.

        Args:
            query: Query sequence (must be same length as library sequences)

        Returns:
            Tuple of (best_match_idx or None, distance)
        """
        if len(query) != self.sequence_length:
            return None, -1

        query_upper = query.upper()
        best_idx = None
        best_dist = self.max_mismatches + 1
        match_count = 0

        for idx, lib_seq in enumerate(self.sequences_upper):
            dist = self._hamming_distance(query_upper, lib_seq)

            if dist == 0:
                return idx, 0

            if dist <= self.max_mismatches:
                match_count += 1
                if dist < best_dist:
                    best_dist = dist
                    best_idx = idx

        # Only return unique matches
        if match_count == 1:
            return best_idx, best_dist

        return None, best_dist if match_count > 0 else -1

    def _extract_grna(self, read: str) -> Optional[str]:
        """Extract gRNA from read using flanking sequences."""
        read_upper = read.upper()

        # Common flanking patterns
        upstream_patterns = ["ACCG", "CACCG"]
        downstream_patterns = ["GTTT", "GTTTT"]

        for up in upstream_patterns:
            pos = read_upper.find(up)
            if pos >= 0:
                start = pos + len(up)
                end = start + self.sequence_length
                if end <= len(read):
                    return read[start:end]

        for down in downstream_patterns:
            pos = read_upper.find(down)
            if pos >= self.sequence_length:
                start = pos - self.sequence_length
                return read[start:pos]

        return None

    def count_reads(
        self,
        reads: List[str],
        count_multiple: bool = False,
        search_reverse_complement: bool = False,
        use_mismatch_fallback: bool = True
    ) -> Dict[int, int]:
        """
        Count gRNA occurrences with optional mismatch tolerance.

        Args:
            reads: List of sequencing reads
            count_multiple: If True, count multiple matches per read
            search_reverse_complement: Also search reverse complement
            use_mismatch_fallback: Try mismatch matching for unmatched reads

        Returns:
            Dictionary mapping library index to count
        """
        print(f"Counting gRNAs in {len(reads)} reads (max mismatches: {self.max_mismatches})...")
        start_time = time.time()

        counts = defaultdict(int)
        exact_matches = 0
        mismatch_matches = 0
        unmatched = 0

        for read in reads:
            matched = False

            # First try exact match
            for end_idx, library_idx in self.exact_counter.automaton.iter(read):
                counts[library_idx] += 1
                exact_matches += 1
                matched = True
                if not count_multiple:
                    break

            # Fallback to mismatch matching
            if not matched and use_mismatch_fallback and self.max_mismatches > 0:
                grna = self._extract_grna(read)
                if grna:
                    idx, dist = self._find_mismatch_match(grna)
                    if idx is not None:
                        counts[idx] += 1
                        mismatch_matches += 1
                        matched = True

            if not matched:
                unmatched += 1

        count_time = time.time() - start_time
        print(f"Counting completed in {count_time:.3f} seconds")
        print(f"  - Exact matches: {exact_matches}")
        print(f"  - Mismatch matches: {mismatch_matches}")
        print(f"  - Unmatched reads: {unmatched}")

        return dict(counts)

    def get_count_table(self, counts: Dict[int, int]) -> pd.DataFrame:
        """Generate a count table from counting results."""
        return self.exact_counter.get_count_table(counts)


class AhoCorasickExplainer:
    """
    Educational class to explain how the Aho-Corasick algorithm works
    for gRNA counting.
    """

    @staticmethod
    def explain_algorithm():
        """Print explanation of the Aho-Corasick algorithm."""
        explanation = """
        =====================================================
        AHO-CORASICK ALGORITHM FOR gRNA COUNTING
        =====================================================

        The Aho-Corasick algorithm is a string matching algorithm that
        can find all occurrences of multiple patterns in a text in linear time.

        KEY CONCEPTS:
        -------------

        1. TRIE (Prefix Tree):
           - A tree where each node represents a character
           - Paths from root to leaves spell out the patterns
           - Shared prefixes share nodes (saves memory)

           Example for patterns: ['ACG', 'AC', 'CG', 'GC']

                    root
                   / | \\
                  A  C  G
                 /   |   \\
                C    G    C
               /
              G

        2. FAILURE LINKS:
           - Point to the longest proper suffix that is also a prefix
           - Allow the algorithm to continue searching without backtracking
           - Built using BFS after the trie is constructed

        3. OUTPUT LINKS:
           - Connect states to other patterns that are suffixes
           - Ensure we don't miss patterns that are substrings of others

        ALGORITHM PHASES:
        ----------------

        Phase 1: BUILD (Preprocessing)
        - Insert all gRNA sequences into the trie
        - Compute failure links using BFS
        - Time: O(m) where m = total length of all patterns

        Phase 2: SEARCH (Matching)
        - Process text character by character
        - At each position, follow success or failure links
        - Report matches at accepting states
        - Time: O(n + z) where n = text length, z = matches

        TOTAL TIME: O(m + n + z)

        WHY IT'S PERFECT FOR CRISPR SCREENING:
        -------------------------------------
        - Libraries have 50k-100k guides (patterns)
        - NGS produces millions of reads (text)
        - Need to count all occurrences efficiently
        - Linear time regardless of pattern count!

        COMPARISON:
        ----------
        Naive approach: O(n * m * k) for k patterns
        Aho-Corasick:   O(n + m + z)

        For 100M reads × 100k guides: orders of magnitude faster!
        """
        print(explanation)


if __name__ == "__main__":
    # Print algorithm explanation
    AhoCorasickExplainer.explain_algorithm()

    # Demo with small example
    from grna_library import generate_grna_library
    from sequence_generator import generate_simulated_reads

    print("\n" + "="*60)
    print("DEMONSTRATION")
    print("="*60 + "\n")

    # Generate small library
    library = generate_grna_library(num_guides=50, seed=42)
    print(f"Generated library with {len(library)} guides\n")

    # Generate reads
    reads, true_counts = generate_simulated_reads(
        library,
        num_reads=5000,
        error_rate=0,  # No errors for demo
        noise_fraction=0.1,
        seed=42
    )

    # Count using Aho-Corasick
    counter = GRNACounter(library)
    counted = counter.count_reads(reads)
    result_table = counter.get_count_table(counted)

    print("\nTop 10 counted gRNAs:")
    print(result_table[['guide_id', 'gene_name', 'sequence', 'count']].head(10).to_string())
