"""
Question Generation Module for CRISPR gRNA Counting Tool

This module implements a question-generation approach to provide comprehensive
answers about CRISPR screening, gRNA counting, and the Aho-Corasick algorithm.

Methodology:
-----------
When given a question, this module:
1. Generates relevant sub-questions to explore the topic more deeply
2. Provides answers to each sub-question
3. Combines the individual answers into a comprehensive final response

This approach ensures thorough coverage of complex topics and helps users
understand the full context of their queries.
"""

from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import re


@dataclass
class SubQuestion:
    """Represents a sub-question with its answer."""
    question: str
    answer: str
    relevance: str  # Why this sub-question is relevant


@dataclass
class QuestionResponse:
    """Complete response to a main question."""
    main_question: str
    sub_questions: List[SubQuestion]
    combined_answer: str


class QuestionGenerator:
    """
    Generates sub-questions and comprehensive answers for CRISPR-related topics.

    This class follows the methodology of:
    1. Taking a main question
    2. Generating relevant sub-questions to explore the topic
    3. Answering each sub-question
    4. Combining answers for a comprehensive response
    """

    def __init__(self):
        """Initialize the question generator with knowledge base."""
        self.knowledge_base = self._build_knowledge_base()
        self.question_patterns = self._build_question_patterns()

    def _build_knowledge_base(self) -> Dict[str, Dict]:
        """
        Build knowledge base about CRISPR screening and gRNA counting.

        Returns:
            Dictionary mapping topics to their information.
        """
        return {
            "crispr_basics": {
                "title": "CRISPR Screening Basics",
                "content": """CRISPR (Clustered Regularly Interspaced Short Palindromic Repeats)
                screening is a powerful genetic technique used to identify genes involved in
                specific biological processes. In a CRISPR screen, thousands of guide RNAs (gRNAs)
                are used to systematically knock out genes across the genome. After selection
                (e.g., drug treatment, growth conditions), the abundance of each gRNA is measured
                by sequencing. Guides targeting essential genes will be depleted, while guides
                targeting genes that confer resistance will be enriched."""
            },
            "grna_library": {
                "title": "gRNA Libraries",
                "content": """A gRNA library is a collection of guide RNA sequences designed to
                target specific genes. Common libraries include Brunello (human), GeCKO, and TKOv3.
                These libraries typically contain 50,000-100,000 unique guide sequences, with
                multiple guides targeting each gene (usually 4-10 guides per gene). Each guide
                is a 20-nucleotide sequence that directs Cas9 to cut at a specific genomic location."""
            },
            "aho_corasick": {
                "title": "Aho-Corasick Algorithm",
                "content": """The Aho-Corasick algorithm is a string matching algorithm that
                efficiently finds all occurrences of multiple patterns in a text. It works by:
                1. Building a trie (prefix tree) from all patterns
                2. Adding failure links for efficient backtracking
                3. Processing text in a single pass, following the automaton

                Time complexity: O(n + m + z) where n=text length, m=pattern length, z=matches.
                This is ideal for gRNA counting where we search millions of reads for thousands
                of patterns simultaneously."""
            },
            "sequencing": {
                "title": "NGS Sequencing",
                "content": """Next-Generation Sequencing (NGS) produces millions of short reads
                (typically 50-300bp) from a DNA sample. In CRISPR screening, reads contain the
                guide RNA sequences flanked by vector sequences. The sequencing process introduces
                small errors (0.1-0.5% per base), which the counting algorithm must handle.
                Common formats include FASTQ (with quality scores) and FASTA (sequences only)."""
            },
            "counting_results": {
                "title": "Interpreting Count Results",
                "content": """gRNA counts represent how many times each guide was detected in the
                sequencing data. Key metrics include:
                - Total reads: Number of sequences processed
                - Mapped reads: Reads containing a recognized guide
                - Mapping rate: Percentage of reads with matches
                - Count distribution: How counts are spread across guides

                Highly enriched or depleted guides suggest their target genes are important for
                the phenotype being studied."""
            },
            "algorithm_efficiency": {
                "title": "Algorithm Efficiency",
                "content": """The Aho-Corasick algorithm provides significant efficiency gains:
                - Naive approach: O(n * k) where k = number of patterns (100k guides)
                - Aho-Corasick: O(n + m + z) - linear in input size

                For typical experiments with 100M reads and 100k guides, Aho-Corasick can be
                orders of magnitude faster than searching for each pattern individually."""
            }
        }

    def _build_question_patterns(self) -> Dict[str, List[str]]:
        """
        Build patterns to identify question topics and generate sub-questions.

        Returns:
            Dictionary mapping topic patterns to related sub-questions.
        """
        return {
            "what.*crispr|crispr.*what": [
                "What is the biological mechanism of CRISPR?",
                "How are gRNA libraries designed?",
                "What makes CRISPR screening useful for research?"
            ],
            "how.*count|count.*how|algorithm": [
                "What algorithm is used for counting?",
                "Why is Aho-Corasick efficient for this task?",
                "How does the trie structure work?",
                "What are failure links and why are they important?"
            ],
            "library|guide|grna": [
                "What is a gRNA library?",
                "How many guides are typically in a library?",
                "How are guides designed for each gene?"
            ],
            "sequencing|read|fastq|fasta": [
                "What sequencing format is used?",
                "How do sequencing errors affect counting?",
                "What is a typical read length?"
            ],
            "result|count|output|interpret": [
                "What do the counts represent?",
                "How do I interpret enrichment vs depletion?",
                "What is a good mapping rate?"
            ],
            "efficien|fast|performance|speed": [
                "What is the time complexity of the algorithm?",
                "How does it compare to naive approaches?",
                "Why is Aho-Corasick optimal for this use case?"
            ]
        }

    def identify_topic(self, question: str) -> List[str]:
        """
        Identify the topic(s) of a question.

        Args:
            question: The input question string.

        Returns:
            List of identified topics.
        """
        question_lower = question.lower()
        topics = []

        topic_keywords = {
            "crispr_basics": ["crispr", "screen", "knockout", "gene"],
            "grna_library": ["library", "guide", "grna", "brunello", "gecko"],
            "aho_corasick": ["algorithm", "aho", "corasick", "trie", "automaton"],
            "sequencing": ["sequenc", "read", "fastq", "fasta", "ngs"],
            "counting_results": ["result", "count", "output", "interpret", "enrich", "deplet"],
            "algorithm_efficiency": ["efficien", "fast", "performance", "complex", "speed"]
        }

        for topic, keywords in topic_keywords.items():
            if any(kw in question_lower for kw in keywords):
                topics.append(topic)

        # Default to basics if no specific topic identified
        if not topics:
            topics = ["crispr_basics", "aho_corasick"]

        return topics

    def generate_sub_questions(self, main_question: str) -> List[Tuple[str, str]]:
        """
        Generate relevant sub-questions for a main question.

        Args:
            main_question: The main question to analyze.

        Returns:
            List of (sub_question, relevance) tuples.
        """
        sub_questions = []
        question_lower = main_question.lower()

        # Find matching patterns and add their sub-questions
        for pattern, questions in self.question_patterns.items():
            if re.search(pattern, question_lower):
                for q in questions:
                    relevance = f"Helps explore aspects of '{pattern.split('|')[0]}'"
                    sub_questions.append((q, relevance))

        # Add topic-specific questions
        topics = self.identify_topic(main_question)
        for topic in topics:
            topic_qs = self._get_topic_sub_questions(topic)
            for q, r in topic_qs:
                if (q, r) not in sub_questions:
                    sub_questions.append((q, r))

        # Ensure we always have some sub-questions
        if not sub_questions:
            sub_questions = [
                ("What is the purpose of gRNA counting?", "Establishes context"),
                ("How does the algorithm work?", "Explains the mechanism"),
                ("What are the key outputs?", "Describes results")
            ]

        return sub_questions[:6]  # Limit to 6 sub-questions

    def _get_topic_sub_questions(self, topic: str) -> List[Tuple[str, str]]:
        """Get sub-questions for a specific topic."""
        topic_questions = {
            "crispr_basics": [
                ("What is a CRISPR screen?", "Provides biological context"),
                ("Why count gRNAs?", "Explains the purpose")
            ],
            "grna_library": [
                ("What is a gRNA library?", "Defines key input data"),
                ("How are guides structured?", "Explains sequence format")
            ],
            "aho_corasick": [
                ("How does pattern matching work?", "Explains core algorithm"),
                ("What is a trie?", "Describes data structure")
            ],
            "sequencing": [
                ("What is NGS?", "Explains sequencing technology"),
                ("How are reads formatted?", "Describes input format")
            ],
            "counting_results": [
                ("What do counts mean?", "Interprets output"),
                ("How to identify hits?", "Guides analysis")
            ],
            "algorithm_efficiency": [
                ("What is the time complexity?", "Quantifies performance"),
                ("How does it scale?", "Addresses scalability")
            ]
        }
        return topic_questions.get(topic, [])

    def answer_sub_question(self, sub_question: str) -> str:
        """
        Generate an answer for a sub-question.

        Args:
            sub_question: The sub-question to answer.

        Returns:
            Answer string.
        """
        question_lower = sub_question.lower()

        # Match to knowledge base topics and generate answers
        for topic, info in self.knowledge_base.items():
            topic_keywords = {
                "crispr_basics": ["crispr", "screen", "biological", "knockout"],
                "grna_library": ["library", "guide", "grna", "design"],
                "aho_corasick": ["algorithm", "pattern", "trie", "matching", "aho", "automaton", "failure"],
                "sequencing": ["sequenc", "read", "ngs", "format", "fastq"],
                "counting_results": ["count", "result", "interpret", "output", "hit"],
                "algorithm_efficiency": ["complex", "efficien", "time", "scale", "performance"]
            }

            if topic in topic_keywords:
                if any(kw in question_lower for kw in topic_keywords[topic]):
                    return info["content"].strip().replace("                ", " ")

        # Default answer
        return "This topic relates to the gRNA counting process using the Aho-Corasick algorithm."

    def generate_response(self, main_question: str) -> QuestionResponse:
        """
        Generate a complete response with sub-questions and combined answer.

        This method implements the full question-answering methodology:
        1. Generate relevant sub-questions
        2. Answer each sub-question
        3. Combine into comprehensive response

        Args:
            main_question: The main question to answer.

        Returns:
            QuestionResponse with sub-questions and combined answer.
        """
        # Generate sub-questions
        sub_q_pairs = self.generate_sub_questions(main_question)

        # Answer each sub-question
        sub_questions = []
        for question, relevance in sub_q_pairs:
            answer = self.answer_sub_question(question)
            sub_questions.append(SubQuestion(
                question=question,
                answer=answer,
                relevance=relevance
            ))

        # Combine answers
        combined = self._combine_answers(main_question, sub_questions)

        return QuestionResponse(
            main_question=main_question,
            sub_questions=sub_questions,
            combined_answer=combined
        )

    def _combine_answers(
        self,
        main_question: str,
        sub_questions: List[SubQuestion]
    ) -> str:
        """
        Combine sub-question answers into a comprehensive response.

        Args:
            main_question: The original question.
            sub_questions: List of answered sub-questions.

        Returns:
            Combined comprehensive answer.
        """
        # Build combined answer
        parts = []

        # Introduction
        topics = self.identify_topic(main_question)
        if topics:
            intro_topic = topics[0]
            if intro_topic in self.knowledge_base:
                parts.append(f"To answer your question about {self.knowledge_base[intro_topic]['title']}:\n")

        # Add unique insights from sub-questions
        seen_content = set()
        for sq in sub_questions:
            # Avoid duplicate content
            content_key = sq.answer[:50]
            if content_key not in seen_content:
                parts.append(f"- {sq.answer[:200]}...")
                seen_content.add(content_key)

        # Summary
        parts.append("\nIn summary: The gRNA counting tool uses the Aho-Corasick algorithm " +
                    "to efficiently count guide RNA sequences in CRISPR screening experiments, " +
                    "processing millions of reads against thousands of patterns in linear time.")

        return "\n".join(parts)

    def format_response(self, response: QuestionResponse) -> str:
        """
        Format a QuestionResponse for display.

        Args:
            response: The QuestionResponse to format.

        Returns:
            Formatted string for display.
        """
        lines = []
        lines.append("=" * 70)
        lines.append("QUESTION ANALYSIS")
        lines.append("=" * 70)
        lines.append(f"\nMain Question: {response.main_question}\n")

        lines.append("-" * 70)
        lines.append("SUB-QUESTIONS GENERATED")
        lines.append("-" * 70)

        for i, sq in enumerate(response.sub_questions, 1):
            lines.append(f"\n{i}. {sq.question}")
            lines.append(f"   Relevance: {sq.relevance}")
            lines.append(f"   Answer: {sq.answer[:150]}...")

        lines.append("\n" + "-" * 70)
        lines.append("COMBINED ANSWER")
        lines.append("-" * 70)
        lines.append(f"\n{response.combined_answer}")
        lines.append("\n" + "=" * 70)

        return "\n".join(lines)


def ask_question(question: str, verbose: bool = True) -> QuestionResponse:
    """
    Convenience function to ask a question and get a comprehensive answer.

    Args:
        question: The question to answer.
        verbose: If True, print formatted response.

    Returns:
        QuestionResponse object.
    """
    generator = QuestionGenerator()
    response = generator.generate_response(question)

    if verbose:
        print(generator.format_response(response))

    return response


# Example questions for demonstration
EXAMPLE_QUESTIONS = [
    "How does the gRNA counting algorithm work?",
    "What is CRISPR screening and why is it useful?",
    "Why is Aho-Corasick efficient for pattern matching?",
    "How do I interpret the counting results?",
    "What format should my sequencing data be in?",
]


if __name__ == "__main__":
    print("=" * 70)
    print("QUESTION GENERATOR DEMONSTRATION")
    print("=" * 70)
    print("\nThis module generates sub-questions to provide comprehensive answers")
    print("about CRISPR screening and gRNA counting.\n")

    # Demonstrate with an example question
    example = "How does the gRNA counting algorithm work?"
    print(f"Example question: '{example}'\n")

    response = ask_question(example)

    print("\n" + "=" * 70)
    print("OTHER EXAMPLE QUESTIONS YOU CAN ASK:")
    print("=" * 70)
    for q in EXAMPLE_QUESTIONS[1:]:
        print(f"  - {q}")
