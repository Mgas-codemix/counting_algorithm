# CRISPR Guide RNA Counter

A tool to count guide RNAs in your sequencing data from CRISPR screening experiments.

---

## What Does This Tool Do?

When you run a CRISPR screening experiment, you need to know how many times each guide RNA (gRNA) appears in your sequencing data. This tells you which genes were "hit" more or less often, helping you identify important genes.

This tool:
1. Takes your sequencing reads (the DNA sequences from your sequencer)
2. Searches for all your guide RNA sequences in those reads
3. Counts how many times each guide appears
4. Gives you a simple table with the results

**Why is it fast?** It uses a clever algorithm called "Aho-Corasick" that can search for thousands of guide sequences at the same time, instead of one by one. Think of it like searching for many words in a book simultaneously, rather than reading the book once for each word.

---

## Before You Start

### What You Need

1. **Python 3.8 or newer** installed on your computer
2. **Your data files:**
   - A guide RNA library file (CSV format) with your guide sequences
   - Your sequencing reads file (FASTQ or FASTA format)

### Installing the Tool

Open your terminal (command line) and follow these steps:

```bash
# Step 1: Go to the folder where you downloaded this tool
cd counting_algorithm

# Step 2: Install the required packages
pip install -r requirements.txt
```

That's it! The tool is ready to use.

---

## How to Use It

### Option 1: Try the Demo First (Recommended for Beginners)

If you want to see how the tool works before using your own data, run the demo:

```bash
python main.py
```

This will:
- Create a fake guide RNA library (for testing)
- Generate fake sequencing reads (simulating what comes from a sequencer)
- Count all the guides
- Save the results in a folder called `output/`

You'll see something like this on your screen:

```
STEP 1: Generating gRNA Reference Library
STEP 2: Generating Simulated NovaSeq Reads
STEP 3: Counting gRNAs with Aho-Corasick Algorithm
STEP 4: Results Summary

TOP 20 COUNTED gRNAs (OUTPUT TABLE)
   guide_id     gene_name             sequence  count
   sgTP53_4          TP53 CGGGCCAATTACCTGTCTTA    524
   sgBRCA1_3        BRCA1 CGCGATACCCTACCATACCA    589
   ...
```

### Option 2: Use Your Own Data

If you have your own guide library and sequencing data:

```bash
python main.py --library your_library.csv --reads your_reads.fastq --output results.csv
```

**Replace:**
- `your_library.csv` with the path to your guide RNA library file
- `your_reads.fastq` with the path to your sequencing file
- `results.csv` with whatever name you want for the output file

---

## Preparing Your Files

### Your Guide RNA Library File

This should be a CSV (comma-separated) file with at least these columns:

| guide_id | gene_name | sequence |
|----------|-----------|----------|
| sg_001 | TP53 | ACGTACGTACGTACGTACGT |
| sg_002 | TP53 | TGCATGCATGCATGCATGCA |
| sg_003 | BRCA1 | GGCCAATTGGCCAATTGGCC |

- **guide_id**: A unique name for each guide (can be anything)
- **gene_name**: The gene this guide targets
- **sequence**: The 20-letter DNA sequence of the guide (using A, T, G, C)

Save this as a `.csv` file (you can create it in Excel and "Save As" CSV).

### Your Sequencing Reads File

This is the file you get from your sequencing facility. It's usually in FASTQ format (file ending in `.fastq` or `.fq`).

The tool also accepts FASTA format. If your file is FASTA, add `--format fasta`:

```bash
python main.py --library your_library.csv --reads your_reads.fasta --format fasta --output results.csv
```

---

## Understanding the Results

After running, you'll get a CSV file with your results. Open it in Excel or any spreadsheet program:

| guide_id | gene_name | sequence | count |
|----------|-----------|----------|-------|
| sg_001 | TP53 | ACGTACGT... | 1523 |
| sg_002 | TP53 | TGCATGCA... | 892 |
| sg_003 | BRCA1 | GGCCAATT... | 2341 |

**What the columns mean:**
- **guide_id**: The name of each guide
- **gene_name**: Which gene this guide targets
- **sequence**: The DNA sequence of the guide
- **count**: **How many times this guide was found in your reads** (this is the important number!)

Guides with higher counts appeared more in your sample. Depending on your experiment:
- In a **positive selection screen**: High counts = genes that help cells survive
- In a **negative selection screen**: Low counts = genes essential for survival

---

## Common Tasks

### Change the Number of Guides or Reads in Demo Mode

```bash
# Use 1000 guides and 200,000 reads
python main.py --num-guides 1000 --num-reads 200000
```

### Save Output to a Different Folder

```bash
python main.py --output-dir my_results
```

### Learn About the Algorithm

If you're curious about how the counting works:

```bash
python main.py --explain
```

---

## Troubleshooting

### "Command not found: python"

Try using `python3` instead:
```bash
python3 main.py
```

### "No module named 'ahocorasick'"

You need to install the required packages:
```bash
pip install -r requirements.txt
```

Or if that doesn't work:
```bash
pip3 install -r requirements.txt
```

### "File not found" Error

Make sure you're in the right folder and your file paths are correct. Use the full path if needed:
```bash
python main.py --library /full/path/to/your_library.csv --reads /full/path/to/your_reads.fastq
```

### The Counts Seem Too Low

This could happen if:
1. Your guide sequences don't match what's in the reads (check for typos)
2. The guides are in reverse complement orientation - the tool searches forward by default
3. There are too many sequencing errors - **try using mismatch tolerance** (see below)

---

## Advanced Features

### Mismatch-Tolerant Counting

Real sequencing data has errors (~0.1-1% per base). By default, the tool only counts **exact matches**. If you're missing counts due to sequencing errors, enable mismatch tolerance:

```bash
# Allow 1 mismatch (recommended for most cases)
python main.py --max-mismatches 1

# Allow up to 2 mismatches (use with caution - may cause false matches)
python main.py --max-mismatches 2
```

**How it works:**
1. First tries exact matching (fast, using Aho-Corasick)
2. For unmatched reads, extracts the gRNA region using flanking sequences
3. Finds the best match within the allowed mismatch threshold
4. Only counts **unique matches** to avoid ambiguity

**Example output with mismatch tolerance:**
```
Counting gRNAs in 100000 reads (max mismatches: 1)...
  - Exact matches: 85234
  - Mismatch matches: 8921
  - Unmatched reads: 5845
```

### Quality Control Reports

The tool automatically runs quality control checks on your library and counting results:

```bash
# Run with QC (default)
python main.py

# Skip QC for faster runs
python main.py --no-qc
```

**QC checks include:**
- **Library validation**: Invalid characters, duplicates, GC content distribution
- **Counting metrics**: Mapping rate, zero-count guides, count distribution
- **Bias detection**: GC content bias analysis

A QC report is saved to `output/qc_report.txt` with details like:
```
CRISPR gRNA COUNTING - QUALITY CONTROL REPORT
----------------------------------------------
1. LIBRARY QUALITY
   Total guides: 500
   Valid guides: 500
   GC distribution: {'low (<30%)': 12, 'normal (30-70%)': 476, 'high (>70%)': 12}

2. COUNTING QUALITY
   Total reads: 100,000
   Mapped reads: 94,155 (94.2%)
   Guides with zero counts: 23 (4.6%)
   Gini coefficient: 0.456

3. BIAS ANALYSIS
   GC bias correlation: 0.023
   Interpretation: No significant GC bias detected
```

### Understanding QC Metrics

| Metric | Good Value | Warning Signs |
|--------|------------|---------------|
| Mapping rate | >80% | <50% suggests library/read mismatch |
| Zero-count guides | <20% | >50% indicates high dropout |
| Gini coefficient | 0.3-0.7 | >0.9 = very unequal distribution |
| GC bias correlation | <0.1 | >0.3 = significant GC bias |

### GC Content Normalization

If GC bias is detected, you can apply normalization:

```bash
python main.py --gc-normalize
```

This adjusts counts based on GC content to correct for PCR amplification bias.

### R Visualization and QC for DEG Analysis

An R script is provided for detailed visualization and QC specifically for downstream differential expression (DEG) analysis:

```bash
# Run from command line
Rscript qc_visualization.R output/grna_counts.csv output/qc/

# Or in R
source("qc_visualization.R")
results <- run_qc_analysis("output/grna_counts.csv")
```

**Generated outputs:**
- `qc_report_main.png` - Combined QC visualization
- `qc_warning_guides.csv` - Guides flagged for potential issues
- `counts_with_qc_flags.csv` - Full data with QC annotations
- Individual plots: distribution, outliers, GC bias, etc.

**QC flags for DEG analysis:**
| Flag | Description | Impact on DEG |
|------|-------------|---------------|
| `flag_zero_count` | No reads detected | Will cause errors in DEG tools |
| `flag_low_count` | Bottom 10% or <10 reads | Unreliable fold-change estimates |
| `flag_high_outlier` | IQR outlier (high) | May dominate analysis |
| `flag_extreme_gc` | GC <25% or >75% | Potential amplification bias |

**Visualizations include:**
- Count distribution (histogram + log-transformed)
- Cumulative distribution (Lorenz curve for inequality)
- Outlier detection plot
- Gene-level boxplots
- GC content vs count scatter plot
- QC status summary bar chart

---

## Handling Common CRISPR Experiment Issues

### Issue: Low Mapping Rate

**Symptoms:** <50% of reads match guides

**Possible causes & solutions:**
1. **Wrong library file** - Verify your library matches the experiment
2. **Adapter contamination** - Trim adapters before counting
3. **High error rate** - Use `--max-mismatches 1`

### Issue: Many Guides with Zero Counts

**Symptoms:** >30% of guides have count=0

**Possible causes & solutions:**
1. **Strong selection** - Expected in some screens (not an error)
2. **Library representation issues** - Check original library QC
3. **Sequencing depth too low** - More reads needed

### Issue: Extreme Count Distribution

**Symptoms:** A few guides dominate all counts (Gini >0.9)

**Possible causes & solutions:**
1. **Strong selection pressure** - Expected in drug screens
2. **PCR bias** - Consider GC normalization
3. **Contamination** - Check for specific guide dominance

---

## Quick Reference

| Command | What It Does |
|---------|--------------|
| `python main.py` | Run demo with default settings |
| `python main.py --help` | Show all available options |
| `python main.py --explain` | Explain how the algorithm works |
| `python main.py -l LIB -r READS -o OUT` | Count guides from your files |
| `python main.py --max-mismatches 1` | Enable mismatch-tolerant counting |
| `python main.py --no-qc` | Skip quality control checks |
| `python main.py --gc-normalize` | Apply GC content normalization |

---

## Getting Help

If you run into problems:
1. Check the Troubleshooting section above
2. Make sure your files are in the correct format
3. Try running the demo first to make sure the tool works on your computer

---

## Summary

1. **Install**: `pip install -r requirements.txt`
2. **Try demo**: `python main.py`
3. **Use your data**: `python main.py -l library.csv -r reads.fastq -o counts.csv`
4. **Handle errors**: Add `--max-mismatches 1` if counts are low due to sequencing errors
5. **Check quality**: Review `output/qc_report.txt` for potential issues
6. **Open results**: Look at the CSV file in Excel

That's all you need to count guide RNAs in your CRISPR screening data!
