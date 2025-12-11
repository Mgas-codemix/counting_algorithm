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
3. There are too many sequencing errors

---

## Quick Reference

| Command | What It Does |
|---------|--------------|
| `python main.py` | Run demo with default settings |
| `python main.py --help` | Show all available options |
| `python main.py --explain` | Explain how the algorithm works |
| `python main.py -l LIB -r READS -o OUT` | Count guides from your files |

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
4. **Open results**: Look at the CSV file in Excel

That's all you need to count guide RNAs in your CRISPR screening data!
