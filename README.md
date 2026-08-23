# mtDNA Mutation Analyser

A command-line tool for analysing point mutations in human mitochondrial DNA (mtDNA). Built in Python using only the standard library.

## Overview

This tool loads a mitochondrial DNA sequence from a FASTA file, detects Open Reading Frames (ORFs), and allows the user to analyse point mutations — both manually and in bulk. It classifies mutations as **Silent**, **Missense**, or **Nonsense**, and performs statistical analysis to determine whether a given ORF's mutation profile differs significantly from the full mtDNA baseline.

> Uses NCBI Translation Table 1 - Standard Genetic code.

## Features

- **FASTA loading & validation** — reads and validates mtDNA sequence
- **ORF detection** — finds all ORFs ≥ 90 bp on the forward strand
- **Manual mutation analysis** — enter a mutation (e.g. `A25G`), get codon-level breakdown and mutation type
- **Batch analysis** — automatically tests all possible single-base substitutions across a selected ORF
- **Chi-square test** — compares ORF mutation profile against full mtDNA distribution (df=2, α=0.05)
- **Hotspot detection** — identifies the position most susceptible to missense mutations
- **CSV export** — saves mutation history to file

## Project Structure

```
├── main.py           # Entry point, menu system
├── analysis.py       # Batch analysis, chi-square, p-value, hotspot detection
├── fasta_utils.py    # FASTA loading and DNA validation
├── orf_utils.py      # ORF detection algorithm
├── mutation_utils.py # Mutation parsing, codon operations, classification
├── constants.py      # Standard genetic codon table
├── models.py         # MutationRecord class
├── storage.py        # CSV export
├── decorators.py     # Timing decorator
├── exceptions.py     # Custom exceptions
└── mtDNA.fasta       # Human mitochondrial DNA sequence
```

## How to Run

```bash
python main.py
```

No external dependencies — uses Python standard library only.

## Statistical Method

Batch analysis tests every possible single-base substitution at every position of a selected ORF. Results are compared against the full mtDNA mutation profile using a chi-square goodness-of-fit test:

```
χ² = Σ (observed - expected)² / expected
```

Expected values are normalised to match the ORF's total mutation count. P-value is computed analytically for df=2:

```
p = e^(-χ²/2)
```

## Example Output

```
ORF_19 — Batch Analysis
Silent: 205
Nonsense: 25
Missense: 670

Chi-square result: 6.78
P-value: 0.034
Statistically significant difference from full DNA distribution!

Hotspot position: 47, Missense count: 3
```

## Author

Recep Berke Küçükyağcı  
BSc Computer Science Student
Polish-Japanese Academy of Information Technology (PJATK), Warsaw
