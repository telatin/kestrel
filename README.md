# Kestrel: a simple K-mer Taxonomy Classifier

![Kestrel Banner](docs/kestrel.svg)

A fast k-mer-based taxonomic classifier written in Nim, inspired by Sepia and Kraken2. 
Klass reads taxonomy information directly from FASTA sequence headers and uses efficient k-mer/minimizer algorithms for accurate read classification.

## Features

- **Header-based taxonomy**: Reads GTDB-style taxonomy directly from FASTA headers
- **Efficient k-mer processing**: Canonical k-mers with optional minimizer compression
- **Spaced k-mers**: Supports shaped k-mers (e.g., OOOO--OOOO) for improved specificity
- **LCA resolution**: Handles k-mer conflicts using Lowest Common Ancestor algorithms
- **Quality filtering**: Phred score-based quality filtering for FASTQ reads
- **Compressed formats**: Supports gzipped FASTA/FASTQ files via readfx
- **Fast binary database**: Compact k-mer database with MurmurHash3 indexing

## Installation

```bash
# Install dependencies and build binaries
nimble build
```

This creates two executables in the `bin/` directory:
- `kestrel-build` - Database building tool
- `kestrel-classify` - Read classification tool

## Usage

### 1. Database Building

Build a k-mer database from reference genomes with taxonomy in headers:

```bash
./bin/kestrel-build [options] -o <database> <FASTA_FILES>...
```

**Options:**
- `-o, --output DB` - Output database directory (required)
- `-k, --kmer-size INT` - K-mer size [default: 25]
- `-s, --kmer-shape STRING` - K-mer shape pattern (e.g., OOOO--OOOO)
- `-m, --minimizer-size INT` - Minimizer size (0 = no minimizers) [default: 0]
- `--verbose` - Print verbose logging
- `--debug` - Print debug information

**Example:**
```bash
# Build database with k-mer size 31
./bin/kestrel-build --verbose -o mydb -k 31 genome1.fasta genome2.fasta

# Build with spaced k-mers (k=8, window=10)
./bin/kestrel-build --verbose -o mydb -s "OOOO--OOOO" genome1.fasta genome2.fasta

# Build with minimizers for smaller database
./bin/kestrel-build -o mydb -k 31 -m 15 *.fasta
```

:warning: **FASTA Header Format:**
Sequences must have GTDB-style taxonomy in the header comment:
```
>sequence_name d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
ATCGATCGATCG...
```

### 2. Read Classification

Classify FASTQ/FASTA reads against a built database:

```bash
./bin/kestrel-classify [options] -d <database> -o <output> <FASTQ_FILES>...
```

**Options:**
- `-d, --db DIR` - K-mer database directory (required)
- `-o, --output PREFIX` - Output prefix for results (required)
- `-q, --quality INT` - Minimum quality score (Phred) [default: 15]
- `-m, --min-hits INT` - Minimum k-mer hits for classification [default: 3]
- `--verbose` - Print verbose logging

**Example:**
```bash
# Classify single-end reads
./bin/kestrel-classify --verbose -d mydb -o sample1 reads.fastq.gz

# Classify with quality filtering
./bin/kestrel-classify -d mydb -o sample1 -q 20 reads.fastq.gz

# Classify multiple files
./bin/kestrel-classify -d mydb -o batch sample1.fq.gz sample2.fq.gz
```

## Output Files

### Classification Results

**`{prefix}_classification.txt`** - Per-read classifications:
```
READ_ID    TAXONOMY    HIT_COUNT    TOTAL_KMERS
READ:1     d__Bacteria;p__Proteobacteria;...;s__Escherichia coli    25    30
READ:2     d__Bacteria;p__Proteobacteria;...;g__Escherichia         18    28  
READ:3     no hits                                                   0     25
```

**`{prefix}_summary.txt`** - Summary statistics:
```
Taxonomy                                                     Reads  Avg_Score  Total_Bases
d__Bacteria;p__Proteobacteria;...;s__Escherichia coli       15420    0.823     4635600
d__Bacteria;p__Proteobacteria;...;g__Escherichia             2100    0.645      630000
no hits                                                       1250    0.000      375000
```

### Database Files

**`params.json`** - Database parameters and taxonomy mappings
**`kmers.bin`** - Binary k-mer to taxonomy mapping

## Algorithm Details

### K-mer Processing
- Uses 2-bit nucleotide encoding (A=0, C=1, G=2, T=3)
- Computes canonical k-mers (lexicographically smaller of forward/reverse complement)
- Supports spaced k-mers with custom patterns (O = include base, - = skip base)
- Optional minimizer compression using sliding window with deque optimization
- Skips k-mers containing ambiguous bases (N, Y, R, etc.)

### Spaced K-mers
Kestrel supports spaced k-mers (also called shaped k-mers) which can improve classification specificity by focusing on conserved positions while skipping variable regions:

- **Pattern format**: String of 'O' (include) and '-' (skip) characters
- **Examples**:
  - `OOOOOOOO` = traditional 8-mer (equivalent to `--kmer-size 8`)
  - `OOOO--OOOO` = 8-mer with 2-base gap (window size = 10)
  - `OOO-O-OOO` = 7-mer with single-base gaps (window size = 9)
- **Benefits**: Can reduce false positives by ignoring hypervariable positions
- **Usage**: Use `--kmer-shape` instead of `--kmer-size` (mutually exclusive)

### Taxonomy Resolution
- Builds hierarchical taxonomy graph from GTDB lineages
- Maps full taxonomy strings to most specific taxonomic level
- Uses LCA (Lowest Common Ancestor) to resolve k-mer conflicts
- Handles ties in classification using hierarchical scoring

### Quality Control
- Phred+33 quality score filtering for FASTQ reads
- Replaces low-quality bases with 'N' before k-mer extraction
- Configurable minimum hit thresholds for classification confidence

## Test Data

The `data/` directory contains example files:

**Reference genomes** (`data/ref/`):
- `ecoli.fasta` - E. coli reference with taxonomy header
- `senterica.fasta` - Salmonella enterica reference
- `efergusonii.fasta` - E. fergusonii reference

**Test reads** (`data/reads/`):
- `Ecoli1_1000K_*.fastq.gz` - E. coli test reads (paired-end)
- `Efaecium_1000K_*.fastq.gz` - E. faecium test reads

**Example workflow:**
```bash
# Build test database
./bin/kestrel-build --verbose -o data/testdb -k 25 data/ref/{eco,sen,efe}*.fasta

# Classify test reads  
./bin/kestrel-classify --verbose -d data/testdb -o test_results data/reads/Ecoli1_1000K_1.fastq.gz
```

 
## Requirements

- Nim >= 2.0
- Dependencies (installed via nimble):
  - docopt - Command line parsing
  - readfx - FASTA/FASTQ parsing with compression support

## License

(C) Quadram Institute Bioscience

MIT License
