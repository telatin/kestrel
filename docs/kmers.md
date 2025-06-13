# K-mers in Klass

This document explains how k-mers are implemented and used in the Klass taxonomic classification system.

## What are K-mers?

K-mers are subsequences of length k from a biological sequence. For example, the DNA sequence "ATCGATCG" contains the following 3-mers (k=3):
- ATC
- TCG
- CGA
- GAT
- ATC
- TCG

## K-mer Encoding

Klass uses an efficient 2-bit encoding scheme to represent nucleotides:

- A = 0 (00 in binary)
- C = 1 (01 in binary) 
- G = 2 (10 in binary)
- T = 3 (11 in binary)

Each k-mer is stored as a single 64-bit unsigned integer (`uint64`), allowing for fast bitwise operations and compact memory usage.

## Maximum K-mer Size

**The maximum k-mer size in Klass is 31.**

### Why 31?

This limit exists due to the 64-bit integer storage constraint:

1. **Storage calculation**: Each nucleotide requires 2 bits
2. **Maximum capacity**: 64 bits ÷ 2 bits per nucleotide = 32 nucleotides
3. **Safe limit**: 31 is used to avoid potential overflow in bit manipulation operations

### Code Location

The k-mer size validation occurs in `src/kbuild.nim`:

```nim
if config.kmerSize < 1 or config.kmerSize > 31:
  stderr.writeLine("Error: K-mer size must be between 1 and 31")
  return 1
```

## Canonical K-mers

Klass uses canonical k-mers to handle the double-stranded nature of DNA. For each k-mer, both the forward and reverse complement sequences are considered, and the lexicographically smaller one is stored.

For example:
- Forward: ATCG
- Reverse complement: CGAT
- Canonical: ATCG (lexicographically smaller)

This approach reduces memory usage by half and ensures that k-mers from both DNA strands are treated as equivalent.

## K-mer Extraction Process

1. **Sequence scanning**: The sequence is scanned with a sliding window of size k
2. **Nucleotide validation**: Invalid nucleotides (N, ambiguous bases) reset the k-mer counter
3. **Bit shifting**: New nucleotides are added using bit shifts: `currentKmer = (currentKmer << 2) | newNucleotide`
4. **Masking**: A bit mask ensures only the last k nucleotides are retained
5. **Canonicalization**: Each k-mer is converted to its canonical form before storage

## Default K-mer Size

The default k-mer size in Klass is **25**, which provides a good balance between:
- **Specificity**: Longer k-mers are more specific but may miss matches due to sequencing errors
- **Sensitivity**: Shorter k-mers are more sensitive but less specific
- **Performance**: 25-mers fit comfortably within the 64-bit constraint while maintaining good discriminatory power

## Impact on Classification

K-mer size affects classification accuracy:

- **Smaller k-mers** (k < 20): Higher sensitivity but lower specificity, may lead to false positives
- **Larger k-mers** (k > 30): Higher specificity but lower sensitivity, may miss true matches due to sequencing errors
- **Optimal range** (k = 20-30): Good balance for most microbial classification tasks

## Implementation Details

Key functions for k-mer handling are located in `src/utils.nim`:

- `extractKmers()`: Extracts all valid k-mers from a sequence
- `canonical()`: Converts a k-mer to its canonical form
- `reverseComplement()`: Computes the reverse complement of a k-mer
- `nucToNumber()`: Converts nucleotide characters to 2-bit encoding

The k-mer extraction is optimized for performance using bit manipulation operations rather than string operations, making it suitable for processing large genomic datasets.