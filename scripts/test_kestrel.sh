#!/bin/bash

# Test script for Klass taxonomy classifier
# This script builds a database and classifies test reads

set -e  # Exit on any error

if [[ ! -e ./bin/kestrel-build ]]; then
  echo "This script is meant to be run from the repository root."
  exit 1
fi

echo "=== Klass Test Script ==="
echo "Testing database building and read classification"
echo

# Clean up any previous test results
echo "Cleaning up previous test results..."
rm -rf data/testdb test_results_*
echo

# Step 1: Build the binaries
echo "Step 1: Building Klass binaries..."
nimble build
if [ $? -eq 0 ]; then
    echo "✓ Build successful"
else
    echo "✗ Build failed"
    exit 1
fi
echo

# Step 2: Build test database
echo "Step 2: Building test database..."
echo "Using reference genomes: data/ref/*.fasta"
echo "K-mer size: 21 (small for faster testing)"

./bin/kestrel-build --verbose -o data/testdb -k 21 data/ref/*.fasta

if [ $? -eq 0 ]; then
    echo "✓ Database build successful"
    echo "Database files:"
    ls -lh data/testdb/
else
    echo "✗ Database build failed"
    exit 1
fi
echo

# Step 3: Test classification on single file
echo "Step 3: Testing classification on single FASTQ file..."
echo "Classifying: data/reads/Ecoli1_1000K_1.fastq.gz"

./bin/kestrel-classify --verbose -d data/testdb -o test_results_single data/reads/Ecoli1_1000K_1.fastq.gz

if [ $? -eq 0 ]; then
    echo "✓ Single file classification successful"
    echo
    echo "Results preview:"
    echo "--- Classification file (first 5 lines) ---"
    head -5 test_results_single_classification.txt
    echo
    echo "--- Summary file ---"
    cat test_results_single_summary.txt
else
    echo "✗ Single file classification failed"
    exit 1
fi
echo

# Step 4: Test classification on multiple files (if they exist)
echo "Step 4: Testing classification on multiple files..."
if [ -f "data/reads/Ecoli2_1000K_1.fastq.gz" ]; then
    echo "Classifying multiple E. coli files..."
    
    ./bin/kestrel-classify --verbose -d data/testdb -o test_results_multi \
        data/reads/Ecoli1_1000K_1.fastq.gz \
        data/reads/Ecoli2_1000K_1.fastq.gz
    
    if [ $? -eq 0 ]; then
        echo "✓ Multiple file classification successful"
        echo
        echo "Multi-file results summary:"
        head -10 test_results_multi_summary.txt
    else
        echo "✗ Multiple file classification failed"
        exit 1
    fi
else
    echo "Skipping multiple file test (additional files not found)"
fi
echo

# Step 5: Performance summary
echo "Step 5: Performance Summary"
echo "=========================="

# Database stats
if [ -f "data/testdb/params.json" ]; then
    echo "Database statistics:"
    grep -E '"kmer_size"|"num_kmers"' data/testdb/params.json | sed 's/^/  /'
    echo "  Database size: $(du -h data/testdb/kmers.bin | cut -f1)"
fi

# Classification stats
if [ -f "test_results_single_summary.txt" ]; then
    echo
    echo "Classification statistics:"
    total_reads=$(tail -n +2 test_results_single_summary.txt | awk '{sum += $2} END {print sum}')
    classified_reads=$(tail -n +2 test_results_single_summary.txt | awk '$1 != "no hits" {sum += $2} END {print sum}')
    
    if [ ! -z "$total_reads" ] && [ ! -z "$classified_reads" ]; then
        classification_rate=$(echo "scale=1; $classified_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "N/A")
        echo "  Total reads processed: $total_reads"
        echo "  Successfully classified: $classified_reads"
        echo "  Classification rate: ${classification_rate}%"
    fi
fi

echo
echo "=== Test Completed Successfully! ==="
echo
echo "Generated files:"
echo "  Database: data/testdb/"
echo "  Single file results: test_results_single_*"
if [ -f "test_results_multi_summary.txt" ]; then
    echo "  Multi file results: test_results_multi_*"
fi
echo
echo "To test with your own data:"
echo "  1. Build database: ./bin/kestrel-build -o mydb -k 25 your_genomes.fasta"
echo "  2. Classify reads: ./bin/kestrel-classify -d mydb -o results your_reads.fastq.gz"
