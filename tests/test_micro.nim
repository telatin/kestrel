import unittest, os, strutils, osproc, json

suite "Micro database classification test":
  const testDataDir = "tests/data"
  const microDbPath = testDataDir / "micro_db.fa"
  const samplePath = testDataDir / "sample.fq.gz"
  const outputDbDir = "tests/temp_db"
  const outputPrefix = "tests/temp_results"
  
  setup:
    # Ensure test data exists
    doAssert fileExists(microDbPath), "Micro database file not found: " & microDbPath
    doAssert fileExists(samplePath), "Sample FASTQ file not found: " & samplePath
  
  teardown:
    # Clean up intermediate files
    if dirExists(outputDbDir):
      removeDir(outputDbDir)
    for ext in ["_classification.txt", "_summary.txt"]:
      let file = outputPrefix & ext
      if fileExists(file):
        removeFile(file)
  
  test "Build database with kmer size 29 and classify all reads":
    # Step 1: Build database with kmer size 29
    let buildCmd = "bin/kestrel-build -k 29 -o " & outputDbDir & " " & microDbPath & " --verbose"
    let buildResult = execCmdEx(buildCmd)
    
    check buildResult.exitCode == 0
    check dirExists(outputDbDir)
    check fileExists(outputDbDir / "kmers.bin")
    check fileExists(outputDbDir / "params.json")
    
    # Verify kmer size in params
    let paramsContent = readFile(outputDbDir / "params.json")
    let params = parseJson(paramsContent)
    check params["kmer_size"].getInt() == 29
    
    # Step 2: Classify reads
    let classifyCmd = "bin/kestrel-classify -d " & outputDbDir & " -o " & outputPrefix & " " & samplePath & " --verbose"
    let classifyResult = execCmdEx(classifyCmd)
    
    check classifyResult.exitCode == 0
    check fileExists(outputPrefix & "_classification.txt")
    check fileExists(outputPrefix & "_summary.txt")
    
    # Step 3: Verify all reads are classified
    let classificationContent = readFile(outputPrefix & "_classification.txt")
    let lines = classificationContent.strip().split('\n')
    
    # Skip header line, count classified reads
    var classifiedCount = 0
    var totalReads = 0
    
    for line in lines[1..^1]:  # Skip header
      if line.strip().len > 0:
        totalReads += 1
        let fields = line.split('\t')
        # Check if read has a taxonomic classification (not "no hits")
        if fields.len > 1 and fields[1] != "no hits":
          classifiedCount += 1
    
    # Verify all reads are classified
    check totalReads > 0
    check classifiedCount == totalReads
    
    echo "Successfully classified ", classifiedCount, " out of ", totalReads, " reads"
