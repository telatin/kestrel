import readfx
import docopt, strutils
import os, json, streams, times
import ./utils
import sequtils, tables, sets, options

# Helper function to format numbers with thousand separators
proc formatWithCommas(num: int): string =
  let numStr = $num
  result = ""
  let len = numStr.len
  for i, c in numStr:
    if i > 0 and (len - i) mod 3 == 0:
      result.add(",")
    result.add(c)

const NimblePkgVersion {.strdefine.} = "undef"
 
const version = if NimblePkgVersion == "undef": "<prerelease>"
                else: NimblePkgVersion

# Configuration type
type
  BuildDbConfig = object
    kmerSize: int
    minimizerSize: int
    kmerShape: Option[KmerShape]
    verbose: bool
    debug: bool
    outputDir: string

 
proc collectTaxonomies(fastaFiles: seq[string], config: BuildDbConfig): seq[string] =
  ## Collect all unique taxonomies from FASTA headers
  result = @[]
  var taxonomySet = initHashSet[string]()
  
  for fastaFile in fastaFiles:
    if config.verbose:
      stderr.writeLine("Processing file: ", fastaFile)
    
    for record in readfq(fastaFile):
      let taxonomy = record.comment.strip()
      
      if not isValidTaxonomy(taxonomy):
        stderr.writeLine("Warning: Invalid taxonomy in ", fastaFile, " for ", record.name, ": ", taxonomy)
        continue
      
      if taxonomy notin taxonomySet:
        taxonomySet.incl(taxonomy)
        result.add(taxonomy)
        
        if config.debug:
          stderr.writeLine("Found taxonomy: ", taxonomy)

proc buildKmerDatabase(fastaFiles: seq[string], config: BuildDbConfig, dbParams: DatabaseParams): Table[uint64, TaxonomyId] =
  ## Build k-mer to taxonomy mapping
  result = initTable[uint64, TaxonomyId]()
  var kmerConflicts = initTable[uint64, HashSet[TaxonomyId]]()
  
  for fastaFile in fastaFiles:
    if config.verbose:
      stderr.writeLine("Extracting k-mers from: ", fastaFile)
    
    var sequenceCount = 0
    var processedSequences = 0
    var skippedSequences = 0
    var lastProgressTime = cpuTime()
    var lastProgressPercent = -1  # Track last printed percentage
    
    # First pass: count total sequences for progress tracking
    if config.verbose:
      stderr.write("Counting sequences... ")
      for record in readfq(fastaFile):
        sequenceCount += 1
      stderr.writeLine(sequenceCount, " sequences found")
    
    # Second pass: process sequences with progress updates
    for record in readfq(fastaFile):
      let taxonomy = record.comment.strip()
      processedSequences += 1
      
      # Progress update only when percentage increases by at least 1% or every 10 seconds
      if config.verbose:
        let currentTime = cpuTime()
        let progress = if sequenceCount > 0: (processedSequences * 100) div sequenceCount else: 0
        
        if (progress >= lastProgressPercent + 10) or (currentTime - lastProgressTime) > 30.0:
          stderr.writeLine("  ", progress, "% - ", formatWithCommas(result.len), " unique k-mers so far (", 
                          formatWithCommas(processedSequences), "/", formatWithCommas(sequenceCount), ")")
          lastProgressTime = currentTime
          lastProgressPercent = progress
      
      if not dbParams.taxonomyMap.hasKey(taxonomy):
        if config.debug:
          stderr.writeLine("Warning: Unknown taxonomy: ", taxonomy)
        skippedSequences += 1
        continue
      
      let taxonomyId = dbParams.taxonomyMap[taxonomy]
      let sequence = record.sequence
      
      # Extract k-mers, minimizers, or shaped k-mers
      let kmers = if dbParams.kmerShape.isSome:
                    extractKmersWithShape(sequence, dbParams.kmerShape.get())
                  elif config.minimizerSize > 0:
                    extractMinimizers(sequence, config.kmerSize, config.minimizerSize)
                  else:
                    extractKmers(sequence, config.kmerSize)
      
      # Process each k-mer
      for kmer in kmers:
        if result.hasKey(kmer):
          # Handle conflicts using LCA
          let existingTaxonId = result[kmer]
          if existingTaxonId != taxonomyId:
            let lcaTaxonId = findLCA(existingTaxonId, taxonomyId, dbParams.lineageGraph)
            result[kmer] = lcaTaxonId
            
            # Track conflicts for reporting
            if not kmerConflicts.hasKey(kmer):
              kmerConflicts[kmer] = initHashSet[TaxonomyId]()
            kmerConflicts[kmer].incl(existingTaxonId)
            kmerConflicts[kmer].incl(taxonomyId)
        else:
          result[kmer] = taxonomyId
    
    if config.verbose:
      stderr.writeLine("File completed: ", fastaFile)
      stderr.writeLine("  Processed: ", processedSequences, " sequences")
      stderr.writeLine("  Skipped: ", skippedSequences, " sequences (invalid taxonomy)")
      stderr.writeLine("  Current k-mers: ", result.len)
  
  if config.verbose:
    stderr.writeLine("Total k-mers: ", result.len)
    stderr.writeLine("K-mer conflicts resolved: ", kmerConflicts.len)

proc saveDatabase(kmerDb: Table[uint64, TaxonomyId], dbParams: DatabaseParams, config: BuildDbConfig) =
  ## Save database to files
  let paramsFile = config.outputDir / "params.json"
  let dbFile = config.outputDir / "kmers.bin"
  let lineageFile = config.outputDir / "lineage.bin"
  
  # Save parameters as JSON
  var paramsJson = %* {
    "kmer_size": config.kmerSize,
    "minimizer_size": config.minimizerSize,
    "value_bits": 24,  # Use 24 bits for taxonomy IDs (supports ~16M taxa)
    "num_kmers": kmerDb.len,
    "taxonomies": {}
  }
  
  # Add kmer shape if present
  if dbParams.kmerShape.isSome:
    let shape = dbParams.kmerShape.get()
    paramsJson["kmer_shape"] = %* {
      "pattern": shape.pattern,
      "window_size": shape.windowSize
    }
  
  # Add taxonomy mappings
  for taxon, id in dbParams.taxonomyMap:
    paramsJson["taxonomies"][taxon] = %id
  
  # Write parameters
  let paramsStream = newFileStream(paramsFile, fmWrite)
  paramsStream.write(paramsJson.pretty())
  paramsStream.close()
  
  # Save lineage graph in binary format
  let lineageStream = newFileStream(lineageFile, fmWrite)
  
  # Write number of lineage entries
  lineageStream.write(dbParams.lineageGraph.len.uint64)
  
  # Write child -> parent mappings
  for child, parent in dbParams.lineageGraph:
    lineageStream.write(child)
    lineageStream.write(parent)
  
  lineageStream.close()
  
  # Save k-mer database in binary format
  let dbStream = newFileStream(dbFile, fmWrite)
  
  # Write header: number of k-mers
  dbStream.write(kmerDb.len.uint64)
  
  # Write k-mer -> taxonomy pairs
  for kmer, taxonomyId in kmerDb:
    dbStream.write(kmer)
    dbStream.write(taxonomyId)
  
  dbStream.close()
  
  if config.verbose:
    stderr.writeLine("Database saved to: ", config.outputDir)
    stderr.writeLine("  Parameters: ", paramsFile)
    stderr.writeLine("  Lineage graph: ", lineageFile)
    stderr.writeLine("  K-mers: ", dbFile)

proc main(): int =
  let args = docopt("""
  klass-build: generate a taxonomy database from a set of sequences

  Usage: 
    klass-build [options] -o <database> <FASTA_FILES>...

  Options:
    -o, --output DB         Output database directory
    -k, --kmer-size INT     K-mer size [default: 25]
    -s, --kmer-shape STRING K-mer shape pattern (e.g., OOOO--OOOO)
    -m, --minimizer-size INT Minimizer size (0 = no minimizers) [default: 0]
  
  Other options:
    --verbose               Print verbose log
    --debug                 Print debug log  
    --help                  Show help
  """, version="0.1.0", argv=commandLineParams())

  # Parse command line arguments
  var config = BuildDbConfig(
    debug: bool(args["--debug"]),
    verbose: bool(args["--verbose"]),
    outputDir: $args["--output"],
    minimizerSize: if args["--minimizer-size"]: parseInt($args["--minimizer-size"]) else: 0,
    kmerShape: none(KmerShape)
  )
  
  # Handle k-mer size or shape - only one should be specified
  let hasKmerSize = args["--kmer-size"] and $args["--kmer-size"] != "25"  # Check if user provided non-default value
  let hasKmerShape = args["--kmer-shape"]
  
  if hasKmerSize and hasKmerShape:
    stderr.writeLine("Error: Cannot specify both --kmer-size and --kmer-shape")
    return 1
  
  if hasKmerShape:
    try:
      config.kmerShape = some(parseKmerShape($args["--kmer-shape"]))
      config.kmerSize = config.kmerShape.get().kmerSize
    except ValueError as e:
      stderr.writeLine("Error parsing kmer shape: ", e.msg)
      return 1
  else:
    config.kmerSize = parseInt($args["--kmer-size"])

  # Validate k-mer size
  if config.kmerSize < 1 or config.kmerSize > 31:
    stderr.writeLine("Error: K-mer size must be between 1 and 31")
    return 1
  
  if config.minimizerSize > 0 and config.minimizerSize >= config.kmerSize:
    stderr.writeLine("Error: Minimizer size must be smaller than k-mer size")
    return 1

  if config.verbose:
    # Print kmer size, minimiser size
    for file in @(args["<FASTA_FILES>"]):
      stderr.writeLine("Input file:       ", args["<FASTA_FILES>"])
    if config.kmerShape.isSome:
      let shape = config.kmerShape.get()
      stderr.writeLine("K-mer shape:      ", shape.pattern)
      stderr.writeLine("K-mer size:       ", shape.kmerSize)
      stderr.writeLine("Window size:      ", shape.windowSize)
    else:
      stderr.writeLine("K-mer size:       ", config.kmerSize)
      stderr.writeLine("Window size:      ", config.kmerSize)  # For traditional k-mers, window size = k-mer size
    stderr.writeLine("Minimizer size:   ", config.minimizerSize)
    stderr.writeLine("Output directory: ", config.outputDir)
  # Create output directory
  if not dirExists(config.outputDir):
    try:
      createDir(config.outputDir)
      if config.verbose:
        stderr.writeLine("Created output directory: ", config.outputDir)
    except OSError as e:
      stderr.writeLine("Error creating output directory: ", e.msg)
      return 1

  # Validate input files
  let fastaFiles = @(args["<FASTA_FILES>"])
  for fastaFile in fastaFiles:
    if not fileExists(fastaFile):
      stderr.writeLine("File not found: ", fastaFile)
      return 1

  try:
    # Step 1: Collect taxonomies from all files
    if config.verbose:
      stderr.writeLine("Step 1: Collecting taxonomies...")
    
    let taxonomies = collectTaxonomies(fastaFiles, config)
    if taxonomies.len == 0:
      stderr.writeLine("Error: No valid taxonomies found")
      return 1
    
    if config.verbose:
      stderr.writeLine("Found ", taxonomies.len, " unique taxonomies")

    # Step 2: Build taxonomy graph
    if config.verbose:
      stderr.writeLine("Step 2: Building taxonomy graph...")
    
    let (taxonomyMap, taxonomyLookup, lineageGraph) = buildTaxonomyGraph(taxonomies)
    
    if config.debug:
      stderr.writeLine("Taxonomy map has ", taxonomyMap.len, " entries")
      for taxon, id in taxonomyMap:
        stderr.writeLine("  ", taxon, " -> ", id)
    
    let dbParams = DatabaseParams(
      kmerSize: config.kmerSize,
      minimizerSize: config.minimizerSize,
      valueBits: 24,
      taxonomyMap: taxonomyMap,
      taxonomyLookup: taxonomyLookup,
      lineageGraph: lineageGraph,
      kmerShape: config.kmerShape
    )

    # Step 3: Build k-mer database
    if config.verbose:
      stderr.writeLine("Step 3: Building k-mer database...")
    
    let kmerDb = buildKmerDatabase(fastaFiles, config, dbParams)
    
    if kmerDb.len == 0:
      stderr.writeLine("Error: No k-mers extracted")
      return 1

    # Step 4: Save database
    if config.verbose:
      stderr.writeLine("Step 4: Saving database...")
    
    saveDatabase(kmerDb, dbParams, config)
    
    if config.verbose:
      stderr.writeLine("Database build completed successfully!")
    
    return 0
    
  except Exception as e:
    stderr.writeLine("Error during database building: ", e.msg)
    return 1

when isMainModule:
  quit(main())