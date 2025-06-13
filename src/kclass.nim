import readfx
import docopt, strutils
import os, json, streams, tables
import ./utils
import sequtils, math, algorithm

const NimblePkgVersion {.strdefine.} = "undef"
 
const version = if NimblePkgVersion == "undef": "<prerelease>"
                else: NimblePkgVersion

# Configuration type
type
  ClassifyConfig = object
    verbose: bool
    debug: bool
    dbDir: string
    outputPrefix: string
    qualityThreshold: int
    minKmerHits: int

  Database = object
    params: DatabaseParams
    kmerTable: Table[uint64, TaxonomyId]

  ClassificationResult = object
    readId: string
    taxonomy: string
    hitCount: int
    totalKmers: int
    confidence: float

proc loadDatabase(dbDir: string, verbose: bool): Database =
  ## Load database from files
  let paramsFile = dbDir / "params.json"
  let dbFile = dbDir / "kmers.bin"
  let lineageFile = dbDir / "lineage.bin"
  
  if not fileExists(paramsFile):
    raise newException(IOError, "Parameters file not found: " & paramsFile)
  
  if not fileExists(dbFile):
    raise newException(IOError, "Database file not found: " & dbFile)
    
  if not fileExists(lineageFile):
    raise newException(IOError, "Lineage file not found: " & lineageFile)
  
  # Load parameters
  if verbose:
    stderr.writeLine("Loading parameters from: ", paramsFile)
  
  let paramsJson = parseJson(readFile(paramsFile))
  
  let kmerSize = paramsJson["kmer_size"].getInt()
  let minimizerSize = paramsJson["minimizer_size"].getInt()
  let valueBits = paramsJson["value_bits"].getInt().uint32
  
  # Build taxonomy maps
  var taxonomyMap = initTable[string, TaxonomyId]()
  var taxonomyLookup = initTable[TaxonomyId, string]()
  
  for taxon, idNode in paramsJson["taxonomies"]:
    let id = idNode.getInt().TaxonomyId
    taxonomyMap[taxon] = id
    taxonomyLookup[id] = taxon
  
  # Load pre-built lineage graph
  if verbose:
    stderr.writeLine("Loading lineage graph from: ", lineageFile)
  
  let lineageStream = newFileStream(lineageFile, fmRead)
  let numLineageEntries = lineageStream.readUint64()
  
  var lineageGraph = initTable[TaxonomyId, TaxonomyId]()
  
  for i in 0..<numLineageEntries:
    let child = lineageStream.readUint32()
    let parent = lineageStream.readUint32()
    lineageGraph[child] = parent
  
  lineageStream.close()
  
  let params = DatabaseParams(
    kmerSize: kmerSize,
    minimizerSize: minimizerSize,
    valueBits: valueBits,
    taxonomyMap: taxonomyMap,
    taxonomyLookup: taxonomyLookup,
    lineageGraph: lineageGraph
  )
  
  # Load k-mer table
  if verbose:
    stderr.writeLine("Loading k-mer database from: ", dbFile)
  
  let dbStream = newFileStream(dbFile, fmRead)
  let numKmers = dbStream.readUint64()
  
  var kmerTable = initTable[uint64, TaxonomyId]()
  
  for i in 0..<numKmers:
    let kmer = dbStream.readUint64()
    let taxonomyId = dbStream.readUint32()
    kmerTable[kmer] = taxonomyId
  
  dbStream.close()
  
  if verbose:
    stderr.writeLine("Loaded ", kmerTable.len, " k-mers")
    stderr.writeLine("Database loaded successfully")
  
  result = Database(params: params, kmerTable: kmerTable)

proc filterSequenceByQuality(sequence: string, quality: string, threshold: int): string =
  ## Replace low-quality bases with 'N'
  if quality.len != sequence.len:
    return sequence  # No quality filtering if lengths don't match
  
  result = ""
  for i in 0..<sequence.len:
    let phredScore = ord(quality[i]) - 33  # Convert from Phred+33
    if phredScore >= threshold:
      result.add(sequence[i])
    else:
      result.add('N')

proc classifySequence(sequence: string, database: Database, config: ClassifyConfig): ClassificationResult =
  ## Classify a single sequence
  let params = database.params
  
  # Extract k-mers or minimizers
  let kmers = if params.minimizerSize > 0:
                extractMinimizers(sequence, params.kmerSize, params.minimizerSize)
              else:
                extractKmers(sequence, params.kmerSize)
  
  if kmers.len == 0:
    return ClassificationResult(
      taxonomy: "no hits",
      hitCount: 0,
      totalKmers: 0,
      confidence: 0.0
    )
  
  # Count hits per taxonomy
  var hitCounts = initTable[TaxonomyId, int]()
  var totalHits = 0
  
  for kmer in kmers:
    if database.kmerTable.hasKey(kmer):
      let taxonomyId = database.kmerTable[kmer]
      if not hitCounts.hasKey(taxonomyId):
        hitCounts[taxonomyId] = 0
      hitCounts[taxonomyId] += 1
      totalHits += 1
  
  if totalHits < config.minKmerHits:
    return ClassificationResult(
      taxonomy: "no hits",
      hitCount: totalHits,
      totalKmers: kmers.len,
      confidence: 0.0
    )
  
  # Find best classification
  var bestTaxonomyId: TaxonomyId = 0
  var bestCount = 0
  
  for taxonomyId, count in hitCounts:
    if count > bestCount:
      bestCount = count
      bestTaxonomyId = taxonomyId
  
  # Handle ties using LCA
  var tiedTaxa: seq[TaxonomyId] = @[]
  for taxonomyId, count in hitCounts:
    if count == bestCount:
      tiedTaxa.add(taxonomyId)
  
  if tiedTaxa.len > 1:
    # Find LCA of all tied taxa
    bestTaxonomyId = tiedTaxa[0]
    for i in 1..<tiedTaxa.len:
      bestTaxonomyId = findLCA(bestTaxonomyId, tiedTaxa[i], params.lineageGraph)
  
  # Get taxonomy string
  let taxonomy = if params.taxonomyLookup.hasKey(bestTaxonomyId):
                   params.taxonomyLookup[bestTaxonomyId]
                 else:
                   "unclassified"
  
  let confidence = bestCount.float / kmers.len.float
  
  result = ClassificationResult(
    taxonomy: taxonomy,
    hitCount: bestCount,
    totalKmers: kmers.len,
    confidence: confidence
  )

proc processReads(inputFiles: seq[string], database: Database, config: ClassifyConfig) =
  ## Process input reads and classify them
  let classificationFile = config.outputPrefix & "_classification.txt"
  let summaryFile = config.outputPrefix & "_summary.txt"
  
  var classificationStream = newFileStream(classificationFile, fmWrite)
  var results: seq[ClassificationResult] = @[]
  var taxonomyCounts = initTable[string, int]()
  var taxonomyBases = initTable[string, int]()
  var taxonomyScores = initTable[string, float]()
  var processedReads = 0
  
  for inputFile in inputFiles:
    if config.verbose:
      stderr.writeLine("Processing: ", inputFile)
    
    for record in readfq(inputFile):
      processedReads += 1
      
      # Progress update every 50,000 reads
      if config.verbose and processedReads mod 50000 == 0:
        stderr.write("\rProcessed ", processedReads, " reads...")
        stderr.flushFile()
      var sequence = record.sequence
      
      # Apply quality filtering if quality scores are available
      if record.quality.len > 0:
        sequence = filterSequenceByQuality(sequence, record.quality, config.qualityThreshold)
      
      # Classify sequence
      var result = classifySequence(sequence, database, config)
      result.readId = record.name
      
      # Write classification result
      classificationStream.writeLine(result.readId, "\t", result.taxonomy, "\t", 
                                   result.hitCount, "\t", result.totalKmers)
      
      # Update summary statistics
      if not taxonomyCounts.hasKey(result.taxonomy):
        taxonomyCounts[result.taxonomy] = 0
        taxonomyBases[result.taxonomy] = 0
        taxonomyScores[result.taxonomy] = 0.0
      
      taxonomyCounts[result.taxonomy] += 1
      taxonomyBases[result.taxonomy] += sequence.len
      taxonomyScores[result.taxonomy] += result.confidence
      
      results.add(result)
  
  # Clear progress line if we printed any
  if config.verbose and processedReads >= 50000:
    stderr.write("\r")
    stderr.flushFile()
  
  classificationStream.close()
  
  # Write summary file
  if config.verbose:
    stderr.writeLine("Writing summary to: ", summaryFile)
  
  var summaryStream = newFileStream(summaryFile, fmWrite)
  summaryStream.writeLine("Taxonomy\tReads\tAvg_Score\tTotal_Bases")
  
  # Sort by read count (descending)
  var sortedTaxa: seq[(string, int)] = @[]
  for taxonomy, count in taxonomyCounts:
    sortedTaxa.add((taxonomy, count))
  
  sortedTaxa.sort(proc(a, b: (string, int)): int = cmp(b[1], a[1]), SortOrder.Ascending)
  
  for (taxonomy, count) in sortedTaxa:
    let avgScore = if count > 0: taxonomyScores[taxonomy] / count.float else: 0.0
    let totalBases = taxonomyBases[taxonomy]
    summaryStream.writeLine(taxonomy, "\t", count, "\t", avgScore.formatFloat(ffDecimal, 3), "\t", totalBases)
  
  summaryStream.close()
  
  if config.verbose:
    stderr.writeLine("Classification completed!")
    stderr.writeLine("  Total reads processed: ", results.len)
    stderr.writeLine("  Classification file: ", classificationFile)
    stderr.writeLine("  Summary file: ", summaryFile)

proc main(): int =
  let args = docopt("""
  kestrel-classify: classify sequences using a k-mer database made with kestrel-build

  Usage: 
    kestrel-classify [options] -d <database> -o <output> <FASTQ_FILES>...

  Parameters:
    -d, --db DIR            K-mer database directory
    -o, --output PREFIX     Output prefix for results
    
  Options:
    -q, --quality INT       Minimum quality score (Phred) [default: 15]
    -m, --min-hits INT      Minimum k-mer hits for classification [default: 3]
    
  Other options:
    --verbose               Print verbose log
    --debug                 Print debug log  
    --help                  Show help
  """, version="0.1.0", argv=commandLineParams())

  # Parse command line arguments
  let config = ClassifyConfig(
    verbose: bool(args["--verbose"]),
    debug: bool(args["--debug"]),
    dbDir: $args["--db"],
    outputPrefix: $args["--output"],
    qualityThreshold: parseInt($args["--quality"]),
    minKmerHits: parseInt($args["--min-hits"])
  )

  # Validate database directory
  if not dirExists(config.dbDir):
    stderr.writeLine("Database directory not found: ", config.dbDir)
    return 1

  # Validate input files
  let inputFiles = @(args["<FASTQ_FILES>"])
  for inputFile in inputFiles:
    if not fileExists(inputFile):
      stderr.writeLine("File not found: ", inputFile)
      return 1

  try:
    # Load database
    if config.verbose:
      stderr.writeLine("Loading database from: ", config.dbDir)
    
    let database = loadDatabase(config.dbDir, config.verbose)
    
    # Process reads
    processReads(inputFiles, database, config)
    
    return 0
    
  except Exception as e:
    stderr.writeLine("Error during classification: ", e.msg)
    return 1

when isMainModule:
  quit(main())