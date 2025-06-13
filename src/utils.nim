import tables, sequtils, strutils, deques, sets, hashes

# Nucleotide encoding constants
const
  TOGGLE* = 0xe37e28c4271b5a2d'u64
  NUC_A* = 0'u64
  NUC_C* = 1'u64
  NUC_G* = 2'u64
  NUC_T* = 3'u64
  NUC_INVALID* = 4'u64

# Taxonomy types
type
  TaxonomyId* = uint32
  KmerHash* = uint64
  
  TaxonomyMap* = Table[string, TaxonomyId]
  TaxonomyLookup* = Table[TaxonomyId, string]
  TaxonomyGraph* = Table[TaxonomyId, TaxonomyId]  # child -> parent mapping
  
  DatabaseParams* = object
    kmerSize*: int
    minimizerSize*: int
    valueBits*: uint32
    taxonomyMap*: TaxonomyMap
    taxonomyLookup*: TaxonomyLookup
    lineageGraph*: TaxonomyGraph




proc isValidTaxonomy*(taxonomy: string): bool =
  ## Validates taxonomy strings for both GTDB and SILVA formats
  ## GTDB format: d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...
  ## SILVA format: k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...
  
  if taxonomy.len == 0:
    return false
  
  let levels = taxonomy.split(';')
  if levels.len == 0:
    return false
  
  # Valid prefixes for each taxonomic level
  const validPrefixes = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
  
  for i, level in levels:
    let trimmedLevel = level.strip()
    
    if trimmedLevel.len == 0:
      return false
    
    # Check if level starts with a valid prefix
    var hasValidPrefix = false
    
    # For first level, accept both d__ (GTDB) and k__ (SILVA)
    if i == 0:
      if trimmedLevel.startsWith("d__") or trimmedLevel.startsWith("k__"):
        hasValidPrefix = true
    else:
      # For other levels, check standard prefixes
      let expectedPrefix = case i:
        of 1: "p__"
        of 2: "c__"
        of 3: "o__"
        of 4: "f__"
        of 5: "g__"
        of 6: "s__"
        else: ""
      
      if expectedPrefix != "" and trimmedLevel.startsWith(expectedPrefix):
        hasValidPrefix = true
    
    if not hasValidPrefix:
      return false
    
    if trimmedLevel.len <= 3:
      return false
    
    let taxonName = trimmedLevel[3..^1].strip()
    
    if taxonName.len == 0:
      return false
    
    # More permissive character validation for SILVA compatibility
    # Allow: letters, numbers, spaces, underscores, hyphens, dots, parentheses, forward slashes, colons
    for c in taxonName:
      if not (c.isAlphaNumeric() or c in [' ', '_', '-', '.', '(', ')', '/', ':']):
        return false
  
  return true

# Nucleotide encoding lookup table
proc nucToNumber*(nuc: char): uint64 =
  case nuc.toUpperAscii()
  of 'A': NUC_A
  of 'C': NUC_C
  of 'G': NUC_G
  of 'T': NUC_T
  else: NUC_INVALID

# Reverse complement of a k-mer using bit manipulation
proc reverseComplement*(kmer: uint64, kmerSize: int): uint64 =
  var result = kmer
  let totalBits = kmerSize * 2
  
  # Complement all bits
  result = not result
  
  # Reverse the order of nucleotides (2-bit pairs)
  # Swap consecutive pairs
  result = ((result and 0x3333333333333333'u64) shl 2) or ((result and 0xCCCCCCCCCCCCCCCC'u64) shr 2)
  # Swap consecutive nibbles
  result = ((result and 0x0F0F0F0F0F0F0F0F'u64) shl 4) or ((result and 0xF0F0F0F0F0F0F0F0'u64) shr 4)
  # Swap consecutive bytes
  result = ((result and 0x00FF00FF00FF00FF'u64) shl 8) or ((result and 0xFF00FF00FF00FF00'u64) shr 8)
  # Swap consecutive 16-bit words
  result = ((result and 0x0000FFFF0000FFFF'u64) shl 16) or ((result and 0xFFFF0000FFFF0000'u64) shr 16)
  # Swap consecutive 32-bit words
  result = ((result and 0x00000000FFFFFFFF'u64) shl 32) or ((result and 0xFFFFFFFF00000000'u64) shr 32)
  
  # Shift to align with k-mer size
  result shr (64 - totalBits)

# Get canonical k-mer (lexicographically smaller of forward and reverse complement)
proc canonical*(kmer: uint64, kmerSize: int): uint64 =
  let revComp = reverseComplement(kmer, kmerSize)
  min(kmer, revComp)

# MurmurHash3 implementation
proc murmurHash3*(key: uint64): uint64 =
  var k = key
  k = k xor (k shr 33)
  k = k * 0xff51afd7ed558ccd'u64
  k = k xor (k shr 33)
  k = k * 0xc4ceb9fe1a85ec53'u64
  k = k xor (k shr 33)
  k

# Extract k-mers from a sequence using sliding window
proc extractKmers*(sequence: string, kmerSize: int): seq[uint64] =
  result = @[]
  if sequence.len < kmerSize:
    return
  
  var currentKmer = 0'u64
  let mask = (1'u64 shl (kmerSize * 2)) - 1
  
  # Initialize first k-mer
  var validKmerLen = 0
  for i in 0..<sequence.len:
    let nuc = nucToNumber(sequence[i])
    
    if nuc == NUC_INVALID:
      # Reset on invalid nucleotide
      currentKmer = 0
      validKmerLen = 0
      continue
    
    # Shift and add new nucleotide
    currentKmer = (currentKmer shl 2) or nuc
    currentKmer = currentKmer and mask
    validKmerLen += 1
    
    # When we have a complete k-mer
    if validKmerLen >= kmerSize:
      result.add(canonical(currentKmer, kmerSize))

# Sliding window minimizer extraction
proc extractMinimizers*(sequence: string, kmerSize: int, minimizerSize: int): seq[uint64] =
  result = @[]
  if sequence.len < kmerSize or minimizerSize == 0:
    return extractKmers(sequence, kmerSize)
  
  var window = initDeque[(uint64, int)]()
  var candidate = 0'u64
  let mask = (1'u64 shl (minimizerSize * 2)) - 1
  let toggle = TOGGLE and mask
  
  var validPos = 0
  
  for i in 0..<sequence.len:
    let nuc = nucToNumber(sequence[i])
    
    if nuc == NUC_INVALID:
      # Reset on invalid nucleotide
      candidate = 0
      validPos = 0
      window.clear()
      continue
    
    # Update candidate minimizer
    candidate = (candidate shl 2) or nuc
    candidate = candidate and mask
    validPos += 1
    
    if validPos >= minimizerSize:
      let minimizer = canonical(candidate, minimizerSize) xor toggle
      
      # Remove larger minimizers from back
      while window.len > 0 and window.peekLast()[0] >= minimizer:
        discard window.popLast()
      
      # Add current minimizer
      window.addLast((minimizer, i))
      
      # Remove expired minimizers from front
      let windowStart = i - kmerSize + minimizerSize + 1
      while window.len > 0 and window.peekFirst()[1] < windowStart:
        discard window.popFirst()
      
      # If we have a complete k-mer, record the minimizer
      if validPos >= kmerSize:
        if window.len > 0:
          result.add(window.peekFirst()[0] xor toggle)

# Parse taxonomy string and return individual levels
proc parseTaxonomyLevels*(taxonomy: string): seq[string] =
  taxonomy.split(';').mapIt(it.strip())

# Build taxonomy graph from taxonomy strings
proc buildTaxonomyGraph*(taxonomies: seq[string]): (TaxonomyMap, TaxonomyLookup, TaxonomyGraph) =
  var taxonomyMap = initTable[string, TaxonomyId]()
  var taxonomyLookup = initTable[TaxonomyId, string]()
  var lineageGraph = initTable[TaxonomyId, TaxonomyId]()
  var nextId = 1'u32
  
  # Add root
  taxonomyMap["root"] = 0'u32
  taxonomyLookup[0'u32] = "root"
  
  for taxonomy in taxonomies:
    let levels = parseTaxonomyLevels(taxonomy)
    var parentId = 0'u32  # root
    
    # First pass: build the individual taxonomy levels
    for level in levels:
      if level.len == 0:
        continue
        
      if not taxonomyMap.hasKey(level):
        taxonomyMap[level] = nextId
        taxonomyLookup[nextId] = level
        nextId += 1
      
      let currentId = taxonomyMap[level]
      if currentId != parentId:  # Don't link to self
        lineageGraph[currentId] = parentId
      parentId = currentId
    
    # Second pass: map the full taxonomy string to the most specific level (last level)
    if levels.len > 0:
      let mostSpecificLevel = levels[levels.len - 1]
      if taxonomyMap.hasKey(mostSpecificLevel):
        taxonomyMap[taxonomy] = taxonomyMap[mostSpecificLevel]
  
  (taxonomyMap, taxonomyLookup, lineageGraph)

# Find Lowest Common Ancestor of two taxonomy IDs
proc findLCA*(taxon1: TaxonomyId, taxon2: TaxonomyId, graph: TaxonomyGraph): TaxonomyId =
  if taxon1 == taxon2:
    return taxon1
  
  # Trace path from taxon1 to root
  var path1 = initHashSet[TaxonomyId]()
  var current = taxon1
  path1.incl(current)
  
  while graph.hasKey(current):
    current = graph[current]
    path1.incl(current)
  
  # Trace path from taxon2 until we find common ancestor
  current = taxon2
  if current in path1:
    return current
    
  while graph.hasKey(current):
    current = graph[current]
    if current in path1:
      return current
  
  # Return root if no common ancestor found
  return 0'u32

# Hash table operations for sepia mode
proc populateSepia*(key: uint64, value: TaxonomyId, valueBits: uint32): uint32 =
  var cell = (key shr 32).uint32
  cell = cell and (not ((1'u32 shl valueBits) - 1))
  cell = cell or value
  cell

proc getSepia*(key: uint64, table: seq[uint32], capacity: int, valueBits: uint32): TaxonomyId =
  let hc = murmurHash3(key)
  let compactedKey = (hc shr (32 + valueBits)).uint32
  
  var idx = (hc mod capacity.uint64).int
  let secondHash = ((hc shr 32) or 1).uint64
  
  for i in 0..<capacity:
    let cell = table[idx]
    let storedKey = cell shr valueBits
    
    if storedKey == compactedKey:
      return cell and ((1'u32 shl valueBits) - 1)
    elif cell == 0:
      return 0  # Not found
    
    # Quadratic probing
    idx = ((idx.uint64 + secondHash) mod capacity.uint64).int
  
  return 0  # Not found after full probe

