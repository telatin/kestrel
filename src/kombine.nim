import readfx
import docopt, strutils
import os, json, streams, times
import ./utils
import sequtils, tables, sets
#[
  A program to generate a FASTA file starting from one or more FASTA file, 
  merging all the sequences that  have the same "comment" (in our case, taxonomy)
  separating them with a stretch of Ns.

]# 
const NimblePkgVersion {.strdefine.} = "undef"
 
const version = if NimblePkgVersion == "undef": "<prerelease>"
                else: NimblePkgVersion

 


proc main(): int =
  let args = docopt("""
  kestrel-combine-seqs: generate a taxonomy database from a set of sequences

  Usage: 
    kestrel-combine-seqs [options] [-o <FASTA_OUT>] <FASTA_FILES>...

  Options:
    -o, --output FASTA_OUT     Output file, STDOUT if not specified
    -p, --prefix STRING        Output sequence prefix [default: "seq."]
    -n, --n-length INT         Length of Ns to separate sequences [default: 32]
    -i, --ignore-taxonomy      Combine sequences by comment even if not a valid taxonomy
  
  Other options:
    --verbose               Print verbose log
    --debug                 Print debug log  
    --help                  Show help
  """, version="0.1.0", argv=commandLineParams())
 
  # Parse arguments
  let verbose = bool(args["--verbose"])
  let debug = bool(args["--debug"])
  let ignoreTaxonomy = bool(args["--ignore-taxonomy"])
  let outputFile = if args["--output"]: $args["--output"] else: ""
  let seqPrefix = ($args["--prefix"]).strip(chars = {'"'})
  let nLength = parseInt($args["--n-length"])

  # Validate input files
  let fastaFiles = @(args["<FASTA_FILES>"])
  for fastaFile in fastaFiles:
    if not fileExists(fastaFile):
      stderr.writeLine("File not found: ", fastaFile)
      return 1

  try:
    # Prepare a table to collect sequences by taxonomy [taxonomy -> [sequences]]
    var taxonomySequences = initTable[string, seq[string]]()
    var processedSeqs = 0
    var skippedSeqs = 0
    
    # Read files and collect sequences by taxonomy
    for fastaFile in fastaFiles:
      if verbose:
        stderr.writeLine("Processing file: ", fastaFile)
      
      for record in readfq(fastaFile):
        let taxonomy = record.comment.strip()
        processedSeqs += 1
        
        # Validate taxonomy if required
        if not ignoreTaxonomy:
          if not isValidTaxonomy(taxonomy):
            if debug:
              stderr.writeLine("Invalid taxonomy in record: ", record.name, " - ", taxonomy)
            skippedSeqs += 1
            continue
        
        # Add sequence to taxonomy group
        if not taxonomySequences.hasKey(taxonomy):
          taxonomySequences[taxonomy] = @[]
        
        taxonomySequences[taxonomy].add(record.sequence)
        
        if debug and processedSeqs mod 10000 == 0:
          stderr.writeLine("Processed ", processedSeqs, " sequences, found ", taxonomySequences.len, " taxonomies")
    
    if verbose:
      stderr.writeLine("Total sequences processed: ", processedSeqs)
      stderr.writeLine("Sequences skipped (invalid taxonomy): ", skippedSeqs)
      stderr.writeLine("Unique taxonomies found: ", taxonomySequences.len)
    
    # Prepare output stream
    let outputStream = if outputFile != "":
                        newFileStream(outputFile, fmWrite)
                       else:
                        newFileStream(stdout)
    
    # Generate N-string separator
    let nSeparator = "N".repeat(nLength)
    var seqCounter = 1
    
    # Write combined sequences
    for taxonomy, sequences in taxonomySequences:
      if sequences.len == 0:
        continue
      
      # Create sequence ID with count of sequences
      let seqId = seqPrefix & $seqCounter & "_" & $sequences.len
      seqCounter += 1
      
      # Combine all sequences for this taxonomy with N separators
      let combinedSequence = sequences.join(nSeparator)
      
      # Write FASTA record
      outputStream.writeLine(">" & seqId & " " & taxonomy)
      outputStream.writeLine(combinedSequence)
      
      if debug:
        stderr.writeLine("Combined ", sequences.len, " sequences for taxonomy: ", taxonomy, " (length: ", combinedSequence.len, ")")
    
    if outputFile != "":
      outputStream.close()
    
    if verbose:
      stderr.writeLine("Output written to: ", if outputFile != "": outputFile else: "STDOUT")
      stderr.writeLine("Combined sequences created: ", seqCounter - 1)
    
    return 0
    
  except Exception as e:
    stderr.writeLine("Error during sequence combination: ", e.msg)
    return 1

when isMainModule:
  quit(main())