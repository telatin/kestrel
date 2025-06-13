# Package
version       = "0.1.0"
author        = "Andrea Telatin"
description   = "Klass k-mer classifier"
license       = "MIT"

# Dependencies
requires "nim >= 2.0",  "docopt#v0.7.1", "readfx", "iterutils", "argparse",  "colorize"

srcDir = "src"
binDir = "bin" 

namedBin = {
    "kbuild": "kestrel-build",
    "kclass": "kestrel-classify",
    "kombine": "kestrel-combine-seqs",
}.toTable()
