#!/usr/bin/env python3
"""
Convert FASTA headers from taxonomic names to sequential names with taxonomic comments.
"""

def convert_fasta_headers(input_file, output_file):
    """
    Convert FASTA headers from taxonomic format to sequential naming.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str): Path to output FASTA file
    """
    # Taxonomic rank prefixes in order
    rank_prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    
    sequence_counter = 1
    processed_lines = 0
    
    try:
        with open(input_file, 'r', buffering=8192) as infile, \
             open(output_file, 'w', buffering=8192) as outfile:
            
            for line in infile:
                line = line.strip()
                
                if line.startswith('>'):
                    # Remove the '>' and split by semicolon
                    taxonomy = line[1:].split(';')
                    
                    # Filter out empty strings
                    taxonomy = [tax.strip() for tax in taxonomy if tax.strip()]
                    
                    # Add rank prefixes to taxonomic levels
                    formatted_taxonomy = []
                    for i, taxon in enumerate(taxonomy):
                        if i < len(rank_prefixes):
                            formatted_taxonomy.append(f"{rank_prefixes[i]}{taxon}")
                        else:
                            # If more levels than prefixes, just add the taxon
                            formatted_taxonomy.append(taxon)
                    
                    # Create new header with sequential name and taxonomic comment
                    new_header = f">seq_{sequence_counter} {';'.join(formatted_taxonomy)}"
                    outfile.write(new_header + '\n')
                    sequence_counter += 1
                else:
                    # Write sequence lines as-is
                    outfile.write(line + '\n')
                
                processed_lines += 1
                
                # Progress indicator for large files
                if processed_lines % 10000 == 0:
                    print(f"Processed {processed_lines} lines, {sequence_counter-1} sequences...")
                    outfile.flush()  # Force write to disk
    
    except KeyboardInterrupt:
        print("\nOperation interrupted by user")
        raise
    except OSError as e:
        print(f"OS Error: {e}")
        print("This might be due to:")
        print("- Large file size causing timeout")
        print("- Network storage issues")
        print("- Insufficient disk space")
        print("- File permissions")
        raise

def main():
    """
    Main function to run the FASTA header conversion.
    """
    import argparse
    
    # Set up command line argument parser
    parser = argparse.ArgumentParser(
        description="Convert FASTA headers from taxonomic format to sequential naming"
    )
    parser.add_argument(
        "input_file", 
        help="Input FASTA file path"
    )
    parser.add_argument(
        "output_file", 
        help="Output FASTA file path"
    )
    
    # Parse command line arguments
    args = parser.parse_args()
    
    try:
        convert_fasta_headers(args.input_file, args.output_file)
        print(f"Successfully converted {args.input_file} to {args.output_file}")
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
