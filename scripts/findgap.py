#!/usr/bin/env python3
"""
Find gap positions in each chromosome of a FASTA file
Usage: python find_gaps.py input.fasta > gaps.bed
"""

import re
import sys
from Bio import SeqIO

def find_gaps_in_fasta(fasta_file, gap_char='N', min_gap_size=1):
    """
    Find gap positions in FASTA file
    
    Parameters:
    fasta_file: input FASTA file path
    gap_char: character representing gap (default 'N')
    min_gap_size: minimum gap size (default 1)
    """
    
    gaps_found = []
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id
            sequence = str(record.seq).upper()
            
            print(f"Analyzing sequence: {seq_id} (length: {len(sequence)})", file=sys.stderr)
            
            # use regular expression to find consecutive gap regions
            gap_pattern = re.compile(f'{re.escape(gap_char)}+')
            
            for match in gap_pattern.finditer(sequence):
                start = match.start()
                end = match.end()
                gap_length = end - start
                
                if gap_length >= min_gap_size:
                    gaps_found.append({
                        'chrom': seq_id,
                        'start': start,
                        'end': end,
                        'length': gap_length
                    })
                    
                    # output BED format
                    print(f"{seq_id}\t{start}\t{end}\tgap_{len(gaps_found)}\t{gap_length}")
    
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Exception occurred while processing file - {e}", file=sys.stderr)
        sys.exit(1)
    
    return gaps_found

def main():
    if len(sys.argv) != 2:
        print("Usage: python find_gaps.py input.fasta")
        print("Output: BED format gap position information")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    
    print(f"# Gap position report - Input file: {fasta_file}", file=sys.stderr)
    print(f"# Output format: chromosome\\tstart position\\tend position\\tgap name\\tgap length", file=sys.stderr)
    print("# BED format output:", file=sys.stderr)
    
    # output BED header
    print("track name=Gaps description=\"Genome gaps\" visibility=2")
    
    gaps = find_gaps_in_fasta(fasta_file)
    
    # statistics output to stderr
    print(f"\n# Statistics:", file=sys.stderr)
    print(f"# Total gaps found: {len(gaps)}", file=sys.stderr)
    
    if gaps:
        total_gap_length = sum(gap['length'] for gap in gaps)
        max_gap = max(gaps, key=lambda x: x['length'])
        min_gap = min(gaps, key=lambda x: x['length'])
        
        print(f"# Total gap length: {total_gap_length} bp", file=sys.stderr)
        print(f"# Maximum gap: {max_gap['chrom']}:{max_gap['start']}-{max_gap['end']} ({max_gap['length']} bp)", file=sys.stderr)
        print(f"# Minimum gap: {min_gap['chrom']}:{min_gap['start']}-{min_gap['end']} ({min_gap['length']} bp)", file=sys.stderr)
        print(f"# Average gap length: {total_gap_length/len(gaps):.1f} bp", file=sys.stderr)

if __name__ == "__main__":
    main()