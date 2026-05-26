#!/usr/bin/env python3
"""
Find gap positions in each chromosome of a FASTA file
Usage: findgap.py input.fasta [--min-size SIZE] [--gap-char CHAR] [--output OUTPUT]
"""

import re
import sys
import argparse
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
    
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Exception occurred while processing file - {e}", file=sys.stderr)
        sys.exit(1)
    
    return gaps_found

def main():
    parser = argparse.ArgumentParser(
        description='Find gap positions in FASTA file',
        epilog='Example: findgap.py genome.fasta --min-size 10 --output gaps.bed'
    )
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('--min-size', '-m', type=int, default=1,
                        help='Minimum gap size (default: 1)')
    parser.add_argument('--gap-char', '-c', default='N',
                        help='Gap character (default: N)')
    parser.add_argument('--output', '-o', help='Output BED file (default: stdout)')
    
    args = parser.parse_args()
    
    # Print header to stderr
    print(f"# Gap position report - Input file: {args.input}", file=sys.stderr)
    print(f"# Gap character: {args.gap_char}, Minimum size: {args.min_size}", file=sys.stderr)
    print("# Output format: chromosome\\tstart\\tend\\tgap_name\\tgap_length", file=sys.stderr)
    
    # Find gaps
    gaps = find_gaps_in_fasta(args.input, args.gap_char, args.min_size)
    
    # Prepare output (修改后的逻辑)
    output_lines = []
    output_lines.append("track name=Gaps description=\"Genome gaps\" visibility=2")
    
    # 按染色体分组gap
    gaps_by_chrom = {}
    for gap in gaps:
        chrom = gap['chrom']
        if chrom not in gaps_by_chrom:
            gaps_by_chrom[chrom] = []
        gaps_by_chrom[chrom].append(gap)
    
    # 按染色体名称排序（保证输出有序）
    sorted_chroms = sorted(gaps_by_chrom.keys())
    
    # 按染色体内部计数生成命名
    for chrom in sorted_chroms:
        chrom_gaps = gaps_by_chrom[chrom]
        # 按gap的起始位置排序，保证序号按物理位置递增
        chrom_gaps_sorted = sorted(chrom_gaps, key=lambda x: x['start'])
        for idx, gap in enumerate(chrom_gaps_sorted, 1):
            # 提取染色体编号（适配chr1/1/Chr02等格式）
            chrom_num = re.findall(r'\d+', chrom)[0] if re.findall(r'\d+', chrom) else chrom
            gap_name = f"gap{chrom_num}-{idx}"
            # 构建BED行
            line = f"{gap['chrom']}\t{gap['start']}\t{gap['end']}\t{gap_name}\t{gap['length']}"
            output_lines.append(line)
    
    # Write output
    if args.output:
        with open(args.output, 'w') as f:
            f.write('\n'.join(output_lines))
        print(f"\n# Output written to: {args.output}", file=sys.stderr)
    else:
        print('\n'.join(output_lines))
    
    # Statistics
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
