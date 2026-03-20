#!/usr/bin/env python3
"""
Chromosome Analysis Tool - Simplified Command Line Version
Automatically detects telomere patterns from sequence data
"""

import sys
import re
import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
from collections import OrderedDict, Counter
from Bio import SeqIO

# ============================================================================
# Telomere pattern database
# ============================================================================

TELOMERE_PATTERNS = {
    'plant': {
        'name': 'Plant (Arabidopsis-type)',
        'forward': 'TTTAGGG',
        'reverse': 'CCCTAAA',
        'description': 'Most plants',
        'weight': 1.0
    },
    'human': {
        'name': 'Human/Vertebrate',
        'forward': 'TTAGGG',
        'reverse': 'CCCTAA',
        'description': 'Humans, vertebrates, most mammals',
        'weight': 1.0
    },
    'insect': {
        'name': 'Insect',
        'forward': 'TTAGG',
        'reverse': 'CCTAA',
        'description': 'Most insects',
        'weight': 0.9
    },
    'nematode': {
        'name': 'Nematode',
        'forward': 'TTAGGC',
        'reverse': 'GCCTAA',
        'description': 'C. elegans and other nematodes',
        'weight': 0.9
    },
    'tetrahymena': {
        'name': 'Tetrahymena',
        'forward': 'TTGGGG',
        'reverse': 'CCCCAA',
        'description': 'Some ciliates',
        'weight': 0.8
    },
    'oxytricha': {
        'name': 'Oxytricha',
        'forward': 'TTTTGGGG',
        'reverse': 'CCCCAAAA',
        'description': 'Some ciliates',
        'weight': 0.8
    },
    'yeast': {
        'name': 'Yeast',
        'forward': 'TGTGGGTGTGGTG',
        'reverse': None,
        'description': 'S. cerevisiae (variable length)',
        'weight': 0.7
    },
    'trypanosome': {
        'name': 'Trypanosome',
        'forward': 'TTAGGG',
        'reverse': 'CCCTAA',
        'description': 'Trypanosoma brucei (same as human)',
        'weight': 1.0
    },
    'aspergillus': {
        'name': 'Aspergillus',
        'forward': 'TTAGGGTCAACA',
        'reverse': None,
        'description': 'Some fungi',
        'weight': 0.7
    }
}

def reverse_complement(seq):
    """Reverse complement of DNA sequence"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp.get(b, b) for b in seq[::-1])

def detect_telomere_species(sequences, sample_size=5, window=5000):
    """
    Automatically detect telomere type from sequences
    
    Samples from sequence ends and counts occurrences of various telomere patterns
    """
    print("\n[Detection] Automatically detecting telomere type...")
    
    # Prepare sequence samples
    samples = []
    for i, record in enumerate(sequences):
        if i >= sample_size:
            break
        seq = str(record.seq).upper()
        if len(seq) > window * 2:
            samples.append(seq[:window])  # 5' end
            samples.append(seq[-window:]) # 3' end
        else:
            samples.append(seq)
    
    if not samples:
        print("  Warning: Cannot sample sequences, using default plant telomere type")
        return 'plant'
    
    # Count occurrences of each telomere pattern
    pattern_counts = Counter()
    
    for sample in samples:
        for species, info in TELOMERE_PATTERNS.items():
            forward = info['forward']
            if forward:
                # Count forward repeats (at least 3 repeats)
                f_pattern = forward * 3
                if f_pattern in sample:
                    count = sample.count(forward)
                    pattern_counts[species] += count * info['weight']
            
            reverse = info.get('reverse')
            if reverse:
                r_pattern = reverse * 3
                if r_pattern in sample:
                    count = sample.count(reverse)
                    pattern_counts[species] += count * info['weight']
    
    if not pattern_counts:
        print("  No known telomere patterns detected, using default plant telomere type")
        return 'plant'
    
    # Find the highest scoring species
    top_species = pattern_counts.most_common(1)[0][0]
    
    # Show detection results
    print(f"  Detected telomere type: {TELOMERE_PATTERNS[top_species]['name']}")
    print(f"  Pattern: {TELOMERE_PATTERNS[top_species]['forward']} / "
          f"{TELOMERE_PATTERNS[top_species].get('reverse', 'N/A')}")
    
    # Show other possible matches
    if len(pattern_counts) > 1:
        others = pattern_counts.most_common(3)[1:]
        if others:
            other_strs = []
            for s in others:
                other_strs.append(f"{TELOMERE_PATTERNS[s[0]]['name']}({s[1]:.0f})")
            print(f"  Other possibilities: {', '.join(other_strs)}")
    
    return top_species

def compile_patterns(species):
    """Compile regex patterns for specified species"""
    info = TELOMERE_PATTERNS[species]
    patterns = {}
    
    if info['forward']:
        patterns[f'{species}_forward'] = re.compile(r'(' + re.escape(info['forward']) + r'){4,}')
    
    if info.get('reverse'):
        patterns[f'{species}_reverse'] = re.compile(r'(' + re.escape(info['reverse']) + r'){4,}')
    
    # Also add human patterns as fallback (many species share these)
    if species != 'human':
        patterns['human_fallback'] = re.compile(r'(TTAGGG){4,}')
        patterns['human_fallback_rev'] = re.compile(r'(CCCTAA){4,}')
    
    return patterns

# ============================================================================
# Core analysis functions
# ============================================================================

def detect_internal_telomeres(sequence, patterns, min_length=50, flank_size=5000):
    """Detect internal telomeres (not at chromosome ends)"""
    seq_len = len(sequence)
    internal_telomeres = []
    
    internal_start = flank_size
    internal_end = seq_len - flank_size
    
    for pattern_name, pattern_regex in patterns.items():
        for match in pattern_regex.finditer(sequence):
            start, end = match.start(), match.end()
            length = end - start
            
            if length >= min_length:
                if start >= internal_start and end <= internal_end:
                    internal_telomeres.append({
                        'start': start,
                        'end': end,
                        'length': length,
                        'pattern': pattern_name
                    })
    
    # Merge adjacent telomeres
    if internal_telomeres:
        internal_telomeres.sort(key=lambda x: x['start'])
        merged = []
        current = internal_telomeres[0]
        
        for tel in internal_telomeres[1:]:
            if tel['start'] <= current['end'] + 500:
                current['end'] = max(current['end'], tel['end'])
                current['length'] = current['end'] - current['start']
            else:
                merged.append(current)
                current = tel
        merged.append(current)
        return merged
    
    return []

def analyze_chromosome(args):
    """Analyze telomeres for a single chromosome"""
    record, species, min_repeats, search_window, internal_flank = args
    
    seq_id = record.id
    seq_length = len(record.seq)
    sequence = str(record.seq).upper()
    
    # Compile patterns for this species
    patterns = compile_patterns(species)
    
    # Terminal telomere detection
    left_seq = sequence[:search_window]
    right_seq = sequence[-search_window:] if seq_length > search_window else sequence
    
    # Get telomere sequences for this species
    forward = TELOMERE_PATTERNS[species]['forward']
    reverse = TELOMERE_PATTERNS[species].get('reverse')
    
    left_repeats = left_seq.count(forward)
    right_repeats = right_seq.count(forward)
    
    if reverse:
        left_repeats += left_seq.count(reverse)
        right_repeats += right_seq.count(reverse)
    
    # Internal telomere detection
    internal_telomeres = detect_internal_telomeres(
        sequence, patterns, min_length=50, flank_size=internal_flank
    )
    
    return {
        'seq_id': seq_id,
        'length': seq_length,
        'left': left_repeats >= min_repeats,
        'right': right_repeats >= min_repeats,
        'left_repeats': left_repeats,
        'right_repeats': right_repeats,
        'internal': internal_telomeres,
        'internal_count': len(internal_telomeres)
    }

def find_gaps(sequence, min_gap=1):
    """Find gaps (N's) in sequence"""
    gaps = []
    start = None
    
    for i, base in enumerate(sequence):
        if base == 'N':
            if start is None:
                start = i
        else:
            if start is not None:
                length = i - start
                if length >= min_gap:
                    gaps.append((start, i, length))
                start = None
    
    if start is not None:
        length = len(sequence) - start
        if length >= min_gap:
            gaps.append((start, len(sequence), length))
    
    return gaps

def analyze_fragments(sequence, gaps, min_fragment=None):
    """Analyze fragments between gaps"""
    fragments = []
    seq_len = len(sequence)
    
    if not gaps:
        fragments.append((0, seq_len, seq_len))
    else:
        if gaps[0][0] > 0:
            fragments.append((0, gaps[0][0], gaps[0][0]))
        
        for i in range(len(gaps) - 1):
            gap_end = gaps[i][1]
            next_start = gaps[i + 1][0]
            if next_start > gap_end:
                fragments.append((gap_end, next_start, next_start - gap_end))
        
        last_end = gaps[-1][1]
        if last_end < seq_len:
            fragments.append((last_end, seq_len, seq_len - last_end))
    
    if min_fragment:
        kept = [f for f in fragments if f[2] >= min_fragment]
        discarded = [f for f in fragments if f[2] < min_fragment]
    else:
        kept = fragments
        discarded = []
    
    return {
        'total': len(fragments),
        'kept': kept,
        'discarded': discarded,
        'kept_count': len(kept),
        'discarded_count': len(discarded)
    }

def draw_diagram(seq_id, length, gaps, tel_info, internal, width=60):
    """Draw chromosome diagram"""
    if length == 0:
        return ""
    
    scale = length / width
    diagram = ['█'] * width
    
    # Draw gaps
    for start, end, _ in gaps:
        s = max(0, min(int(start / scale), width - 1))
        e = max(0, min(int(end / scale), width - 1))
        for i in range(s, e + 1):
            diagram[i] = '○'
    
    # Draw internal telomeres
    for tel in internal:
        pos = int(((tel['start'] + tel['end']) // 2) / scale)
        pos = max(0, min(pos, width - 1))
        diagram[pos] = '★'
    
    # Draw terminal telomeres
    if tel_info.get('left'):
        diagram[0] = '▶'
    if tel_info.get('right'):
        diagram[-1] = '◀'
    
    # Statistics
    gap_len = sum(g[2] for g in gaps)
    gap_pct = (gap_len / length * 100) if length > 0 else 0
    
    info = f"{seq_id} | {length:,}bp | gaps:{len(gaps)}({gap_pct:.1f}%)"
    if internal:
        info += f" | int-tel:{len(internal)}"
    info += f" | 5':{'Y' if tel_info.get('left') else 'N'}({tel_info.get('left_repeats', 0)})"
    info += f" 3':{'Y' if tel_info.get('right') else 'N'}({tel_info.get('right_repeats', 0)})"
    
    return f"{info}\n5' {''.join(diagram)} 3'"

# ============================================================================
# Main function
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Chromosome Analysis Tool - Automatically detects telomere types',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s genome.fasta                    # Auto-detect telomere type and analyze
  %(prog)s genome.fasta --min-fragment 1m  # Keep only fragments >1Mb
  %(prog)s genome.fasta --species human    # Manually specify species
  %(prog)s genome.fasta --no-diagram       # Don't show diagrams
        """
    )
    
    parser.add_argument('query', help='Input FASTA file')
    parser.add_argument('--species', choices=list(TELOMERE_PATTERNS.keys()),
                       help='Manually specify species (default: auto-detect)')
    parser.add_argument('--min-repeats', type=int, default=5,
                       help='Minimum telomere repeats (default: 5)')
    parser.add_argument('--window', type=int, default=5000,
                       help='Telomere search window size (default: 5000bp)')
    parser.add_argument('--min-gap', type=int, default=1,
                       help='Minimum gap size (default: 1bp)')
    parser.add_argument('--min-fragment', 
                       help='Minimum fragment length (e.g., 100k, 1m, 500k)')
    parser.add_argument('--threads', type=int, default=8,
                       help='Number of threads (default: 8)')
    parser.add_argument('--no-diagram', action='store_true',
                       help='Do not show chromosome diagrams')
    parser.add_argument('--output', '-o', default='analysis',
                       help='Output file prefix (default: analysis)')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Parse minimum fragment length
    min_frag = None
    if args.min_fragment:
        min_frag = args.min_fragment.lower()
        if min_frag.endswith('k'):
            min_frag = int(float(min_frag[:-1]) * 1000)
        elif min_frag.endswith('m'):
            min_frag = int(float(min_frag[:-1]) * 1000000)
        else:
            min_frag = int(min_frag)
        print(f"[Parameter] Minimum fragment length: {min_frag:,} bp")
    
    # Read sequences
    print(f"[Input] Reading file: {args.query}")
    try:
        records = list(SeqIO.parse(args.query, "fasta"))
    except Exception as e:
        print(f"[Error] Failed to read file: {e}")
        return 1
        
    print(f"[Input] Loaded {len(records):,} sequences")
    
    if not records:
        print("[Error] No sequences found")
        return 1
    
    # Determine species
    if args.species:
        species = args.species
        print(f"[Species] Using specified species: {TELOMERE_PATTERNS[species]['name']}")
    else:
        species = detect_telomere_species(records)
    
    # Parallel chromosome analysis
    print(f"\n[Analysis] Analyzing {len(records)} chromosomes with {args.threads} threads...")
    
    analyze_args = [(r, species, args.min_repeats, args.window, args.window) 
                    for r in records]
    
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(analyze_chromosome, arg) for arg in analyze_args]
        
        for i, future in enumerate(as_completed(futures), 1):
            if i % max(1, len(records)//10) == 0:
                print(f"  Progress: {i}/{len(records)} ({i/len(records)*100:.0f}%)")
            try:
                results.append(future.result())
            except Exception as e:
                print(f"  Warning: Analysis failed: {e}")
    
    # Sort by original order
    id_to_result = {r['seq_id']: r for r in results}
    results = [id_to_result[r.id] for r in records if r.id in id_to_result]
    
    # Open output files
    try:
        gap_file = open(f"{args.output}_gaps.txt", 'w')
        telo_file = open(f"{args.output}_telomeres.txt", 'w')
        internal_file = open(f"{args.output}_internal.txt", 'w')
    except Exception as e:
        print(f"[Error] Cannot create output files: {e}")
        return 1
    
    gap_file.write("#Chromosome\tStart\tEnd\tLength\n")
    telo_file.write("#Chromosome\t5'_telomere\t5'_repeats\t3'_telomere\t3'_repeats\n")
    internal_file.write("#Chromosome\tStart\tEnd\tLength\n")
    
    # Statistics
    total_len = 0
    total_gaps = 0
    total_internal = 0
    telo_left = 0
    telo_right = 0
    telo_both = 0
    
    print(f"\n{'='*80}")
    print(f"Results - Species: {TELOMERE_PATTERNS[species]['name']}")
    print(f"{'='*80}")
    
    for res in results:
        seq_id = res['seq_id']
        # Get original sequence
        record = [r for r in records if r.id == seq_id][0]
        sequence = str(record.seq).upper()
        
        # Find gaps
        gaps = find_gaps(sequence, args.min_gap)
        
        # Analyze fragments
        fragments = analyze_fragments(sequence, gaps, min_frag)
        
        # Write to files
        for g in gaps:
            gap_file.write(f"{seq_id}\t{g[0]}\t{g[1]}\t{g[2]}\n")
        
        telo_file.write(f"{seq_id}\t{'yes' if res['left'] else 'no'}\t{res['left_repeats']}\t"
                       f"{'yes' if res['right'] else 'no'}\t{res['right_repeats']}\n")
        
        for tel in res['internal']:
            internal_file.write(f"{seq_id}\t{tel['start']}\t{tel['end']}\t{tel['length']}\n")
        
        # Update statistics
        total_len += res['length']
        total_gaps += len(gaps)
        total_internal += res['internal_count']
        if res['left']:
            telo_left += 1
        if res['right']:
            telo_right += 1
        if res['left'] and res['right']:
            telo_both += 1
        
        # Show diagram
        if not args.no_diagram:
            print(draw_diagram(seq_id, res['length'], gaps, res, res['internal']))
            
            if fragments['discarded_count'] > 0:
                print(f"  Discarded small fragments: {fragments['discarded_count']} "
                      f"(total length: {sum(f[2] for f in fragments['discarded']):,} bp)")
            print()
    
    # Close files
    gap_file.close()
    telo_file.close()
    internal_file.close()
    
    # Generate filtered FASTA
    if min_frag:
        try:
            filtered_file = open(f"{args.output}_filtered.fa", 'w')
            filtered_count = 0
            
            for res in results:
                seq_id = res['seq_id']
                record = [r for r in records if r.id == seq_id][0]
                sequence = str(record.seq).upper()
                gaps = find_gaps(sequence, args.min_gap)
                fragments = analyze_fragments(sequence, gaps, min_frag)
                
                if fragments['kept']:
                    kept_parts = []
                    for start, end, _ in fragments['kept']:
                        kept_parts.append(sequence[start:end])
                    
                    if kept_parts:
                        new_seq = 'N'.join(kept_parts)
                        filtered_file.write(f">{seq_id}_filtered\n")
                        for i in range(0, len(new_seq), 80):
                            filtered_file.write(new_seq[i:i+80] + '\n')
                        filtered_count += 1
            
            filtered_file.close()
            print(f"[Output] Saved filtered sequences: {filtered_count} chromosomes")
        except Exception as e:
            print(f"[Error] Cannot create filtered file: {e}")
    
    # Final statistics
    elapsed = time.time() - start_time
    
    print(f"\n{'='*80}")
    print("Summary Statistics")
    print(f"{'='*80}")
    print(f"Processing time: {elapsed:.1f} seconds")
    print(f"Total chromosomes: {len(results):,}")
    print(f"Total length: {total_len:,} bp ({total_len/1e6:.2f} Mb)")
    print(f"Total gaps: {total_gaps:,}")
    print(f"Internal telomeres: {total_internal:,}")
    print(f"5' telomeres: {telo_left}/{len(results)} ({telo_left/len(results)*100:.1f}%)")
    print(f"3' telomeres: {telo_right}/{len(results)} ({telo_right/len(results)*100:.1f}%)")
    print(f"Both ends: {telo_both}/{len(results)} ({telo_both/len(results)*100:.1f}%)")
    
    print(f"\nOutput files:")
    print(f"  - {args.output}_gaps.txt (gap coordinates)")
    print(f"  - {args.output}_telomeres.txt (telomere information)")
    print(f"  - {args.output}_internal.txt (internal telomeres)")
    if min_frag:
        print(f"  - {args.output}_filtered.fa (filtered sequences)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())