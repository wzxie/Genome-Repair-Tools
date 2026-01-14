#!/usr/bin/env python3

import sys
import re
import os
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
from collections import OrderedDict
import itertools

def detect_internal_telomeres(sequence, telomere_patterns, min_length=100, flank_size=5000):
    seq_len = len(sequence)
    internal_telomeres = []
    
    internal_start = flank_size
    internal_end = seq_len - flank_size
    
    for pattern_name, pattern_regex in telomere_patterns.items():
        for match in pattern_regex.finditer(sequence):
            start, end = match.start(), match.end()
            length = end - start
            
            if length >= min_length:
                if start >= internal_start and end <= internal_end:
                    internal_telomeres.append({
                        'start': start,
                        'end': end,
                        'length': length,
                        'pattern': pattern_name,
                        'type': 'internal'
                    })
    
    if internal_telomeres:
        internal_telomeres.sort(key=lambda x: x['start'])
        merged = []
        current = internal_telomeres[0]
        
        for tel in internal_telomeres[1:]:
            if tel['start'] <= current['end'] + 1000:
                current['end'] = max(current['end'], tel['end'])
                current['length'] = current['end'] - current['start']
            else:
                merged.append(current)
                current = tel
        
        merged.append(current)
        return merged
    
    return []

def compile_telomere_patterns_for_plant():
    patterns = {
        'Plant_TTTAGGG': re.compile(r'(TTTAGGG){10,}'),
        'Plant_CCCTAAA': re.compile(r'(CCCTAAA){10,}'),
        'TTAGGG': re.compile(r'(TTAGGG){10,}'),
        'CCCTAA': re.compile(r'(CCCTAA){10,}'),
    }
    return patterns

def analyze_single_chromosome_telomeres(args):
    record, min_repeats, search_window, internal_flank, internal_patterns = args
    
    seq_id = record.id
    seq_length = len(record.seq)
    sequence = str(record.seq).upper()
    
    plant_telomere_patterns = ['TTTAGGG', 'CCCTAAA']
    left_seq = sequence[:search_window]
    right_seq = sequence[-search_window:] if seq_length > search_window else sequence
    
    left_total_repeats = sum(left_seq.count(pattern) for pattern in plant_telomere_patterns)
    right_total_repeats = sum(right_seq.count(pattern) for pattern in plant_telomere_patterns)
    
    left_has_telomere = "yes" if left_total_repeats >= min_repeats else "no"
    right_has_telomere = "yes" if right_total_repeats >= min_repeats else "no"
    
    internal_telomeres = detect_internal_telomeres(
        sequence, 
        internal_patterns, 
        min_length=100,
        flank_size=internal_flank
    )
    
    return {
        'seq_id': seq_id,
        'left': left_has_telomere == "yes",
        'right': right_has_telomere == "yes",
        'left_repeats': left_total_repeats,
        'right_repeats': right_total_repeats,
        'internal_telomeres': internal_telomeres,
        'internal_count': len(internal_telomeres),
        'seq_length': seq_length
    }

def find_gaps_in_sequence_optimized(sequence, min_gap_size=1):
    gaps = []
    gap_start = None
    
    for i, base in enumerate(sequence):
        if base == 'N':
            if gap_start is None:
                gap_start = i
        else:
            if gap_start is not None:
                gap_length = i - gap_start
                if gap_length >= min_gap_size:
                    gaps.append((gap_start, i, gap_length))
                gap_start = None
    
    if gap_start is not None:
        gap_length = len(sequence) - gap_start
        if gap_length >= min_gap_size:
            gaps.append((gap_start, len(sequence), gap_length))
    
    return gaps

def analyze_fragments_optimized(sequence, gaps, min_fragment_length=None):
    sequence_length = len(sequence)
    fragments = []
    
    if not gaps:
        fragments.append((0, sequence_length, sequence_length, sequence))
    else:
        if gaps[0][0] > 0:
            frag_length = gaps[0][0]
            fragments.append((0, gaps[0][0], frag_length, sequence[0:gaps[0][0]]))
        
        for i in range(len(gaps) - 1):
            gap_end = gaps[i][1]
            next_gap_start = gaps[i + 1][0]
            
            if next_gap_start > gap_end:
                frag_length = next_gap_start - gap_end
                fragments.append((gap_end, next_gap_start, frag_length, sequence[gap_end:next_gap_start]))
        
        last_gap_end = gaps[-1][1]
        if last_gap_end < sequence_length:
            frag_length = sequence_length - last_gap_end
            fragments.append((last_gap_end, sequence_length, frag_length, sequence[last_gap_end:]))
    
    if not fragments:
        return {
            'total_fragments': 0,
            'fragments': [],
            'largest_fragment': 0,
            'smallest_fragment': 0,
            'small_fragments': []
        }
    
    fragment_lengths = [frag[2] for frag in fragments]
    small_fragments = [frag for frag in fragments if frag[2] < 1000000]
    
    result = {
        'total_fragments': len(fragments),
        'fragments': fragments,
        'largest_fragment': max(fragment_lengths) if fragment_lengths else 0,
        'smallest_fragment': min(fragment_lengths) if fragment_lengths else 0,
        'small_fragments': small_fragments
    }
    
    if min_fragment_length is not None:
        kept_fragments = [frag for frag in fragments if frag[2] >= min_fragment_length]
        discarded_fragments = [frag for frag in fragments if frag[2] < min_fragment_length]
        
        total_kept_length = sum(frag[2] for frag in kept_fragments)
        
        result.update({
            'kept_fragments': len(kept_fragments),
            'discarded_fragments': len(discarded_fragments),
            'total_kept_length': total_kept_length,
            'total_discarded_length': sum(frag[2] for frag in discarded_fragments),
            'kept_fragments_list': kept_fragments,
            'discarded_fragments_list': discarded_fragments,
            'original_sequence': sequence
        })
    
    return result

def create_spaced_diagram_to_file(f, seq_id, sequence_length, gaps, telomere_info, 
                                 internal_telomeres, width=50, spacing=2):
    if sequence_length == 0:
        f.write(" " * width + "\n")
        return " " * width
    
    scale = sequence_length / width if sequence_length > 0 else 1
    diagram = ['█'] * width
    
    for start, end, gap_length in gaps:
        start_idx = int(start / scale)
        end_idx = int(end / scale)
        
        start_idx = max(0, min(start_idx, width-1))
        end_idx = max(0, min(end_idx, width-1))
        
        for pos in range(start_idx, end_idx + 1):
            diagram[pos] = '○'
    
    if internal_telomeres:
        for tel in internal_telomeres:
            tel_mid = (tel['start'] + tel['end']) // 2
            tel_idx = int(tel_mid / scale)
            tel_idx = max(0, min(tel_idx, width-1))
            diagram[tel_idx] = '★'
    
    if telomere_info:
        if telomere_info.get('left', False):
            diagram[0] = '▶'
        if telomere_info.get('right', False):
            diagram[-1] = '◀'
    
    diagram_str = ''.join(diagram)
    
    info_parts = []
    info_parts.append(f"{sequence_length:,}bp")
    
    if gaps:
        total_gap_length = sum(gap[2] for gap in gaps)
        gap_percent = (total_gap_length / sequence_length * 100) if sequence_length > 0 else 0
        info_parts.append(f"gaps:{len(gaps)}({gap_percent:.1f}%)")
    
    if internal_telomeres:
        info_parts.append(f"intTel:{len(internal_telomeres)}")
    
    if telomere_info:
        left = "Y" if telomere_info.get('left') else "N"
        right = "Y" if telomere_info.get('right') else "N"
        left_repeats = telomere_info.get('left_repeats', 0)
        right_repeats = telomere_info.get('right_repeats', 0)
        info_parts.append(f"5'({left}[{left_repeats}]) 3'({right}[{right_repeats}])")
    
    info_str = " | ".join(info_parts)
    
    f.write(f"{'='*60}\n")
    f.write(f"Chromosome: {seq_id}\n")
    f.write(f"{'-'*60}\n")
    f.write(f"Stats: {info_str}\n")
    f.write(f"Diagram: 5' {diagram_str} 3'\n")
    f.write(f"{'='*60}\n\n")
    
    return diagram_str

def process_single_sequence(record, telomere_data_dict, min_gap_size, analyze_fragments, 
                           min_fragment_length, gap_output, fragment_output, internal_patterns):
    seq_id = record.id
    sequence_str = str(record.seq).upper()
    
    gaps = find_gaps_in_sequence_optimized(sequence_str, min_gap_size)
    
    if gaps:
        with open(gap_output, 'a') as f:
            for i, gap in enumerate(gaps, 1):
                start, end, length = gap
                f.write(f"{seq_id}\t{start}\t{end}\tgap_{i}\n")
    
    internal_telomeres = []
    if internal_patterns:
        internal_telomeres = detect_internal_telomeres(
            sequence_str, 
            internal_patterns, 
            min_length=100,
            flank_size=5000
        )
    else:
        seq_telomere_info = telomere_data_dict.get(seq_id, {}) if telomere_data_dict else {}
        internal_telomeres = seq_telomere_info.get('internal_telomeres', [])
    
    fragment_info = None
    if analyze_fragments:
        fragment_info = analyze_fragments_optimized(sequence_str, gaps, min_fragment_length)
        
        if fragment_info and 'fragments' in fragment_info:
            with open(fragment_output, 'a') as f:
                for start, end, length, seq in fragment_info['fragments']:
                    if min_fragment_length is None:
                        status = "small" if length < 1000000 else "large"
                    else:
                        status = "kept" if length >= min_fragment_length else "discarded"
                    f.write(f"{seq_id}\t{start}\t{end}\t{status}\t{length}\n")
    
    seq_telomere_info = telomere_data_dict.get(seq_id, {}) if telomere_data_dict else {}
    
    return_data = {
        'seq_id': seq_id,
        'sequence_length': len(sequence_str),
        'gaps': gaps,
        'telomere_info': {
            'left': seq_telomere_info.get('left', False),
            'right': seq_telomere_info.get('right', False),
            'left_repeats': seq_telomere_info.get('left_repeats', 0),
            'right_repeats': seq_telomere_info.get('right_repeats', 0)
        } if seq_telomere_info else {},
        'fragment_info': {
            'total_fragments': fragment_info.get('total_fragments', 0) if fragment_info else 0,
            'kept_fragments': fragment_info.get('kept_fragments', 0) if fragment_info else 0,
            'discarded_fragments': fragment_info.get('discarded_fragments', 0) if fragment_info else 0,
            'kept_fragments_list': fragment_info.get('kept_fragments_list', []) if fragment_info else [],
            'discarded_fragments_list': fragment_info.get('discarded_fragments_list', []) if fragment_info else [],
            'total_kept_length': fragment_info.get('total_kept_length', 0) if fragment_info else 0,
            'original_sequence': fragment_info.get('original_sequence', '') if fragment_info else ''
        } if fragment_info else None,
        'internal_telomeres': internal_telomeres,
        'gap_count': len(gaps),
        'internal_telomere_count': len(internal_telomeres),
        'sequence_str': sequence_str
    }
    
    return return_data

def process_single_sequence_wrapper(args):
    record, telomere_data_dict, min_gap_size, analyze_fragments, min_fragment_length, gap_output, fragment_output, internal_patterns = args
    return process_single_sequence(record, telomere_data_dict, min_gap_size, analyze_fragments, 
                                  min_fragment_length, gap_output, fragment_output, internal_patterns)

def analyze_chromosomes(
    input_file: str,
    detect_telomeres: bool = True,
    analyze_fragments: bool = True,
    save_diagrams: bool = True,
    min_gap_size: int = 1,
    min_fragment_length: int = None,
    threads: int = 16,
    min_repeats: int = 5,
    search_window: int = 2000,
    internal_flank: int = 5000,
    diagram_width: int = 50,
    output_prefix: str = "analysis",
    memory_optimized: bool = True,
    verbose: bool = False
) -> dict:
    
    if verbose:
        print("="*70)
        print("CHROMOSOME ANALYSIS TOOL - MODULE VERSION")
        print("="*70)
    
    total_start_time = time.time()
    
    gap_output = f"{output_prefix}_gap_coordinates.txt"
    fragment_output = f"{output_prefix}_fragment_info.txt"
    diagram_output = f"{output_prefix}_chromosome_diagrams.txt"
    filtered_fasta = f"{output_prefix}_filtered_sequences.fa"
    telomere_output = f"{output_prefix}_telomere_info.txt"
    internal_telomere_output = f"{output_prefix}_internal_telomeres.txt"
    
    open(gap_output, 'w').close()
    if analyze_fragments:
        open(fragment_output, 'w').close()
    
    if verbose:
        print(f"\n[Step 1] Loading sequences from {input_file}...")
    
    load_start = time.time()
    
    from Bio import SeqIO
    records = []
    with open(input_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            records.append(record)
    
    load_time = time.time() - load_start
    if verbose:
        print(f"[Load] Loaded {len(records):,} chromosomes in {load_time:.1f} seconds")
    
    if not records:
        raise ValueError(f"No sequences found in input file: {input_file}")
    
    original_order = [record.id for record in records]
    
    telomere_data = OrderedDict()
    all_internal_telomeres = []
    internal_patterns = None
    
    if detect_telomeres:
        telomere_start = time.time()
        
        internal_patterns = compile_telomere_patterns_for_plant()
        total_records = len(records)
        
        if verbose:
            print(f"[Telomere] Processing {total_records} chromosomes with {threads} threads...")
        
        telomere_args = [(record, min_repeats, search_window, internal_flank, internal_patterns) 
                         for record in records]
        
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(analyze_single_chromosome_telomeres, arg): i 
                       for i, arg in enumerate(telomere_args)}
            
            completed = 0
            for future in as_completed(futures):
                try:
                    result = future.result()
                    telomere_data[result['seq_id']] = {
                        'left': result['left'],
                        'right': result['right'],
                        'left_repeats': result['left_repeats'],
                        'right_repeats': result['right_repeats'],
                        'internal_telomeres': result['internal_telomeres'],
                        'internal_count': result['internal_count'],
                        'seq_length': result['seq_length']
                    }
                    
                    completed += 1
                    if verbose and completed % max(1, total_records // 10) == 0:
                        print(f"[Telomere] Progress: {completed}/{total_records} "
                              f"({completed/total_records*100:.0f}%)")
                        
                except Exception as e:
                    print(f"[Error] Telomere detection failed: {e}")
                    continue
        
        with open(telomere_output, 'w') as telo_file, open(internal_telomere_output, 'w') as internal_file:
            telo_file.write("#Chromosome\t5prime_telomere\t5prime_repeats\t3prime_telomere\t3prime_repeats\n")
            internal_file.write("#Chromosome\tStart\tEnd\tType\tPattern\tLength\n")
            
            for seq_id in original_order:
                if seq_id in telomere_data:
                    data = telomere_data[seq_id]
                    
                    left_str = "yes" if data['left'] else "no"
                    right_str = "yes" if data['right'] else "no"
                    
                    telo_file.write(f"{seq_id}\t{left_str}\t{data['left_repeats']}\t"
                                   f"{right_str}\t{data['right_repeats']}\n")
                    
                    for tel in data['internal_telomeres']:
                        internal_file.write(f"{seq_id}\t{tel['start']}\t{tel['end']}\t"
                                           f"internal\t{tel['pattern']}\t{tel['length']}\n")
                        all_internal_telomeres.append({
                            'chr': seq_id,
                            'start': tel['start'],
                            'end': tel['end'],
                            'pattern': tel['pattern'],
                            'length': tel['length']
                        })
        
        telomere_time = time.time() - telomere_start
        if verbose:
            print(f"[Telomere] Detection completed in {telomere_time:.1f} seconds")
    else:
        internal_patterns = compile_telomere_patterns_for_plant()
    
    if verbose:
        print(f"\n[Step 2] Analyzing gaps and fragments with {threads} threads...")
    
    analysis_start = time.time()
    
    analysis_args = [(record, telomere_data, min_gap_size, analyze_fragments, 
                      min_fragment_length, gap_output, fragment_output, internal_patterns) 
                     for record in records]
    
    analysis_results = []
    fragment_data_dict = {}
    
    results_by_seq_id = {}
    
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_single_sequence_wrapper, arg): i 
                   for i, arg in enumerate(analysis_args)}
        
        completed = 0
        for future in as_completed(futures):
            try:
                result = future.result()
                results_by_seq_id[result['seq_id']] = result
                
                if result.get('fragment_info'):
                    fragment_data_dict[result['seq_id']] = result['fragment_info']
                
                completed += 1
                if verbose and completed % max(1, len(records)//10) == 0:
                    progress = completed / len(records) * 100
                    print(f"[Progress] {completed:,}/{len(records):,} ({progress:.0f}%)")
                    
            except Exception as e:
                print(f"[Error] Analysis failed: {e}")
                continue
    
    for seq_id in original_order:
        if seq_id in results_by_seq_id:
            analysis_results.append(results_by_seq_id[seq_id])
    
    analysis_time = time.time() - analysis_start
    if verbose:
        print(f"[Analysis] Gap and fragment analysis completed in {analysis_time:.1f} seconds")
    
    if save_diagrams:
        if verbose:
            print(f"\n[Step 3] Saving chromosome diagrams to file...")
        
        diagram_start = time.time()
        
        with open(diagram_output, 'w') as f:
            f.write("="*70 + "\n")
            f.write("CHROMOSOME STRUCTURE DIAGRAMS\n")
            f.write("="*70 + "\n\n")
            f.write("Legend: █ = sequence, ○ = gap, ★ = internal telomere, ▶ = 5' telomere, ◀ = 3' telomere\n\n")
            
            for result in analysis_results:
                create_spaced_diagram_to_file(
                    f,
                    result['seq_id'],
                    result['sequence_length'],
                    result['gaps'],
                    result['telomere_info'],
                    result['internal_telomeres'],
                    width=diagram_width
                )
            
            f.write("="*70 + "\n")
            f.write("END OF DIAGRAMS\n")
            f.write("="*70 + "\n")
        
        diagram_time = time.time() - diagram_start
        if verbose:
            print(f"[Diagram] Diagrams saved in {diagram_time:.1f} seconds")
    
    if analyze_fragments and min_fragment_length is not None:
        if verbose:
            print(f"\n[Step 4] Creating filtered FASTA...")
        
        filter_start = time.time()
        
        sequences_written = 0
        
        with open(filtered_fasta, 'w') as f:
            for result in analysis_results:
                seq_id = result['seq_id']
                fragment_info = result.get('fragment_info')
                
                if fragment_info and 'kept_fragments_list' in fragment_info:
                    kept_fragments = fragment_info['kept_fragments_list']
                    if not kept_fragments:
                        continue
                    
                    new_sequence_parts = []
                    for i, (start, end, length, seq) in enumerate(kept_fragments):
                        if i > 0:
                            new_sequence_parts.append('N' * 100)
                        new_sequence_parts.append(seq)
                    
                    new_sequence = ''.join(new_sequence_parts)
                    
                    if new_sequence:
                        total_kept_fragments = len(kept_fragments)
                        total_kept_length = fragment_info.get('total_kept_length', 0)
                        original_length = result['sequence_length']
                        
                        header = (f">{seq_id}_filtered_"
                                 f"frags{total_kept_fragments}_"
                                 f"len{total_kept_length}_"
                                 f"orig{original_length}")
                        
                        f.write(f'{header}\n')
                        for i in range(0, len(new_sequence), 80):
                            f.write(new_sequence[i:i+80] + '\n')
                        
                        sequences_written += 1
        
        filter_time = time.time() - filter_start
        if verbose:
            print(f"[Filter] Filtered FASTA created in {filter_time:.1f} seconds")
            print(f"[Filter] Wrote {sequences_written} sequences to {filtered_fasta}")
    
    total_sequences = len(analysis_results)
    total_length = sum(r['sequence_length'] for r in analysis_results)
    total_gaps = sum(r['gap_count'] for r in analysis_results)
    total_internal_telomeres = sum(r['internal_telomere_count'] for r in analysis_results)
    
    telomere_stats = {}
    if telomere_data:
        sequences_with_terminal = sum(1 for data in telomere_data.values() if data['left'] or data['right'])
        sequences_with_both = sum(1 for data in telomere_data.values() if data['left'] and data['right'])
        sequences_with_internal = sum(1 for data in telomere_data.values() if data['internal_count'] > 0)
        total_internal_telo = sum(data['internal_count'] for data in telomere_data.values())
        
        telomere_stats = {
            'sequences_with_terminal': sequences_with_terminal,
            'sequences_with_both': sequences_with_both,
            'sequences_with_internal': sequences_with_internal,
            'total_internal_telomeres': total_internal_telo
        }
    
    fragment_stats = {}
    if analyze_fragments:
        fragment_info_list = [r for r in analysis_results if r['fragment_info']]
        if fragment_info_list:
            total_fragments = sum(info['fragment_info']['total_fragments'] for info in fragment_info_list)
            
            fragment_stats['total_fragments'] = total_fragments
            
            if min_fragment_length is not None:
                total_kept = sum(info['fragment_info'].get('kept_fragments', 0) for info in fragment_info_list)
                total_discarded = sum(info['fragment_info'].get('discarded_fragments', 0) for info in fragment_info_list)
                total_kept_bp = sum(info['fragment_info'].get('total_kept_length', 0) for info in fragment_info_list)
                
                fragment_stats.update({
                    'kept_fragments': total_kept,
                    'discarded_fragments': total_discarded,
                    'kept_bases': total_kept_bp,
                    'kept_percentage': (total_kept_bp / total_length * 100) if total_length > 0 else 0
                })
    
    total_time = time.time() - total_start_time
    processing_speed = total_length / total_time / 1e6 if total_time > 0 else 0
    
    results = {
        'summary': {
            'total_sequences': total_sequences,
            'total_length': total_length,
            'total_gaps': total_gaps,
            'total_internal_telomeres': total_internal_telomeres,
            'processing_time': total_time,
            'processing_speed_mbp_sec': processing_speed,
            'threads_used': threads
        },
        'telomere_stats': telomere_stats,
        'fragment_stats': fragment_stats,
        'output_files': {
            'gap_coordinates': gap_output,
            'fragment_info': fragment_output if analyze_fragments else None,
            'telomere_info': telomere_output if detect_telomeres else None,
            'internal_telomeres': internal_telomere_output if detect_telomeres else None,
            'chromosome_diagrams': diagram_output if save_diagrams else None,
            'filtered_fasta': filtered_fasta if (analyze_fragments and min_fragment_length is not None) else None
        },
        'analysis_results': analysis_results
    }
    
    if verbose:
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE")
        print("="*70)
        print(f"Total analysis time: {total_time:.1f} seconds")
    else:
        print("\nANALYSIS COMPLETE")
        print(f"Sequences: {total_sequences:,}")
        print(f"Total length: {total_length:,} bp")
        print(f"Total gaps: {total_gaps:,}")
        print(f"Internal telomeres: {total_internal_telomeres:,}")
        print(f"Processing speed: {processing_speed:.1f} Mbp/sec")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Chromosome Analysis Tool - Can be used as module or standalone',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standalone usage with quiet mode
  python chromosome_analyzer.py -q genome.fa --detect-telomeres --analyze-fragments --save-diagrams
  
  # Basic usage with minimal output
  python chromosome_analyzer.py -q genome.fa --detect-telomeres --analyze-fragments
  
  # With all options
  python chromosome_analyzer.py -q genome.fa --detect-telomeres --analyze-fragments --min-fragment-length 1000000 --threads 32
  
  # Module usage example (in another Python script):
  from chromosome_analyzer import analyze_chromosomes
  
  results = analyze_chromosomes(
      input_file="genome.fa",
      detect_telomeres=True,
      analyze_fragments=True,
      threads=16,
      output_prefix="my_analysis",
      verbose=False
  )
  
  print(f"Analyzed {results['summary']['total_sequences']} chromosomes")
  print(f"Results saved to: {results['output_files']['gap_coordinates']}")
        """
    )
    
    parser.add_argument('-q', '--query', required=True,
                       help='Input FASTA file (required)')
    parser.add_argument('-w', '--width', type=int, default=50,
                       help='Diagram width (default: 50)')
    parser.add_argument('-m', '--min-gap-size', type=int, default=1,
                       help='Minimum gap size (default: 1)')
    parser.add_argument('-o', '--output-prefix', default='analysis',
                       help='Output file prefix (default: analysis)')
    parser.add_argument('-t', '--threads', type=int, default=16,
                       help='Number of threads (default: 16)')
    parser.add_argument('--save-diagrams', action='store_true',
                       help='Save chromosome diagrams to file')
    parser.add_argument('--detect-telomeres', action='store_true',
                       help='Multi-threaded telomere detection')
    parser.add_argument('--analyze-fragments', action='store_true',
                       help='Analyze fragments between gaps')
    parser.add_argument('--min-fragment-length', type=int,
                       help='Minimum fragment length to keep (bp)')
    parser.add_argument('--min-repeats', type=int, default=5,
                       help='Minimum telomere repeats (default: 5)')
    parser.add_argument('--search-window', type=int, default=2000,
                       help='Telomere search window (bp, default: 2000)')
    parser.add_argument('--internal-flank', type=int, default=5000,
                       help='Internal telomere flank size (bp, default: 5000)')
    parser.add_argument('--no-diagrams', action='store_true',
                       help='Do not save diagrams (overrides --save-diagrams)')
    parser.add_argument('--quiet', action='store_true',
                       help='Reduce screen output (only show errors and final summary)')
    parser.add_argument('--verbose', action='store_true',
                       help='Show detailed progress information')
    parser.add_argument('--version', action='version', version='Chromosome Analyzer v1.0')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args()
    
    try:
        results = analyze_chromosomes(
            input_file=args.query,
            detect_telomeres=args.detect_telomeres,
            analyze_fragments=args.analyze_fragments,
            save_diagrams=args.save_diagrams and not args.no_diagrams,
            min_gap_size=args.min_gap_size,
            min_fragment_length=args.min_fragment_length,
            threads=args.threads,
            min_repeats=args.min_repeats,
            search_window=args.search_window,
            internal_flank=args.internal_flank,
            diagram_width=args.width,
            output_prefix=args.output_prefix,
            memory_optimized=True,
            verbose=args.verbose or not args.quiet
        )
        
        if not args.verbose and args.quiet:
            summary = results['summary']
            telomere_stats = results.get('telomere_stats', {})
            
            print(f"\nSUMMARY:")
            print(f"  Chromosomes: {summary['total_sequences']:,}")
            print(f"  Total length: {summary['total_length']:,} bp")
            print(f"  Total gaps: {summary['total_gaps']:,}")
            print(f"  Internal telomeres: {summary['total_internal_telomeres']:,}")
            
            if args.detect_telomeres and telomere_stats:
                print(f"  With terminal telomeres: {telomere_stats.get('sequences_with_terminal', 0):,}")
                print(f"  With internal telomeres: {telomere_stats.get('sequences_with_internal', 0):,}")
            
            if args.min_fragment_length:
                fragment_stats = results.get('fragment_stats', {})
                if 'kept_fragments' in fragment_stats:
                    print(f"  Kept fragments: {fragment_stats['kept_fragments']:,}")
                    print(f"  Discarded fragments: {fragment_stats['discarded_fragments']:,}")
                    print(f"  Bases kept: {fragment_stats['kept_bases']:,} ({fragment_stats['kept_percentage']:.1f}%)")
        
        print(f"\nOUTPUT FILES:")
        for name, path in results['output_files'].items():
            if path and os.path.exists(path):
                print(f"  ✓ {name}: {path}")
        
        if args.detect_telomeres:
            print("\nNote: Internal telomeres are marked with ★ in diagrams")
            print("      See internal_telomeres.txt for coordinates")
        
        print("\n" + "="*70)
        print("All tasks completed! ✓")
        
    except Exception as e:
        print(f"\n[Error] Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

def run_as_subscript():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-q', '--query', required=True)
    parser.add_argument('-w', '--width', type=int, default=50)
    parser.add_argument('-m', '--min-gap-size', type=int, default=1)
    parser.add_argument('-o', '--output-prefix', default='analysis')
    parser.add_argument('-t', '--threads', type=int, default=16)
    parser.add_argument('--save-diagrams', action='store_true')
    parser.add_argument('--detect-telomeres', action='store_true')
    parser.add_argument('--analyze-fragments', action='store_true')
    parser.add_argument('--min-fragment-length', type=int)
    parser.add_argument('--min-repeats', type=int, default=5)
    parser.add_argument('--search-window', type=int, default=2000)
    parser.add_argument('--internal-flank', type=int, default=5000)
    parser.add_argument('--no-diagrams', action='store_true')
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    
    args = parser.parse_args()
    
    return analyze_chromosomes(
        input_file=args.query,
        detect_telomeres=args.detect_telomeres,
        analyze_fragments=args.analyze_fragments,
        save_diagrams=args.save_diagrams and not args.no_diagrams,
        min_gap_size=args.min_gap_size,
        min_fragment_length=args.min_fragment_length,
        threads=args.threads,
        min_repeats=args.min_repeats,
        search_window=args.search_window,
        internal_flank=args.internal_flank,
        diagram_width=args.width,
        output_prefix=args.output_prefix,
        memory_optimized=True,
        verbose=args.verbose or not args.quiet
    )