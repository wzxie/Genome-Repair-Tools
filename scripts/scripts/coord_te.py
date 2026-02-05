#!/usr/bin/env python3
"""
Telomere Alignment Tool - Complex alignment, generate coords files
Part 1: Extract telomere contigs, run nucmer alignment
"""

import sys
import os
import re
import subprocess
import concurrent.futures
import argparse
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import time
from typing import Dict, List, Tuple, Optional, Any, Set

VERBOSE = False

def printv(message: str):
    if VERBOSE:
        print(message)

def set_environment_limits():
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'

def find_telomere_sequences(sequence, min_repeats=20, min_telomere_length=500):
    telomere_patterns = [
        (re.compile(r'(TTTAGGG){%d,}' % min_repeats), "Plant telomere forward"),
        (re.compile(r'(CCCTAAA){%d,}' % min_repeats), "Plant telomere reverse"),
        (re.compile(r'(TTAGGG){%d,}' % min_repeats), "Universal telomere forward"),
        (re.compile(r'(CCCTAA){%d,}' % min_repeats), "Universal telomere reverse"),
    ]
    
    matches = []
    sequence_upper = sequence.upper()
    
    for pattern, pattern_name in telomere_patterns:
        for match in pattern.finditer(sequence_upper):
            if len(match.group()) >= min_telomere_length:
                matches.append({
                    'pattern': pattern_name,
                    'sequence': match.group(),
                    'start': match.start(),
                    'end': match.end(),
                    'length': len(match.group())
                })
    
    return matches

def extract_telomere_contigs(input_file, output_file, min_repeats=20, min_telomere_length=500):
    printv(f"Extracting telomere-containing contigs (repeats {min_repeats} times, length {min_telomere_length}bp)...")
    
    telomere_contigs = []
    total_contigs = 0
    
    try:
        for record in SeqIO.parse(input_file, "fasta"):
            total_contigs += 1
            sequence = str(record.seq)
            
            matches = find_telomere_sequences(sequence, min_repeats, min_telomere_length)
            
            if matches:
                telomere_info = []
                for i, match in enumerate(matches, 1):
                    telomere_info.append(
                        f"telomere_{i}:pos={match['start']}-{match['end']}:len={match['length']}"
                    )
                
                if telomere_info:
                    record.description = f"{record.description} [Telomere: {', '.join(telomere_info)}]"
                    telomere_contigs.append(record)
    
    except Exception as e:
        print(f"Error: {e}")
        raise
    
    if telomere_contigs:
        SeqIO.write(telomere_contigs, output_file, "fasta")
        printv(f"Extracted {len(telomere_contigs)} telomere-containing contigs")
    else:
        printv("No telomere-containing contigs found")
        return 0, total_contigs
    
    return len(telomere_contigs), total_contigs

def analyze_chromosome_telomeres_simple(fasta_file: str, min_repeats: int = 5, search_window: int = 2000) -> Dict:
    printv(f"Analyzing chromosome telomeres (only detecting ends)...")
    
    chromosome_status = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq_len = len(record.seq)
        sequence = str(record.seq).upper()
        
        left_seq = sequence[:search_window]
        right_seq = sequence[-search_window:] if seq_len > search_window else sequence
        
        plant_telomere_patterns = ['TTTAGGG', 'CCCTAAA']
        left_total_repeats = sum(left_seq.count(pattern) for pattern in plant_telomere_patterns)
        right_total_repeats = sum(right_seq.count(pattern) for pattern in plant_telomere_patterns)
        
        has_5prime = left_total_repeats >= min_repeats
        has_3prime = right_total_repeats >= min_repeats
        
        status = "complete"
        if not has_5prime and not has_3prime:
            status = "missing_both"
        elif not has_5prime:
            status = "missing_5prime"
        elif not has_3prime:
            status = "missing_3prime"
        
        chromosome_status[seq_id] = {
            'length': seq_len,
            'has_5prime': has_5prime,
            'has_3prime': has_3prime,
            'status': status,
            'needs_repair': status != "complete",
            'left_repeats': left_total_repeats,
            'right_repeats': right_total_repeats
        }
        
        status_symbol = "✓" if status == "complete" else "⚠️"
        printv(f"  {status_symbol} {seq_id}: {seq_len:,} bp - {status} (5':{left_total_repeats} repeats, 3':{right_total_repeats} repeats)")
    
    needs_repair = sum(1 for info in chromosome_status.values() if info['needs_repair'])
    printv(f"Chromosomes needing repair: {needs_repair}")
    
    return chromosome_status

def prepare_chromosome_ends(fasta_file: str, chromosome_status: Dict, output_dir: str) -> Dict:
    printv(f"Preparing chromosome end sequences...")
    
    query_data = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        
        if seq_id not in chromosome_status or not chromosome_status[seq_id]['needs_repair']:
            continue
        
        seq_len = len(record.seq)
        status = chromosome_status[seq_id]['status']
        
        if status in ['missing_5prime', 'missing_both']:
            end_seq = record.seq[:5000000]
            end_record = SeqRecord(end_seq, id=f"{seq_id}_5prime", description="")
            end_file = os.path.join(output_dir, f"{seq_id}_5prime.fa")
            SeqIO.write([end_record], end_file, "fasta")
            
            if seq_id not in query_data:
                query_data[seq_id] = {}
            query_data[seq_id]['5prime_file'] = end_file
        
        if status in ['missing_3prime', 'missing_both']:
            end_seq = record.seq[-5000000:] if seq_len > 5000000 else record.seq
            end_record = SeqRecord(end_seq, id=f"{seq_id}_3prime", description="")
            end_file = os.path.join(output_dir, f"{seq_id}_3prime.fa")
            SeqIO.write([end_record], end_file, "fasta")
            
            if seq_id not in query_data:
                query_data[seq_id] = {}
            query_data[seq_id]['3prime_file'] = end_file
        
        query_data[seq_id]['sequence'] = record
    
    printv(f"Prepared end sequences for {len(query_data)} chromosomes")
    return query_data

def run_nucmer_with_filter(ref_fa: str, query_fa: str, output_prefix: str, threads: int = 32) -> Optional[str]:
    try:
        env = os.environ.copy()
        env['OPENBLAS_NUM_THREADS'] = '1'
        
        nucmer_cmd = [
            "nucmer",
            "-c", "1000",
            "-l", "100",
            "--batch=500000000",
            "-t", str(threads),
            "-p", output_prefix,
            ref_fa, query_fa
        ]
        
        subprocess.run(nucmer_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        
        delta_file = f"{output_prefix}.delta"
        if not os.path.exists(delta_file):
            return None
        
        filter_cmd = ["delta-filter", "-i", "-r", "-l", "10000", delta_file]
        filtered_delta = f"{output_prefix}.filtered.delta"
        
        with open(filtered_delta, 'w') as outfile:
            subprocess.run(filter_cmd, stdout=outfile, stderr=subprocess.DEVNULL, env=env)
        
        coords_cmd = ["show-coords", "-r", "-l", filtered_delta]
        coords_file = f"{output_prefix}.coords"
        
        with open(coords_file, 'w') as outfile:
            subprocess.run(coords_cmd, stdout=outfile, stderr=subprocess.DEVNULL, env=env)
        
        has_alignments = False
        if os.path.exists(coords_file):
            with open(coords_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('[') and not line.startswith('=') and not line.startswith('NUCMER'):
                        has_alignments = True
                        break
        
        for f in [delta_file, filtered_delta]:
            if os.path.exists(f):
                os.remove(f)
        
        if not has_alignments:
            if os.path.exists(coords_file):
                os.remove(coords_file)
            return None
        
        return coords_file
        
    except subprocess.CalledProcessError as e:
        for ext in ['.delta', '.filtered.delta', '.coords']:
            f = f"{output_prefix}{ext}"
            if os.path.exists(f):
                os.remove(f)
        printv(f"nucmer alignment failed: {e}")
        return None
    except Exception as e:
        printv(f"Error during alignment: {e}")
        return None

def process_single_alignment(contig_file: str, query_file: str, output_prefix: str, threads: int = 1) -> Optional[str]:
    return run_nucmer_with_filter(contig_file, query_file, output_prefix, threads)

def process_chromosome_telomeres_directional(
    seq_id: str, 
    query_data: Dict, 
    telomere_contigs: List[str], 
    output_dir: str, 
    threads: int
) -> Dict:
    results = {
        '5prime': {'coords_files': [], 'matched_contigs': set()},
        '3prime': {'coords_files': [], 'matched_contigs': set()}
    }
    
    if '5prime_file' in query_data[seq_id]:
        query_file = query_data[seq_id]['5prime_file']
        
        printv(f"  Processing 5' end alignment...")
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(threads, len(telomere_contigs))) as executor:
            future_to_contig = {}
            for i, contig_file in enumerate(telomere_contigs, 1):
                contig_name = os.path.basename(contig_file).replace('.fa', '')
                prefix = os.path.join(output_dir, f"{seq_id}_5prime_{contig_name}")
                future = executor.submit(process_single_alignment, contig_file, query_file, prefix, 1)
                future_to_contig[future] = (contig_file, contig_name)
            
            for future in concurrent.futures.as_completed(future_to_contig):
                contig_file, contig_name = future_to_contig[future]
                try:
                    coords_file = future.result()
                    if coords_file and os.path.exists(coords_file):
                        results['5prime']['coords_files'].append(coords_file)
                        results['5prime']['matched_contigs'].add(contig_file)
                except Exception as e:
                    continue
    
    if '3prime_file' in query_data[seq_id]:
        query_file = query_data[seq_id]['3prime_file']
        
        printv(f"  Processing 3' end alignment...")
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(threads, len(telomere_contigs))) as executor:
            future_to_contig = {}
            for i, contig_file in enumerate(telomere_contigs, 1):
                contig_name = os.path.basename(contig_file).replace('.fa', '')
                prefix = os.path.join(output_dir, f"{seq_id}_3prime_{contig_name}")
                future = executor.submit(process_single_alignment, contig_file, query_file, prefix, 1)
                future_to_contig[future] = (contig_file, contig_name)
            
            for future in concurrent.futures.as_completed(future_to_contig):
                contig_file, contig_name = future_to_contig[future]
                try:
                    coords_file = future.result()
                    if coords_file and os.path.exists(coords_file):
                        results['3prime']['coords_files'].append(coords_file)
                        results['3prime']['matched_contigs'].add(contig_file)
                except Exception as e:
                    continue
    
    all_successful_coords = results['5prime']['coords_files'] + results['3prime']['coords_files']
    if all_successful_coords:
        merged_coords = os.path.join(output_dir, f"{seq_id}_successful.coords")
        if merge_coords_files_clean(all_successful_coords, merged_coords):
            results['merged_coords'] = merged_coords
    
    return results

def merge_coords_files_clean(coords_files: List[str], output_file: str) -> bool:
    if not coords_files:
        return False
    
    all_data_lines = []
    headers_written = False
    
    for coords_file in coords_files:
        if not os.path.exists(coords_file):
            continue
            
        with open(coords_file, 'r') as f:
            lines = f.readlines()
        
        data_started = False
        for line in lines:
            line = line.strip()
            
            if line.startswith('[') or line.startswith('='):
                data_started = True
                if not headers_written:
                    all_data_lines.append("NUCMER\n")
                    all_data_lines.append("[S1]\t[E1]\t|[S2]\t[E2]\t|[LEN 1]\t[LEN 2]\t|[%% IDY]\t|[TAGS]\n")
                    all_data_lines.append("=====================================================================================\n")
                    headers_written = True
                continue
            
            if data_started and line and not line.startswith('NUCMER'):
                all_data_lines.append(line + '\n')
    
    if all_data_lines:
        with open(output_file, 'w') as f:
            f.writelines(all_data_lines)
        
        for coords_file in coords_files:
            if os.path.exists(coords_file) and coords_file != output_file:
                os.remove(coords_file)
        
        return True
    
    return False

def cleanup_unused_files(output_dir: str, all_results: Dict, telomere_contigs: List[str]):
    printv(f"Cleaning up unused files...")
    
    all_matched_contigs = set()
    for results in all_results.values():
        all_matched_contigs.update(results.get('5prime', {}).get('matched_contigs', set()))
        all_matched_contigs.update(results.get('3prime', {}).get('matched_contigs', set()))
    
    for contig_file in telomere_contigs:
        if contig_file not in all_matched_contigs and os.path.exists(contig_file):
            os.remove(contig_file)
    
    printv(f"Cleanup completed")

def run_alignment_pipeline(
    contigs: str,
    query: str,
    output: str,
    threads: int = 32,
    min_repeats: int = 5,
    contig_min_repeats: int = 20,
    contig_min_length: int = 500,
    skip_extract: bool = False,
    verbose: bool = False
) -> Dict[str, Any]:
    
    global VERBOSE
    VERBOSE = verbose
    
    start_time = time.time()
    
    try:
        set_environment_limits()
        
        for f in [contigs, query]:
            if not os.path.exists(f):
                return {
                    'success': False,
                    'error': f"File does not exist: {f}"
                }
        
        os.makedirs(output, exist_ok=True)
        
        if verbose:
            print("="*70)
            print("Telomere Alignment Tool - Complex alignment, generate coords files (Part 1)")
            print("="*70)
            print(f"Original contigs: {contigs}")
            print(f"Query genome: {query}")
            print(f"Output directory: {output}")
            print(f"nucmer threads: {threads}")
            print(f"Chromosome telomere detection: minimum {min_repeats} repeats")
            print(f"Contig telomere detection: minimum {contig_min_repeats} repeats, minimum length {contig_min_length}bp")
            print("="*70)
        
        telomere_contigs_file = None
        telomere_count = 0
        total_count = 0
        
        if not skip_extract:
            if verbose:
                print(f"\nStep 1: Extracting telomere-containing contigs (using contig telomere detection parameters)")
            telomere_contigs_file = os.path.join(output, "telomere_contigs.fa")
            telomere_count, total_count = extract_telomere_contigs(
                contigs, telomere_contigs_file, 
                contig_min_repeats, contig_min_length
            )
            if telomere_count == 0:
                return {
                    'success': False,
                    'error': "No telomere-containing contigs found"
                }
        else:
            telomere_contigs_file = contigs
            if verbose:
                print(f"\nStep 1: Using provided contigs file")
            total_count = len(list(SeqIO.parse(telomere_contigs_file, "fasta")))
            telomere_count = total_count
        
        if verbose:
            print(f"\nStep 2: Analyzing chromosome telomere status (only detecting end telomeres)")
        chromosome_status = analyze_chromosome_telomeres_simple(
            query, 
            min_repeats=min_repeats,
            search_window=2000
        )
        
        analysis_file = os.path.join(output, "telomere_analysis.json")
        with open(analysis_file, 'w') as f:
            json.dump(chromosome_status, f, indent=2)
        
        needs_repair = [seq_id for seq_id, info in chromosome_status.items() if info['needs_repair']]
        if not needs_repair:
            if verbose:
                print(f"\nAll chromosomes have complete telomeres, no further analysis needed")
            return {
                'success': True,
                'telomere_contigs': telomere_contigs_file,
                'chromosome_status': chromosome_status,
                'query_data': {},
                'coords_file': None,
                'summary': {
                    'total_chromosomes': len(chromosome_status),
                    'needs_repair': 0,
                    'all_complete': True
                },
                'report_file': None,
                'output_dir': output
            }
        
        if verbose:
            print(f"\nChromosomes needing telomere analysis: {len(needs_repair)}")
            for seq_id in needs_repair:
                print(f"  - {seq_id}: {chromosome_status[seq_id]['status']}")
        
        if verbose:
            print(f"\nStep 3: Preparing chromosome end sequences")
        query_data = prepare_chromosome_ends(query, chromosome_status, output)
        
        if not query_data:
            return {
                'success': False,
                'error': "No chromosomes found for analysis"
            }
        
        query_data_file = os.path.join(output, "query_data.json")
        query_data_simple = {}
        for seq_id, data in query_data.items():
            query_data_simple[seq_id] = {
                '5prime_file': data.get('5prime_file', ''),
                '3prime_file': data.get('3prime_file', ''),
                'sequence_id': data['sequence'].id,
                'sequence_length': len(data['sequence'].seq),
                'sequence_description': data['sequence'].description
            }
        
        with open(query_data_file, 'w') as f:
            json.dump(query_data_simple, f, indent=2)
        
        if verbose:
            print(f"\nStep 4: Preparing telomere contigs")
        temp_contig_dir = os.path.join(output, "temp_contigs")
        os.makedirs(temp_contig_dir, exist_ok=True)
        
        telomere_contigs_list = []
        for i, record in enumerate(SeqIO.parse(telomere_contigs_file, "fasta"), 1):
            contig_file = os.path.join(temp_contig_dir, f"contig_{i:04d}.fa")
            SeqIO.write([record], contig_file, "fasta")
            telomere_contigs_list.append(contig_file)
        
        if verbose:
            print(f"Prepared {len(telomere_contigs_list)} telomere contigs for alignment")
        
        if verbose:
            print(f"\nStep 5: Starting telomere alignment")
            print("="*70)
        
        all_results = {}
        
        threads_per_chromosome = max(1, threads // len(query_data))
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(threads, len(query_data))) as executor:
            future_to_seq = {}
            for seq_id in query_data:
                future = executor.submit(
                    process_chromosome_telomeres_directional,
                    seq_id, query_data, telomere_contigs_list,
                    output, threads_per_chromosome
                )
                future_to_seq[future] = seq_id
            
            for future in concurrent.futures.as_completed(future_to_seq):
                seq_id = future_to_seq[future]
                try:
                    results = future.result()
                    all_results[seq_id] = results
                    
                    has_5prime = len(results.get('5prime', {}).get('coords_files', [])) > 0
                    has_3prime = len(results.get('3prime', {}).get('coords_files', [])) > 0
                    
                    if verbose:
                        if has_5prime or has_3prime:
                            print(f"✓ {seq_id}: Successful alignment")
                            if has_5prime:
                                print(f"  5' end: {len(results['5prime']['coords_files'])} alignments")
                            if has_3prime:
                                print(f"  3' end: {len(results['3prime']['coords_files'])} alignments")
                        else:
                            print(f"✗ {seq_id}: No successful alignments")
                            
                except Exception as e:
                    if verbose:
                        print(f"✗ {seq_id}: Alignment exception - {e}")
        
        cleanup_unused_files(output, all_results, telomere_contigs_list)
        
        if os.path.exists(temp_contig_dir):
            shutil.rmtree(temp_contig_dir, ignore_errors=True)
        
        if verbose:
            print(f"\nStep 6: Merging all chromosome successful alignment results")
        
        all_successful_coords = []
        for seq_id, results in all_results.items():
            if 'merged_coords' in results and os.path.exists(results['merged_coords']):
                all_successful_coords.append(results['merged_coords'])
        
        final_coords = None
        if all_successful_coords:
            final_coords = os.path.join(output, "all_successful_matches.coords")
            if merge_coords_files_clean(all_successful_coords, final_coords):
                if verbose:
                    print(f"Merged {len(all_successful_coords)} chromosome coords files to: {final_coords}")
        
        if verbose:
            print(f"\nStep 7: Saving alignment result information")
        
        results_summary = {}
        successful_chromosomes = 0
        total_successful_matches = 0
        
        for seq_id, results in all_results.items():
            has_5prime = len(results.get('5prime', {}).get('coords_files', [])) > 0
            has_3prime = len(results.get('3prime', {}).get('coords_files', [])) > 0
            
            results_summary[seq_id] = {
                'has_5prime_matches': has_5prime,
                'has_3prime_matches': has_3prime,
                '5prime_match_count': len(results.get('5prime', {}).get('coords_files', [])),
                '3prime_match_count': len(results.get('3prime', {}).get('coords_files', [])),
                'merged_coords': results.get('merged_coords', '')
            }
            
            if has_5prime or has_3prime:
                successful_chromosomes += 1
                total_successful_matches += len(results.get('5prime', {}).get('coords_files', [])) + len(results.get('3prime', {}).get('coords_files', []))
        
        summary_file = os.path.join(output, "alignment_summary.json")
        with open(summary_file, 'w') as f:
            json.dump(results_summary, f, indent=2)
        
        if verbose:
            print(f"\nStep 8: Generating final report")
            print("="*70)
        
        report_file = os.path.join(output, "alignment_report.txt")
        with open(report_file, 'w') as f:
            f.write("Telomere Alignment Tool - Alignment Result Report (Part 1)\n")
            f.write("="*60 + "\n")
            f.write(f"Original contigs: {contigs}\n")
            f.write(f"Query genome: {query}\n")
            f.write(f"Output directory: {output}\n")
            f.write(f"nucmer threads: {threads}\n")
            f.write(f"Chromosome telomere detection: {min_repeats} repeats\n")
            f.write(f"Contig telomere detection: {contig_min_repeats} repeats, {contig_min_length}bp\n\n")
            
            f.write("Alignment statistics:\n")
            f.write(f"- Total chromosomes: {len(chromosome_status)}\n")
            f.write(f"- Chromosomes needing repair: {len(needs_repair)}\n")
            f.write(f"- Successfully aligned chromosomes: {successful_chromosomes}/{len(query_data)}\n")
            f.write(f"- Total successful alignments: {total_successful_matches}\n")
            f.write(f"- Telomere contigs: {telomere_count}/{total_count}\n\n")
            
            f.write("Detailed alignment by chromosome:\n")
            for seq_id, summary in results_summary.items():
                f.write(f"{seq_id}:\n")
                f.write(f"  5' end alignments: {summary['5prime_match_count']}\n")
                f.write(f"  3' end alignments: {summary['3prime_match_count']}\n")
                if summary['merged_coords']:
                    f.write(f"  Merged coords file: {os.path.basename(summary['merged_coords'])}\n")
                f.write("\n")
        
        total_time = time.time() - start_time
        
        result = {
            'success': True,
            'telomere_contigs': telomere_contigs_file,
            'chromosome_status': chromosome_status,
            'query_data': query_data_simple,
            'coords_file': final_coords,
            'summary': {
                'total_chromosomes': len(chromosome_status),
                'needs_repair': len(needs_repair),
                'successful_chromosomes': successful_chromosomes,
                'total_successful_matches': total_successful_matches,
                'telomere_contigs_count': telomere_count,
                'total_contigs_count': total_count,
                'processing_time': total_time
            },
            'report_file': report_file,
            'output_dir': output,
            'files': {
                'telomere_analysis': analysis_file,
                'query_data': query_data_file,
                'alignment_summary': summary_file,
                'alignment_report': report_file,
                'coords_file': final_coords
            }
        }
        
        if verbose:
            print(f"\nAlignment completed!")
            print(f"Processing time: {total_time:.2f} seconds")
            print(f"Main output files:")
            for file_desc, file_path in [
                ("Telomere contigs", telomere_contigs_file),
                ("Chromosome telomere analysis", analysis_file),
                ("Chromosome end information", query_data_file),
                ("All successful alignments coords file", final_coords),
                ("Alignment result summary", summary_file),
                ("Alignment report", report_file)
            ]:
                if file_path and os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    print(f"  {file_desc}: {os.path.basename(file_path)} ({file_size:,} bytes)")
        
        return result
        
    except Exception as e:
        import traceback
        error_trace = traceback.format_exc()
        if verbose:
            print(f"Error during alignment process: {e}")
            print(error_trace)
        
        return {
            'success': False,
            'error': str(e),
            'traceback': error_trace
        }

def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description="Telomere Alignment Tool - Complex alignment, generate coords files (Part 1)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Complete alignment pipeline (chromosomes: 5 repeat detection, contigs: 20 repeat detection)
  python telomere_alignment.py -c contigs.fasta -q genome.fasta -o results
  
  # Using more threads
  python telomere_alignment.py -c contigs.fasta -q genome.fasta -o results -t 64
  
  # Custom parameters
  python telomere_alignment.py -c contigs.fasta -q genome.fasta -o results \
    --min-repeats 5 --contig-min-repeats 20 --contig-min-length 500
  
  # Skip telomere extraction step (use existing telomere contigs)
  python telomere_alignment.py -c telomere_contigs.fasta -q genome.fasta -o results --skip-extract
        """
    )
    
    parser.add_argument("-c", "--contigs", required=True, help="Original contigs file")
    parser.add_argument("-q", "--query", required=True, help="Genome file to analyze")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    
    parser.add_argument("-t", "--threads", type=int, default=32, help="nucmer threads, default 32")
    parser.add_argument("--min-repeats", type=int, default=5, help="Minimum chromosome telomere repeats, default 5")
    parser.add_argument("--contig-min-repeats", type=int, default=20, help="Minimum contig telomere repeats, default 20")
    parser.add_argument("--contig-min-length", type=int, default=500, help="Minimum contig telomere length, default 500bp")
    parser.add_argument("--skip-extract", action="store_true", help="Skip telomere extraction step")
    parser.add_argument("--verbose", action="store_true", help="Show detailed output")
    
    return parser.parse_args()

def main():
    args = parse_command_line_args()
    
    result = run_alignment_pipeline(
        contigs=args.contigs,
        query=args.query,
        output=args.output,
        threads=args.threads,
        min_repeats=args.min_repeats,
        contig_min_repeats=args.contig_min_repeats,
        contig_min_length=args.contig_min_length,
        skip_extract=args.skip_extract,
        verbose=args.verbose
    )
    
    if result['success']:
        summary = result['summary']
        print("\n" + "="*70)
        print("Alignment completed!")
        print("="*70)
        print(f"Total chromosomes: {summary['total_chromosomes']}")
        print(f"Chromosomes needing repair: {summary['needs_repair']}")
        print(f"Successfully aligned chromosomes: {summary['successful_chromosomes']}")
        print(f"Total successful alignments: {summary['total_successful_matches']}")
        print(f"Telomere contigs: {summary['telomere_contigs_count']}/{summary['total_contigs_count']}")
        print(f"Processing time: {summary['processing_time']:.2f} seconds")
        
        if result['coords_file']:
            print(f"\nMerged coords file: {result['coords_file']}")
        
        print(f"Detailed report: {result['report_file']}")
        
        print(f"\nMain output files:")
        files_to_check = [
            (os.path.join(args.output, "telomere_contigs.fa"), "Telomere contigs"),
            (os.path.join(args.output, "telomere_analysis.json"), "Chromosome telomere analysis"),
            (os.path.join(args.output, "query_data.json"), "Chromosome end information"),
            (os.path.join(args.output, "all_successful_matches.coords"), "All successful alignments coords file"),
            (os.path.join(args.output, "alignment_summary.json"), "Alignment result summary"),
            (os.path.join(args.output, "alignment_report.txt"), "Alignment report")
        ]
        
        for file_path, description in files_to_check:
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path)
                print(f"  {description}: {os.path.basename(file_path)} ({file_size:,} bytes)")
        
        print("="*70)
        print("\nNext step: Run part 2 script for sequence extraction and merging")
        print("Example command: python telomere_extraction_and_merge.py -o results -c telomere_contigs.fa -q genome.fasta")
    else:
        print(f"\nError: {result['error']}")
        if args.verbose and 'traceback' in result:
            print(f"Detailed error information:\n{result['traceback']}")
        sys.exit(1)

if __name__ == "__main__":
    main()