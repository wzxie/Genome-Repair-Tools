#!/usr/bin/env python3
"""
Genome Assembly and Gap Filling Integrated Pipeline - Simplified Version
Automated pipeline for assembly, gap detection, and filling
Version: 1.2 (Added TGS mode control)
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from collections import defaultdict
import re


try:
    sys.path.append('.')  # Add current directory to path
    from chromosome_analyzer import analyze_chromosomes
    
    from gap_filler import gap_filler as run_gap_filler
    from gap_filler import GapFillerAPI  # Import API class

    from assembly_script import (
        AssemblyPipeline, 
        AssemblyConfig, 
        run_pipeline_from_dict
    )
    

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
except ImportError as e:
    print(f"Failed to import modules: {e}")
    print("Please ensure all scripts and dependencies are installed")
    print("Dependencies: pip install biopython")
    sys.exit(1)

class AssemblyGapFiller:
    """Integrated Genome Assembly and Gap Filling Pipeline"""
    
    def __init__(self):
        self.work_dir = Path.cwd()
        self.timestamp = time.strftime("%Y%m%d_%H%M%S")
        self.verbose = False
        
        self.temp_files = []
        
    def log(self, message: str, level: str = "INFO"):
        """Log output"""
        prefix = {
            "INFO": "[INFO]",
            "WARN": "[WARNING]",
            "ERROR": "[ERROR]",
            "STEP": "[STEP]",
            "SUCCESS": "[SUCCESS]"
        }.get(level, "[INFO]")
        
        if level == "STEP":
            print(f"\n{'='*70}")
            print(f"{prefix} {message}")
            print(f"{'='*70}")
        elif level == "SUCCESS":
            print(f"\n{'✓'*3} {message} {'✓'*3}")
        else:
            print(f"{prefix} {message}")
    
    def cleanup_temp_files(self):
        """Clean up temporary files"""
        if not self.verbose:
            for temp_file in self.temp_files:
                if os.path.exists(temp_file):
                    try:
                        os.remove(temp_file)
                    except:
                        pass
        self.temp_files.clear()
    
    def register_temp_file(self, filepath: str):
        """Register a temporary file"""
        self.temp_files.append(filepath)
    
    def analyze_genome_with_module(self, input_fasta: str, output_prefix: str, 
                                  threads: int = 16, filter_fragments: bool = False,
                                  min_fragment_length: int = 100000) -> Optional[Dict]:
        """Analyze genome using chromosome_analyzer module"""
        self.log(f"Analyzing genome: {os.path.basename(input_fasta)}", "INFO")
        
        try:
            results = analyze_chromosomes(
                input_file=input_fasta,
                detect_telomeres=True,
                analyze_fragments=True,
                save_diagrams=True,
                min_gap_size=1,
                min_fragment_length=min_fragment_length if filter_fragments else None,
                threads=threads,
                min_repeats=5,
                search_window=2000,
                internal_flank=5000,
                diagram_width=50,
                output_prefix=output_prefix,
                memory_optimized=True,
                verbose=self.verbose
            )
            
            self.log(f"Analysis completed: {results['summary']['total_sequences']} chromosomes", "INFO")
            
            for suffix in ['_gap_coordinates.txt', '_chromosome_diagrams.txt', 
                          '_telomere_info.txt', '_internal_telomeres.txt', 
                          '_fragment_info.txt']:
                temp_file = f"{output_prefix}{suffix}"
                if os.path.exists(temp_file):
                    self.register_temp_file(temp_file)
            
            return results
            
        except Exception as e:
            self.log(f"Genome analysis failed: {e}", "ERROR")
            return None
    
    def get_filtered_fasta_from_analysis(self, analysis_prefix: str) -> Optional[str]:
        """Get filtered FASTA file from analysis results"""
        filtered_file = f"{analysis_prefix}_filtered_sequences.fa"
        if os.path.exists(filtered_file):
            self.register_temp_file(filtered_file)
            return filtered_file
        return None
    
    def find_gaps_in_sequence(self, sequence: str, gap_char: str = 'N', 
                             min_gap_size: int = 1) -> List[Tuple[int, int, int]]:
        """Find gaps in a single sequence"""
        gaps = []
        gap_pattern = re.compile(f'{re.escape(gap_char)}+')
        
        for match in gap_pattern.finditer(sequence):
            start = match.start()
            end = match.end()
            gap_length = end - start
            
            if gap_length >= min_gap_size:
                gaps.append((start, end, gap_length))
        
        return gaps
    
    def analyze_and_sort_chromosomes(self, input_fasta: str, 
                                    no_gap_output: str = "0_gap.fa",
                                    with_gap_output: str = "with_gap.fa") -> Tuple[str, str, Dict]:
        """Analyze FASTA file and sort chromosomes"""
        self.log(f"Detecting gaps and sorting chromosomes: {os.path.basename(input_fasta)}", "INFO")
        
        try:
            for output_file in [no_gap_output, with_gap_output]:
                output_dir = os.path.dirname(output_file)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir, exist_ok=True)
            
            records = list(SeqIO.parse(input_fasta, "fasta"))
            
            no_gap_records = []
            with_gap_records = []
            gap_stats = defaultdict(lambda: {'count': 0, 'total_length': 0})
            
            for record in records:
                sequence = str(record.seq)
                gaps = self.find_gaps_in_sequence(sequence, 'N', 1)
                
                if not gaps:
                    no_gap_records.append(record)
                else:
                    with_gap_records.append(record)
                    
                    gap_stats[record.id]['count'] = len(gaps)
                    gap_stats[record.id]['total_length'] = sum(g[2] for g in gaps)
                    gap_stats[record.id]['positions'] = [(g[0], g[1]) for g in gaps]

            if no_gap_records:
                SeqIO.write(no_gap_records, no_gap_output, "fasta")
                self.log(f"Chromosomes without gaps: {len(no_gap_records)} → {os.path.basename(no_gap_output)}", "INFO")
            
            if with_gap_records:
                SeqIO.write(with_gap_records, with_gap_output, "fasta")
                self.log(f"Chromosomes with gaps: {len(with_gap_records)} → {os.path.basename(with_gap_output)}", "INFO")

            total_chromosomes = len(records)
            total_with_gaps = len(with_gap_records)
            total_gap_count = sum(stats['count'] for stats in gap_stats.values())
            total_gap_length = sum(stats['total_length'] for stats in gap_stats.values())
            
            summary = {
                'total_chromosomes': total_chromosomes,
                'chromosomes_without_gaps': len(no_gap_records),
                'chromosomes_with_gaps': total_with_gaps,
                'total_gap_count': total_gap_count,
                'total_gap_length': total_gap_length,
                'gap_stats': dict(gap_stats)
            }
            
            self.log(f"Sorting completed: {len(no_gap_records)} without gaps, {total_with_gaps} with gaps, {total_gap_count} total gaps", "SUCCESS")
            
            return no_gap_output, with_gap_output, summary
            
        except Exception as e:
            self.log(f"Chromosome sorting failed: {e}", "ERROR")
            return "", "", {}
    
    def simple_split_contigs(self, input_file: str, output_file: str) -> bool:
        """Simple contig splitting method, avoiding complex processing"""
        self.log(f"Splitting contigs: {os.path.basename(input_file)}", "INFO")
        
        try:
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                in_sequence = False
                current_seq = []
                seq_count = 0
                
                for line in infile:
                    line = line.strip()
                    
                    if line.startswith('>'):
                        if current_seq:
                            sequence = ''.join(current_seq)
                            parts = re.split(r'N{10,}', sequence.upper())  # At least 10 Ns
                            for i, part in enumerate(parts):
                                if part and len(part) > 0:
                                    seq_count += 1
                                    outfile.write(f">split_{seq_count}\n")
                                    for j in range(0, len(part), 80):
                                        outfile.write(part[j:j+80] + "\n")
                        
                        # Reset current sequence
                        current_seq = []
                    
                    else:
                        current_seq.append(line)
                
                if current_seq:
                    sequence = ''.join(current_seq)
                    parts = re.split(r'N{10,}', sequence.upper())
                    for i, part in enumerate(parts):
                        if part and len(part) > 0:
                            seq_count += 1
                            outfile.write(f">split_{seq_count}\n")
                            for j in range(0, len(part), 80):
                                outfile.write(part[j:j+80] + "\n")
            
            if seq_count > 0:
                self.log(f"Splitting completed: {seq_count} sequences total", "SUCCESS")
                return True
            else:
                output_dir = os.path.dirname(output_file)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir, exist_ok=True)
                import shutil
                shutil.copy2(input_file, output_file)
                self.log(f"Splitting failed, copying file directly", "WARN")
                return True
                
        except Exception as e:
            self.log(f"Simple splitting failed: {e}", "ERROR")
            return False
    
    def assemble_with_module(self, hifi_files: List[str] = None,
                            ont_files: List[str] = None,
                            clr_files: List[str] = None,
                            threads: int = 32,
                            output_prefix: str = "assembly",
                            assembly_tools: List[str] = None) -> Optional[str]:
        """Perform assembly using assembly_script module"""
        self.log("Assembling contigs", "STEP")
        self.log(f"Assembly tools: {assembly_tools or ['all']}", "INFO")
        
        try:
            config_dict = {
                'hifi_files': hifi_files or [],
                'ont_ul_files': ont_files or [],
                'clr_files': clr_files or [],
                'threads': threads,
                'output_prefix': output_prefix,
                'tools_to_run': assembly_tools or ['all']
            }
            
            if not (config_dict['hifi_files'] or config_dict['ont_ul_files'] or config_dict['clr_files']):
                self.log("Warning: No sequencing data provided, skipping assembly step", "WARN")
                return None
            
            self.log(f"Input data: HiFi={len(config_dict['hifi_files'])}, ONT={len(config_dict['ont_ul_files'])}, CLR={len(config_dict['clr_files'])}", "INFO")
            
            results = run_pipeline_from_dict(config_dict)
            
            if not results or 'results' not in results:
                self.log("Assembly failed or no results generated", "ERROR")
                return None
            
            if 'merged' in results['results']:
                merged_result = results['results'].get('merged', {})
                if merged_result.get('status') == 'success' and merged_result.get('file'):
                    assembly_file = str(merged_result['file'])
                    
                    if os.path.exists(assembly_file):
                        self.log(f"Assembly completed: {assembly_file}", "SUCCESS")
                        return assembly_file
            
            for tool_name, tool_result in results.get('results', {}).items():
                if tool_name != 'merged' and tool_result.get('status') == 'success':
                    files = tool_result.get('files', [])
                    if files:
                        assembly_file = str(files[0])
                        self.log(f"Using {tool_name} assembly result: {assembly_file}", "INFO")
                        return assembly_file
            
            self.log("No valid assembly results found", "WARN")
            return None
            
        except Exception as e:
            self.log(f"Assembly failed: {e}", "ERROR")
            return None
    
    def merge_contigs_files(self, contig_files: List[str], output_file: str) -> bool:
        """Merge multiple contigs files"""
        self.log(f"Merging {len(contig_files)} contigs files", "INFO")
        
        try:
            all_records = []
            seen_ids = set()
            
            for contig_file in contig_files:
                if os.path.exists(contig_file):
                    self.log(f"Processing file: {os.path.basename(contig_file)}", "DEBUG")
                    for record in SeqIO.parse(contig_file, "fasta"):
                        original_id = record.id
                        if original_id in seen_ids:
                            base_id = original_id
                            counter = 1
                            while f"{base_id}_{counter}" in seen_ids:
                                counter += 1
                            record.id = f"{base_id}_{counter}"
                            record.description = ""
                        
                        seen_ids.add(record.id)
                        all_records.append(record)
            
            if all_records:
                output_dir = os.path.dirname(output_file)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir, exist_ok=True)
                    
                SeqIO.write(all_records, output_file, "fasta")
                self.log(f"Merging completed: {output_file} ({len(all_records)} sequences)", "SUCCESS")
                return True
            else:
                self.log("No contig sequences found", "ERROR")
                return False
                
        except Exception as e:
            self.log(f"Merging contigs failed: {e}", "ERROR")
            return False
    
    def fill_gaps_with_api(self, draft_genome: str, contigs: str, 
                          output_prefix: str, threads: int = 16,
                          use_tgs: bool = False) -> Optional[str]:
        """Fill gaps using GapFillerAPI interface"""
        self.log(f"Filling gaps (TGS mode: {'ON' if use_tgs else 'OFF'})", "STEP")
        
        try:
            gap_patches = [contigs] if isinstance(contigs, str) else contigs
            
            if '/' in output_prefix:
                output_dir = os.path.dirname(output_prefix)
                base_prefix = os.path.basename(output_prefix)
            else:
                output_dir = "."
                base_prefix = output_prefix
            
            os.makedirs(output_dir, exist_ok=True)
            
            filler = GapFillerAPI(verbose=self.verbose)
            
            result = filler.fill_gaps(
                draft_fasta=draft_genome,
                gap_patches=gap_patches,
                output_dir=output_dir,
                output_prefix=base_prefix,
                flank_size=5000,
                min_alignment_length=1000,
                min_identity=40.0,
                max_fill_length=1000000,  # Max 1Mb fill length
                aligner='minimap2',
                threads=threads,
                overwrite=True,
                keep_temp=False,
                return_dict=True,
                use_tgs=use_tgs,  # Use passed parameter to control TGS
                tgs_minmap="-x asm5",  # Fixed as asm5
                tgs_threads=threads,
                tgs_type="ont",
                keep_tgs_temp=self.verbose,
                preprocess_contigs=True,
                min_gap_length=100,  # Minimum 100bp gap to fill
                min_contig_length=1000  # Minimum 1kb contig to use
            )
            
            if result and 'filled_fasta' in result:
                filled_genome = result['filled_fasta']
                if os.path.exists(filled_genome):
                    gaps_closed = result.get('gaps_closed', 0)
                    gaps_remaining = result.get('gaps_remaining', 0)
                    
                    if use_tgs and 'tgs_stage' in result and result['tgs_stage']:
                        tgs_closed = result['tgs_stage'].get('gaps_closed_by_tgs', 0)
                        self.log(f"TGS stage filled: {tgs_closed} gaps", "INFO")
                    
                    self.log(f"Gap filling completed: {gaps_closed} gaps filled, {gaps_remaining} gaps remaining", "SUCCESS")
                    
                    for file_key in ['detail_file', 'modified_agp', 'stat_file', 'combined_report_file']:
                        if file_key in result and os.path.exists(result[file_key]):
                            self.register_temp_file(result[file_key])
                    
                    try:
                        filler.cleanup()
                    except Exception as cleanup_error:
                        self.log(f"API cleanup warning: {cleanup_error}", "WARN")
                    
                    return filled_genome
            
            return None
            
        except Exception as e:
            self.log(f"API gap filling failed: {e}", "ERROR")
            return None
    
    def process_contigs_files(self, contigs_files: List[str]) -> List[str]:
        """Process all contigs files, return processed file list"""
        processed_files = []
        
        for i, contig_file in enumerate(contigs_files):
            if not os.path.exists(contig_file):
                self.log(f"Contigs file not found: {contig_file}", "ERROR")
                continue
            
            self.log(f"Processing contigs file {i+1}: {os.path.basename(contig_file)}", "INFO")
            
            file_size = os.path.getsize(contig_file)
            if file_size == 0:
                self.log(f"File is empty: {contig_file}", "WARN")
                continue
            
            base_name = os.path.basename(contig_file)
            name_without_ext = os.path.splitext(base_name)[0]
            processed_file = f"{name_without_ext}_processed_{self.timestamp}.fa"
            
            if self.simple_split_contigs(contig_file, processed_file):
                # Check processed file
                if os.path.exists(processed_file) and os.path.getsize(processed_file) > 0:
                    processed_files.append(processed_file)
                    self.register_temp_file(processed_file)
                else:
                    self.log(f"Processed file invalid: {processed_file}", "WARN")
            else:
                self.log(f"Processing failed: {contig_file}", "WARN")
        
        return processed_files
    
    def filter_small_fragments(self, input_fasta: str, 
                              min_fragment_length: int = 100000,
                              output_prefix: str = "filtered") -> str:
        """Filter out fragments smaller than specified length"""
        self.log(f"Filtering fragments smaller than {min_fragment_length}bp", "INFO")
        
        try:
            analysis_results = self.analyze_genome_with_module(
                input_fasta=input_fasta,
                output_prefix=output_prefix,
                threads=16,
                filter_fragments=True,
                min_fragment_length=min_fragment_length
            )
            
            if analysis_results:
                # Get filtered file
                filtered_file = f"{output_prefix}_filtered_sequences.fa"
                
                if os.path.exists(filtered_file):
                    original_length = self.get_fasta_length(input_fasta)
                    filtered_length = self.get_fasta_length(filtered_file)
                    
                    if original_length > 0:
                        filtered_percent = (original_length - filtered_length) / original_length * 100
                        self.log(f"Filtering results: original {original_length:,}bp → filtered {filtered_length:,}bp "
                                f"({filtered_percent:.1f}% reduction)", "SUCCESS")
                    else:
                        self.log(f"Filtering results: filtered {filtered_length:,}bp", "SUCCESS")
                    
                    return filtered_file

            self.log("Filtering failed, returning original file", "WARN")
            return input_fasta
            
        except Exception as e:
            self.log(f"Fragment filtering failed: {e}", "ERROR")
            return input_fasta
    
    def get_fasta_length(self, fasta_file: str) -> int:
        """Calculate total length of FASTA file"""
        total_length = 0
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                total_length += len(record.seq)
        except Exception as e:
            self.log(f"Failed to calculate FASTA length {fasta_file}: {e}", "WARN")
        return total_length
    
    def process_optimized_filling(self, draft_genome: str, contigs_file: str,
                                 threads: int = 16,
                                 output_dir: str = "optimized_filling",
                                 use_tgs: bool = False) -> Tuple[str, List[Dict]]:
        """Optimized filling process: first round filling, filter fragments, second round filling"""
        self.log(f"Starting optimized filling process (TGS mode: {'ON' if use_tgs else 'OFF'})", "STEP")
        
        os.makedirs(output_dir, exist_ok=True)
        
        results_log = []
        
        self.log("Stage 1: Initial analysis (no fragment filtering)", "STEP")

        no_gap_file = os.path.join(output_dir, "initial_no_gap.fa")
        with_gap_file = os.path.join(output_dir, "initial_with_gap.fa")
        
        no_gap_file, with_gap_file, initial_gap_summary = self.analyze_and_sort_chromosomes(
            draft_genome,
            no_gap_file,
            with_gap_file
        )
        
        if not os.path.exists(with_gap_file) or os.path.getsize(with_gap_file) == 0:
            self.log("Initial chromosomes have no gaps, no filling needed", "SUCCESS")
            return no_gap_file, results_log
        
        results_log.append({
            'stage': 'initial_analysis',
            'no_gap_chromosomes': initial_gap_summary.get('chromosomes_without_gaps', 0),
            'with_gap_chromosomes': initial_gap_summary.get('chromosomes_with_gaps', 0),
            'total_gaps': initial_gap_summary.get('total_gap_count', 0)
        })
        
        self.log("Stage 2: First round filling", "STEP")
        
        round1_prefix = os.path.join(output_dir, "round1_filled")
        filled_round1 = self.fill_gaps_with_api(
            with_gap_file,
            contigs_file,
            round1_prefix,
            threads,
            use_tgs=use_tgs  # Pass use_tgs parameter
        )
        
        if not filled_round1:
            self.log("First round filling failed", "ERROR")
            return draft_genome, results_log
        
        round1_no_gap = os.path.join(output_dir, "round1_no_gap.fa")
        round1_with_gap = os.path.join(output_dir, "round1_with_gap.fa")
        
        _, _, round1_gap_summary = self.analyze_and_sort_chromosomes(
            filled_round1,
            round1_no_gap,
            round1_with_gap
        )
        
        results_log.append({
            'stage': 'round1_filling',
            'filled_file': filled_round1,
            'gaps_remaining': round1_gap_summary.get('total_gap_count', 0),
            'chromosomes_with_gaps': round1_gap_summary.get('chromosomes_with_gaps', 0),
            'use_tgs': use_tgs  # Record TGS usage status
        })

        if round1_gap_summary.get('total_gap_count', 0) == 0:
            self.log("No gaps after first round filling, process complete", "SUCCESS")
            return filled_round1, results_log

        self.log("Stage 3: Filtering small fragments after filling (<100kb)", "STEP")
        
        filtered_file = self.filter_small_fragments(
            filled_round1,
            min_fragment_length=100000,  # Fixed 100kb
            output_prefix=os.path.join(output_dir, "post_fill_filtered")
        )
        
        if filtered_file == filled_round1:
            self.log("No fragments need filtering", "INFO")
            filtered_genome = filled_round1
        else:
            filtered_no_gap = os.path.join(output_dir, "filtered_no_gap.fa")
            filtered_with_gap = os.path.join(output_dir, "filtered_with_gap.fa")
            
            _, _, filtered_gap_summary = self.analyze_and_sort_chromosomes(
                filtered_file,
                filtered_no_gap,
                filtered_with_gap
            )
            
            results_log.append({
                'stage': 'post_fill_filtering',
                'filtered_file': filtered_file,
                'gaps_after_filter': filtered_gap_summary.get('total_gap_count', 0),
                'chromosomes_after_filter': filtered_gap_summary.get('chromosomes_with_gaps', 0)
            })
            
            filtered_genome = filtered_file
        
        if not os.path.exists(filtered_with_gap) or filtered_gap_summary.get('total_gap_count', 0) == 0:
            self.log("No gaps after filtering, process complete", "SUCCESS")
            return filtered_genome, results_log
        
        self.log("Stage 4: Second round filling", "STEP")
        
        round2_prefix = os.path.join(output_dir, "round2_filled")
        filled_round2 = self.fill_gaps_with_api(
            filtered_genome,
            contigs_file,
            round2_prefix,
            threads,
            use_tgs=use_tgs  # Pass use_tgs parameter
        )
        
        if not filled_round2:
            self.log("Second round filling failed, returning first round result", "WARN")
            return filled_round1, results_log
        
        round2_no_gap = os.path.join(output_dir, "round2_no_gap.fa")
        round2_with_gap = os.path.join(output_dir, "round2_with_gap.fa")
        
        _, _, round2_gap_summary = self.analyze_and_sort_chromosomes(
            filled_round2,
            round2_no_gap,
            round2_with_gap
        )
        
        results_log.append({
            'stage': 'round2_filling',
            'filled_file': filled_round2,
            'gaps_remaining': round2_gap_summary.get('total_gap_count', 0),
            'chromosomes_with_gaps': round2_gap_summary.get('chromosomes_with_gaps', 0),
            'use_tgs': use_tgs  # Record TGS usage status
        })
        
        self.log("Optimized filling process completed", "SUCCESS")
        
        return filled_round2, results_log
    
    def create_final_genome(self, original_no_gap: str, filled_genome: str, 
                           output_file: str, prefix: str = "final") -> Dict:
        """Create final genome file"""
        self.log("Creating final genome", "STEP")
        
        try:
            output_dir = os.path.dirname(output_file) if '/' in output_file else "."
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            
            final_records = []
            
            if os.path.exists(original_no_gap):
                for record in SeqIO.parse(original_no_gap, "fasta"):
                    final_records.append(record)
            
            filled_no_gap = os.path.join(output_dir, f"{prefix}_filled_no_gap.fa")
            filled_with_gap = os.path.join(output_dir, f"{prefix}_filled_with_gap.fa")
            
            filled_no_gap, filled_with_gap, filled_summary = self.analyze_and_sort_chromosomes(
                filled_genome,
                filled_no_gap,
                filled_with_gap
            )
            
            if os.path.exists(filled_no_gap):
                for record in SeqIO.parse(filled_no_gap, "fasta"):
                    # Modify ID to avoid duplication
                    record.id = f"filled_{record.id}"
                    final_records.append(record)
            
            if os.path.exists(filled_with_gap):
                for record in SeqIO.parse(filled_with_gap, "fasta"):
                    record.id = f"unfilled_{record.id}"
                    final_records.append(record)
            
            SeqIO.write(final_records, output_file, "fasta")
            
            total_bases = sum(len(record.seq) for record in final_records)
            filled_chromosomes = len([r for r in final_records if r.id.startswith('filled_')])
            unfilled_chromosomes = len([r for r in final_records if r.id.startswith('unfilled_')])
            original_chromosomes = len(final_records) - filled_chromosomes - unfilled_chromosomes
            
            summary = {
                'output_file': output_file,
                'total_chromosomes': len(final_records),
                'original_no_gap': original_chromosomes,
                'filled_no_gap': filled_chromosomes,
                'still_with_gap': unfilled_chromosomes,
                'total_bases': total_bases,
                'filled_summary': filled_summary
            }
            
            self.log(f"Final genome created: {output_file}", "SUCCESS")
            self.log(f"Chromosome statistics: original {original_chromosomes}, filled without gaps {filled_chromosomes}, still with gaps {unfilled_chromosomes}", "INFO")
            self.log(f"Total length: {total_bases:,} bp", "INFO")
            
            return summary
            
        except Exception as e:
            self.log(f"Final genome creation failed: {e}", "ERROR")
            return {}
    
    def run_complete_pipeline(self, query_fasta: str, contigs_files: List[str] = None,
                             hifi_files: List[str] = None, ont_files: List[str] = None,
                             clr_files: List[str] = None, threads: int = 16,
                             output_prefix: str = "filled_genome",
                             verbose: bool = False,
                             assembly_tools: List[str] = None,
                             use_tgs: bool = False) -> Dict:
        """Run complete assembly and gap filling pipeline"""
        self.log(f"Starting integrated genome assembly and gap filling pipeline (TGS mode: {'ON' if use_tgs else 'OFF'})", "STEP")
        self.verbose = verbose
        
        if assembly_tools:
            self.log(f"Specified assembly tools: {', '.join(assembly_tools)}", "INFO")
        else:
            self.log("Using default assembly tool configuration", "INFO")
        
        final_output = f"{output_prefix}_final.fa"
        summary_file = f"{output_prefix}_summary.json"
        
        output_dir = os.path.dirname(output_prefix) if '/' in output_prefix else "."
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            self.log(f"Created output directory: {output_dir}", "INFO")
        
        self.log("Step 1: Initial genome analysis", "STEP")
        analysis_prefix = os.path.join(output_dir, f"initial_analysis_{self.timestamp}")
        initial_analysis = self.analyze_genome_with_module(
            query_fasta, 
            analysis_prefix, 
            threads,
            filter_fragments=False  # No fragment filtering
        )
        
        if not initial_analysis:
            return {"status": "error", "message": "Initial genome analysis failed"}

        self.log("Step 2: Detect gaps and sort chromosomes", "STEP")
        original_no_gap = os.path.join(output_dir, f"{os.path.basename(output_prefix)}_original_no_gap.fa")
        original_with_gap = os.path.join(output_dir, f"{os.path.basename(output_prefix)}_original_with_gap.fa")
        
        original_no_gap, original_with_gap, gap_summary = self.analyze_and_sort_chromosomes(
            query_fasta,
            original_no_gap,
            original_with_gap
        )
        
        if not os.path.exists(original_with_gap) or os.path.getsize(original_with_gap) == 0:
            self.log("All chromosomes have no gaps, no filling needed", "SUCCESS")
            if os.path.exists(original_no_gap):
                import shutil
                shutil.copy2(original_no_gap, final_output)
            
            result = {
                "status": "success",
                "message": "All chromosomes have no gaps",
                "final_genome": final_output,
                "no_gap_chromosomes": gap_summary.get('chromosomes_without_gaps', 0)
            }
            return result

        self.log("Step 3: Prepare contigs files", "STEP")
        contig_files_to_use = []

        if hifi_files or ont_files or clr_files:
            self.log("Raw data detected, performing assembly", "INFO")
            assembly_prefix = os.path.join(output_dir, f"assembled_{self.timestamp}")
            assembled_contigs = self.assemble_with_module(
                hifi_files=hifi_files,
                ont_files=ont_files,
                clr_files=clr_files,
                threads=threads,
                output_prefix=assembly_prefix,
                assembly_tools=assembly_tools
            )
            
            if assembled_contigs:
                # Process assembled contigs
                processed_assembled = self.process_contigs_files([assembled_contigs])
                if processed_assembled:
                    contig_files_to_use.extend(processed_assembled)
                    self.log(f"Using assembly result: {os.path.basename(assembled_contigs)}", "INFO")
        else:
            self.log("No raw sequencing data provided, skipping assembly step", "INFO")

        if contigs_files:
            processed_user_contigs = self.process_contigs_files(contigs_files)
            if processed_user_contigs:
                contig_files_to_use.extend(processed_user_contigs)
                self.log(f"Using user-provided contigs files: {len(contigs_files)} files", "INFO")
        
        if not contig_files_to_use:
            self.log("Error: No usable contigs files", "ERROR")
            # Try to use original contigs files (unprocessed)
            if contigs_files:
                self.log("Trying to use original contigs files", "WARN")
                contig_files_to_use = contigs_files
            elif 'assembled_contigs' in locals() and assembled_contigs:
                self.log("Trying to use original assembly file", "WARN")
                contig_files_to_use = [assembled_contigs]
            else:
                return {"status": "error", "message": "No usable contigs files"}
        
        if len(contig_files_to_use) == 1:
            merged_contigs = contig_files_to_use[0]
        else:
            merged_contigs = os.path.join(output_dir, f"merged_contigs_{self.timestamp}.fa")
            if not self.merge_contigs_files(contig_files_to_use, merged_contigs):
                # If merging fails, use first file
                self.log("Merging failed, using first contigs file", "WARN")
                merged_contigs = contig_files_to_use[0] if contig_files_to_use else None
        
        if not merged_contigs or not os.path.exists(merged_contigs):
            self.log("Error: No valid contigs files", "ERROR")
            return {"status": "error", "message": "No valid contigs files"}

        self.log("Step 4: Optimized filling process", "STEP")
        self.log(f"Process: First round filling → Filter <100kb fragments → Second round filling", "INFO")
        
        optimized_dir = os.path.join(output_dir, "optimized_filling")
        final_filled_genome, filling_results = self.process_optimized_filling(
            original_with_gap,
            merged_contigs,
            threads=threads,
            output_dir=optimized_dir,
            use_tgs=use_tgs  # Pass use_tgs parameter
        )

        self.log("Step 5: Create final genome", "STEP")
        final_summary = self.create_final_genome(
            original_no_gap,
            final_filled_genome,
            final_output,
            os.path.basename(output_prefix)
        )

        self.log("Step 6: Final analysis", "STEP")
        final_analysis_prefix = os.path.join(output_dir, f"final_analysis_{self.timestamp}")
        final_analysis = self.analyze_genome_with_module(
            final_output,
            final_analysis_prefix,
            threads,
            filter_fragments=False
        )

        complete_summary = {
            "status": "success",
            "timestamp": self.timestamp,
            "version": "1.2 (Integrated Pipeline)",
            "input_file": query_fasta,
            "output_file": final_output,
            "threads_used": threads,
            "gap_analysis": gap_summary,
            "initial_analysis": initial_analysis.get('summary', {}) if initial_analysis else {},
            "final_analysis": final_analysis.get('summary', {}) if final_analysis else {},
            "contigs_used": contig_files_to_use,
            "assembly_tools": assembly_tools or ["all"],
            "gap_filler_mode": f"{'TGS mode' if use_tgs else 'Flanking alignment mode'} (asm5)",
            "use_tgs": use_tgs,
            "min_fragment_length": 100000,
            "min_gap_size": 1,
            "max_gap_size": 1000000,
            "filling_stages": filling_results,
            "final_summary": final_summary,
            "temp_files": self.temp_files
        }
        
        with open(summary_file, 'w') as f:
            json.dump(complete_summary, f, indent=2, ensure_ascii=False)
        
        self.log(f"Complete report saved: {summary_file}", "INFO")

        if not self.verbose:
            self.cleanup_temp_files()
        
        self.log("Genome assembly and gap filling pipeline completed!", "SUCCESS")
        
        return complete_summary


# ============================================================================
# Standardized interface functions
# ============================================================================

def assembly_gap_fill(query_fasta: str, 
                     contigs_files: List[str] = None,
                     hifi_files: List[str] = None,
                     ont_files: List[str] = None,
                     clr_files: List[str] = None,
                     threads: int = 16,
                     output_prefix: str = "filled_genome",
                     verbose: bool = False,
                     assembly_tools: List[str] = None,
                     use_tgs: bool = False) -> Dict:
    """
    Main function for integrated genome assembly and gap filling pipeline
    
    Parameters:
        query_fasta: Genome FASTA file that needs gap filling
        contigs_files: List of reference contigs files (optional)
        hifi_files: List of HiFi sequencing data files (optional)
        ont_files: List of ONT sequencing data files (optional)
        clr_files: List of CLR sequencing data files (optional)
        threads: Number of threads to use (default: 16)
        output_prefix: Output file prefix (default: filled_genome)
        verbose: Show detailed output (default: False)
        assembly_tools: List of assembly tools, e.g., ['hifiasm', 'flye'] (default: None, use all available)
        use_tgs: Whether to enable TGS-GapCloser for pre-filling (default: False)
    
    Returns:
        Dictionary containing filling results with keys:
        - status: Status ("success" or "error")
        - output_file: Final output file path
        - summary: Detailed statistics
        - message: Status message (success or error information)
    """
    try:
        if not os.path.exists(query_fasta):
            return {"status": "error", "message": f"Input file not found: {query_fasta}"}

        if not contigs_files and not hifi_files and not ont_files and not clr_files:
            return {"status": "error", "message": "Please provide at least one of: contigs files, HiFi data, ONT data, or CLR data"}

        if contigs_files:
            for contig_file in contigs_files:
                if not os.path.exists(contig_file):
                    return {"status": "error", "message": f"Contigs file not found: {contig_file}"}

        controller = AssemblyGapFiller()
        controller.verbose = verbose
        
        # Run complete pipeline
        results = controller.run_complete_pipeline(
            query_fasta=query_fasta,
            contigs_files=contigs_files,
            hifi_files=hifi_files,
            ont_files=ont_files,
            clr_files=clr_files,
            threads=threads,
            output_prefix=output_prefix,
            verbose=verbose,
            assembly_tools=assembly_tools,
            use_tgs=use_tgs
        )
        
        # Standardize output format
        if results.get('status') == 'success':
            return {
                "status": "success",
                "output_file": results.get('output_file', ''),
                "summary": results,
                "message": "Genome assembly and gap filling pipeline completed successfully"
            }
        else:
            return {
                "status": "error",
                "output_file": "",
                "summary": results,
                "message": results.get('message', 'Unknown error')
            }
            
    except Exception as e:
        return {
            "status": "error",
            "output_file": "",
            "summary": {},
            "message": f"Pipeline execution failed: {str(e)}"
        }


def assembly_gap_fill_from_dict(config: Dict) -> Dict:
    """
    Run genome assembly and gap filling pipeline from configuration dictionary
    
    Parameters:
        config: Configuration dictionary with keys:
            - query_fasta: Genome FASTA file that needs gap filling
            - contigs_files: List of reference contigs files (optional)
            - hifi_files: List of HiFi sequencing data files (optional)
            - ont_files: List of ONT sequencing data files (optional)
            - clr_files: List of CLR sequencing data files (optional)
            - threads: Number of threads to use (default: 16)
            - output_prefix: Output file prefix (default: filled_genome)
            - verbose: Show detailed output (default: False)
            - assembly_tools: List of assembly tools (optional)
            - use_tgs: Whether to enable TGS-GapCloser (default: False)
    
    Returns:
        Dictionary containing filling results
    """
    return assembly_gap_fill(
        query_fasta=config.get('query_fasta', ''),
        contigs_files=config.get('contigs_files'),
        hifi_files=config.get('hifi_files'),
        ont_files=config.get('ont_files'),
        clr_files=config.get('clr_files'),
        threads=config.get('threads', 16),
        output_prefix=config.get('output_prefix', 'filled_genome'),
        verbose=config.get('verbose', False),
        assembly_tools=config.get('assembly_tools'),
        use_tgs=config.get('use_tgs', False)
    )


# ============================================================================
# Main function
# ============================================================================

def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description="Integrated Genome Assembly and Gap Filling Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Usage examples:
  # Basic usage: Fill gaps using reference contigs (TGS OFF by default)
  python assembly_gap_fill.py -q draft_genome.fa -c reference_contigs.fa
  
  # Enable TGS-GapCloser for pre-filling
  python assembly_gap_fill.py -q draft_genome.fa -c reference_contigs.fa --use-tgs
  
  # Use HiFi data, enable TGS mode
  python assembly_gap_fill.py -q draft_genome.fa --hifi reads.hifi.fastq.gz --use-tgs
  
  # Use ONT data, only run flye assembly, enable TGS mode
  python assembly_gap_fill.py -q draft_genome.fa --ont reads.ont.fastq.gz --run-flye --use-tgs
  
  # Run multiple assembly tools, TGS OFF (default)
  python assembly_gap_fill.py -q draft.fa --hifi hifi.fastq.gz --run-hifiasm --run-flye
  
  # Use default configuration (all available tools), enable TGS
  python assembly_gap_fill.py -q draft.fa --hifi hifi.fastq.gz --use-tgs
  
  # Provide both contigs and sequencing data, TGS OFF
  python assembly_gap_fill.py -q draft_genome.fa -c ref.fa --hifi hifi.fastq.gz
  
  # Use multiple threads, enable TGS
  python assembly_gap_fill.py -q draft.fa -c ref.fa -t 32 --use-tgs
  
  # Specify output directory and prefix, TGS OFF
  python assembly_gap_fill.py -q draft.fa -c ref.fa -o results/ -p my_genome
  
  # Show verbose output, enable TGS
  python assembly_gap_fill.py -q draft.fa -c ref.fa -v --use-tgs

Process description:
  1. Initial genome analysis
  2. Gap detection and chromosome sorting
  3. Contigs preparation (assemble if sequencing data provided, otherwise use provided contigs)
  4. Optimized filling process (first round → filter <100kb fragments → second round)
  5. Gap filling mode: TGS OFF by default, can be enabled with --use-tgs
  6. Final genome creation
  7. Results analysis and report generation
  
Default parameters:
  - TGS mode: OFF
  - Minimum gap to fill: 1bp
  - Maximum gap to fill: 1Mb
  - Fragment filtering threshold: 100kb
  - Assembly tools: All available (if not specified)
        """
    )
    
    # Required parameters
    parser.add_argument('-q', '--query', required=True,
                       help='Genome FASTA file that needs gap filling')
    
    # Optional parameters
    parser.add_argument('-c', '--contigs', action='append', default=[],
                       help='Reference contigs files (FASTA format), can be used multiple times')
    parser.add_argument('--hifi', '--pacbio-hifi', action='append', default=[],
                       help='HiFi sequencing data files (FASTQ/GZ format)')
    parser.add_argument('--ont', '--ont-ul', action='append', default=[],
                       help='ONT sequencing data files (FASTQ/GZ format)')
    parser.add_argument('--clr', '--pacbio-clr', action='append', default=[],
                       help='CLR sequencing data files (FASTQ/GZ format)')
    
    # Assembly tool selection parameters
    parser.add_argument('--run-hifiasm', action='store_true',
                       help='Run hifiasm assembly tool')
    parser.add_argument('--run-verkko', action='store_true',
                       help='Run verkko assembly tool')
    parser.add_argument('--run-nextdenovo', action='store_true',
                       help='Run nextDenovo assembly tool')
    parser.add_argument('--run-flye', action='store_true',
                       help='Run flye assembly tool')
    parser.add_argument('--run-shasta', action='store_true',
                       help='Run shasta assembly tool')
    
    # TGS control parameters
    parser.add_argument('--use-tgs', action='store_true',
                       help='Enable TGS-GapCloser for pre-filling (default: OFF)')
    
    # General parameters
    parser.add_argument('-t', '--threads', type=int, default=16,
                       help='Number of threads (default: 16)')
    parser.add_argument('-o', '--output-dir', default='.',
                       help='Output directory (default: current directory)')
    parser.add_argument('-p', '--output-prefix', default='filled_genome',
                       help='Output file prefix (default: filled_genome)')
    
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed output information')
    parser.add_argument('--version', action='version', version='Integrated Genome Assembly and Gap Filling Pipeline v1.2')
    
    # If no arguments, show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.query):
        print(f"[ERROR] Input file not found: {args.query}")
        sys.exit(1)
    
    if args.contigs:
        for contig_file in args.contigs:
            if not os.path.exists(contig_file):
                print(f"[ERROR] Contigs file not found: {contig_file}")
                sys.exit(1)
    
    # Check if at least contigs or sequencing data is provided
    if not args.contigs and not args.hifi and not args.ont and not args.clr:
        print("[ERROR] Please provide at least one of the following:")
        print("        -c reference contigs files")
        print("        --hifi HiFi sequencing data")
        print("        --ont ONT sequencing data")
        print("        --clr CLR sequencing data")
        sys.exit(1)
    
    # Build complete output prefix
    if args.output_dir != '.':
        full_output_prefix = os.path.join(args.output_dir, args.output_prefix)
    else:
        full_output_prefix = args.output_prefix
    
    if args.output_dir != '.':
        summary_file = os.path.join(args.output_dir, f"{args.output_prefix}_summary.json")
    else:
        summary_file = f"{args.output_prefix}_summary.json"
    
    # Determine which assembly tools to use
    assembly_tools_list = []
    
    # If user specified tools, use those
    if args.run_hifiasm:
        assembly_tools_list.append('hifiasm')
    if args.run_verkko:
        assembly_tools_list.append('verkko')
    if args.run_nextdenovo:
        assembly_tools_list.append('nextdenovo')
    if args.run_flye:
        assembly_tools_list.append('flye')
    if args.run_shasta:
        assembly_tools_list.append('shasta')
    
    # If no tools specified, use default (all)
    if not assembly_tools_list:
        assembly_tools_list = None  # Use default configuration
        print("[INFO] No assembly tools specified, will use all available tools")
    else:
        print(f"[INFO] Specified assembly tools: {', '.join(assembly_tools_list)}")
 
    try:
        # Use standardized interface
        results = assembly_gap_fill(
            query_fasta=args.query,
            contigs_files=args.contigs,
            hifi_files=args.hifi,
            ont_files=args.ont,
            clr_files=args.clr,
            threads=args.threads,
            output_prefix=full_output_prefix,
            verbose=args.verbose,
            assembly_tools=assembly_tools_list,
            use_tgs=args.use_tgs  # Pass command line parameter
        )
        
        # Show brief results
        if results.get('status') == 'success':
            print(f"\n{'='*70}")
            print("Genome assembly and gap filling pipeline completed successfully!")
            print(f"{'='*70}")
            print(f"Input file: {args.query}")
            print(f"Output file: {results.get('output_file', 'Unknown')}")
            print(f"Threads: {args.threads}")
            print(f"TGS mode: {'ON' if args.use_tgs else 'OFF'}")  # Show TGS status
            
            summary = results.get('summary', {})
            print(f"Assembly tools: {', '.join(summary.get('assembly_tools', ['All available tools']))}")
            
            if 'final_summary' in summary:
                final_summary = summary['final_summary']
                print(f"\nChromosome statistics:")
                print(f"  Total: {final_summary.get('total_chromosomes', 'Unknown')}")
                print(f"  Total length: {final_summary.get('total_bases', 0):,} bp")
                print(f"  Original without gaps: {final_summary.get('original_no_gap', 0)}")
                print(f"  Filled without gaps: {final_summary.get('filled_no_gap', 0)}")
                print(f"  Still with gaps: {final_summary.get('still_with_gap', 0)}")
            
            # Show filling stage results
            if 'filling_stages' in summary:
                print(f"\nFilling stage results:")
                for stage in summary['filling_stages']:
                    stage_name = stage.get('stage', 'unknown')
                    gaps = stage.get('gaps_remaining', 0)
                    if stage_name == 'initial_analysis':
                        print(f"  Initial gap count: {stage.get('total_gaps', 0):,}")
                    elif stage_name == 'round1_filling':
                        print(f"  Gaps remaining after first round: {gaps:,}")
                    elif stage_name == 'post_fill_filtering':
                        print(f"  Gaps remaining after filtering: {gaps:,}")
                    elif stage_name == 'round2_filling':
                        print(f"  Gaps remaining after second round: {gaps:,}")
            
            print(f"\nDetailed report: {summary_file}")
            
        else:
            print(f"\nPipeline failed: {results.get('message', 'Unknown error')}")
            sys.exit(1)
        
    except KeyboardInterrupt:
        print("\n\nOperation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] Pipeline execution failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()