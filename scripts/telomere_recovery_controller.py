#!/usr/bin/env python3
"""
Telomere Recovery Controller Script
Purpose: Unified controller for three telomere processing sub-modules to recover chromosome telomeres
Process:
1. Use coord_te.py to analyze chromosome telomere status and generate alignments
2. Use extract_te.py to extract repair sequences from alignment results
3. Use add_te.py to merge repair sequences into chromosomes
Adapter pattern support added
Updated for enhanced extract_te.py and add_te.py
Optimized file organization and cleanup logic
Fixed path handling issues
"""

import sys
import os
import argparse
import textwrap
import json
import time
import shutil
import glob
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

# Import functions from three sub-modules
# Note: Ensure all three scripts are in the same directory or Python path
try:
    from coord_te import run_alignment_pipeline
    from extract_te import extract_telomere_regions
    from add_te import merge_telomere_sequences
    
    MODULES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Cannot import sub-modules - {e}")
    print("Telomere recovery functionality may be limited")
    MODULES_AVAILABLE = False

def print_header(step: int, title: str):
    """Print step header"""
    print(f"\n{'='*80}")
    print(f"Step {step}: {title}")
    print('='*80)

def cleanup_temp_files(output_dir: str = ".", cleanup_all: bool = True):
    """
    Clean up temporary files
    
    Parameters:
    ----------
    output_dir : str
        Output directory
    cleanup_all : bool
        Whether to clean all temporary files (including intermediate results)
    """
    print(f"\nCleaning temporary files (directory: {output_dir})...")
    
    # Clean possible temporary directories
    temp_dirs = [
        "temp_contigs",
        "alignment_results",
        "extraction_results",
        "merge_results",
        "temp",
        "tmp"
    ]
    
    for dir_name in temp_dirs:
        dir_path = os.path.join(output_dir, dir_name)
        if os.path.exists(dir_path):
            try:
                shutil.rmtree(dir_path)
                print(f"  Cleaned directory: {dir_path}")
            except Exception as e:
                print(f"  Warning: Cannot clean directory {dir_path}: {e}")
    
    # Clean possible temporary files
    temp_files = [
        "telomere_contigs.fa",
        "telomere_analysis.json",
        "query_data.json",
        "alignment_summary.json",
        "alignment_report.txt",
        "all_successful_matches.coords",
        "all_extracted_regions.fa",
        "extraction_report.txt",
        "repair_summary.fasta",
        "repair_report.txt",
        "repaired_genome.fasta",
        "telomere_matches.coords",
        "telomere_contigs.fasta"
    ]
    
    for file_name in temp_files:
        file_path = os.path.join(output_dir, file_name)
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"  Cleaned file: {file_path}")
            except Exception as e:
                print(f"  Warning: Cannot clean file {file_path}: {e}")
    
    # Clean pattern-matched temporary files (like ptg*.fa, Chr*.fa, etc.)
    if cleanup_all:
        patterns_to_clean = [
            "ptg*.fa",
            "Chr*.fa",
            "*_5prime.fa",
            "*_3prime.fa",
            "*_*.fa",  # Matches pattern_id_chromosome_start_end.fa format
            "temp_*.fa",
            "tmp_*.fa",
            "query_*.fa",
            "contig_*.fa"
        ]
        
        for pattern in patterns_to_clean:
            for file_path in glob.glob(os.path.join(output_dir, pattern)):
                try:
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                        print(f"  Cleaned pattern file: {os.path.basename(file_path)}")
                except Exception as e:
                    print(f"  Warning: Cannot clean file {file_path}: {e}")

def find_coords_file_in_output(output_dir: str) -> Optional[str]:
    """
    Find coords file in output directory
    
    Returns the first found coords file path, or None
    """
    coords_patterns = [
        "all_successful_matches.coords",
        "*.coords",
        "*matches*.coords",
        "telomere_matches.coords"
    ]
    
    for pattern in coords_patterns:
        for file_path in glob.glob(os.path.join(output_dir, pattern)):
            if os.path.isfile(file_path):
                return file_path
    
    # Recursively search subdirectories
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith('.coords') and 'match' in file.lower():
                return os.path.join(root, file)
    
    return None

def find_telomere_contigs_in_output(output_dir: str) -> Optional[str]:
    """
    Find telomere contigs file in output directory
    
    Returns the first found telomere contigs file path, or None
    """
    contigs_patterns = [
        "telomere_contigs.fa",
        "telomere_contigs.fasta",
        "*telomere*.fa",
        "*telomere*.fasta"
    ]
    
    for pattern in contigs_patterns:
        for file_path in glob.glob(os.path.join(output_dir, pattern)):
            if os.path.isfile(file_path):
                return file_path
    
    # Recursively search subdirectories
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if 'telomere' in file.lower() and (file.endswith('.fa') or file.endswith('.fasta')):
                return os.path.join(root, file)
    
    return None

def run_telomere_recovery(
    genome_file: str,
    contigs_file: str,
    output_dir: str = "telomere_recovery_output",
    output_file: str = "recovered_genome.fasta",
    threads: int = 32,
    verbose: bool = False,
    keep_intermediate: bool = False,  # New: whether to keep intermediate files
    # extract_te.py new parameters
    max_total_length: int = 5000000,
    max_alignment_keep: int = 10000,
    # add_te.py new parameters
    minimap2_mode: str = 'asm5',
    min_match_length: int = 100,
    min_mapq: int = 20,
    selection_strategy: str = 'minimal_extension'
) -> Dict[str, Any]:
    """
    Main telomere recovery function
    
    Parameters:
    ----------
    genome_file : str
        Genome FASTA file path
    contigs_file : str
        Contigs file for gap filling
    output_dir : str, optional
        Output directory (default: "telomere_recovery_output")
    output_file : str, optional
        Output filename (default: "recovered_genome.fasta")
    threads : int, optional
        Number of threads (default: 32)
    verbose : bool, optional
        Whether to show detailed output (default: False)
    keep_intermediate : bool, optional
        Whether to keep intermediate files (default: False)
    max_total_length : int, optional
        Maximum extraction total length(bp) (default: 5000000)
    max_alignment_keep : int, optional
        Maximum alignment keep length(bp) (default: 10000)
    minimap2_mode : str, optional
        minimap2 alignment mode (default: 'asm5')
    min_match_length : int, optional
        Minimum match length(bp) (default: 100)
    min_mapq : int, optional
        Minimum mapping quality (default: 20)
    selection_strategy : str, optional
        Selection strategy: 'first_success', 'minimal_extension', 'balanced' (default: 'minimal_extension')
    
    Returns:
    -------
    Dict[str, Any]
        Results containing the entire process
    """
    
    start_time = time.time()
    all_results = {
        'success': False,
        'output_dir': output_dir,
        'output_file': os.path.join(output_dir, output_file),
        'steps': {},
        'statistics': {},
        'chromosomes_repaired': 0,
        'total_extension': 0,
        'processing_time': 0
    }
    
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        print("Telomere Recovery Controller")
        print("="*80)
        print(f"Input genome: {genome_file}")
        print(f"Input contigs: {contigs_file}")
        print(f"Output directory: {output_dir}")
        print(f"Output file: {output_file}")
        print(f"Threads: {threads}")
        print(f"Keep intermediate files: {'Yes' if keep_intermediate else 'No'}")
        print(f"Extraction length limit: ≤{max_total_length:,}bp (alignment ≤{max_alignment_keep:,}bp)")
        print(f"Merge strategy: {selection_strategy} (minimap2 mode: {minimap2_mode})")
        print("="*80)
        
        # Check module availability
        if not MODULES_AVAILABLE:
            error_msg = "Telomere recovery modules unavailable, ensure coord_te.py, extract_te.py, add_te.py are in Python path"
            print(f"Error: {error_msg}")
            all_results['error'] = error_msg
            return all_results
        
        # Step 1: Analyze chromosome telomere status and generate alignments
        print_header(1, "Analyze chromosome telomere status and generate alignments")
        
        # Clean up any old files
        if os.path.exists(output_dir):
            cleanup_temp_files(output_dir, cleanup_all=False)
        
        # Use absolute paths to ensure coord_te.py handles correctly
        abs_output_dir = os.path.abspath(output_dir)
        
        alignment_result = run_alignment_pipeline(
            contigs=contigs_file,
            query=genome_file,
            output=abs_output_dir,  # Use absolute path
            threads=threads,
            min_repeats=5,            # Chromosome telomere detection parameter
            contig_min_repeats=20,    # Contig telomere detection parameter
            contig_min_length=500,    # Contig telomere minimum length
            skip_extract=False,       # Don't skip telomere extraction
            verbose=verbose
        )
        
        if not alignment_result['success']:
            print(f"Step 1 failed: {alignment_result.get('error', 'Unknown error')}")
            return {
                **all_results,
                'error': f"Step 1 failed: {alignment_result.get('error', 'Unknown error')}"
            }
        
        all_results['steps']['alignment'] = alignment_result
        
        # Check if any chromosomes need repair
        summary = alignment_result.get('summary', {})
        needs_repair = summary.get('needs_repair', 0)
        
        if needs_repair == 0:
            print("All chromosomes have complete telomeres, no repair needed")
            # Directly copy original genome to output file
            final_output = os.path.join(output_dir, output_file)
            shutil.copy2(genome_file, final_output)
            all_results.update({
                'success': True,
                'chromosomes_repaired': 0,
                'total_extension': 0,
                'processing_time': time.time() - start_time
            })
            
            # Clean temporary files (if not keeping intermediate files)
            if not keep_intermediate:
                cleanup_temp_files(output_dir, cleanup_all=True)
            
            return all_results
        
        print(f"Found {needs_repair} chromosomes needing telomere repair")
        
        # Get alignment result files
        coords_file = alignment_result.get('coords_file')
        telomere_contigs = alignment_result.get('telomere_contigs')
        
        # Debug information
        if verbose:
            print(f"Debug - File paths from alignment_result:")
            print(f"  coords_file: {coords_file}")
            print(f"  telomere_contigs: {telomere_contigs}")
        
        # If paths from alignment_result are invalid, try to find in output directory
        if not coords_file or not os.path.exists(coords_file):
            if verbose:
                print(f"Warning: coords file from result doesn't exist, searching in output directory...")
            coords_file = find_coords_file_in_output(output_dir)
        
        if not telomere_contigs or not os.path.exists(telomere_contigs):
            if verbose:
                print(f"Warning: telomere contigs file from result doesn't exist, searching in output directory...")
            telomere_contigs = find_telomere_contigs_in_output(output_dir)
        
        # Verify files exist
        if not coords_file or not os.path.exists(coords_file):
            print(f"Error: Alignment result file not found")
            print(f"Searched in directory {output_dir}, but no .coords file found")
            
            # List directory contents for debugging
            if verbose:
                print(f"Directory contents:")
                for item in os.listdir(output_dir):
                    item_path = os.path.join(output_dir, item)
                    if os.path.isfile(item_path):
                        print(f"  File: {item}")
                    else:
                        print(f"  Directory: {item}")
            
            return {
                **all_results,
                'error': f"No valid alignment result file generated, please check coord_te.py module output"
            }
        
        if not telomere_contigs or not os.path.exists(telomere_contigs):
            print(f"Error: Telomere contigs file not found")
            return {
                **all_results,
                'error': "No valid telomere contigs file found"
            }
        
        print(f"Alignment result file: {os.path.basename(coords_file)}")
        print(f"Telomere contigs file: {os.path.basename(telomere_contigs)}")
        
        # Step 2: Extract repair sequences from alignment results
        print_header(2, "Extract repair sequences from alignment results")
        
        extraction_result = extract_telomere_regions(
            coords_file=coords_file,
            contig_file=telomere_contigs,
            output_dir=output_dir,     # In output directory
            fragment_type="both",      # Extract both 5' and 3' ends
            min_repeats=20,            # Minimum telomere repeats
            min_telomere_length=500,   # Minimum telomere length
            similarity_threshold=95.0, # Alignment similarity threshold
            min_alignment_length=100,  # Minimum alignment length
            extend_before=0,           # Don't extend before alignment
            extend_after=0,            # Don't extend after alignment
            max_total_length=max_total_length,     # New: maximum total length
            max_alignment_keep=max_alignment_keep, # New: maximum alignment keep length
            verbose=verbose
        )
        
        if not extraction_result['success']:
            print(f"Step 2 failed: {extraction_result.get('error', 'Unknown error')}")
            return {
                **all_results,
                'error': f"Step 2 failed: {extraction_result.get('error', 'Unknown error')}"
            }
        
        all_results['steps']['extraction'] = extraction_result
        
        # Get extracted repair sequence file
        extracted_file = extraction_result.get('all_extracted_file')
        
        if not extracted_file or not os.path.exists(extracted_file):
            # Try to find in output directory
            extracted_file_candidate = os.path.join(output_dir, "all_extracted_regions.fa")
            if os.path.exists(extracted_file_candidate):
                extracted_file = extracted_file_candidate
            else:
                print(f"Error: Extracted repair sequence file not found")
                return {
                    **all_results,
                    'error': "No valid repair sequence file generated"
                }
        
        stats = extraction_result.get('statistics', {})
        total_regions = stats.get('total_regions', 0)
        
        if total_regions == 0:
            print("No valid repair sequences extracted")
            return {
                **all_results,
                'error': "No valid repair sequences extracted"
            }
        
        print(f"Successfully extracted {total_regions} repair sequence regions")
        print(f"Extracted file: {extracted_file}")
        
        # Step 3: Merge repair sequences into chromosomes
        print_header(3, "Merge repair sequences into chromosomes")
        
        # Build parameters for merge_telomere_sequences
        merge_params = {
            'genome_file': genome_file,
            'repair_file': extracted_file,
            'output_dir': output_dir,      # In output directory
            'search_range': 5000000,       # Search range
            'min_similarity': 0.6,         # Minimum similarity
            # New parameters
            'minimap2_mode': minimap2_mode,
            'min_match_length': min_match_length,
            'min_mapq': min_mapq,
            'selection_strategy': selection_strategy,
            'verbose': verbose
        }
        
        merge_result = merge_telomere_sequences(**merge_params)
        
        if not merge_result['success']:
            print(f"Step 3 failed: {merge_result.get('error', 'Unknown error')}")
            return {
                **all_results,
                'error': f"Step 3 failed: {merge_result.get('error', 'Unknown error')}"
            }
        
        all_results['steps']['merge'] = merge_result
        
        # Get repaired genome file
        repaired_genome = merge_result.get('repaired_genome')
        
        if not repaired_genome or not os.path.exists(repaired_genome):
            # Try to find in output directory
            repaired_genome_candidate = os.path.join(output_dir, "repaired_genome.fasta")
            if os.path.exists(repaired_genome_candidate):
                repaired_genome = repaired_genome_candidate
            else:
                print(f"Error: Repaired genome file not found")
                return {
                    **all_results,
                    'error': "No repaired genome file generated"
                }
        
        # Rename or copy to final output file
        final_output = os.path.join(output_dir, output_file)
        os.makedirs(os.path.dirname(final_output), exist_ok=True)
        if repaired_genome != final_output:
             shutil.copy2(repaired_genome, final_output)
             print(f"Repaired genome copied to: {final_output}")

        
        # Collect statistics
        merge_stats = merge_result.get('statistics', {})
        chromosomes_repaired = merge_stats.get('success_chromosomes', 0)
        total_extension = merge_stats.get('total_extension', 0)
        
        # Clean temporary files (if not keeping intermediate files)
        if not keep_intermediate:
            print(f"\nCleaning intermediate files...")
            # Only clean temporary files, keep final results
            cleanup_temp_files(output_dir, cleanup_all=True)
            
            # Ensure final output file is not deleted
            if os.path.exists(final_output):
                print(f"Preserving final output file: {final_output}")
            
            # Preserve important report files
            important_files = [
                os.path.join(output_dir, "repair_report.txt"),
                os.path.join(output_dir, "extraction_report.txt"),
                os.path.join(output_dir, "alignment_report.txt")
            ]
            
            for file_path in important_files:
                if os.path.exists(file_path):
                    print(f"Preserving report file: {os.path.basename(file_path)}")
        else:
            print(f"\nPreserving all intermediate files in directory: {output_dir}")
            
            # List intermediate files
            intermediate_patterns = [
                "ptg*.fa",
                "Chr*.fa",
                "*_5prime.fa",
                "*_3prime.fa",
                "*_*.fa",
                "telomere_contigs.fa",
                "all_successful_matches.coords",
                "all_extracted_regions.fa"
            ]
            
            intermediate_count = 0
            for pattern in intermediate_patterns:
                for file_path in glob.glob(os.path.join(output_dir, pattern)):
                    if os.path.isfile(file_path):
                        intermediate_count += 1
            
            if intermediate_count > 0:
                print(f"  Preserved {intermediate_count} intermediate files")
        
        # Ensure essential files exist
        essential_files_to_check = [
            ("Repaired genome", final_output),
            ("Repair report", os.path.join(output_dir, "repair_report.txt")),
            ("Extraction report", os.path.join(output_dir, "extraction_report.txt")),
            ("Alignment report", os.path.join(output_dir, "alignment_report.txt"))
        ]
        
        print(f"\nFinal output files:")
        for desc, file_path in essential_files_to_check:
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path)
                print(f"  {desc}: {os.path.basename(file_path)} ({file_size:,} bytes)")
            else:
                print(f"  {desc}: Not found")
        
        processing_time = time.time() - start_time
        
        all_results.update({
            'success': True,
            'chromosomes_repaired': chromosomes_repaired,
            'total_extension': total_extension,
            'processing_time': processing_time,
            'statistics': {
                'chromosomes_needing_repair': needs_repair,
                'chromosomes_repaired': chromosomes_repaired,
                'extraction_regions': total_regions,
                'total_extension_bp': total_extension,
                'selection_strategy': selection_strategy,
                'minimap2_mode': minimap2_mode,
                'max_total_length': max_total_length,
                'total_time_seconds': processing_time,
                'keep_intermediate': keep_intermediate,
                'output_directory': output_dir
            }
        })
        
        return all_results
        
    except Exception as e:
        import traceback
        error_trace = traceback.format_exc()
        print(f"Error during telomere recovery: {e}")
        if verbose:
            print(f"Detailed error information:\n{error_trace}")
        
        return {
            **all_results,
            'success': False,
            'error': str(e),
            'traceback': error_trace
        }


# ==================== ADAPTER CLASS ====================

class TelomereRecoveryAdapter:
    """
    Adapter class for telomere_recovery_controller.py
    Provides unified controller interface
    Updated for enhanced extract_te.py and add_te.py
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        Initialize adapter
        
        Args:
            config: Configuration dictionary
        """
        self.config = config or {}
        self.result = None
        self.status = "initialized"
        
        # Check module availability
        if not MODULES_AVAILABLE:
            self.modules_available = False
            print(f"Warning: Telomere recovery modules unavailable, adapter functionality limited")
        else:
            self.modules_available = True
    
    def validate_config(self, config: Dict[str, Any]) -> Tuple[bool, str]:
        """
        Validate configuration parameters
        
        Args:
            config: Configuration dictionary
            
        Returns:
            (is_valid, error_message)
        """
        # Check required parameters
        required_params = ['query_fasta', 'contigs_file']
        
        for param in required_params:
            if param not in config:
                return False, f"Missing required parameter: {param}"
            
            # Check file existence
            if param in ['query_fasta', 'contigs_file']:
                if not os.path.exists(config[param]):
                    return False, f"File does not exist: {config[param]}"
        
        # Check thread count
        if 'threads' in config and config['threads'] <= 0:
            return False, "Thread count must be greater than 0"
        
        # Check new parameters
        if 'max_total_length' in config and config['max_total_length'] <= 0:
            return False, "max_total_length must be greater than 0"
        
        if 'max_alignment_keep' in config and config['max_alignment_keep'] <= 0:
            return False, "max_alignment_keep must be greater than 0"
        
        if 'min_match_length' in config and config['min_match_length'] <= 0:
            return False, "min_match_length must be greater than 0"
        
        if 'min_mapq' in config and config['min_mapq'] < 0:
            return False, "min_mapq must be greater than or equal to 0"
        
        if 'selection_strategy' in config:
            valid_strategies = ['first_success', 'minimal_extension', 'balanced']
            if config['selection_strategy'] not in valid_strategies:
                return False, f"selection_strategy must be one of {', '.join(valid_strategies)}"
        
        if 'minimap2_mode' in config:
            valid_modes = ['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb']
            if config['minimap2_mode'] not in valid_modes:
                return False, f"minimap2_mode must be one of {', '.join(valid_modes)}"
        
        # Check module availability
        if not self.modules_available:
            return False, "Telomere recovery modules unavailable, ensure coord_te.py, extract_te.py, add_te.py are in Python path"
        
        return True, ""
    
    def run(self, config: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Run pipeline
        
        Args:
            config: Configuration dictionary, if None use initialization configuration
            
        Returns:
            Result dictionary
        """
        if config is not None:
            self.config = config
        
        # Validate configuration
        is_valid, error_msg = self.validate_config(self.config)
        if not is_valid:
            self.status = "config_error"
            return {
                "status": "error",
                "error": error_msg,
                "module": "telomere_recovery",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
        
        try:
            self.status = "running"
            
            # Prepare parameters (including new parameters)
            params = {
                'genome_file': self.config['query_fasta'],
                'contigs_file': self.config['contigs_file'],
                'output_dir': self.config.get('output_dir', 'telomere_recovery_output'),
                'output_file': self.config.get('output_file', 'recovered_genome.fasta'),
                'threads': self.config.get('threads', 32),
                'verbose': self.config.get('verbose', False),
                'keep_intermediate': self.config.get('keep_intermediate', False),
                # extract_te.py new parameters
                'max_total_length': self.config.get('max_total_length', 5000000),
                'max_alignment_keep': self.config.get('max_alignment_keep', 10000),
                # add_te.py new parameters
                'minimap2_mode': self.config.get('minimap2_mode', 'asm5'),
                'min_match_length': self.config.get('min_match_length', 100),
                'min_mapq': self.config.get('min_mapq', 20),
                'selection_strategy': self.config.get('selection_strategy', 'minimal_extension')
            }
            
            # Run telomere recovery process
            raw_result = run_telomere_recovery(**params)
            
            # Standardize result format
            self.result = self._standardize_result(raw_result)
            
            if self.result.get('status') == 'success':
                self.status = "completed"
            else:
                self.status = "failed"
            
            return self.result
            
        except Exception as e:
            self.status = "failed"
            error_result = {
                "status": "error",
                "error": str(e),
                "module": "telomere_recovery",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
            
            if hasattr(e, '__traceback__'):
                import traceback
                error_result["traceback"] = traceback.format_exc()
            
            self.result = error_result
            return error_result
    
    def _standardize_result(self, raw_result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Standardize result format
        
        Args:
            raw_result: Raw result
            
        Returns:
            Standardized result
        """
        if not raw_result.get('success', False):
            error_msg = raw_result.get('error', 'Unknown error')
            return {
                "status": "error",
                "error": error_msg,
                "module": "telomere_recovery",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
        
        # Success result
        output_dir = raw_result.get('output_dir', 'telomere_recovery_output')
        output_file = raw_result.get('output_file', os.path.join(output_dir, 'recovered_genome.fasta'))
        
        standardized = {
            "status": "success",
            "module": "telomere_recovery",
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "pipeline_version": "telomere_recovery_v2",
            "input_files": {
                "query_fasta": self.config.get('query_fasta', ''),
                "contigs_file": self.config.get('contigs_file', '')
            },
            "output_files": {
                "final_genome": output_file,
                "output_directory": output_dir,
                "repair_report": os.path.join(output_dir, "repair_report.txt"),
                "extraction_report": os.path.join(output_dir, "extraction_report.txt"),
                "alignment_report": os.path.join(output_dir, "alignment_report.txt")
            },
            "statistics": {
                "chromosomes_needing_repair": raw_result.get('statistics', {}).get('chromosomes_needing_repair', 0),
                "chromosomes_repaired": raw_result.get('chromosomes_repaired', 0),
                "extraction_regions": raw_result.get('statistics', {}).get('extraction_regions', 0),
                "total_extension_bp": raw_result.get('total_extension', 0),
                "selection_strategy": raw_result.get('statistics', {}).get('selection_strategy', 'minimal_extension'),
                "minimap2_mode": raw_result.get('statistics', {}).get('minimap2_mode', 'asm5'),
                "max_total_length": raw_result.get('statistics', {}).get('max_total_length', 5000000),
                "keep_intermediate": raw_result.get('statistics', {}).get('keep_intermediate', False),
                "processing_time": raw_result.get('processing_time', 0)
            },
            "processing_info": {
                "threads_used": self.config.get('threads', 32),
                "timestamp": time.strftime("%Y%m%d_%H%M%S"),
                "steps_completed": list(raw_result.get('steps', {}).keys()),
                "parameters_used": {
                    "max_total_length": self.config.get('max_total_length', 5000000),
                    "max_alignment_keep": self.config.get('max_alignment_keep', 10000),
                    "minimap2_mode": self.config.get('minimap2_mode', 'asm5'),
                    "min_match_length": self.config.get('min_match_length', 100),
                    "min_mapq": self.config.get('min_mapq', 20),
                    "selection_strategy": self.config.get('selection_strategy', 'minimal_extension'),
                    "keep_intermediate": self.config.get('keep_intermediate', False)
                }
            }
        }
        
        # Add repair details
        steps = raw_result.get('steps', {})
        if steps:
            standardized["step_details"] = {
                step: {
                    "success": step_result.get('success', False),
                    "summary": step_result.get('summary', {})
                }
                for step, step_result in steps.items()
            }
        
        # Add raw result (optional)
        if self.config.get('include_raw_result', False):
            standardized["raw_result"] = raw_result
        
        return standardized
    
    def get_status(self) -> str:
        """
        Get current status
        
        Returns:
            Status string
        """
        return self.status
    
    def get_result(self) -> Dict[str, Any]:
        """
        Get result
        
        Returns:
            Result dictionary
        """
        return self.result
    
    def cleanup(self, output_dir: str = None) -> bool:
        """
        Clean up temporary files
        
        Args:
            output_dir: Output directory, if None use configuration directory
            
        Returns:
            Whether successful
        """
        try:
            if output_dir is None:
                output_dir = self.config.get('output_dir', 'telomere_recovery_output')
            
            # Only clean temporary files, preserve final results
            cleanup_temp_files(output_dir, cleanup_all=True)
            
            # Ensure final output file is not deleted
            final_output = os.path.join(output_dir, 
                                      self.config.get('output_file', 'recovered_genome.fasta'))
            if os.path.exists(final_output):
                print(f"Preserving final output file: {final_output}")
            
            return True
        except Exception as e:
            print(f"Cleanup failed: {e}")
            return False
    
    def get_default_config(self) -> Dict[str, Any]:
        """
        Get default configuration
        
        Returns:
            Default configuration dictionary
        """
        return {
            "query_fasta": "",  # Required
            "contigs_file": "",  # Required
            "output_dir": "telomere_recovery_output",
            "output_file": "recovered_genome.fasta",
            "threads": 32,
            "verbose": False,
            "keep_intermediate": False,
            "include_raw_result": False,
            "cleanup_temp_files": True,
            # extract_te.py new parameters
            "max_total_length": 5000000,
            "max_alignment_keep": 10000,
            # add_te.py new parameters
            "minimap2_mode": "asm5",
            "min_match_length": 100,
            "min_mapq": 20,
            "selection_strategy": "minimal_extension"
        }
    
    def get_module_info(self) -> Dict[str, Any]:
        """
        Get module information
        
        Returns:
            Module information dictionary
        """
        return {
            "name": "TelomereRecoveryAdapter",
            "version": "2.0",
            "description": "Telomere recovery controller adapter - Recover and repair chromosome telomeres (adapted for enhanced sub-modules)",
            "capabilities": [
                "Chromosome telomere status analysis",
                "Telomere contig alignment",
                "Intelligent length-limited repair sequence extraction",
                "minimap2 single-sided alignment merging",
                "Multiple selection strategies (first_success/minimal_extension/balanced)",
                "Telomere integrity verification",
                "Organized file output",
                "Optional intermediate file preservation",
                "Automatic file finding and path repair"
            ],
            "dependencies": ["coord_te", "extract_te", "add_te"],
            "modules_available": self.modules_available,
            "parameters": {
                "required": ["query_fasta", "contigs_file"],
                "optional": {
                    "output_dir": "Output directory path",
                    "output_file": "Output filename",
                    "threads": "Number of threads",
                    "verbose": "Verbose output mode",
                    "keep_intermediate": "Whether to keep intermediate files",
                    # extract_te.py parameters
                    "max_total_length": "Maximum extraction total length(bp)",
                    "max_alignment_keep": "Maximum alignment keep length(bp)",
                    # add_te.py parameters
                    "minimap2_mode": "minimap2 alignment mode (asm5/asm10/asm20/map-ont/map-pb)",
                    "min_match_length": "Minimum match length(bp)",
                    "min_mapq": "Minimum mapping quality(0-60)",
                    "selection_strategy": "Selection strategy (first_success/minimal_extension/balanced)"
                }
            }
        }
    
    def check_dependencies(self) -> Dict[str, bool]:
        """
        Check dependency modules
        
        Returns:
            Dependency module status dictionary
        """
        dependencies = {
            "coord_te": False,
            "extract_te": False,
            "add_te": False
        }
        
        try:
            import importlib
            for module in dependencies.keys():
                try:
                    importlib.import_module(module)
                    dependencies[module] = True
                except ImportError:
                    pass
        
        except Exception:
            pass
        
        return dependencies


def get_controller_api(config: Dict[str, Any] = None):
    """
    Get controller API instance (factory function)
    
    Args:
        config: Configuration dictionary
        
    Returns:
        TelomereRecoveryAdapter instance
    """
    return TelomereRecoveryAdapter(config)


def main():
    """Command line main function"""
    parser = argparse.ArgumentParser(
        description="Telomere Recovery Controller Script - Unified controller for three telomere processing sub-modules (adapted for enhanced version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Complete process:
          1. Analyze chromosome telomere status, identify chromosomes lacking telomeres
          2. Use contigs alignment and extract repair sequences (supports intelligent length limits)
          3. Merge repair sequences into corresponding chromosome ends (supports multiple selection strategies)
          
        New features:
          - extract_te.py: Intelligent length limits, prevent extracting overly long sequences
          - add_te.py: minimap2 single-sided alignment, supports multiple selection strategies
          - Organized file output: All results saved in specified directory
          - Optional intermediate file preservation
          - Automatic file finding and path repair
          
        Usage examples:
          # Basic usage (using default output directory)
          python telomere_recovery_controller.py -q genome.fasta -c contigs.fasta
          
          # Specify output directory and filename
          python telomere_recovery_controller.py -q genome.fasta -c contigs.fasta \
            -o recovered_genome.fasta -d recovery_results
          
          # Keep intermediate files (for debugging)
          python telomere_recovery_controller.py -q genome.fasta -c contigs.fasta \
            --keep-intermediate -v
          
          # Custom extraction length limits
          python telomere_recovery_controller.py -q genome.fasta -c contigs.fasta \
            --max-total-length 1000000 --max-alignment-keep 5000
          
          # Custom merge strategy
          python telomere_recovery_controller.py -q genome.fasta -c contigs.fasta \
            --strategy balanced --minimap2-mode asm10 --min-match-length 200
          
          # Adapter mode
          python telomere_recovery_controller.py --adapter-mode --config config.json
          
          # Show module information
          python telomere_recovery_controller.py --module-info
          
          # Check dependencies
          python telomere_recovery_controller.py --check-deps
        """)
    )
    
    parser.add_argument("-q", "--query", required=False, 
                       help="Genome FASTA file")
    parser.add_argument("-c", "--contigs", required=False,
                       help="Contigs file for gap filling")
    parser.add_argument("-d", "--output-dir", default="telomere_recovery_output",
                       help="Output directory (default: telomere_recovery_output)")
    parser.add_argument("-o", "--output", default="recovered_genome.fasta",
                       help="Output filename (default: recovered_genome.fasta)")
    parser.add_argument("-t", "--threads", type=int, default=32,
                       help="Number of threads (default: 32)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Show detailed output")
    parser.add_argument("--keep-intermediate", action="store_true",
                       help="Keep intermediate files (for debugging)")
    
    # extract_te.py new parameters
    parser.add_argument("--max-total-length", type=int, default=5000000,
                       help="Maximum extraction total length(bp) (default: 5,000,000)")
    parser.add_argument("--max-alignment-keep", type=int, default=10000,
                       help="Maximum alignment keep length(bp) (default: 10000)")
    
    # add_te.py new parameters
    parser.add_argument("--minimap2-mode", default='asm5',
                       choices=['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb'],
                       help="minimap2 alignment mode (default: asm5)")
    parser.add_argument("--min-match-length", type=int, default=100,
                       help="Minimum match length(bp) (default: 100)")
    parser.add_argument("--min-mapq", type=int, default=20,
                       help="Minimum mapping quality (default: 20)")
    parser.add_argument("--strategy", default='minimal_extension',
                       choices=['first_success', 'minimal_extension', 'balanced'],
                       help="Selection strategy (default: minimal_extension)")
    
    # Adapter mode parameters
    parser.add_argument("--adapter-mode", action="store_true",
                       help="Use adapter mode")
    parser.add_argument("--config", type=str,
                       help="JSON configuration file path (adapter mode use)")
    parser.add_argument("--module-info", action="store_true",
                       help="Show module information")
    parser.add_argument("--check-deps", action="store_true",
                       help="Check dependency modules")
    parser.add_argument("--skip-cleanup", action="store_true",
                       help="Skip temporary file cleanup (traditional mode only)")
    
    args = parser.parse_args()
    
    # Show module information
    if args.module_info:
        adapter = TelomereRecoveryAdapter()
        info = adapter.get_module_info()
        print(json.dumps(info, indent=2, ensure_ascii=False))
        sys.exit(0)
    
    # Check dependencies
    if args.check_deps:
        adapter = TelomereRecoveryAdapter()
        deps = adapter.check_dependencies()
        print("Dependency module check:")
        for module, available in deps.items():
            status = "✓ Available" if available else "✗ Not available"
            print(f"  {module}: {status}")
        sys.exit(0)
    
    # Adapter mode
    if args.adapter_mode:
        if args.config:
            # Load from configuration file
            try:
                with open(args.config, 'r') as f:
                    config = json.load(f)
            except Exception as e:
                print(f"Cannot load configuration file: {e}")
                sys.exit(1)
        else:
            # Build configuration from command line arguments
            if not args.query or not args.contigs:
                print("Adapter mode requires providing -q/--query and -c/--contigs parameters")
                sys.exit(1)
                
            config = {
                'query_fasta': args.query,
                'contigs_file': args.contigs,
                'output_dir': args.output_dir,
                'output_file': args.output,
                'threads': args.threads,
                'verbose': args.verbose,
                'keep_intermediate': args.keep_intermediate,
                # New parameters
                'max_total_length': args.max_total_length,
                'max_alignment_keep': args.max_alignment_keep,
                'minimap2_mode': args.minimap2_mode,
                'min_match_length': args.min_match_length,
                'min_mapq': args.min_mapq,
                'selection_strategy': args.strategy
            }
        
        # Create adapter and run
        adapter = TelomereRecoveryAdapter(config)
        result = adapter.run()
        
        # Output result
        print(json.dumps(result, indent=2, ensure_ascii=False))
        
        if result.get('status') == 'success':
            sys.exit(0)
        else:
            sys.exit(1)
    
    # Traditional command line mode
    else:
        # Check input files
        if not args.query:
            print("Error: Please provide genome file (-q/--query)")
            sys.exit(1)
        
        if not os.path.exists(args.query):
            print(f"Error: Genome file does not exist: {args.query}")
            sys.exit(1)
        
        if not args.contigs:
            print("Error: Please provide contigs file (-c/--contigs)")
            sys.exit(1)
        
        if not os.path.exists(args.contigs):
            print(f"Error: Contigs file does not exist: {args.contigs}")
            sys.exit(1)
        
        # Run telomere recovery process
        result = run_telomere_recovery(
            genome_file=args.query,
            contigs_file=args.contigs,
            output_dir=args.output_dir,
            output_file=args.output,
            threads=args.threads,
            verbose=args.verbose,
            keep_intermediate=args.keep_intermediate,
            # New parameters
            max_total_length=args.max_total_length,
            max_alignment_keep=args.max_alignment_keep,
            minimap2_mode=args.minimap2_mode,
            min_match_length=args.min_match_length,
            min_mapq=args.min_mapq,
            selection_strategy=args.strategy
        )
        
        # Output result
        print("\n" + "="*80)
        print("Telomere recovery process completed!")
        print("="*80)
        
        if result['success']:
            stats = result['statistics']
            print(f"Processing results:")
            print(f"  Chromosomes needing repair: {stats.get('chromosomes_needing_repair', 0)}")
            print(f"  Chromosomes successfully repaired: {stats.get('chromosomes_repaired', 0)}")
            print(f"  Extraction regions: {stats.get('extraction_regions', 0)}")
            print(f"  Total extension length: {stats.get('total_extension_bp', 0):,} bp")
            print(f"  Selection strategy: {stats.get('selection_strategy', 'minimal_extension')}")
            print(f"  Output directory: {stats.get('output_directory', args.output_dir)}")
            print(f"  Keep intermediate files: {'Yes' if stats.get('keep_intermediate', False) else 'No'}")
            print(f"  Total processing time: {stats.get('total_time_seconds', 0):.2f} seconds")
            
            final_output = result.get('output_file', os.path.join(args.output_dir, args.output))
            print(f"\nFinal output:")
            print(f"  Repaired genome: {final_output}")
            print(f"  Repair report: {os.path.join(args.output_dir, 'repair_report.txt')}")
            print(f"  Extraction report: {os.path.join(args.output_dir, 'extraction_report.txt')}")
            print(f"  Alignment report: {os.path.join(args.output_dir, 'alignment_report.txt')}")
            
            # If keeping intermediate files, list some
            if args.keep_intermediate:
                import glob
                intermediate_files = glob.glob(os.path.join(args.output_dir, "*.fa"))
                if intermediate_files:
                    print(f"\nIntermediate files ({len(intermediate_files)}):")
                    for i, file_path in enumerate(intermediate_files[:5]):  # Show first 5
                        file_name = os.path.basename(file_path)
                        file_size = os.path.getsize(file_path)
                        print(f"  {file_name} ({file_size:,} bytes)")
                    if len(intermediate_files) > 5:
                        print(f"  ... {len(intermediate_files)-5} more files")
            
            # Clean temporary files (unless specified to skip and not keeping intermediate files)
            if not args.skip_cleanup and not args.keep_intermediate:
                print(f"\nCleaning temporary files...")
                cleanup_temp_files(args.output_dir, cleanup_all=True)
                print(f"Temporary files cleaned")
        else:
            print(f"Telomere recovery failed!")
            print(f"Error: {result.get('error', 'Unknown error')}")
            
            if args.verbose and 'traceback' in result:
                print(f"\nDetailed error information:")
                print(result['traceback'])
            
            sys.exit(1)
        
        print("="*80)


if __name__ == "__main__":
    main()
