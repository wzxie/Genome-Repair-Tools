#!/usr/bin/env python3
"""
Fast Parallel Nucmer Analyzer - Unified API Version
Chromosome folder organization, support resume from breakpoint, can be called as submodule
Support analyzing specified chromosomes only
"""

import sys
import os
import subprocess
import concurrent.futures
import shutil
import argparse
import multiprocessing
import json
import time
import threading
import traceback
import logging
from typing import Dict, List, Tuple, Optional, Any, Union
from pathlib import Path

class NucmerAnalyzerAPI:
    """Unified API for Nucmer analysis"""
    
    def __init__(self, verbose: bool = True, log_file: Optional[str] = None):
        """
        Initialize NucmerAnalyzer API
        
        Args:
            verbose: Show detailed output
            log_file: Log file path (optional)
        """
        self.verbose = verbose
        self.log_file = log_file
        self.temp_dirs = []
        self._setup_logging()
    
    def _setup_logging(self):
        """Configure logging system"""
        self.logger = logging.getLogger('NucmerAnalyzer')
        self.logger.setLevel(logging.INFO if self.verbose else logging.WARNING)
        
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO if self.verbose else logging.WARNING)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.INFO)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
    
    def log(self, message: str, level: str = "info"):
        """Unified logging"""
        if level == "info":
            self.logger.info(message)
        elif level == "warning":
            self.logger.warning(message)
        elif level == "error":
            self.logger.error(message)
        elif level == "debug":
            self.logger.debug(message)
    
    def run_analysis(
        self,
        reference_fasta: str,
        query_fasta: str,
        output_dir: str = ".",
        output_prefix: str = "nucmer_analysis",
        threads: int = 0,
        size_threshold: int = 1000000,
        chunk_size: int = 10000000,
        batch_size: int = 10,
        min_alignment_length: int = 10000,
        overwrite: bool = False,
        keep_temp: bool = False,
        chromosomes: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Run nucmer alignment analysis
        
        Args:
            reference_fasta: Reference genome file path
            query_fasta: Query genome file path
            output_dir: Output directory
            output_prefix: Output file prefix
            threads: Number of threads (0=auto-detect)
            size_threshold: Large contig threshold (bp)
            chunk_size: Small contigs combination size (bp)
            batch_size: Batch processing size
            min_alignment_length: Minimum alignment length (bp)
            overwrite: Overwrite existing results
            keep_temp: Keep temporary files
            chromosomes: List of chromosomes to analyze (None=all)
            
        Returns:
            Dictionary containing analysis results
        """
        self.log(f"Starting nucmer alignment analysis", "info")
        self.log(f"Reference genome: {reference_fasta}", "info")
        self.log(f"Query genome: {query_fasta}", "info")
        self.log(f"Output directory: {output_dir}", "info")
        
        if chromosomes:
            self.log(f"Analyzing specified chromosomes: {', '.join(chromosomes)}", "info")
        
        start_time = time.time()
        
        os.makedirs(output_dir, exist_ok=True)
        
        report_file = os.path.join(output_dir, f"{output_prefix}_report.json")
        if os.path.exists(report_file) and not overwrite:
            self.log(f"Previous analysis results detected, loading report: {report_file}", "info")
            try:
                with open(report_file, 'r') as f:
                    result = json.load(f)
                result["api_info"] = {
                    "execution_time": 0,
                    "note": "Loaded from existing report"
                }
                return result
            except Exception as e:
                self.log(f"Failed to load report, re-analyzing: {e}", "warning")
        
        self._check_dependencies()
        
        self._set_environment_limits()
        
        if threads == 0:
            threads = self._get_optimal_threads()
        
        self.log(f"Using threads: {threads}", "info")
        
        analyzer = _FastNucmerParallelAnalyzer(verbose=self.verbose, specified_chromosomes=chromosomes)
        
        try:
            result = analyzer.analyze(
                ref_fasta=reference_fasta,
                query_fasta=query_fasta,
                output_dir=output_dir,
                threads=threads,
                size_threshold=size_threshold,
                chunk_size=chunk_size,
                batch_size=batch_size,
                min_alignment_length=min_alignment_length
            )
            
            result["api_info"] = {
                "execution_time": time.time() - start_time,
                "input_files": {
                    "reference_fasta": reference_fasta,
                    "query_fasta": query_fasta
                },
                "parameters": {
                    "threads": threads,
                    "size_threshold": size_threshold,
                    "chunk_size": chunk_size,
                    "batch_size": batch_size,
                    "min_alignment_length": min_alignment_length,
                    "specified_chromosomes": chromosomes if chromosomes else "all"
                },
                "output_structure": result.get("output_structure", {}),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
            }
            
            if not keep_temp:
                self._cleanup_temp_files(result)
            
            api_report_file = os.path.join(output_dir, f"{output_prefix}_api_report.json")
            with open(api_report_file, 'w') as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
            
            self.log(f"API format report saved: {api_report_file}", "info")
            
            self._print_summary(result, chromosomes)
            
            return result
            
        except Exception as e:
            self.log(f"Nucmer analysis failed: {str(e)}", "error")
            raise
    
    def _check_dependencies(self):
        """Check required dependencies"""
        for cmd in ["nucmer", "delta-filter", "show-coords"]:
            try:
                subprocess.run([cmd, "--version"], 
                             stdout=subprocess.DEVNULL, 
                             stderr=subprocess.DEVNULL,
                             check=False)
            except FileNotFoundError:
                raise RuntimeError(f"{cmd} not installed or not in PATH. Install MUMmer: conda install -c bioconda mummer")
        
        try:
            from Bio import SeqIO
        except ImportError:
            raise RuntimeError("BioPython not installed. Install: conda install -c bioconda biopython")
    
    def _set_environment_limits(self):
        """Set environment variable limits for threads"""
        env_vars = {
            'OPENBLAS_NUM_THREADS': '1',
            'GOTO_NUM_THREADS': '1',
            'OMP_NUM_THREADS': '1',
            'MKL_NUM_THREADS': '1',
            'VECLIB_MAXIMUM_THREADS': '1',
            'NUMEXPR_NUM_THREADS': '1'
        }
        os.environ.update(env_vars)
    
    def _get_optimal_threads(self) -> int:
        """Get optimal number of threads"""
        try:
            return max(1, multiprocessing.cpu_count())
        except:
            return 8
    
    def _cleanup_temp_files(self, result: Dict[str, Any]):
        """Clean up temporary files"""
        temp_dirs_to_clean = []
        
        if "temp_directory" in result:
            temp_dirs_to_clean.append(result["temp_directory"])
        
        if "output_structure" in result:
            base_dir = result["output_structure"].get("base_dir", "")
            if base_dir:
                for item in os.listdir(base_dir):
                    item_path = os.path.join(base_dir, item)
                    if os.path.isdir(item_path) and item.startswith("tmp_"):
                        temp_dirs_to_clean.append(item_path)
        
        for temp_dir in temp_dirs_to_clean:
            if os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    self.log(f"Cleaned temporary directory: {temp_dir}", "info")
                except Exception as e:
                    self.log(f"Failed to clean temporary directory {temp_dir}: {e}", "warning")
    
    def _print_summary(self, result: Dict[str, Any], chromosomes: Optional[List[str]] = None):
        """Print analysis summary"""
        self.log("\n" + "="*60, "info")
        self.log("Nucmer alignment analysis completed!", "info")
        self.log("="*60, "info")
        
        if chromosomes:
            self.log(f"Analyzed specified chromosomes: {', '.join(chromosomes)}", "info")
        
        stats = result.get("statistics", {})
        self.log(f"Total runtime: {stats.get('total_time', 0):.1f} seconds", "info")
        self.log(f"Alignment time: {stats.get('alignment_time', 0):.1f} seconds", "info")
        self.log(f"Tasks completed: {stats.get('completed_tasks', 0)}/{stats.get('total_tasks', 0)}", "info")
        self.log(f"Success rate: {stats.get('success_rate', 0):.1f}%", "info")
        
        processed_chromosomes = result.get("chromosome_results", {})
        self.log(f"\nProcessed {len(processed_chromosomes)} chromosomes:", "info")
        for seq_id, chr_data in processed_chromosomes.items():
            syntenic_info = result.get("syntenic_results", {}).get(seq_id, {})
            contig_count = syntenic_info.get('contig_count', 0)
            self.log(f"  {seq_id}: {chr_data.get('seq_length', 0):,} bp, syntenic contigs: {contig_count}", "info")
        
        output_dir = result.get("output_structure", {}).get("base_dir", ".")
        self.log(f"\nOutput directory: {output_dir}", "info")
    
    def batch_analysis(
        self,
        reference_fasta_list: List[str],
        query_fasta_list: List[str],
        output_dir: str = ".",
        batch_prefix: str = "batch",
        **kwargs
    ) -> Dict[str, Any]:
        """
        Run multiple nucmer analyses in batch
        
        Args:
            reference_fasta_list: List of reference genome files
            query_fasta_list: List of query genome files
            output_dir: Output directory
            batch_prefix: Batch processing prefix
            **kwargs: Other parameters passed to run_analysis
            
        Returns:
            Batch processing summary
        """
        if len(reference_fasta_list) != len(query_fasta_list):
            raise ValueError("reference_fasta_list and query_fasta_list must have same length")
        
        self.log(f"Starting batch analysis for {len(reference_fasta_list)} alignments", "info")
        
        results = {}
        for i, (ref_fasta, query_fasta) in enumerate(zip(reference_fasta_list, query_fasta_list)):
            self.log(f"Processing alignment {i+1}/{len(reference_fasta_list)}: {os.path.basename(ref_fasta)} vs {os.path.basename(query_fasta)}", "info")
            
            output_subdir = os.path.join(output_dir, f"{batch_prefix}_{i+1:03d}")
            output_prefix = f"analysis_{i+1:03d}"
            
            try:
                result = self.run_analysis(
                    reference_fasta=ref_fasta,
                    query_fasta=query_fasta,
                    output_dir=output_subdir,
                    output_prefix=output_prefix,
                    **kwargs
                )
                results[f"{os.path.basename(ref_fasta)}_{os.path.basename(query_fasta)}"] = {
                    "success": True,
                    "result": result,
                    "output_dir": output_subdir,
                    "output_prefix": output_prefix
                }
            except Exception as e:
                results[f"{os.path.basename(ref_fasta)}_{os.path.basename(query_fasta)}"] = {
                    "success": False,
                    "error": str(e),
                    "output_dir": output_subdir,
                    "output_prefix": output_prefix
                }
                self.log(f"Analysis failed: {e}", "warning")
        
        summary = {
            "total": len(reference_fasta_list),
            "successful": sum(1 for r in results.values() if r["success"]),
            "failed": sum(1 for r in results.values() if not r["success"]),
            "results": results,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        summary_file = os.path.join(output_dir, f"{batch_prefix}_summary.json")
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        self.log(f"Batch analysis completed, successful: {summary['successful']}/{summary['total']}", "info")
        self.log(f"Summary report saved: {summary_file}", "info")
        
        return summary
    
    def extract_syntenic_contigs(
        self,
        analysis_result: Dict[str, Any],
        chromosome: Optional[str] = None,
        output_file: Optional[str] = None,
        min_identity: float = 95.0,
        min_coverage: float = 0.5
    ) -> Dict[str, Any]:
        """
        Extract syntenic contigs from nucmer analysis results
        
        Args:
            analysis_result: Nucmer analysis result dictionary
            chromosome: Specify chromosome (optional, None=all)
            output_file: Output file path (optional)
            min_identity: Minimum alignment identity (%)
            min_coverage: Minimum coverage ratio
            
        Returns:
            Extraction result dictionary
        """
        self.log(f"Starting syntenic contig extraction", "info")
        
        output_structure = analysis_result.get("output_structure", {})
        base_dir = output_structure.get("base_dir", ".")
        
        if not os.path.exists(base_dir):
            raise FileNotFoundError(f"Output directory not found: {base_dir}")
        
        chromosomes_to_process = []
        if chromosome:
            if chromosome in analysis_result.get("chromosome_results", {}):
                chromosomes_to_process = [chromosome]
            else:
                raise ValueError(f"Chromosome {chromosome} not found in analysis results")
        else:
            chromosomes_to_process = list(analysis_result.get("chromosome_results", {}).keys())
        
        self.log(f"Processing syntenic contigs for {len(chromosomes_to_process)} chromosomes", "info")
        
        extraction_results = {}
        all_extracted_contigs = []
        
        for chr_name in chromosomes_to_process:
            chr_data = analysis_result["chromosome_results"].get(chr_name, {})
            chr_dir = chr_data.get("dir", "")
            
            if not os.path.exists(chr_dir):
                self.log(f"Chromosome directory not found: {chr_dir}", "warning")
                continue
            
            coords_files = []
            for file in os.listdir(chr_dir):
                if file.endswith('_merged.coords'):
                    coords_files.append(os.path.join(chr_dir, file))
            
            if not coords_files:
                self.log(f"No coords file found: {chr_dir}", "warning")
                continue
            
            coords_file = coords_files[0]
            ref_fasta = analysis_result.get("parameters", {}).get("ref_fasta", "")
            
            if not os.path.exists(ref_fasta):
                self.log(f"Reference genome file not found: {ref_fasta}", "warning")
                continue
            
            contig_count, contig_file = _extract_syntenic_contigs_from_coords(
                coords_file=coords_file,
                ref_fasta=ref_fasta,
                seq_id=chr_name,
                output_dir=chr_dir,
                min_alignment_length=10000,
                min_identity=min_identity,
                min_coverage=min_coverage
            )
            
            extraction_results[chr_name] = {
                "contig_count": contig_count,
                "contig_file": contig_file,
                "coords_file": coords_file
            }
            
            if contig_file and os.path.exists(contig_file):
                try:
                    from Bio import SeqIO
                    for record in SeqIO.parse(contig_file, "fasta"):
                        all_extracted_contigs.append(record)
                except Exception as e:
                    self.log(f"Failed to read contigs file: {e}", "warning")
        
        if output_file and all_extracted_contigs:
            try:
                from Bio import SeqIO
                SeqIO.write(all_extracted_contigs, output_file, "fasta")
                self.log(f"Merged syntenic contigs saved: {output_file}", "info")
            except Exception as e:
                self.log(f"Failed to save merged contigs: {e}", "warning")
                output_file = None
        
        result = {
            "extraction_results": extraction_results,
            "total_contigs": len(all_extracted_contigs),
            "merged_file": output_file,
            "chromosomes_processed": len(chromosomes_to_process)
        }
        
        self.log(f"Total extracted {len(all_extracted_contigs)} syntenic contigs", "info")
        
        return result
    
    def cleanup(self):
        """Clean up all temporary files"""
        for temp_dir in self.temp_dirs:
            if os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    self.log(f"Cleaned temporary directory: {temp_dir}", "info")
                except Exception as e:
                    self.log(f"Failed to clean temporary directory {temp_dir}: {e}", "warning")
        self.temp_dirs.clear()
    
    def __del__(self):
        """Destructor, automatically clean temporary files"""
        self.cleanup()

class _FastNucmerParallelAnalyzer:
    """Internal fast parallel nucmer analyzer"""
    
    def __init__(self, verbose: bool = True, specified_chromosomes: Optional[List[str]] = None):
        self.verbose = verbose
        self.specified_chromosomes = specified_chromosomes
        
        self.start_time = time.time()
        self.completed_tasks = 0
        self.failed_tasks = 0
        self.total_tasks = 0
        self.last_progress_update = 0
        self.progress_interval = 1.0
        
        self.lock = threading.Lock()
    
    def check_dependencies(self):
        """Check required dependencies"""
        for cmd in ["nucmer", "delta-filter", "show-coords"]:
            try:
                subprocess.run([cmd, "--version"], 
                             stdout=subprocess.DEVNULL, 
                             stderr=subprocess.DEVNULL,
                             check=False)
            except FileNotFoundError:
                print(f"Error: {cmd} not installed or not in PATH")
                print("Install MUMmer: conda install -c bioconda mummer")
                sys.exit(1)
        
        try:
            from Bio import SeqIO
        except ImportError:
            print("Error: BioPython not installed")
            print("Install: conda install -c bioconda biopython")
            sys.exit(1)
    
    def set_environment_limits(self):
        """Set environment variable limits for threads"""
        env_vars = {
            'OPENBLAS_NUM_THREADS': '1',
            'GOTO_NUM_THREADS': '1',
            'OMP_NUM_THREADS': '1',
            'MKL_NUM_THREADS': '1',
            'VECLIB_MAXIMUM_THREADS': '1',
            'NUMEXPR_NUM_THREADS': '1'
        }
        os.environ.update(env_vars)
    
    def get_optimal_threads(self) -> int:
        """Get optimal number of threads"""
        try:
            return max(1, multiprocessing.cpu_count())
        except:
            return 8
    
    def update_progress(self, force: bool = False):
        """Update progress display"""
        current_time = time.time()
        
        if not force and current_time - self.last_progress_update < self.progress_interval:
            return
        
        with self.lock:
            if self.total_tasks > 0:
                progress = (self.completed_tasks + self.failed_tasks) / self.total_tasks * 100
                elapsed = time.time() - self.start_time
                
                if self.verbose:
                    sys.stdout.write(f"\rProgress: {progress:.1f}% ({self.completed_tasks + self.failed_tasks}/{self.total_tasks}) "
                                   f"[{elapsed:.0f}s]")
                    sys.stdout.flush()
            
            self.last_progress_update = current_time
    
    def create_chromosome_folders(self, query_fasta: str, output_base_dir: str) -> Dict:
        """Create independent folders for each chromosome"""
        if not os.path.exists(query_fasta):
            raise FileNotFoundError(f"Query file not found: {query_fasta}")
        
        query_data = {}
        os.makedirs(output_base_dir, exist_ok=True)
        
        if self.verbose:
            print(f"Reading query genome: {query_fasta}")
        
        try:
            from Bio import SeqIO
            for record in SeqIO.parse(query_fasta, "fasta"):
                seq_id = record.id
                
                if self.specified_chromosomes and seq_id not in self.specified_chromosomes:
                    continue
                
                safe_seq_id = self.sanitize_filename(seq_id)
                chr_dir = os.path.join(output_base_dir, safe_seq_id)
                os.makedirs(chr_dir, exist_ok=True)
                
                chr_fasta = os.path.join(chr_dir, f"{safe_seq_id}.fa")
                SeqIO.write([record], chr_fasta, "fasta")
                
                temp_dir = os.path.join(chr_dir, "temp")
                os.makedirs(temp_dir, exist_ok=True)
                
                query_data[safe_seq_id] = {
                    'original_id': seq_id,
                    'chr_dir': chr_dir,
                    'chr_fasta': chr_fasta,
                    'temp_dir': temp_dir,
                    'seq_length': len(record.seq)
                }
        except Exception as e:
            raise RuntimeError(f"Failed to read query file: {e}")
        
        if not query_data:
            if self.specified_chromosomes:
                if self.verbose:
                    print(f"Warning: Specified chromosomes not found, trying sanitized name matching...")
                for record in SeqIO.parse(query_fasta, "fasta"):
                    seq_id = record.id
                    safe_seq_id = self.sanitize_filename(seq_id)
                    
                    if safe_seq_id in self.specified_chromosomes:
                        chr_dir = os.path.join(output_base_dir, safe_seq_id)
                        os.makedirs(chr_dir, exist_ok=True)
                        
                        chr_fasta = os.path.join(chr_dir, f"{safe_seq_id}.fa")
                        SeqIO.write([record], chr_fasta, "fasta")
                        
                        temp_dir = os.path.join(chr_dir, "temp")
                        os.makedirs(temp_dir, exist_ok=True)
                        
                        query_data[safe_seq_id] = {
                            'original_id': seq_id,
                            'chr_dir': chr_dir,
                            'chr_fasta': chr_fasta,
                            'temp_dir': temp_dir,
                            'seq_length': len(record.seq)
                        }
        
        if not query_data:
            if self.specified_chromosomes:
                raise ValueError(f"Specified chromosomes not found in query file: {self.specified_chromosomes}")
            else:
                raise ValueError("No valid sequences found in query file")
        
        if self.verbose:
            if self.specified_chromosomes:
                print(f"  Created independent folders for {len(query_data)} specified chromosomes")
            else:
                print(f"  Created independent folders for {len(query_data)} chromosomes")
        
        return query_data
    
    def sanitize_filename(self, filename: str) -> str:
        """Sanitize filename to avoid path issues"""
        unsafe_chars = ['<', '>', ':', '"', '/', '\\', '|', '?', '*', ' ']
        for char in unsafe_chars:
            filename = filename.replace(char, '_')
        if len(filename) > 100:
            filename = filename[:50] + "..." + filename[-47:]
        return filename
    
    def load_checkpoint(self, checkpoint_file: str) -> Dict:
        """Load checkpoint data"""
        if os.path.exists(checkpoint_file):
            try:
                with open(checkpoint_file, 'r') as f:
                    return json.load(f)
            except:
                return {}
        return {}
    
    def save_checkpoint(self, checkpoint_file: str, data: Dict):
        """Save checkpoint data"""
        try:
            with open(checkpoint_file, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            if self.verbose:
                print(f"Warning: Failed to save checkpoint: {e}")
    
    def split_contigs_optimized(self, ref_fasta: str, output_base_dir: str, **kwargs):
        """Optimized contig splitting"""
        if not os.path.exists(ref_fasta):
            raise FileNotFoundError(f"Reference file not found: {ref_fasta}")
        
        if self.verbose:
            print(f"Reading reference genome: {ref_fasta}")
        
        from Bio import SeqIO
        
        split_dir = os.path.join(output_base_dir, "shared_split_contigs")
        checkpoint_file = os.path.join(split_dir, "split_checkpoint.json")
        
        checkpoint = self.load_checkpoint(checkpoint_file)
        
        if checkpoint.get('split_completed', False):
            if self.verbose:
                print("  Loading existing contig split results...")
            
            large_files = []
            small_files = []
            
            large_dir = os.path.join(split_dir, "large")
            if os.path.exists(large_dir):
                large_files = sorted([
                    os.path.join(large_dir, f) 
                    for f in os.listdir(large_dir) 
                    if f.endswith('.fa')
                ])
            
            small_dir = os.path.join(split_dir, "small")
            if os.path.exists(small_dir):
                small_files = sorted([
                    os.path.join(small_dir, f) 
                    for f in os.listdir(small_dir) 
                    if f.endswith('.fa')
                ])
            
            if self.verbose:
                print(f"  Resume from breakpoint: {len(large_files)} large contigs, {len(small_files)} small contig combinations")
            
            return large_files, small_files, split_dir
        
        os.makedirs(split_dir, exist_ok=True)
        
        try:
            records = list(SeqIO.parse(ref_fasta, "fasta"))
        except Exception as e:
            raise RuntimeError(f"Failed to read reference file: {e}")
        
        if not records:
            raise ValueError("No valid sequences found in reference file")
        
        size_threshold = kwargs.get('size_threshold', 1000000)
        chunk_size = kwargs.get('chunk_size', 10000000)
        
        large_records = [r for r in records if len(r.seq) >= size_threshold]
        small_records = [r for r in records if len(r.seq) < size_threshold]
        
        large_dir = os.path.join(split_dir, "large")
        small_dir = os.path.join(split_dir, "small")
        os.makedirs(large_dir, exist_ok=True)
        os.makedirs(small_dir, exist_ok=True)
        
        large_files = []
        if large_records:
            for i, record in enumerate(large_records):
                file_path = os.path.join(large_dir, f"large_{i+1}.fa")
                SeqIO.write([record], file_path, "fasta")
                large_files.append(file_path)
        
        small_files = []
        if small_records:
            small_records_sorted = sorted(small_records, key=lambda x: len(x.seq), reverse=True)
            
            chunks = []
            current_chunk = []
            current_size = 0
            
            for record in small_records_sorted:
                record_size = len(record.seq)
                if current_size + record_size <= chunk_size or not current_chunk:
                    current_chunk.append(record)
                    current_size += record_size
                else:
                    chunks.append(current_chunk)
                    current_chunk = [record]
                    current_size = record_size
            
            if current_chunk:
                chunks.append(current_chunk)
            
            for i, chunk in enumerate(chunks):
                file_path = os.path.join(small_dir, f"small_chunk_{i+1}.fa")
                SeqIO.write(chunk, file_path, "fasta")
                small_files.append(file_path)
        
        checkpoint = {
            'split_completed': True,
            'large_files_count': len(large_files),
            'small_files_count': len(small_files),
            'size_threshold': size_threshold,
            'chunk_size': chunk_size,
            'timestamp': time.time()
        }
        self.save_checkpoint(checkpoint_file, checkpoint)
        
        if self.verbose:
            print(f"  Split results: {len(large_files)} large contigs, {len(small_files)} small contig combinations")
        
        return large_files, small_files, split_dir
    
    def run_nucmer_task(self, ref_file: str, query_file: str, output_prefix: str, 
                        contig_size: int, seq_id: str, checkpoint_data: Dict) -> Optional[str]:
        """Run single nucmer task"""
        
        coords_file = f"{output_prefix}.coords"
        task_key = f"{os.path.basename(ref_file)}_{seq_id}"
        
        if checkpoint_data.get(task_key, {}).get('completed', False):
            if os.path.exists(coords_file):
                with self.lock:
                    self.completed_tasks += 1
                self.update_progress(force=False)
                return coords_file
        
        try:
            threads = self.calculate_threads_for_contig(contig_size)
            
            output_dir = os.path.dirname(output_prefix)
            os.makedirs(output_dir, exist_ok=True)
            
            nucmer_cmd = [
                "nucmer",
                "-c", "1000",
                "-l", "100",
                "--batch=500000000",
                "-t", str(threads),
                "-p", output_prefix,
                ref_file,
                query_file
            ]
            
            result = subprocess.run(
                nucmer_cmd,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                env=os.environ.copy()
            )
            
            delta_file = f"{output_prefix}.delta"
            if not os.path.exists(delta_file):
                return None
            
            final_coords_file = self.process_delta_file(delta_file, output_prefix)
            
            if final_coords_file and os.path.exists(final_coords_file):
                checkpoint_data[task_key] = {
                    'completed': True,
                    'timestamp': time.time(),
                    'output_file': final_coords_file
                }
                
                with self.lock:
                    self.completed_tasks += 1
                self.update_progress(force=False)
                
                return final_coords_file
            
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            error_file = f"{output_prefix}_error.log"
            with open(error_file, 'w') as f:
                f.write(f"Nucmer failed: {error_msg}")
        except Exception as e:
            error_file = f"{output_prefix}_error.log"
            with open(error_file, 'w') as f:
                f.write(f"Unknown error: {str(e)}")
        
        with self.lock:
            self.failed_tasks += 1
        self.update_progress(force=False)
        
        return None
    
    def calculate_threads_for_contig(self, contig_size: int) -> int:
        """Calculate threads based on contig size"""
        if contig_size > 10000000:
            return 8
        elif contig_size > 1000000:
            return 4
        elif contig_size > 100000:
            return 2
        else:
            return 1
    
    def process_delta_file(self, delta_file: str, output_prefix: str) -> Optional[str]:
        """Process delta file"""
        try:
            filter_cmd = ["delta-filter", "-i", "-r", "-l", "10000", delta_file]
            filtered_delta = delta_file.replace(".delta", ".filtered.delta")
            
            with open(filtered_delta, 'w') as outfile:
                subprocess.run(filter_cmd, stdout=outfile, stderr=subprocess.PIPE)
            
            coords_cmd = ["show-coords", "-r", "-l", filtered_delta]
            coords_file = filtered_delta.replace(".filtered.delta", ".coords")
            
            with open(coords_file, 'w') as outfile:
                subprocess.run(coords_cmd, stdout=outfile, stderr=subprocess.PIPE)
            
            for f in [delta_file, filtered_delta]:
                if os.path.exists(f):
                    os.remove(f)
            
            return coords_file
            
        except Exception:
            return None
    
    def create_all_tasks(self, large_files: List[str], small_files: List[str], 
                        query_data: Dict, checkpoint_data: Dict,
                        batch_size: int = 10) -> List[Dict]:
        """Create all task list"""
        tasks = []
        
        for i, ref_file in enumerate(large_files):
            contig_size = self.estimate_contig_size(ref_file)
            
            for seq_id, seq_data in query_data.items():
                output_prefix = os.path.join(seq_data['temp_dir'], f"large_{i+1}")
                
                task_key = f"large_{i+1}_{seq_id}"
                if not checkpoint_data.get(task_key, {}).get('completed', False):
                    tasks.append({
                        'ref_file': ref_file,
                        'query_file': seq_data['chr_fasta'],
                        'output_prefix': output_prefix,
                        'contig_size': contig_size,
                        'seq_id': seq_id,
                        'task_key': task_key
                    })
                else:
                    self.completed_tasks += 1
        
        if small_files:
            small_batches = []
            for batch_start in range(0, len(small_files), batch_size):
                batch_end = min(batch_start + batch_size, len(small_files))
                batch_files = small_files[batch_start:batch_end]
                small_batches.append(batch_files)
            
            for batch_idx, batch_files in enumerate(small_batches):
                for seq_id, seq_data in query_data.items():
                    temp_dir = seq_data['temp_dir']
                    merged_file = os.path.join(temp_dir, f"merged_batch_{batch_idx}.fa")
                    
                    if not os.path.exists(merged_file):
                        merged_file = self.merge_contigs_for_batch(batch_files, temp_dir, batch_idx)
                        if not merged_file:
                            continue
                    
                    total_size = sum(self.estimate_contig_size(f) for f in batch_files)
                    avg_size = total_size // len(batch_files) if batch_files else 0
                    
                    output_prefix = os.path.join(temp_dir, f"batch_{batch_idx}")
                    
                    task_key = f"batch_{batch_idx}_{seq_id}"
                    if not checkpoint_data.get(task_key, {}).get('completed', False):
                        tasks.append({
                            'ref_file': merged_file,
                            'query_file': seq_data['chr_fasta'],
                            'output_prefix': output_prefix,
                            'contig_size': avg_size,
                            'seq_id': seq_id,
                            'task_key': task_key
                        })
                    else:
                        self.completed_tasks += 1
        
        self.total_tasks = len(tasks) + self.completed_tasks
        
        if self.verbose:
            print(f"  {len(tasks)} tasks remaining (skipped {self.completed_tasks} completed tasks)")
        
        return tasks
    
    def estimate_contig_size(self, fasta_file: str) -> int:
        """Estimate contig size"""
        try:
            with open(fasta_file, 'r') as f:
                total = 0
                for line in f:
                    if line.startswith('>'):
                        continue
                    total += len(line.strip())
                    break
            return total
        except:
            return 100000
    
    def merge_contigs_for_batch(self, contig_files: List[str], output_dir: str, batch_id: int) -> Optional[str]:
        """Merge contigs for batch processing"""
        try:
            from Bio import SeqIO
            
            merged_records = []
            
            for i, contig_file in enumerate(contig_files):
                for record in SeqIO.parse(contig_file, "fasta"):
                    new_id = f"batch{batch_id}_contig{i}_{record.id}"
                    record.id = new_id
                    record.description = ""
                    merged_records.append(record)
            
            if merged_records:
                merged_file = os.path.join(output_dir, f"merged_batch_{batch_id}.fa")
                SeqIO.write(merged_records, merged_file, "fasta")
                return merged_file
            
        except Exception as e:
            if self.verbose:
                print(f"Warning: Failed to merge contigs: {e}")
        
        return None
    
    def execute_tasks_parallel(self, tasks: List[Dict], max_workers: int, 
                             checkpoint_file: str, checkpoint_data: Dict) -> List[str]:
        """Execute all tasks in parallel"""
        if not tasks:
            return []
        
        results = []
        
        if self.verbose:
            print(f"Starting parallel alignment ({max_workers} threads)...")
        
        checkpoint_lock = threading.Lock()
        last_checkpoint_time = time.time()
        checkpoint_interval = 30
        
        def task_wrapper(task):
            """Task wrapper"""
            result = self.run_nucmer_task(
                task['ref_file'],
                task['query_file'],
                task['output_prefix'],
                task['contig_size'],
                task['seq_id'],
                checkpoint_data
            )
            
            nonlocal last_checkpoint_time
            current_time = time.time()
            if current_time - last_checkpoint_time >= checkpoint_interval:
                with checkpoint_lock:
                    self.save_checkpoint(checkpoint_file, checkpoint_data)
                    last_checkpoint_time = current_time
            
            return result
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {}
            for task in tasks:
                future = executor.submit(task_wrapper, task)
                future_to_task[future] = task
            
            for future in concurrent.futures.as_completed(future_to_task):
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                except Exception:
                    pass
        
        self.save_checkpoint(checkpoint_file, checkpoint_data)
        
        self.update_progress(force=True)
        if self.verbose:
            print()
        
        return results
    
    def merge_results_per_chromosome(self, query_data: Dict):
        """Merge results for each chromosome"""
        if self.verbose:
            print("\nMerging result files for each chromosome...")
        
        for seq_id, seq_data in query_data.items():
            chr_dir = seq_data['chr_dir']
            temp_dir = seq_data['temp_dir']
            
            coords_files = []
            if os.path.exists(temp_dir):
                for file in os.listdir(temp_dir):
                    if file.endswith('.coords'):
                        coords_files.append(os.path.join(temp_dir, file))
            
            if coords_files:
                merged_file = os.path.join(chr_dir, f"{seq_id}_merged.coords")
                self.merge_coords_files(coords_files, merged_file)
                
                if self.verbose:
                    print(f"  {seq_id}: Merged {len(coords_files)} coords files")
    
    def merge_coords_files(self, coords_files: List[str], output_file: str):
        """Merge coords files"""
        if not coords_files:
            return
        
        all_lines = []
        for i, coords_file in enumerate(coords_files):
            try:
                with open(coords_file, 'r') as f:
                    lines = f.readlines()
                    
                    if i == 0:
                        all_lines.extend(lines)
                    else:
                        data_started = False
                        for line in lines:
                            if line.startswith('[') or line.startswith('='):
                                data_started = True
                                continue
                            if data_started and line.strip() and not line.startswith('NUCMER'):
                                all_lines.append(line)
            except Exception:
                continue
        
        if all_lines:
            with open(output_file, 'w') as f:
                f.writelines(all_lines)
    
    def cleanup_and_organize(self, query_data: Dict, split_dir: str):
        """Clean intermediate files and organize final results"""
        if self.verbose:
            print("\nCleaning intermediate files and organizing final results...")
        
        for seq_id, seq_data in query_data.items():
            temp_dir = seq_data['temp_dir']
            if os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    if self.verbose:
                        print(f"  Cleaned temporary files for {seq_id}")
                except:
                    pass
        
        if os.path.exists(split_dir):
            try:
                shutil.rmtree(split_dir, ignore_errors=True)
                if self.verbose:
                    print(f"  Cleaned shared split directory")
            except:
                pass
    
    def analyze(self, ref_fasta: str, query_fasta: str, output_dir: str,
                threads: int = 0, size_threshold: int = 1000000, 
                chunk_size: int = 10000000, batch_size: int = 10,
                min_alignment_length: int = 10000) -> Dict[str, Any]:
        """
        Run complete analysis pipeline
        """
        self.set_environment_limits()
        
        if threads == 0:
            threads = self.get_optimal_threads()
        
        if self.verbose:
            print("="*60)
            print("Fast Parallel Nucmer Analyzer")
            print("="*60)
            print(f"Reference genome: {ref_fasta}")
            print(f"Query genome: {query_fasta}")
            print(f"Output directory: {output_dir}")
            print(f"Threads: {threads}")
            print(f"Large contig threshold: {size_threshold:,} bp")
            print(f"Small contigs combination size: {chunk_size:,} bp")
            print(f"Batch size: {batch_size}")
            if self.specified_chromosomes:
                print(f"Specified chromosomes: {', '.join(self.specified_chromosomes)}")
            print("-"*60)
        
        checkpoint_file = os.path.join(output_dir, "checkpoint.json")
        checkpoint_data = self.load_checkpoint(checkpoint_file)
        
        try:
            query_data = self.create_chromosome_folders(query_fasta, output_dir)
            
            large_files, small_files, split_dir = self.split_contigs_optimized(
                ref_fasta=ref_fasta,
                output_base_dir=output_dir,
                size_threshold=size_threshold,
                chunk_size=chunk_size
            )
            
            tasks = self.create_all_tasks(
                large_files=large_files,
                small_files=small_files,
                query_data=query_data,
                checkpoint_data=checkpoint_data,
                batch_size=batch_size
            )
            
            start_time = time.time()
            results = self.execute_tasks_parallel(
                tasks, 
                max_workers=threads,
                checkpoint_file=checkpoint_file,
                checkpoint_data=checkpoint_data
            )
            alignment_time = time.time() - start_time
            
            self.merge_results_per_chromosome(query_data)
            
            syntenic_results = {}
            for seq_id, seq_data in query_data.items():
                merged_file = os.path.join(seq_data['chr_dir'], f"{seq_id}_merged.coords")
                if os.path.exists(merged_file):
                    contig_count, contig_file = _extract_syntenic_contigs_from_coords(
                        coord_file=merged_file,
                        ref_fasta=ref_fasta,
                        seq_id=seq_id,
                        output_dir=seq_data['chr_dir'],
                        min_alignment_length=min_alignment_length
                    )
                    syntenic_results[seq_id] = {
                        'contig_count': contig_count,
                        'contig_file': contig_file,
                        'chr_dir': seq_data['chr_dir']
                    }
            
            self.cleanup_and_organize(query_data, split_dir)
            
            if os.path.exists(checkpoint_file):
                os.remove(checkpoint_file)
            
            total_time = time.time() - self.start_time
            
            report = {
                'parameters': {
                    'ref_fasta': ref_fasta,
                    'query_fasta': query_fasta,
                    'threads': threads,
                    'size_threshold': size_threshold,
                    'chunk_size': chunk_size,
                    'batch_size': batch_size,
                    'min_alignment_length': min_alignment_length,
                    'specified_chromosomes': self.specified_chromosomes if self.specified_chromosomes else "all"
                },
                'statistics': {
                    'total_time': total_time,
                    'alignment_time': alignment_time,
                    'total_tasks': self.total_tasks,
                    'completed_tasks': self.completed_tasks,
                    'failed_tasks': self.failed_tasks,
                    'success_rate': (self.completed_tasks / self.total_tasks * 100) if self.total_tasks > 0 else 0,
                    'chromosomes': len(query_data),
                    'large_contigs': len(large_files),
                    'small_contigs': len(small_files)
                },
                'chromosome_results': {},
                'syntenic_results': syntenic_results,
                'output_structure': {
                    'base_dir': output_dir,
                    'chromosome_folders': []
                }
            }
            
            for seq_id, seq_data in query_data.items():
                chr_dir = seq_data['chr_dir']
                report['chromosome_results'][seq_id] = {
                    'dir': chr_dir,
                    'original_id': seq_data.get('original_id', seq_id),
                    'seq_length': seq_data['seq_length'],
                    'files': []
                }
                
                for file in os.listdir(chr_dir):
                    if os.path.isfile(os.path.join(chr_dir, file)):
                        report['chromosome_results'][seq_id]['files'].append(file)
                
                report['output_structure']['chromosome_folders'].append(f"{seq_id}/")
            
            report_file = os.path.join(output_dir, "analysis_report.json")
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2, ensure_ascii=False)
            
            if self.verbose:
                self._print_summary(report)
            
            return report
            
        except Exception as e:
            error_report = {
                'error': str(e),
                'timestamp': time.time(),
                'checkpoint_file': checkpoint_file if os.path.exists(checkpoint_file) else None
            }
            error_file = os.path.join(output_dir, "error_report.json")
            with open(error_file, 'w') as f:
                json.dump(error_report, f, indent=2)
            
            raise
    
    def _print_summary(self, report: Dict[str, Any]):
        """Print summary"""
        print("\n" + "="*60)
        print("Analysis completed!")
        print("="*60)
        
        stats = report['statistics']
        print(f"Total runtime: {stats['total_time']:.1f} seconds")
        print(f"Alignment time: {stats['alignment_time']:.1f} seconds")
        print(f"Tasks completed: {stats['completed_tasks']}/{stats['total_tasks']}")
        print(f"Success rate: {stats['success_rate']:.1f}%")
        
        print(f"\nProcessed {stats['chromosomes']} chromosomes:")
        for seq_id, chr_data in report['chromosome_results'].items():
            syntenic_info = report['syntenic_results'].get(seq_id, {})
            contig_count = syntenic_info.get('contig_count', 0)
            print(f"  {seq_id}: {chr_data['seq_length']:,} bp, syntenic contigs: {contig_count}")
        
        print(f"\nOutput directory structure:")
        print(f"  {report['output_structure']['base_dir']}/")
        for chr_folder in report['output_structure']['chromosome_folders']:
            print(f"  ├── {chr_folder}")
        
        print("\n" + "="*60)

def _extract_syntenic_contigs_from_coords(
    coord_file: str, 
    ref_fasta: str, 
    seq_id: str, 
    output_dir: str,
    min_alignment_length: int = 10000,
    min_identity: float = 95.0,
    min_coverage: float = 0.5
) -> Tuple[int, Optional[str]]:
    """Extract syntenic contigs from coords file"""
    if not os.path.exists(coord_file):
        return 0, None
    
    from Bio import SeqIO
    
    output_file = os.path.join(output_dir, f"syntenic_contigs_{seq_id}.fa")
    
    alignments = []
    try:
        with open(coord_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('/') or line.startswith('NUCMER') or \
                   line.startswith('[') or line.startswith('='):
                    continue
                parts = line.split()
                if len(parts) >= 12:
                    try:
                        ref_start = int(parts[0])
                        ref_end = int(parts[1])
                        ref_contig = parts[-2]
                        qry_contig = parts[-1]
                        
                        if qry_contig == seq_id:
                            align_length = abs(ref_end - ref_start) + 1
                            if align_length >= min_alignment_length:
                                if len(parts) > 6:
                                    identity = float(parts[6])
                                    if identity >= min_identity:
                                        alignments.append({
                                            'ref_contig': ref_contig,
                                            'align_length': align_length
                                        })
                                else:
                                    alignments.append({
                                        'ref_contig': ref_contig,
                                        'align_length': align_length
                                    })
                    except:
                        continue
    except Exception as e:
        print(f"Warning: Failed to read coords file {coord_file}: {e}")
        return 0, None
    
    if not alignments:
        return 0, None
    
    ref_contig_sequences = {}
    try:
        for record in SeqIO.parse(ref_fasta, "fasta"):
            ref_contig_sequences[record.id] = record
    except Exception as e:
        print(f"Warning: Unable to read reference genome: {e}")
        return 0, None
    
    contig_stats = {}
    for align in alignments:
        contig_stats[align['ref_contig']] = contig_stats.get(align['ref_contig'], 0) + align['align_length']
    
    sorted_contigs = sorted(contig_stats.items(), key=lambda x: x[1], reverse=True)
    
    extracted_contigs = []
    for contig_name, total_align_length in sorted_contigs:
        if contig_name in ref_contig_sequences:
            extracted_contigs.append(ref_contig_sequences[contig_name])
    
    if extracted_contigs:
        try:
            SeqIO.write(extracted_contigs, output_file, "fasta")
            return len(extracted_contigs), output_file
        except Exception as e:
            print(f"Warning: Failed to save syntenic contigs: {e}")
    
    return 0, None

def run_nucmer_analysis(
    reference_fasta: str, 
    query_fasta: str, 
    output_dir: str,
    threads: int = 0, 
    verbose: bool = True, 
    size_threshold: int = 1000000, 
    chunk_size: int = 10000000,
    batch_size: int = 10, 
    min_alignment_length: int = 10000,
    overwrite: bool = False,
    keep_temp: bool = False,
    chromosomes: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Simplified interface for running nucmer analysis
    
    Example:
        >>> from nucmer_analyzer import run_nucmer_analysis
        >>> result = run_nucmer_analysis(
        ...     reference_fasta="reference.fasta",
        ...     query_fasta="query.fasta",
        ...     output_dir="results",
        ...     threads=16,
        ...     verbose=False,
        ...     chromosomes=["chr1", "chr2"]
        ... )
    """
    analyzer = NucmerAnalyzerAPI(verbose=verbose)
    return analyzer.run_analysis(
        reference_fasta=reference_fasta,
        query_fasta=query_fasta,
        output_dir=output_dir,
        output_prefix="nucmer_analysis",
        threads=threads,
        size_threshold=size_threshold,
        chunk_size=chunk_size,
        batch_size=batch_size,
        min_alignment_length=min_alignment_length,
        overwrite=overwrite,
        keep_temp=keep_temp,
        chromosomes=chromosomes
    )

def main():
    """Command line entry function"""
    parser = argparse.ArgumentParser(
        description="Fast Parallel Nucmer Analyzer - Unified API Version",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Basic usage
  python run_nucmer_parallel.py -g ref.fasta -q query.fasta -o results
  
  # Specify threads, use API mode, analyze specified chromosomes only
  python run_nucmer_parallel.py -g ref.fasta -q query.fasta -o results -T 16 --api-mode --chromosomes chr1 chr2 chr3
  
  # Read chromosome list from file
  python run_nucmer_parallel.py -g ref.fasta -q query.fasta -o results --chromosomes-file chromosomes.txt
  
  # Python module call, analyze specified chromosomes only
  from run_nucmer_parallel import NucmerAnalyzerAPI
  analyzer = NucmerAnalyzerAPI(verbose=True)
  result = analyzer.run_analysis(
      reference_fasta="ref.fasta",
      query_fasta="query.fasta", 
      output_dir="results",
      threads=16,
      chromosomes=["chr1", "chr2"]
  )
        """
    )
    
    parser.add_argument("-g", "--genome", dest="reference_fasta", required=True, 
                       help="Reference genome file")
    parser.add_argument("-q", "--query", dest="query_fasta", required=True, 
                       help="Query genome file")
    parser.add_argument("-o", "--output-dir", required=True, 
                       help="Output directory")
    
    parser.add_argument("-T", "--threads", type=int, default=0, 
                       help="Thread count (0=auto-detect)")
    parser.add_argument("--size-threshold", type=int, default=1000000, 
                       help="Large contig threshold, default 1M bp")
    parser.add_argument("--chunk-size", type=int, default=10000000, 
                       help="Small contigs combination size, default 10M bp")
    parser.add_argument("--batch-size", type=int, default=10, 
                       help="Batch size, default 10")
    parser.add_argument("--min-alignment-length", type=int, default=10000, 
                       help="Minimum alignment length, default 10k bp")
    parser.add_argument("--output-prefix", default="nucmer_analysis", 
                       help="Output file prefix")
    parser.add_argument("--overwrite", action="store_true", 
                       help="Overwrite existing results")
    parser.add_argument("--keep-temp", action="store_true", 
                       help="Keep temporary files")
    parser.add_argument("--api-mode", action="store_true", 
                       help="Use API mode")
    parser.add_argument("--quiet", action="store_true", 
                       help="Quiet mode")
    parser.add_argument("--log-file", help="Log file path")
    
    parser.add_argument("--chromosomes", nargs="+", help="List of chromosomes to analyze")
    parser.add_argument("--chromosomes-file", help="File containing chromosome list, one per line")
    
    args = parser.parse_args()
    
    chromosomes = None
    if args.chromosomes:
        chromosomes = args.chromosomes
    elif args.chromosomes_file:
        try:
            with open(args.chromosomes_file, 'r') as f:
                chromosomes = [line.strip() for line in f if line.strip()]
        except Exception as e:
            print(f"Error: Unable to read chromosome file {args.chromosomes_file}: {e}")
            sys.exit(1)
    
    try:
        if args.api_mode:
            analyzer = NucmerAnalyzerAPI(verbose=not args.quiet, log_file=args.log_file)
            
            result = analyzer.run_analysis(
                reference_fasta=args.reference_fasta,
                query_fasta=args.query_fasta,
                output_dir=args.output_dir,
                output_prefix=args.output_prefix,
                threads=args.threads,
                size_threshold=args.size_threshold,
                chunk_size=args.chunk_size,
                batch_size=args.batch_size,
                min_alignment_length=args.min_alignment_length,
                overwrite=args.overwrite,
                keep_temp=args.keep_temp,
                chromosomes=chromosomes
            )
            
            if not args.quiet:
                print(f"\nDetailed report saved to: {args.output_dir}/{args.output_prefix}_api_report.json")
                print(f"Results for each chromosome saved in respective folders")
            
        else:
            analyzer = _FastNucmerParallelAnalyzer(verbose=not args.quiet, specified_chromosomes=chromosomes)
            
            result = analyzer.analyze(
                ref_fasta=args.reference_fasta,
                query_fasta=args.query_fasta,
                output_dir=args.output_dir,
                threads=args.threads,
                size_threshold=args.size_threshold,
                chunk_size=args.chunk_size,
                batch_size=args.batch_size,
                min_alignment_length=args.min_alignment_length
            )
            
            if not args.quiet:
                print(f"\nDetailed report saved to: {args.output_dir}/analysis_report.json")
        
        if result['statistics']['failed_tasks'] > 0:
            print(f"Warning: {result['statistics']['failed_tasks']} tasks failed")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\nAnalysis interrupted, can continue with same command next time")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()