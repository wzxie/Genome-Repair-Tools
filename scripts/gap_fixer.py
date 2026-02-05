#!/usr/bin/env python3
"""
Genome gap analysis and repair controller
Integrates 3 sub-scripts to implement gap detection, analysis and repair pipeline
Supports parallel processing of multiple chromosomes
Uses aggressive repair mode by default
Added adapter pattern support
New: Support -coords parameter to directly use existing coords file (skip nucmer analysis)
New: Support parameter passing to gap_analyzer.py
Modified: -o parameter supports folder, generates fixed_genome.fasta in folder
"""

import sys
import os
import argparse
import logging
import json
import shutil
import time
import concurrent.futures
import multiprocessing
import glob
from functools import partial
from typing import Dict, List, Tuple, Optional, Any, Callable
from pathlib import Path

try:
    from find_gaps import GapFinderAPI
    from run_nucmer_parallel import NucmerAnalyzerAPI
    from gap_analyzer import GapAnalyzerAPI
    
    ALL_APIS_AVAILABLE = True
except ImportError as e:
    print(f"Error: Cannot import sub-script APIs: {e}")
    print("Please ensure all scripts are in the same directory")
    ALL_APIS_AVAILABLE = False
    sys.exit(1)

class GenomeGapFixer:
    
    def __init__(self, total_threads: int = 4, verbose: bool = True, log_file: Optional[str] = None,
                 repair_mode: str = "aggressive",
                 gap_length: int = 100,
                 max_search_distance: int = 500000,
                 search_step: int = 100000):
        
        self.total_threads = total_threads
        self.verbose = verbose
        self.log_file = log_file
        self.repair_mode = repair_mode
        self.gap_length = gap_length
        self.max_search_distance = max_search_distance
        self.search_step = search_step
        self.work_dir = "gap_fix_workdir"
        self._setup_logging()
        self._setup_apis()
        
        self.current_chromosome = None
        self.current_step = 0
        self.total_steps = 6
        self.gap_history = {}
        self.status = {}
    
    def _setup_logging(self):
        self.logger = logging.getLogger('GenomeGapFixer')
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
        if level == "info":
            self.logger.info(message)
        elif level == "warning":
            self.logger.warning(message)
        elif level == "error":
            self.logger.error(message)
        elif level == "debug":
            self.logger.debug(message)
    
    def _setup_apis(self):
        self.log("Initializing sub-script APIs...", "info")
        
        self.gap_finder = GapFinderAPI(verbose=self.verbose)
        self.nucmer_analyzer = NucmerAnalyzerAPI(verbose=self.verbose)
        self.gap_analyzer = GapAnalyzerAPI(verbose=self.verbose)
        
        self.log(f"All APIs initialized (Total threads: {self.total_threads}, Repair mode: {self.repair_mode})", "info")
    
    def _update_progress(self, step_name: str):
        self.current_step += 1
        progress = (self.current_step / self.total_steps) * 100
        self.log(f"[Progress {progress:.1f}%] Step {self.current_step}/{self.total_steps}: {step_name}", "info")
    
    def _reset_progress(self):
        self.current_step = 0
    
    def _ensure_work_dir(self):
        os.makedirs(self.work_dir, exist_ok=True)
    
    def _cleanup_work_dir(self):
        if os.path.exists(self.work_dir):
            try:
                shutil.rmtree(self.work_dir, ignore_errors=True)
                self.log(f"Cleaned work directory: {self.work_dir}", "info")
            except Exception as e:
                self.log(f"Cannot clean work directory: {e}", "warning")
    
    def calculate_thread_allocation(self, num_chromosomes: int) -> Tuple[int, int]:
        if num_chromosomes == 0:
            return 1, 1
        
        if num_chromosomes <= self.total_threads:
            per_chrom_threads = max(2, self.total_threads // num_chromosomes)
            parallel_chromosomes = min(num_chromosomes, max(1, self.total_threads // per_chrom_threads))
        else:
            per_chrom_threads = 2
            parallel_chromosomes = max(1, self.total_threads // per_chrom_threads)
        
        self.log(f"Thread allocation: {num_chromosomes} chromosomes, total threads {self.total_threads} -> "
                f"{per_chrom_threads} threads per chromosome, {parallel_chromosomes} parallel chromosomes", "info")
        
        return per_chrom_threads, parallel_chromosomes
    
    def find_gaps_in_genome(self, query_fasta: str) -> Dict[str, Any]:
        self.log(f"\nStep 1: Detecting gaps in genome", "info")
        self.log(f"Analyzing file: {query_fasta}", "info")
        
        self._update_progress("Detecting gaps")
        
        try:
            result = self.gap_finder.find_gaps(
                fasta_file=query_fasta,
                output_dir=self.work_dir,
                output_prefix="initial_gap_analysis",
                gap_char="N",
                min_gap_size=100,
                parallel_threshold_mb=100,
                save_json=True,
                save_bed=True,
                save_summary=True
            )
            
            self.gap_history["initial"] = {
                "total_gaps": result["statistics"]["total_gaps"],
                "total_length": result["statistics"]["total_gap_length"],
                "gap_coverage": result["statistics"]["gap_coverage_percent"]
            }
            
            chromosomes_with_gaps = {}
            gaps_by_chromosome = {}
            gap_positions_by_chromosome = {}
            
            for gap in result.get("gaps", []):
                chrom = gap["chrom"]
                if chrom not in chromosomes_with_gaps:
                    chromosomes_with_gaps[chrom] = []
                    gaps_by_chromosome[chrom] = []
                    gap_positions_by_chromosome[chrom] = []
                
                chromosomes_with_gaps[chrom].append(gap)
                gaps_by_chromosome[chrom].append({
                    "start": gap["start"],
                    "end": gap["end"],
                    "length": gap["length"],
                    "position": gap.get("position_human", f"{gap['start']+1}-{gap['end']}")
                })
                
                gap_mid = gap["start"] + gap["length"] // 2
                gap_positions_by_chromosome[chrom].append(gap_mid)
            
            self.log(f"Found {len(chromosomes_with_gaps)} chromosomes with gaps", "info")
            for chrom, gaps in chromosomes_with_gaps.items():
                self.log(f"  {chrom}: {len(gaps)} gaps, total length: {sum(g['length'] for g in gaps):,} bp", "info")
            
            return {
                "result": result,
                "chromosomes_with_gaps": chromosomes_with_gaps,
                "gaps_by_chromosome": gaps_by_chromosome,
                "gap_positions_by_chromosome": gap_positions_by_chromosome,
                "total_gaps": result["statistics"]["total_gaps"],
                "total_gap_length": result["statistics"]["total_gap_length"]
            }
            
        except Exception as e:
            self.log(f"Gap detection failed: {e}", "error")
            raise
    
    def extract_chromosomes_with_gaps(self, query_fasta: str, chromosomes: List[str]) -> Dict[str, str]:
        self.log(f"\nStep 2: Extracting chromosomes with gaps", "info")
        self.log(f"Extracting {len(chromosomes)} chromosomes", "info")
        
        self._update_progress("Extracting chromosomes")
        
        from Bio import SeqIO
        
        chromosome_files = {}
        
        try:
            records = list(SeqIO.parse(query_fasta, "fasta"))
            
            for chrom in chromosomes:
                chrom_records = [r for r in records if r.id == chrom]
                if not chrom_records:
                    self.log(f"Warning: Chromosome {chrom} not found in FASTA file", "warning")
                    continue
                
                chrom_file = os.path.join(self.work_dir, f"{chrom}.fa")
                SeqIO.write(chrom_records, chrom_file, "fasta")
                chromosome_files[chrom] = chrom_file
                
                self.log(f"  Extracted: {chrom} -> {chrom_file}", "info")
            
            return chromosome_files
            
        except Exception as e:
            self.log(f"Chromosome extraction failed: {e}", "error")
            raise
    
    def load_single_coords_file(self, coords_file: str, chromosomes: List[str]) -> Dict[str, Any]:
        if not os.path.exists(coords_file):
            self.log(f"Error: Coords file does not exist: {coords_file}", "error")
            return {}
        
        self.log(f"Using existing coords file: {coords_file}", "info")
        
        if len(chromosomes) > 1:
            self.log(f"Warning: Found {len(chromosomes)} chromosomes with gaps, but only one coords file provided", "warning")
            self.log(f"Will use first chromosome: {chromosomes[0]}", "warning")
        
        target_chromosome = chromosomes[0] if chromosomes else "unknown"
        
        nucmer_results = {
            target_chromosome: {
                "chromosome": target_chromosome,
                "success": True,
                "coords_file": coords_file,
                "output_dir": os.path.dirname(coords_file),
                "message": f"Using specified coords file: {os.path.basename(coords_file)}"
            }
        }
        
        return nucmer_results
    
    def run_nucmer_for_chromosome_parallel(self, chrom_info: Tuple[str, str, str], 
                                         per_chrom_threads: int = 2) -> Dict[str, Any]:
        chrom, chrom_file, reference_fasta = chrom_info
        
        self.log(f"  Starting chromosome {chrom} (threads: {per_chrom_threads})", "info")
        
        chrom_output_dir = os.path.join(self.work_dir, f"nucmer_{chrom}")
        
        try:
            result = self.nucmer_analyzer.run_analysis(
                reference_fasta=reference_fasta,
                query_fasta=chrom_file,
                output_dir=chrom_output_dir,
                output_prefix=f"nucmer_{chrom}",
                threads=per_chrom_threads,
                size_threshold=1000000,
                chunk_size=10000000,
                batch_size=10,
                min_alignment_length=10000,
                overwrite=True,
                keep_temp=False
            )
            
            coords_file = None
            for seq_id, chr_data in result.get("chromosome_results", {}).items():
                chr_dir = chr_data.get("dir", "")
                if os.path.exists(chr_dir):
                    for file in os.listdir(chr_dir):
                        if file.endswith("_merged.coords"):
                            coords_file = os.path.join(chr_dir, file)
                            break
                if coords_file:
                    break
            
            self.log(f"  ✓ Chromosome {chrom} analysis completed", "info")
            
            return {
                "chromosome": chrom,
                "success": True,
                "coords_file": coords_file,
                "output_dir": chrom_output_dir,
                "result": result
            }
            
        except Exception as e:
            self.log(f"  ✗ Synteny analysis failed for chromosome {chrom}: {e}", "error")
            return {
                "chromosome": chrom,
                "success": False,
                "error": str(e),
                "output_dir": chrom_output_dir
            }
    
    def run_nucmer_parallel_for_all_chromosomes(self, chromosome_files: Dict[str, str], 
                                              reference_fasta: str) -> Dict[str, Any]:
        self.log(f"\nStep 3: Parallel synteny analysis", "info")
        self.log(f"Total threads: {self.total_threads}", "info")
        self.log(f"Processing chromosomes: {len(chromosome_files)}", "info")
        
        self._update_progress("Parallel synteny analysis")
        
        per_chrom_threads, parallel_chromosomes = self.calculate_thread_allocation(
            len(chromosome_files)
        )
        
        self.log(f"Thread allocation: {per_chrom_threads} threads per chromosome, {parallel_chromosomes} parallel chromosomes", "info")
        
        chrom_infos = []
        for chrom, chrom_file in chromosome_files.items():
            chrom_infos.append((chrom, chrom_file, reference_fasta))
        
        nucmer_results = {}
        
        try:
            self.log(f"Starting thread pool, max parallelism: {parallel_chromosomes}", "info")
            with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_chromosomes) as executor:
                process_func = partial(self.run_nucmer_for_chromosome_parallel, 
                                     per_chrom_threads=per_chrom_threads)
                
                future_to_chrom = {
                    executor.submit(process_func, chrom_info): chrom_info[0] 
                    for chrom_info in chrom_infos
                }
                
                completed = 0
                total = len(chrom_infos)
                start_time = time.time()
                
                for future in concurrent.futures.as_completed(future_to_chrom):
                    chrom = future_to_chrom[future]
                    try:
                        result = future.result(timeout=7200)
                        nucmer_results[chrom] = result
                        
                        if not result["success"]:
                            self.log(f"  ✗ Chromosome {chrom} analysis failed: {result.get('error', 'Unknown error')}", "warning")
                            
                    except concurrent.futures.TimeoutError:
                        self.log(f"  ⚠ Chromosome {chrom} analysis timeout (2 hours)", "warning")
                        nucmer_results[chrom] = {
                            "chromosome": chrom,
                            "success": False,
                            "error": "Processing timeout"
                        }
                    except Exception as e:
                        self.log(f"  ✗ Chromosome {chrom} processing exception: {e}", "error")
                        nucmer_results[chrom] = {
                            "chromosome": chrom,
                            "success": False,
                            "error": str(e)
                        }
                    
                    completed += 1
                    elapsed = time.time() - start_time
                    avg_time = elapsed / completed if completed > 0 else 0
                    remaining = avg_time * (total - completed) if completed > 0 else 0
                    
                    progress = (completed / total) * 100
                    self.log(f"  Progress: {completed}/{total} ({progress:.1f}%) | "
                            f"Elapsed: {elapsed:.1f}s | Estimated remaining: {remaining:.1f}s", "debug")
            
            successful = sum(1 for r in nucmer_results.values() if r.get("success", False))
            failed = len(nucmer_results) - successful
            
            if failed > 0:
                self.log(f"Parallel analysis completed: {successful} successful, {failed} failed", "warning")
                failed_chroms = [chrom for chrom, res in nucmer_results.items() 
                               if not res.get("success", False)]
                self.log(f"Failed chromosomes: {', '.join(failed_chroms)}", "warning")
            else:
                self.log(f"✓ All {successful} chromosomes analyzed successfully", "info")
            
            return nucmer_results
            
        except Exception as e:
            self.log(f"Parallel processing failed: {e}", "error")
            import traceback
            traceback.print_exc()
            return {}
    
    def analyze_and_repair_gaps(self, chromosome: str, chrom_file: str, coords_file: str, 
                              gap_positions: List[int]) -> str:
        self.log(f"\nStep 4: Analyzing and repairing gaps in chromosome {chromosome}", "info")
        self.log(f"  Using repair mode: {self.repair_mode}", "info")
        self.log(f"  Final gap length: {self.gap_length} bp", "info")
        self.log(f"  Maximum search distance: {self.max_search_distance:,} bp", "info")
        self.log(f"  Search step: {self.search_step:,} bp", "info")
        
        self._update_progress(f"Analyzing and repairing ({chromosome})")
        
        repaired_dir = os.path.join(self.work_dir, f"repaired_{chromosome}")
        os.makedirs(repaired_dir, exist_ok=True)
        repaired_file = os.path.join(repaired_dir, f"repaired_{chromosome}.fa")
        
        try:
            result = self.gap_analyzer.analyze_gaps(
                coords_file=coords_file,
                gap_positions=gap_positions,
                query_fasta=chrom_file,
                output_dir=repaired_dir,
                output_prefix=f"repair_{chromosome}",
                max_search_distance=self.max_search_distance,
                search_step=self.search_step,
                final_gap_length=self.gap_length,
                min_confidence="low",
                repair_mode=self.repair_mode
            )
            
            repaired_fasta = result["api_info"]["output_files"].get("filled_fasta")
            if repaired_fasta and os.path.exists(repaired_fasta):
                if repaired_fasta != repaired_file:
                    shutil.copy2(repaired_fasta, repaired_file)
                
                stats = result.get("summary", {})
                self.log(f"  Repair completed: Processed {stats.get('gaps_repaired', 0)} error regions", "info")
                self.log(f"  Repair mode: {self.repair_mode}", "info")
                self.log(f"  Final gap length: {self.gap_length} bp", "info")
                
                if self.repair_mode == "aggressive" and "repair_summary" in stats:
                    repair_summary = stats["repair_summary"]
                    self.log(f"  Additional aggressive repairs: {repair_summary.get('aggressive_repairs', 0)}", "info")
                
                return repaired_file
            else:
                self.log(f"  Warning: Repair failed, using original file", "warning")
                shutil.copy2(chrom_file, repaired_file)
                return repaired_file
                
        except Exception as e:
            self.log(f"Chromosome {chromosome} repair failed: {e}", "error")
            shutil.copy2(chrom_file, repaired_file)
            return repaired_file
    
    def process_single_chromosome(self, chrom_info: Tuple[str, str, str, List[int], Dict[str, Any]]) -> Dict[str, Any]:
        chrom, chrom_file, reference_fasta, gap_positions, nucmer_result = chrom_info
        
        self.log(f"\n{'='*40}", "info")
        self.log(f"Starting chromosome: {chrom}", "info")
        self.log(f"Repair mode: {self.repair_mode}", "info")
        self.log(f"Final gap length: {self.gap_length} bp", "info")
        self.log(f"{'='*40}", "info")
        
        result = {
            "chromosome": chrom,
            "success": False,
            "error": None,
            "final_file": None,
            "stats": {},
            "repair_mode": self.repair_mode,
            "gap_length": self.gap_length,
            "max_search_distance": self.max_search_distance,
            "search_step": self.search_step
        }
        
        try:
            if not nucmer_result or not nucmer_result.get("success", False):
                result["error"] = "nucmer analysis failed"
                result["final_file"] = chrom_file
                self.log(f"Chromosome {chrom}: nucmer analysis failed, skipping further steps", "warning")
                return result
            
            coords_file = nucmer_result.get("coords_file")
            if not coords_file or not os.path.exists(coords_file):
                result["error"] = "coords file does not exist"
                result["final_file"] = chrom_file
                self.log(f"Chromosome {chrom}: coords file does not exist", "warning")
                return result
            
            if gap_positions:
                repaired_file = self.analyze_and_repair_gaps(chrom, chrom_file, coords_file, 
                                                           gap_positions)
                result["final_file"] = repaired_file
                result["stats"]["gaps_processed"] = len(gap_positions)
            else:
                result["final_file"] = chrom_file
                result["stats"]["gaps_processed"] = 0
            
            result["success"] = True
            self.log(f"Chromosome {chrom} processing completed (mode: {self.repair_mode}, gap length: {self.gap_length}bp)", "info")
            
        except Exception as e:
            result["error"] = str(e)
            result["final_file"] = chrom_file
            self.log(f"Chromosome {chrom} processing failed: {e}", "error")
        
        return result
    
    def merge_all_chromosomes(self, chromosome_files: Dict[str, str], output_fasta: str) -> str:
        self.log(f"\nStep 5: Merging all chromosomes into single file", "info")
        
        self._update_progress("Merging chromosomes")
        
        from Bio import SeqIO
        
        try:
            all_records = []
            
            for chrom, chrom_file in chromosome_files.items():
                if os.path.exists(chrom_file):
                    records = list(SeqIO.parse(chrom_file, "fasta"))
                    all_records.extend(records)
                    self.log(f"  Adding chromosome: {chrom}", "info")
                else:
                    self.log(f"  Warning: Chromosome file does not exist: {chrom_file}", "warning")
            
            if all_records:
                SeqIO.write(all_records, output_fasta, "fasta")
                self.log(f"  Merge completed: {len(all_records)} sequences", "info")
                return output_fasta
            else:
                self.log(f"  Error: No chromosome sequences to merge", "error")
                return ""
                
        except Exception as e:
            self.log(f"Chromosome merging failed: {e}", "error")
            return ""
    
    def analyze_final_genome(self, final_fasta: str) -> Dict[str, Any]:
        self.log(f"\nStep 6: Analyzing gaps in final genome", "info")
        
        self._update_progress("Final gap analysis")
        
        try:
            result = self.gap_finder.find_gaps(
                fasta_file=final_fasta,
                output_dir=self.work_dir,
                output_prefix="final_gap_analysis",
                gap_char="N",
                min_gap_size=100,
                parallel_threshold_mb=100,
                save_json=True,
                save_bed=True,
                save_summary=True
            )
            
            self.gap_history["final"] = {
                "total_gaps": result["statistics"]["total_gaps"],
                "total_length": result["statistics"]["total_gap_length"],
                "gap_coverage": result["statistics"]["gap_coverage_percent"]
            }
            
            if "initial" in self.gap_history:
                initial = self.gap_history["initial"]
                final = self.gap_history["final"]
                
                improvement = {
                    "gaps_analyzed": initial.get("total_gaps", 0),
                    "gaps_remaining": final.get("total_gaps", 0),
                    "gaps_fixed": initial.get("total_gaps", 0) - final.get("total_gaps", 0),
                    "gaps_fixed_percent": ((initial.get("total_gaps", 0) - final.get("total_gaps", 0)) / 
                                         max(initial.get("total_gaps", 1), 1) * 100),
                    "length_fixed": initial.get("total_length", 0) - final.get("total_length", 0),
                    "length_fixed_percent": ((initial.get("total_length", 0) - final.get("total_length", 0)) / 
                                           max(initial.get("total_length", 1), 1) * 100),
                    "coverage_reduction": initial.get("gap_coverage", 0) - final.get("gap_coverage", 0),
                    "repair_mode": self.repair_mode,
                    "gap_length": self.gap_length,
                    "max_search_distance": self.max_search_distance,
                    "search_step": self.search_step
                }
                
                result["improvement"] = improvement
            
            return result
            
        except Exception as e:
            self.log(f"Final gap analysis failed: {e}", "error")
            return {}
    
    def run_complete_pipeline(self, query_fasta: str, reference_fasta: str = None,
                            coords_file: str = None,
                            output_fasta: str = "fixed_genome.fasta") -> Dict[str, Any]:
        
        self.log("=" * 60, "info")
        self.log("Starting genome gap repair pipeline (parallel version)", "info")
        self.log(f"Repair mode: {self.repair_mode} (default: aggressive)", "info")
        self.log(f"Final gap length: {self.gap_length} bp", "info")
        self.log(f"Maximum search distance: {self.max_search_distance:,} bp", "info")
        self.log(f"Search step: {self.search_step:,} bp", "info")
        
        if coords_file:
            self.log(f"Using existing coords file, skipping nucmer analysis: {coords_file}", "info")
        elif reference_fasta:
            self.log(f"Will run nucmer analysis to generate new coords file", "info")
        else:
            self.log("Error: Must provide coords_file or reference_fasta", "error")
            return {"status": "error", "error": "Must provide coords_file or reference_fasta"}
        
        self.log("=" * 60, "info")
        self.log(f"Query genome: {query_fasta}", "info")
        if reference_fasta:
            self.log(f"Reference genome: {reference_fasta}", "info")
        if coords_file:
            self.log(f"Alignment file: {coords_file}", "info")
        self.log(f"Total threads: {self.total_threads}", "info")
        
        start_time = time.time()
        
        try:
            self._ensure_work_dir()
            
            gap_analysis = self.find_gaps_in_genome(query_fasta)
            
            if gap_analysis["total_gaps"] == 0:
                self.log("No gaps found, no processing needed", "info")
                shutil.copy2(query_fasta, output_fasta)
                return {"status": "no_gaps", "output_file": output_fasta}
            
            chromosomes = list(gap_analysis["chromosomes_with_gaps"].keys())
            self.log(f"Need to process {len(chromosomes)} chromosomes with gaps", "info")
            
            chromosome_files = self.extract_chromosomes_with_gaps(query_fasta, chromosomes)
            
            if not chromosome_files:
                self.log("Error: Cannot extract chromosome files", "error")
                return {"status": "error", "error": "Cannot extract chromosome files"}
            
            self.log(f"\n{'='*40}", "info")
            self.log(f"Step 3: Synteny analysis", "info")
            self.log(f"{'='*40}", "info")
            
            if coords_file:
                nucmer_results = self.load_single_coords_file(coords_file, chromosomes)
            else:
                if not reference_fasta:
                    self.log("Error: Must provide reference_fasta or coords_file", "error")
                    return {"status": "error", "error": "Must provide reference_fasta or coords_file"}
                
                nucmer_results = self.run_nucmer_parallel_for_all_chromosomes(
                    chromosome_files, reference_fasta
                )
            
            successful_chromosomes = [
                chrom for chrom, res in nucmer_results.items() 
                if res.get("success", False) and res.get("coords_file") and os.path.exists(res["coords_file"])
            ]
            
            if not successful_chromosomes:
                self.log("Warning: All chromosome synteny analyses failed", "warning")
                successful_chromosomes = list(chromosome_files.keys())
            
            success_count = len(successful_chromosomes)
            total_count = len(chromosomes)
            self.log(f"Total {success_count}/{total_count} chromosomes analyzed successfully", 
                    "info" if success_count > 0 else "warning")
            
            if success_count < total_count:
                failed_chroms = [c for c in chromosomes if c not in successful_chromosomes]
                self.log(f"Failed chromosomes: {', '.join(failed_chroms)}", "warning")
                
                for chrom in failed_chroms:
                    res = nucmer_results.get(chrom, {})
                    if "error" in res:
                        self.log(f"  {chrom}: {res['error']}", "warning")
            
            final_chromosome_files = {}
            chromosome_results = {}
            
            for chrom in successful_chromosomes:
                chrom_file = chromosome_files.get(chrom)
                if not chrom_file or not os.path.exists(chrom_file):
                    self.log(f"Skipping chromosome {chrom}: file does not exist", "warning")
                    final_chromosome_files[chrom] = chrom_file
                    continue
                
                gap_positions = gap_analysis.get("gap_positions_by_chromosome", {}).get(chrom, [])
                nucmer_result = nucmer_results.get(chrom, {})
                
                task = (chrom, chrom_file, reference_fasta, gap_positions, nucmer_result)
                result = self.process_single_chromosome(task)
                
                final_chromosome_files[chrom] = result["final_file"] or chrom_file
                chromosome_results[chrom] = result
                self.status[chrom] = result.get("success", False)
            
            if final_chromosome_files:
                merged_file = self.merge_all_chromosomes(final_chromosome_files, output_fasta)
            else:
                self.log("Warning: No chromosomes processed, using original genome", "warning")
                merged_file = query_fasta
            
            final_analysis = self.analyze_final_genome(merged_file)
            
            total_time = time.time() - start_time
            
            report = {
                "status": "success",
                "input_files": {
                    "query_fasta": query_fasta,
                    "reference_fasta": reference_fasta,
                    "coords_file": coords_file
                },
                "output_files": {
                    "fixed_genome": merged_file,
                    "work_dir": self.work_dir
                },
                "gap_history": self.gap_history,
                "chromosomes_processed": {
                    "total": len(chromosomes),
                    "successful": len(successful_chromosomes),
                    "failed": len(chromosomes) - len(successful_chromosomes)
                },
                "chromosome_results": chromosome_results,
                "processing_time": total_time,
                "parallel_config": {
                    "total_threads": self.total_threads,
                    "successful_chromosomes": successful_chromosomes,
                    "use_existing_coords": coords_file is not None
                },
                "repair_mode": self.repair_mode,
                "analysis_parameters": {
                    "gap_length": self.gap_length,
                    "max_search_distance": self.max_search_distance,
                    "search_step": self.search_step
                },
                "final_analysis": final_analysis
            }
            
            report_dir = os.path.dirname(output_fasta)
            report_file = os.path.join(report_dir, "gap_fix_report.json")
            with open(report_file, "w") as f:
                json.dump(report, f, indent=2, ensure_ascii=False)
            
            self.log("\n" + "=" * 60, "info")
            self.log("Genome gap repair pipeline completed!", "info")
            self.log(f"Repair mode: {self.repair_mode}", "info")
            self.log(f"Final gap length: {self.gap_length} bp", "info")
            self.log(f"Using existing coords file: {coords_file is not None}", "info")
            self.log("=" * 60, "info")
            
            if "initial" in self.gap_history and "final" in self.gap_history:
                initial = self.gap_history["initial"]
                final = self.gap_history["final"]
                
                self.log(f"Initial gap count: {initial['total_gaps']}", "info")
                self.log(f"Final gap count: {final['total_gaps']}", "info")
                self.log(f"Repaired gaps: {initial['total_gaps'] - final['total_gaps']}", "info")
                self.log(f"Initial gap length: {initial['total_length']:,} bp", "info")
                self.log(f"Final gap length: {final['total_length']:,} bp", "info")
                self.log(f"Repaired gap length: {initial['total_length'] - final['total_length']:,} bp", "info")
                self.log(f"Initial gap coverage: {initial['gap_coverage']:.4f}%", "info")
                self.log(f"Final gap coverage: {final['gap_coverage']:.4f}%", "info")
                self.log(f"Repair mode: {self.repair_mode}", "info")
                self.log(f"Final gap length: {self.gap_length} bp", "info")
                self.log(f"Maximum search distance: {self.max_search_distance:,} bp", "info")
                self.log(f"Search step: {self.search_step:,} bp", "info")
                self.log(f"Using existing coords file: {coords_file is not None}", "info")
                self.log(f"Processing time: {total_time:.1f}s", "info")
                self.log(f"Output file: {merged_file}", "info")
                self.log(f"Detailed report: {report_file}", "info")
            
            return report
            
        except Exception as e:
            self.log(f"Pipeline execution failed: {e}", "error")
            import traceback
            traceback.print_exc()
            
            return {
                "status": "error",
                "error": str(e),
                "repair_mode": self.repair_mode,
                "processing_time": time.time() - start_time
            }
    
    def cleanup(self):
        self.log("Cleaning controller resources...", "info")
        self._cleanup_work_dir()

class GenomeGapFixerAdapter:
    
    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or {}
        self.fixer = None
        self.result = None
        self.status = "initialized"
        
        if not ALL_APIS_AVAILABLE:
            raise ImportError("Cannot import all sub-script APIs")
    
    def validate_config(self, config: Dict[str, Any]) -> Tuple[bool, str]:
        required_params = ['query_fasta']
        
        for param in required_params:
            if param not in config:
                return False, f"Missing required parameter: {param}"
            
            if param in ['query_fasta']:
                if not os.path.exists(config[param]):
                    return False, f"File does not exist: {config[param]}"
        
        has_coords = 'coords_file' in config and config['coords_file']
        has_reference = 'reference_fasta' in config and config['reference_fasta']
        
        if not has_coords and not has_reference:
            return False, "Must provide coords_file or reference_fasta"
        
        if has_coords and not os.path.exists(config['coords_file']):
            return False, f"Coords file does not exist: {config['coords_file']}"
        
        if has_reference and not os.path.exists(config['reference_fasta']):
            return False, f"Reference genome file does not exist: {config['reference_fasta']}"
        
        if 'threads' in config and config['threads'] <= 0:
            return False, "Threads must be greater than 0"
        
        repair_mode = config.get('repair_mode', 'aggressive')
        if repair_mode not in ['conservative', 'aggressive']:
            return False, f"Repair mode must be 'conservative' or 'aggressive', current: {repair_mode}"
        
        if 'gap_length' in config and config['gap_length'] <= 0:
            return False, "Gap length must be greater than 0"
        
        if 'max_search_distance' in config and config['max_search_distance'] <= 0:
            return False, "Maximum search distance must be greater than 0"
        
        if 'search_step' in config and config['search_step'] <= 0:
            return False, "Search step must be greater than 0"
        
        return True, ""
    
    def determine_output_path(self, user_output: str) -> str:
        default_filename = "fixed_genome.fasta"
        
        if not user_output:
            return default_filename
        
        if user_output.endswith('/') or user_output.endswith('\\') or os.path.isdir(user_output):
            output_dir = user_output.rstrip('/\\')
            os.makedirs(output_dir, exist_ok=True)
            return os.path.join(output_dir, default_filename)
        
        base_name = os.path.basename(user_output)
        if '.' not in base_name or base_name.endswith('.'):
            os.makedirs(user_output, exist_ok=True)
            return os.path.join(user_output, default_filename)
        
        parent_dir = os.path.dirname(user_output)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        return user_output
    
    def run(self, config: Dict[str, Any] = None) -> Dict[str, Any]:
        if config is not None:
            self.config = config
        
        is_valid, error_msg = self.validate_config(self.config)
        if not is_valid:
            self.status = "config_error"
            return {
                "status": "error",
                "error": error_msg,
                "module": "genome_gap_fixer",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
        
        try:
            self.status = "running"
            
            repair_mode = self.config.get('repair_mode', 'aggressive')
            gap_length = self.config.get('gap_length', 100)
            max_search_distance = self.config.get('max_search_distance', 500000)
            search_step = self.config.get('search_step', 100000)
            
            output_fasta = self.config.get('output_fasta', 'fixed_genome.fasta')
            output_fasta = self.determine_output_path(output_fasta)
            
            self.config['output_fasta'] = output_fasta
            
            self.fixer = GenomeGapFixer(
                total_threads=self.config.get('threads', 4),
                verbose=self.config.get('verbose', True),
                log_file=self.config.get('log_file'),
                repair_mode=repair_mode,
                gap_length=gap_length,
                max_search_distance=max_search_distance,
                search_step=search_step
            )
            
            if 'work_dir' in self.config:
                self.fixer.work_dir = self.config['work_dir']
            
            self.result = self.fixer.run_complete_pipeline(
                query_fasta=self.config['query_fasta'],
                reference_fasta=self.config.get('reference_fasta'),
                coords_file=self.config.get('coords_file'),
                output_fasta=output_fasta
            )
            
            standardized_result = self._standardize_result(self.result)
            
            self.status = "completed"
            return standardized_result
            
        except Exception as e:
            self.status = "failed"
            error_result = {
                "status": "error",
                "error": str(e),
                "module": "genome_gap_fixer",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
            
            if hasattr(e, '__traceback__'):
                import traceback
                error_result["traceback"] = traceback.format_exc()
            
            return error_result
    
    def _standardize_result(self, result: Dict[str, Any]) -> Dict[str, Any]:
        if result.get("status") == "error":
            return {
                "status": "error",
                "error": result.get("error", "Unknown error"),
                "module": "genome_gap_fixer",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
        
        if result.get("status") == "no_gaps":
            return {
                "status": "success",
                "message": "Genome has no gaps, no processing needed",
                "module": "genome_gap_fixer",
                "output_file": result.get("output_file", ""),
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
        
        standardized = {
            "status": "success",
            "module": "genome_gap_fixer",
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "pipeline_version": "simplified_fixer",
            "input_files": {
                "query_fasta": result.get("input_files", {}).get("query_fasta", ""),
                "reference_fasta": result.get("input_files", {}).get("reference_fasta", ""),
                "coords_file": result.get("input_files", {}).get("coords_file", "")
            },
            "output_files": {
                "final_genome": result.get("output_files", {}).get("fixed_genome", ""),
                "report_file": "gap_fix_report.json",
                "work_dir": result.get("output_files", {}).get("work_dir", "")
            },
            "statistics": {
                "gap_history": result.get("gap_history", {}),
                "chromosomes_processed": result.get("chromosomes_processed", {}),
                "chromosome_results": result.get("chromosome_results", {}),
                "repair_mode": result.get("repair_mode", "aggressive"),
                "use_existing_coords": result.get("parallel_config", {}).get("use_existing_coords", False),
                "analysis_parameters": result.get("analysis_parameters", {}),
                "final_analysis": result.get("final_analysis", {})
            },
            "processing_info": {
                "processing_time": result.get("processing_time", 0),
                "threads_used": result.get("parallel_config", {}).get("total_threads", 0),
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
        }
        
        gap_history = result.get("gap_history", {})
        if "initial" in gap_history and "final" in gap_history:
            initial = gap_history["initial"]
            final = gap_history["final"]
            
            standardized["statistics"]["repair_results"] = {
                "gaps_analyzed": initial.get("total_gaps", 0),
                "gaps_remaining": final.get("total_gaps", 0),
                "gaps_fixed": initial.get("total_gaps", 0) - final.get("total_gaps", 0),
                "gaps_fixed_percent": ((initial.get("total_gaps", 0) - final.get("total_gaps", 0)) / 
                                     max(initial.get("total_gaps", 1), 1) * 100),
                "length_fixed": initial.get("total_length", 0) - final.get("total_length", 0),
                "length_fixed_percent": ((initial.get("total_length", 0) - final.get("total_length", 0)) / 
                                       max(initial.get("total_length", 1), 1) * 100),
                "coverage_reduction": initial.get("gap_coverage", 0) - final.get("gap_coverage", 0)
            }
        
        if self.config.get('include_raw_result', False):
            standardized["raw_result"] = result
        
        return standardized
    
    def get_status(self) -> str:
        return self.status
    
    def get_result(self) -> Dict[str, Any]:
        return self.result
    
    def cleanup(self) -> bool:
        try:
            if self.fixer:
                self.fixer.cleanup()
            return True
        except Exception as e:
            print(f"Cleanup failed: {e}")
            return False
    
    def get_default_config(self) -> Dict[str, Any]:
        return {
            "query_fasta": "",
            "reference_fasta": "",
            "coords_file": "",
            "threads": 4,
            "output_fasta": "./gap_fix_output",
            "work_dir": "gap_fix_workdir",
            "repair_mode": "aggressive",
            "gap_length": 100,
            "max_search_distance": 500000,
            "search_step": 100000,
            "verbose": True,
            "log_file": None,
            "include_raw_result": False
        }
    
    def get_module_info(self) -> Dict[str, Any]:
        return {
            "name": "GenomeGapFixerAdapter",
            "version": "simplified_fixer",
            "description": "Genome gap repair controller adapter - integrates 3 sub-scripts",
            "capabilities": [
                "Gap detection and chromosome sorting",
                "Parallel nucmer synteny analysis",
                "Error region analysis and repair",
                "Parallel chromosome processing",
                "Result merging and report generation"
            ],
            "dependencies": ["find_gaps", "run_nucmer_parallel", "gap_analyzer"],
            "repair_modes": ["conservative", "aggressive"],
            "removed": ["gap_filler", "extract_gap_patches"],
            "coords_support": "Directly specify existing coords file via -coords parameter",
            "output_support": "-o parameter supports folder (automatically generates fixed_genome.fasta) or specific filename",
            "analysis_parameters": {
                "gap_length": "Final gap length (default 100bp)",
                "max_search_distance": "Maximum search distance (default 500kb)",
                "search_step": "Search step (default 100kb)"
            }
        }

def get_fixer_api(config: Dict[str, Any] = None):
    return GenomeGapFixerAdapter(config)

def main():
    parser = argparse.ArgumentParser(
        description='Genome gap repair controller - Complete pipeline integrating 3 sub-scripts (parallel version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Basic usage, output to directory (automatically generates fixed_genome.fasta)
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -t 8 -o ./results
  
  # Output to specific directory
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o /path/to/output
  
  # Output to specific file
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o custom_name.fasta
  
  # Full pipeline, run nucmer analysis
  python gap_fixer.py -q my_genome.fasta -c reference.fasta -t 8 -o ./output
  
  # Specify gap length and search parameters
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o results --gap-length 200 --search 1000000 --step 200000
  
  # Use conservative mode
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o ./output --repair-mode conservative
  
  # Use aggressive mode (default)
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o ./output --repair-mode aggressive
  
  # Automatically set based on CPU cores
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o results -t $(nproc)
  
  # Quiet mode
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o output -t 8 --quiet
  
  # Keep work directory
  python gap_fixer.py -q my_genome.fasta -coords existing.coords -o output -t 8 --keep-workdir
  
  # Adapter mode
  python gap_fixer.py --adapter-mode --config config.json
  
  # Show module information
  python gap_fixer.py --module-info
        """
    )
    
    parser.add_argument('-q', '--query', dest='query_fasta', required=False,
                       help='Query genome FASTA file (genome needing gap repair)')
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-coords', dest='coords_file', required=False,
                      help='Alignment coordinates file (coords format), if provided skips nucmer analysis')
    group.add_argument('-c', '--contigs', dest='reference_fasta', required=False,
                      help='Reference genome contig file (reference sequences for nucmer analysis)')
    
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Total threads (for parallel processing), default: 4')
    
    parser.add_argument('--repair-mode', choices=['conservative', 'aggressive'], 
                       default='aggressive',
                       help='Repair mode: conservative | aggressive (default)')
    
    parser.add_argument('--gap-length', type=int, default=100,
                       help='Final gap length (default 100bp)')
    parser.add_argument('--search', type=int, default=500000, dest='max_search_distance',
                       help='Maximum search distance (bp), default 500kb')
    parser.add_argument('--step', type=int, default=100000, dest='search_step',
                       help='Search step (bp), default 100kb')
    
    parser.add_argument('-o', '--output', default='./gap_fix_output/fixed_genome.fasta',
                       help='Output path. Can be directory (will generate fixed_genome.fasta inside) or specific filename. Default: ./gap_fix_output/fixed_genome.fasta')
    parser.add_argument('--work-dir', default='gap_fix_workdir',
                       help='Working directory, default: gap_fix_workdir')
    parser.add_argument('--keep-workdir', action='store_true',
                       help='Keep working directory (default will clean up)')
    parser.add_argument('--verbose', action='store_true', help='Show verbose output')
    parser.add_argument('--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--log-file', help='Log file path')
    
    parser.add_argument('--adapter-mode', action='store_true',
                       help='Use adapter mode')
    parser.add_argument('--config', type=str,
                       help='JSON configuration file path (adapter mode)')
    parser.add_argument('--module-info', action='store_true',
                       help='Show module information')
    
    args = parser.parse_args()
    
    if args.module_info:
        adapter = GenomeGapFixerAdapter()
        info = adapter.get_module_info()
        print(json.dumps(info, indent=2, ensure_ascii=False))
        sys.exit(0)
    
    if args.adapter_mode:
        if args.config:
            try:
                with open(args.config, 'r') as f:
                    config = json.load(f)
            except Exception as e:
                print(f"Cannot load configuration file: {e}")
                sys.exit(1)
        else:
            if not args.query_fasta:
                print("Adapter mode requires -q/--query parameter")
                sys.exit(1)
                
            config = {
                'query_fasta': args.query_fasta,
                'reference_fasta': args.reference_fasta,
                'coords_file': args.coords_file,
                'threads': args.threads,
                'output_fasta': args.output,
                'work_dir': args.work_dir,
                'repair_mode': args.repair_mode,
                'gap_length': args.gap_length,
                'max_search_distance': args.max_search_distance,
                'search_step': args.search_step,
                'verbose': args.verbose and not args.quiet,
                'log_file': args.log_file,
                'include_raw_result': False
            }
        
        adapter = GenomeGapFixerAdapter(config)
        result = adapter.run()
        
        print(json.dumps(result, indent=2, ensure_ascii=False))
        
        if result.get('status') == 'success':
            sys.exit(0)
        else:
            sys.exit(1)
    
    else:
        if not args.query_fasta:
            print("Error: Please provide query genome file (-q/--query)", file=sys.stderr)
            sys.exit(1)
        
        if not os.path.exists(args.query_fasta):
            print(f"Error: Query genome file does not exist: {args.query_fasta}", file=sys.stderr)
            sys.exit(1)
        
        if not args.coords_file and not args.reference_fasta:
            print("Error: Must provide -coords or -c/--contigs parameter", file=sys.stderr)
            print("  Use -coords to provide existing coords file (skips nucmer analysis)", file=sys.stderr)
            print("  Use -c to provide reference genome file (runs nucmer analysis)", file=sys.stderr)
            sys.exit(1)
        
        if args.coords_file and not os.path.exists(args.coords_file):
            print(f"Error: Coords file does not exist: {args.coords_file}", file=sys.stderr)
            sys.exit(1)
        
        if args.reference_fasta and not os.path.exists(args.reference_fasta):
            print(f"Error: Reference genome file does not exist: {args.reference_fasta}", file=sys.stderr)
            sys.exit(1)
        
        if args.threads < 1:
            print(f"Warning: Threads cannot be less than 1, using default 4", file=sys.stderr)
            args.threads = 4
        
        if args.gap_length < 1:
            print(f"Warning: Gap length cannot be less than 1, using default 100", file=sys.stderr)
            args.gap_length = 100
        
        if args.max_search_distance < 1:
            print(f"Warning: Maximum search distance cannot be less than 1, using default 500000", file=sys.stderr)
            args.max_search_distance = 500000
        
        if args.search_step < 1:
            print(f"Warning: Search step cannot be less than 1, using default 100000", file=sys.stderr)
            args.search_step = 100000
        
        max_cpu = multiprocessing.cpu_count()
        if args.threads > max_cpu * 2:
            print(f"Warning: Specified threads ({args.threads}) exceeds 2x system CPU cores ({max_cpu})", file=sys.stderr)
            print(f"Suggested using {max_cpu} threads", file=sys.stderr)
        
        if not ALL_APIS_AVAILABLE:
            print("Error: Cannot import all sub-script APIs", file=sys.stderr)
            print("Please ensure these scripts are in the same directory:", file=sys.stderr)
            print("  - find_gaps.py", file=sys.stderr)
            print("  - run_nucmer_parallel.py", file=sys.stderr)
            print("  - gap_analyzer.py", file=sys.stderr)
            sys.exit(1)
        
        verbose = args.verbose and not args.quiet
        
        def determine_output_path(user_output: str, default_filename: str = "fixed_genome.fasta") -> str:
            if not user_output:
                return default_filename
            
            if user_output.endswith('/') or user_output.endswith('\\') or os.path.isdir(user_output):
                output_dir = user_output.rstrip('/\\')
                os.makedirs(output_dir, exist_ok=True)
                return os.path.join(output_dir, default_filename)
            
            base_name = os.path.basename(user_output)
            if '.' not in base_name or base_name.endswith('.'):
                os.makedirs(user_output, exist_ok=True)
                return os.path.join(user_output, default_filename)
            
            parent_dir = os.path.dirname(user_output)
            if parent_dir:
                os.makedirs(parent_dir, exist_ok=True)
            return user_output
        
        args.output = determine_output_path(args.output)
        
        if args.work_dir == 'gap_fix_workdir':
            output_dir = os.path.dirname(args.output)
            args.work_dir = os.path.join(output_dir, "gap_fix_workdir")
        
        fixer = GenomeGapFixer(
            total_threads=args.threads,
            verbose=verbose, 
            log_file=args.log_file,
            repair_mode=args.repair_mode,
            gap_length=args.gap_length,
            max_search_distance=args.max_search_distance,
            search_step=args.search_step
        )
        fixer.work_dir = args.work_dir
        
        try:
            report = fixer.run_complete_pipeline(
                query_fasta=args.query_fasta,
                reference_fasta=args.reference_fasta,
                coords_file=args.coords_file,
                output_fasta=args.output
            )
            
            if not args.keep_workdir:
                fixer._cleanup_work_dir()
            
            status = report.get("status", "unknown")
            repair_mode = report.get("repair_mode", "aggressive")
            use_existing_coords = report.get("parallel_config", {}).get("use_existing_coords", False)
            analysis_params = report.get("analysis_parameters", {})
            
            if status == "success":
                print(f"\n✓ Gap repair completed successfully!")
                print(f"Repair mode: {repair_mode}")
                print(f"Using existing coords file: {use_existing_coords}")
                print(f"Final gap length: {analysis_params.get('gap_length', 100)} bp")
                print(f"Maximum search distance: {analysis_params.get('max_search_distance', 500000):,} bp")
                print(f"Search step: {analysis_params.get('search_step', 100000):,} bp")
                print(f"Output file: {report['output_files']['fixed_genome']}")
                print(f"Threads used: {args.threads}")
                
                if "gap_history" in report:
                    initial = report["gap_history"].get("initial", {})
                    final = report["gap_history"].get("final", {})
                    
                    if initial and final:
                        print(f"Repair results: {initial.get('total_gaps', 0)} → {final.get('total_gaps', 0)} gaps")
                        print(f"Repaired length: {initial.get('total_length', 0):,} → {final.get('total_length', 0):,} bp")
                
                if "chromosomes_processed" in report:
                    stats = report["chromosomes_processed"]
                    print(f"Chromosomes processed: {stats.get('successful', 0)}/{stats.get('total', 0)} successful")
                
                sys.exit(0)
            elif status == "no_gaps":
                print(f"\n✓ Genome has no gaps, no processing needed")
                sys.exit(0)
            else:
                print(f"\n✗ Gap repair failed: {report.get('error', 'Unknown error')}", file=sys.stderr)
                sys.exit(1)
                
        except KeyboardInterrupt:
            print("\nProcess interrupted by user", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Fixer execution failed: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            sys.exit(1)
        finally:
            fixer.cleanup()

if __name__ == "__main__":
    main()