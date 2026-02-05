#!/usr/bin/env python3
"""
GAP Complete Process Controller - Modular Enhanced Version
Can be used as standalone CLI tool or imported as Python module
Strategy:
1. No patches extracted → Run fixer → Re-run patcher → Fill gaps
2. Patches extracted but all failed → Run fixer → Reuse patches
3. Partially successful patches → Use patcher results directly
"""

import os
import sys
import json
import argparse
import subprocess
import shutil
import time
import re
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Set
from datetime import datetime

SCRIPT_DIR = Path(__file__).parent.absolute()
sys.path.insert(0, str(SCRIPT_DIR))

try:
    from gap_fixer import GenomeGapFixerAdapter, GenomeGapFixer
    from gap_patcher_controller import GapPatcherAdapter, GenomeGapPatcherController
    MODULES_AVAILABLE = True
except ImportError as e:
    MODULES_AVAILABLE = False
    print(f"⚠️  Warning: Cannot import modules: {e}")

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("⚠️  Warning: BioPython not available, some features limited")

class GapCompleteControllerModuleEnhanced:
    
    def __init__(self, args=None, config_dict=None):
        if args is not None:
            self.args = args
        elif config_dict is not None:
            self.args = self._dict_to_args(config_dict)
        else:
            raise ValueError("Must provide args or config_dict")
        
        self.script_dir = SCRIPT_DIR
        
        self.query_fasta = Path(self.args.query_fasta).absolute()
        self.reference_fasta = Path(self.args.reference_fasta).absolute()
        self.coords_file = Path(self.args.coords_file).absolute() if hasattr(self.args, 'coords_file') and self.args.coords_file else None
        
        self.threads = getattr(self.args, 'threads', 8)
        self.output_dir = Path(getattr(self.args, 'output_dir', 'gap_complete_module_enhanced_results')).absolute()
        self.repair_mode = getattr(self.args, 'repair_mode', 'aggressive')
        self.keep_temp = getattr(self.args, 'keep_temp', False)
        
        self.gap_length = getattr(self.args, 'gap_length', 100)
        self.max_search_distance = getattr(self.args, 'max_search_distance', 500000)
        self.search_step = getattr(self.args, 'search_step', 100000)
        
        self.min_gap_size = getattr(self.args, 'min_gap_size', 100)
        
        self.quiet = getattr(self.args, 'quiet', False)
        
        self.chromosome_strategies = {}
        self.partially_successful_chromosomes = []
        self.failed_completely_chromosomes = []
        self.no_patch_chromosomes = []
        self.chromosome_info = {}
        self.chromosome_patches = {}
        self.fixed_chromosomes = {}
        self.gap_details = {}
        self.coordinate_mapping = {}
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self._setup_logging()
    
    def _dict_to_args(self, config_dict):
        class Args:
            def __init__(self, config):
                defaults = {
                    'threads': 8,
                    'output_dir': 'gap_complete_module_enhanced_results',
                    'repair_mode': 'aggressive',
                    'gap_length': 100,
                    'max_search_distance': 500000,
                    'search_step': 100000,
                    'min_gap_size': 100,
                    'keep_temp': False,
                    'quiet': False
                }
                defaults.update(config)
                for key, value in defaults.items():
                    setattr(self, key, value)
        
        return Args(config_dict)
    
    def _setup_logging(self):
        self.log_file = self.output_dir / "complete_controller_module_enhanced.log"
        print(f"Log file: {self.log_file}")
    
    def log(self, message, level="info"):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        if level == "error":
            log_msg = f"[{timestamp}] ❌ {message}"
            prefix = "❌"
        elif level == "warning":
            log_msg = f"[{timestamp}] ⚠️  {message}"
            prefix = "⚠️"
        elif level == "success":
            log_msg = f"[{timestamp}] ✅ {message}"
            prefix = "✅"
        else:
            log_msg = f"[{timestamp}] {message}"
            prefix = "ℹ️"
        
        if not self.quiet:
            if level in ["error", "warning", "success"]:
                print(f"{prefix} {message}")
            else:
                print(f"  {message}")
        
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(log_msg + "\n")
    
    def run_patcher_analysis(self):
        self.log("="*60, "info")
        self.log("Running patcher analysis phase", "info")
        self.log("="*60, "info")
        
        if not self._run_patcher_module():
            return {
                "status": "error", 
                "error": "Patcher execution failed",
                "output_dir": str(self.output_dir),
                "log_file": str(self.log_file)
            }
        
        self._collect_all_chromosome_patches()
        
        return {
            "status": "success",
            "chromosome_strategies": self.chromosome_strategies,
            "chromosome_info": self.chromosome_info,
            "categories": {
                "no_patch": self.no_patch_chromosomes,
                "failed_completely": self.failed_completely_chromosomes,
                "partial_success": self.partially_successful_chromosomes
            },
            "patch_counts": {k: len(v) for k, v in self.chromosome_patches.items()},
            "output_dir": str(self.output_dir),
            "log_file": str(self.log_file)
        }
    
    def run_for_chromosome(self, chrom_name, strategy=None):
        if not self.chromosome_strategies and strategy is None:
            analysis = self.run_patcher_analysis()
            if analysis["status"] != "success":
                return analysis
        
        if strategy is None:
            strategy = self.chromosome_strategies.get(chrom_name, "unknown")
        
        self.log(f"Processing chromosome {chrom_name}, strategy: {strategy}", "info")
        
        if strategy == "no_patch_fixer":
            return self._process_no_patch_chromosome(chrom_name)
        elif strategy == "full_fixer_reuse_patches":
            return self._process_failed_completely_chromosome(chrom_name)
        elif strategy == "partial_success_no_fixer":
            return self._process_partial_success_chromosome(chrom_name)
        elif strategy == "no_gaps":
            return {
                "status": "success",
                "chromosome": chrom_name,
                "strategy": strategy,
                "message": "Chromosome has no GAPs, no processing needed"
            }
        else:
            return {
                "status": "error",
                "error": f"Unknown strategy: {strategy}",
                "chromosome": chrom_name
            }
    
    def run_complete_workflow(self):
        return self.run_complete_pipeline()
    
    def get_results(self):
        return {
            "fixed_chromosomes": {k: str(v) for k, v in self.fixed_chromosomes.items()},
            "strategies": self.chromosome_strategies,
            "categories": {
                "no_patch": self.no_patch_chromosomes,
                "failed_completely": self.failed_completely_chromosomes,
                "partial_success": self.partially_successful_chromosomes
            },
            "output_dir": str(self.output_dir),
            "log_file": str(self.log_file),
            "config": {
                "query_fasta": str(self.query_fasta),
                "reference_fasta": str(self.reference_fasta),
                "threads": self.threads,
                "repair_mode": self.repair_mode
            }
        }
    
    def _run_patcher_module(self):
        self.log("Step 1: Running gap_patcher module", "info")
        
        try:
            patcher_config = {
                "query_fasta": str(self.query_fasta),
                "reference_fasta": str(self.reference_fasta),
                "coords_file": str(self.coords_file) if self.coords_file else "",
                "threads": self.threads,
                "output_dir": str(self.output_dir / "patcher_results"),
                "min_gap_size": self.min_gap_size,
                "keep_temp": True,
                "verbose": not self.quiet,
                "include_raw_result": True
            }
            
            self.log(f"Creating patcher adapter...", "info")
            start_time = time.time()
            
            patcher_adapter = GapPatcherAdapter(patcher_config)
            result = patcher_adapter.run()
            elapsed_time = time.time() - start_time
            
            if result.get("status") == "success":
                self.log(f"Patcher execution successful, time: {elapsed_time:.1f}s", "success")
                
                self._parse_patcher_result(result)
                return True
            else:
                error_msg = result.get("error", "Unknown error")
                self.log(f"Patcher execution failed: {error_msg}", "error")
                return False
                
        except Exception as e:
            self.log(f"Patcher module exception: {e}", "error")
            import traceback
            self.log(traceback.format_exc(), "error")
            return False
    
    def _parse_patcher_result(self, result):
        try:
            raw_result = result.get("raw_result", {})
            if not raw_result:
                raw_result = result
            
            stage_results = raw_result.get("stage_results", {})
            gap_patching = stage_results.get("gap_patching", {})
            
            if not gap_patching:
                self.log("No chromosome filling results found", "warning")
                stats = raw_result.get("statistics", {})
                if stats:
                    gap_finding = stats.get("gap_finding", {})
                    chromosomes_with_gaps = gap_finding.get("chromosomes_with_gaps", {})
                    
                    for chrom_name, gap_count in chromosomes_with_gaps.items():
                        if gap_count > 0:
                            gap_patching[chrom_name] = {
                                "total_gaps": gap_count,
                                "successful_patches": 0,
                                "failed_patches": gap_count,
                                "status": "all_failed"
                            }
            
            if not gap_patching:
                self.log("Cannot parse patcher results", "error")
                return
            
            self.no_patch_chromosomes = []
            self.failed_completely_chromosomes = []
            self.partially_successful_chromosomes = []
            
            for chrom_name, chrom_result in gap_patching.items():
                total_gaps = chrom_result.get("total_gaps", 0)
                successful_patches = chrom_result.get("successful_patches", 0)
                failed_patches = chrom_result.get("failed_patches", 0)
                extracted_patches = successful_patches + failed_patches
                
                self.chromosome_info[chrom_name] = {
                    "original_gaps": total_gaps,
                    "extracted_patches": extracted_patches,
                    "patcher_successful": successful_patches,
                    "patcher_failed": failed_patches,
                    "patcher_success_rate": (successful_patches / total_gaps * 100) if total_gaps > 0 else 0,
                    "patcher_result": chrom_result,
                    "status": chrom_result.get("status", "unknown")
                }
                
                self._collect_gap_details(chrom_name, chrom_result)
                
                strategy = self._determine_chromosome_strategy(
                    chrom_name, total_gaps, extracted_patches, successful_patches
                )
                
                self.chromosome_strategies[chrom_name] = strategy
                
                if strategy == "no_patch_fixer":
                    self.no_patch_chromosomes.append(chrom_name)
                    self.log(f"  Chromosome {chrom_name}: No patches → Full fixer + re-patcher", "warning")
                elif strategy == "full_fixer_reuse_patches":
                    self.failed_completely_chromosomes.append(chrom_name)
                    self.log(f"  Chromosome {chrom_name}: Patches all failed → Fixer + reuse patches", "warning")
                elif strategy == "partial_success_no_fixer":
                    self.partially_successful_chromosomes.append(chrom_name)
                    self.log(f"  Chromosome {chrom_name}: Partial success → No fixer needed", "info")
                elif strategy == "no_gaps":
                    self.log(f"  Chromosome {chrom_name}: No gaps, no processing needed", "info")
                else:
                    self.log(f"  Chromosome {chrom_name}: Unknown strategy", "error")
            
            self.log(f"\nChromosome strategy statistics:", "info")
            self.log(f"  No patch chromosomes: {len(self.no_patch_chromosomes)}", "info")
            self.log(f"  All failed chromosomes: {len(self.failed_completely_chromosomes)}", "info")
            self.log(f"  Partial success chromosomes: {len(self.partially_successful_chromosomes)}", "info")
            
            self._save_chromosome_analysis()
            
        except Exception as e:
            self.log(f"Failed to parse patcher results: {e}", "error")
            import traceback
            self.log(traceback.format_exc(), "error")
    
    def _collect_all_chromosome_patches(self):
        for chrom_name in list(self.chromosome_info.keys()):
            self._collect_chromosome_patches(chrom_name)
    
    def _collect_chromosome_patches(self, chrom_name: str):
        chrom_dir = self._find_chromosome_dir(chrom_name)
        if not chrom_dir:
            return
        
        patches_dir = chrom_dir / "gap_patches"
        if not patches_dir.exists():
            patches_dir = chrom_dir / "patches"
        
        if not patches_dir.exists():
            return
        
        patches = []
        for patch_file in patches_dir.glob("*.fasta"):
            match = re.search(r'gap_(\d+)_patch', patch_file.name)
            if match:
                gap_position = int(match.group(1))
                patch_info = {
                    "original_position": gap_position,
                    "file": patch_file,
                    "status": "available"
                }
                patches.append(patch_info)
        
        if patches:
            patches.sort(key=lambda x: x["original_position"])
            self.chromosome_patches[chrom_name] = patches
            self.log(f"  Chromosome {chrom_name}: Collected {len(patches)} patch files", "info")
    
    def _collect_gap_details(self, chrom_name: str, chrom_result: Dict):
        chrom_dir = self._find_chromosome_dir(chrom_name)
        if not chrom_dir:
            return
        
        gap_report_file = chrom_dir / "gap_analysis" / "gap_report.json"
        if not gap_report_file.exists():
            gap_report_file = chrom_dir / "gap_report.json"
        
        if gap_report_file.exists():
            try:
                with open(gap_report_file, 'r', encoding='utf-8') as f:
                    gap_report = json.load(f)
                
                self.gap_details[chrom_name] = gap_report.get("gaps", [])
                if self.gap_details[chrom_name]:
                    self.log(f"  Chromosome {chrom_name}: Collected details for {len(self.gap_details[chrom_name])} gaps", "info")
            except Exception as e:
                self.log(f"  Failed to read gap report: {e}", "warning")
    
    def _save_chromosome_analysis(self):
        analysis_file = self.output_dir / "chromosome_analysis.json"
        analysis_data = {
            "chromosome_strategies": self.chromosome_strategies,
            "chromosome_info": self.chromosome_info,
            "categories": {
                "no_patch": self.no_patch_chromosomes,
                "failed_completely": self.failed_completely_chromosomes,
                "partial_success": self.partially_successful_chromosomes
            },
            "chromosome_patches": {k: len(v) for k, v in self.chromosome_patches.items()},
            "gap_details": {k: len(v) for k, v in self.gap_details.items()}
        }
        
        with open(analysis_file, 'w', encoding='utf-8') as f:
            json.dump(analysis_data, f, indent=2, ensure_ascii=False)
        
        self.log(f"Chromosome analysis saved: {analysis_file}", "info")
    
    def _process_no_patch_chromosome(self, chrom_name: str) -> Dict[str, Any]:
        chrom_dir = self._find_chromosome_dir(chrom_name)
        chrom_fasta = self._find_chromosome_fasta(chrom_dir, chrom_name)
        coords_file = self._find_coords_file(chrom_dir, chrom_name)
        
        if not all([chrom_dir, chrom_fasta, coords_file]):
            return {
                "status": "error",
                "error": f"Cannot find required files for chromosome {chrom_name}",
                "chromosome": chrom_name
            }
        
        fixed_file = self._run_gap_fixer_module(chrom_dir, chrom_fasta, coords_file, chrom_name)
        if not fixed_file:
            return {
                "status": "error",
                "error": f"gap_fixer repair failed",
                "chromosome": chrom_name
            }
        
        self.log(f"  ✓ gap_fixer repair completed", "success")
        
        self.log(f"  Re-running patcher on fixed chromosome to extract patches", "info")
        new_patches_dir = self._rerun_patcher_for_fixed_chromosome(chrom_dir, fixed_file, chrom_name)
        
        if new_patches_dir and new_patches_dir.exists():
            self.log(f"  Applying newly extracted patches", "info")
            patched_file = self._apply_new_patches_to_fixed(chrom_dir, fixed_file, new_patches_dir, chrom_name)
            
            if patched_file:
                self.log(f"  ✓ New patch application completed: {patched_file.name}", "success")
                self.fixed_chromosomes[chrom_name] = patched_file
                result_file = patched_file
                strategy_applied = "fixer_and_new_patches"
            else:
                self.log(f"  ✗ New patch application failed, using fixer result", "warning")
                self.fixed_chromosomes[chrom_name] = fixed_file
                result_file = fixed_file
                strategy_applied = "fixer_only"
        else:
            self.log(f"  ✗ No new patches extracted, using fixer result", "warning")
            self.fixed_chromosomes[chrom_name] = fixed_file
            result_file = fixed_file
            strategy_applied = "fixer_only"
        
        return {
            "status": "success",
            "chromosome": chrom_name,
            "strategy": "no_patch_fixer",
            "strategy_applied": strategy_applied,
            "result_file": str(result_file),
            "original_length": self._get_sequence_length(chrom_fasta),
            "final_length": self._get_sequence_length(result_file)
        }
    
    def _process_failed_completely_chromosome(self, chrom_name: str) -> Dict[str, Any]:
        chrom_dir = self._find_chromosome_dir(chrom_name)
        chrom_fasta = self._find_chromosome_fasta(chrom_dir, chrom_name)
        coords_file = self._find_coords_file(chrom_dir, chrom_name)
        
        if not all([chrom_dir, chrom_fasta, coords_file]):
            return {
                "status": "error",
                "error": f"Cannot find required files for chromosome {chrom_name}",
                "chromosome": chrom_name
            }
        
        if chrom_name not in self.chromosome_patches or not self.chromosome_patches[chrom_name]:
            self.log(f"  Warning: No patch information, running fixer directly", "warning")
            fixed_file = self._run_gap_fixer_module(chrom_dir, chrom_fasta, coords_file, chrom_name)
            if fixed_file:
                self.fixed_chromosomes[chrom_name] = fixed_file
                return {
                    "status": "success",
                    "chromosome": chrom_name,
                    "strategy": "full_fixer_reuse_patches",
                    "strategy_applied": "fixer_only",
                    "result_file": str(fixed_file),
                    "original_length": self._get_sequence_length(chrom_fasta),
                    "final_length": self._get_sequence_length(fixed_file),
                    "note": "no_patches_available"
                }
            else:
                return {
                    "status": "error",
                    "error": f"gap_fixer repair failed",
                    "chromosome": chrom_name
                }
        
        original_patches = self.chromosome_patches[chrom_name]
        original_gap_positions = [p["original_position"] for p in original_patches]
        
        fixed_file = self._run_gap_fixer_module(chrom_dir, chrom_fasta, coords_file, chrom_name)
        if not fixed_file:
            return {
                "status": "error",
                "error": f"gap_fixer repair failed",
                "chromosome": chrom_name
            }
        
        self.log(f"  ✓ gap_fixer repair completed", "success")
        
        self.log(f"  Analyzing coordinate changes and building mapping", "info")
        coordinate_mapping = self._analyze_coordinate_changes_detailed(
            chrom_name, chrom_fasta, fixed_file, original_gap_positions
        )
        
        if coordinate_mapping:
            self.log(f"  Updating patch coordinates and applying", "info")
            patched_file = self._apply_updated_patches_to_fixed(
                chrom_dir, fixed_file, original_patches, coordinate_mapping, chrom_name
            )
            
            if patched_file:
                self.log(f"  ✓ Patch application completed: {patched_file.name}", "success")
                self.fixed_chromosomes[chrom_name] = patched_file
                result_file = patched_file
                strategy_applied = "fixer_and_reused_patches"
            else:
                self.log(f"  ✗ Patch application failed, using fixer result", "warning")
                self.fixed_chromosomes[chrom_name] = fixed_file
                result_file = fixed_file
                strategy_applied = "fixer_only"
        else:
            self.log(f"  ✗ Cannot establish coordinate mapping, using fixer result", "warning")
            self.fixed_chromosomes[chrom_name] = fixed_file
            result_file = fixed_file
            strategy_applied = "fixer_only"
        
        return {
            "status": "success",
            "chromosome": chrom_name,
            "strategy": "full_fixer_reuse_patches",
            "strategy_applied": strategy_applied,
            "result_file": str(result_file),
            "original_length": self._get_sequence_length(chrom_fasta),
            "final_length": self._get_sequence_length(result_file),
            "patches_available": len(original_patches),
            "patches_applied": len(coordinate_mapping) if coordinate_mapping else 0
        }
    
    def _process_partial_success_chromosome(self, chrom_name: str) -> Dict[str, Any]:
        patcher_final = self.output_dir / "patcher_results" / "final_genome_patched.fasta"
        if not patcher_final.exists():
            for file_path in self.output_dir.glob("**/final_genome_patched.fasta"):
                patcher_final = file_path
                break
        
        if not patcher_final.exists():
            return {
                "status": "error",
                "error": f"Cannot find patcher final output file",
                "chromosome": chrom_name
            }
        
        extracted_file = self._extract_chromosome_from_patcher(patcher_final, chrom_name)
        
        if extracted_file:
            self.fixed_chromosomes[chrom_name] = extracted_file
            return {
                "status": "success",
                "chromosome": chrom_name,
                "strategy": "partial_success_no_fixer",
                "strategy_applied": "patcher_result",
                "result_file": str(extracted_file)
            }
        else:
            return {
                "status": "error",
                "error": f"Cannot extract chromosome {chrom_name} from patcher results",
                "chromosome": chrom_name
            }
    
    def _run_gap_fixer_module(self, chrom_dir: Path, chrom_fasta: Path, 
                            coords_file: Path, chrom_name: str) -> Optional[Path]:
        fix_dir = chrom_dir / "module_fixer_repair"
        fix_dir.mkdir(exist_ok=True)
        
        output_file = fix_dir / f"fixed_{chrom_name}.fasta"
        
        try:
            fixer_config = {
                "query_fasta": str(chrom_fasta),
                "coords_file": str(coords_file),
                "threads": max(2, self.threads // max(1, len(self.no_patch_chromosomes) + len(self.failed_completely_chromosomes))),
                "output_fasta": str(output_file),
                "work_dir": str(fix_dir),
                "repair_mode": self.repair_mode,
                "gap_length": self.gap_length,
                "max_search_distance": self.max_search_distance,
                "search_step": self.search_step,
                "verbose": not self.quiet,
                "include_raw_result": True
            }
            
            self.log(f"  Configuring gap_fixer: {fixer_config['repair_mode']} mode", "info")
            
            start_time = time.time()
            fixer_adapter = GenomeGapFixerAdapter(fixer_config)
            result = fixer_adapter.run()
            elapsed_time = time.time() - start_time
            
            if result.get("status") == "success":
                self.log(f"  gap_fixer completed, time: {elapsed_time:.1f}s", "info")
                
                if output_file.exists() and self._is_valid_fasta(output_file):
                    orig_length = self._get_sequence_length(chrom_fasta)
                    fixed_length = self._get_sequence_length(output_file)
                    length_change = fixed_length - orig_length
                    
                    self.log(f"  Length change: {orig_length:,}bp → {fixed_length:,}bp ({length_change:+,}bp)", "info")
                    return output_file
                else:
                    output_path = result.get("output_files", {}).get("final_genome", "")
                    if output_path and Path(output_path).exists():
                        return Path(output_path)
            
            self.log(f"  gap_fixer did not produce valid output: {result.get('error', 'Unknown error')}", "warning")
            return None
            
        except Exception as e:
            self.log(f"  gap_fixer exception: {e}", "error")
            return None
    
    def _rerun_patcher_for_fixed_chromosome(self, chrom_dir: Path, fixed_file: Path, 
                                          chrom_name: str) -> Optional[Path]:
        repatch_dir = chrom_dir / "repatch_after_fixer"
        repatch_dir.mkdir(exist_ok=True)
        
        try:
            patcher_config = {
                "query_fasta": str(fixed_file),
                "reference_fasta": str(self.reference_fasta),
                "coords_file": "",
                "threads": 4,
                "output_dir": str(repatch_dir),
                "chromosome_name": chrom_name,
                "min_gap_size": self.min_gap_size,
                "keep_temp": True,
                "verbose": not self.quiet,
                "include_raw_result": True
            }
            
            self.log(f"  Re-running patcher to extract patches...", "info")
            start_time = time.time()
            
            patcher_adapter = GapPatcherAdapter(patcher_config)
            result = patcher_adapter.run()
            elapsed_time = time.time() - start_time
            
            if result.get("status") == "success":
                self.log(f"  Re-patcher completed, time: {elapsed_time:.1f}s", "info")
                
                patches_dir = repatch_dir / chrom_name / "gap_patches"
                if not patches_dir.exists():
                    patches_dir = repatch_dir / "gap_patches"
                
                if patches_dir.exists():
                    patch_count = len(list(patches_dir.glob("*.fasta")))
                    self.log(f"  Extracted {patch_count} new patches", "info")
                    return patches_dir
                else:
                    self.log(f"  No new patches extracted", "warning")
                    return None
            else:
                self.log(f"  Re-running patcher failed: {result.get('error', 'Unknown error')}", "warning")
                return None
                
        except Exception as e:
            self.log(f"  Re-running patcher exception: {e}", "error")
            return None
    
    def _apply_new_patches_to_fixed(self, chrom_dir: Path, fixed_file: Path, 
                                  patches_dir: Path, chrom_name: str) -> Optional[Path]:
        apply_dir = chrom_dir / "apply_new_patches"
        apply_dir.mkdir(exist_ok=True)
        
        new_patches = []
        for patch_file in patches_dir.glob("*.fasta"):
            match = re.search(r'gap_(\d+)_patch', patch_file.name)
            if match:
                gap_position = int(match.group(1))
                patch_info = {
                    "original_position": gap_position,
                    "file": patch_file,
                    "status": "new"
                }
                new_patches.append(patch_info)
        
        if not new_patches:
            self.log(f"  No new patches to apply", "warning")
            return None
        
        self.log(f"  {len(new_patches)} new patches to apply", "info")
        
        patch_script = self._find_script("patch_gap.py")
        if not patch_script:
            self.log(f"  Error: Cannot find patch_gap.py script", "error")
            return None
        
        new_patches.sort(key=lambda x: x["original_position"], reverse=True)
        
        current_file = fixed_file
        success_count = 0
        
        for i, patch_info in enumerate(new_patches):
            position = patch_info["original_position"]
            patch_file = patch_info["file"]
            
            self.log(f"  Applying patch {i+1}/{len(new_patches)}: position {position:,}", "info")
            
            output_file = apply_dir / f"temp_patched_{i}.fasta"
            
            if self._apply_single_patch(patch_script, current_file, patch_file, 
                                      position, output_file):
                if output_file.exists() and self._is_valid_fasta(output_file):
                    current_file = output_file
                    success_count += 1
                    self.log(f"    ✓ Patch applied successfully", "info")
                else:
                    self.log(f"    ✗ Invalid file after patch application", "warning")
            else:
                self.log(f"    ✗ Patch application failed", "warning")
        
        if success_count > 0:
            final_file = apply_dir / f"{chrom_name}_repatched.fasta"
            shutil.copy2(current_file, final_file)
            
            if final_file.exists():
                self.log(f"  Successfully applied {success_count}/{len(new_patches)} new patches", "success")
                return final_file
        
        return None
    
    def _analyze_coordinate_changes_detailed(self, chrom_name: str, original_file: Path,
                                           fixed_file: Path, original_gap_positions: List[int]) -> Dict[int, int]:
        if not BIOPYTHON_AVAILABLE:
            self.log("  Warning: BioPython not available, using simplified coordinate mapping", "warning")
            return {pos: pos for pos in original_gap_positions}
        
        try:
            original_seq = SeqIO.read(original_file, "fasta").seq
            fixed_seq = SeqIO.read(fixed_file, "fasta").seq
            
            original_gaps = self._find_gap_regions(str(original_seq), min_gap_size=10)
            fixed_gaps = self._find_gap_regions(str(fixed_seq), min_gap_size=10)
            
            if len(original_gaps) == len(fixed_gaps) and len(original_gaps) > 0:
                original_gaps_sorted = sorted(original_gaps, key=lambda x: x["start"])
                fixed_gaps_sorted = sorted(fixed_gaps, key=lambda x: x["start"])
                
                coordinate_mapping = {}
                for i, (orig_gap, fixed_gap) in enumerate(zip(original_gaps_sorted, fixed_gaps_sorted)):
                    orig_mid = orig_gap["start"] + orig_gap["length"] // 2
                    fixed_mid = fixed_gap["start"] + fixed_gap["length"] // 2
                    
                    closest_original_pos = min(original_gap_positions, key=lambda x: abs(x - orig_mid))
                    coordinate_mapping[closest_original_pos] = fixed_mid
                
                self.log(f"  Established {len(coordinate_mapping)} coordinate mappings", "info")
                return coordinate_mapping
            else:
                self.log(f"  Gap count mismatch: original {len(original_gaps)}, fixed {len(fixed_gaps)}", "warning")
                
                if chrom_name in self.gap_details and self.gap_details[chrom_name]:
                    gap_details = self.gap_details[chrom_name]
                    coordinate_mapping = {}
                    
                    for gap_info in gap_details:
                        orig_pos = gap_info.get("start_position", 0)
                        if orig_pos in original_gap_positions:
                            coordinate_mapping[orig_pos] = orig_pos
                    
                    if coordinate_mapping:
                        self.log(f"  Used gap details to establish {len(coordinate_mapping)} coordinate mappings", "info")
                        return coordinate_mapping
                
                self.log(f"  Using unchanged position assumption", "warning")
                return {pos: pos for pos in original_gap_positions}
                
        except Exception as e:
            self.log(f"  Coordinate mapping analysis failed: {e}", "error")
            return {}
    
    def _apply_updated_patches_to_fixed(self, chrom_dir: Path, fixed_file: Path,
                                      original_patches: List[Dict], 
                                      coordinate_mapping: Dict[int, int],
                                      chrom_name: str) -> Optional[Path]:
        apply_dir = chrom_dir / "apply_updated_patches"
        apply_dir.mkdir(exist_ok=True)
        
        updated_patches = []
        for patch_info in original_patches:
            original_pos = patch_info["original_position"]
            new_pos = coordinate_mapping.get(original_pos)
            
            if new_pos:
                temp_patch = self._create_patch_with_correct_header(patch_info["file"], new_pos)
                if temp_patch:
                    patch_info["new_position"] = new_pos
                    patch_info["temp_file"] = temp_patch
                    updated_patches.append(patch_info)
                    self.log(f"  Patch {original_pos:,} → {new_pos:,} coordinate updated", "info")
        
        if not updated_patches:
            self.log(f"  No valid updated patches", "warning")
            return None
        
        self.log(f"  {len(updated_patches)} updated patches to apply", "info")
        
        patch_script = self._find_script("patch_gap.py")
        if not patch_script:
            self.log(f"  Error: Cannot find patch_gap.py script", "error")
            return None
        
        updated_patches.sort(key=lambda x: x["new_position"], reverse=True)
        
        current_file = fixed_file
        success_count = 0
        
        for i, patch_info in enumerate(updated_patches):
            original_pos = patch_info["original_position"]
            new_pos = patch_info["new_position"]
            temp_patch = patch_info["temp_file"]
            
            self.log(f"  Applying patch {i+1}/{len(updated_patches)}: {original_pos:,}→{new_pos:,}", "info")
            
            output_file = apply_dir / f"temp_patched_{i}.fasta"
            
            if self._apply_single_patch_with_temp_file(patch_script, current_file, 
                                                     temp_patch, new_pos, output_file):
                if output_file.exists() and self._is_valid_fasta(output_file):
                    current_file = output_file
                    success_count += 1
                    self.log(f"    ✓ Patch applied successfully", "info")
                else:
                    self.log(f"    ✗ Invalid file after patch application", "warning")
            else:
                self.log(f"    ✗ Patch application failed", "warning")
            
            if temp_patch and os.path.exists(temp_patch):
                try:
                    os.remove(temp_patch)
                except:
                    pass
        
        if success_count > 0:
            final_file = apply_dir / f"{chrom_name}_updated_patches_applied.fasta"
            shutil.copy2(current_file, final_file)
            
            if final_file.exists():
                self.log(f"  Successfully applied {success_count}/{len(updated_patches)} updated patches", "success")
                return final_file
        
        return None
    
    def _extract_chromosome_from_patcher(self, patcher_final: Path, chrom_name: str) -> Optional[Path]:
        if not BIOPYTHON_AVAILABLE:
            return None
        
        try:
            for record in SeqIO.parse(patcher_final, "fasta"):
                record_id = record.id
                
                if (record_id == chrom_name or 
                    record_id.startswith(f"{chrom_name}_") or
                    chrom_name in record_id.split('_')):
                    
                    output_file = self.output_dir / f"{chrom_name}_partial_success.fasta"
                    record.id = chrom_name
                    record.description = f"partially_successful_{chrom_name}"
                    SeqIO.write([record], output_file, "fasta")
                    
                    if output_file.exists():
                        return output_file
            
            return None
            
        except Exception as e:
            self.log(f"  Failed to extract chromosome: {e}", "error")
            return None
    
    def _determine_chromosome_strategy(self, chrom_name, total_gaps, extracted_patches, successful_patches):
        if total_gaps == 0:
            return "no_gaps"
        
        if extracted_patches == 0:
            return "no_patch_fixer"
        elif successful_patches == 0 and extracted_patches > 0:
            return "full_fixer_reuse_patches"
        elif successful_patches > 0:
            return "partial_success_no_fixer"
        else:
            return "unknown"
    
    def _find_chromosome_dir(self, chrom_name: str) -> Optional[Path]:
        base_dirs = [
            self.output_dir / "patcher_results",
            self.output_dir,
        ]
        
        for base_dir in base_dirs:
            if not base_dir.exists():
                continue
            
            possible_dirs = [
                base_dir / chrom_name,
                base_dir / f"chromosome_{chrom_name}",
                base_dir / f"chr_{chrom_name}",
            ]
            
            for dir_path in possible_dirs:
                if dir_path.exists() and dir_path.is_dir():
                    return dir_path
            
            for item in base_dir.iterdir():
                if item.is_dir() and chrom_name in item.name:
                    return item
        
        return None
    
    def _find_chromosome_fasta(self, chrom_dir: Path, chrom_name: str) -> Optional[Path]:
        if not chrom_dir or not chrom_dir.exists():
            return None
        
        possible_files = [
            chrom_dir / f"{chrom_name}.fa",
            chrom_dir / f"{chrom_name}.fasta",
            chrom_dir / f"{chrom_name}_patched.fasta",
            chrom_dir / "chromosome.fa",
            chrom_dir / "sequence.fa",
        ]
        
        for file_path in possible_files:
            if file_path.exists():
                return file_path
        
        for ext in ["*.fa", "*.fasta"]:
            for file_path in chrom_dir.glob(ext):
                if file_path.exists():
                    return file_path
        
        return None
    
    def _find_coords_file(self, chrom_dir: Path, chrom_name: str) -> Optional[Path]:
        if not chrom_dir or not chrom_dir.exists():
            return None
        
        for pattern in ["*.coords", "*merged.coords"]:
            for file_path in chrom_dir.glob(pattern):
                if file_path.exists():
                    return file_path
        
        parent_dir = chrom_dir.parent
        if parent_dir.exists():
            for pattern in ["*.coords", "*merged.coords"]:
                for file_path in parent_dir.glob(pattern):
                    if file_path.exists():
                        return file_path
        
        if self.coords_file and self.coords_file.exists():
            return self.coords_file
        
        return None
    
    def _is_valid_fasta(self, filepath: Path) -> bool:
        if not filepath.exists():
            return False
        
        if filepath.stat().st_size < 10:
            return False
        
        try:
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>')
        except:
            return False
    
    def _get_sequence_length(self, fasta_file: Path) -> int:
        try:
            length = 0
            with open(fasta_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('>'):
                        length += len(line)
            return length
        except:
            return 0
    
    def _find_gap_regions(self, sequence: str, min_gap_size: int = 50) -> List[Dict]:
        gaps = []
        
        i = 0
        while i < len(sequence):
            if sequence[i].upper() == 'N':
                start = i
                while i < len(sequence) and sequence[i].upper() == 'N':
                    i += 1
                end = i
                gap_length = end - start
                
                if gap_length >= min_gap_size:
                    gaps.append({
                        "start": start,
                        "end": end,
                        "length": gap_length,
                        "mid": start + gap_length // 2
                    })
            else:
                i += 1
        
        return gaps
    
    def _create_patch_with_correct_header(self, patch_file: Path, position: int) -> Optional[str]:
        try:
            with open(patch_file, 'r') as f:
                content = f.read()
            
            temp_fd, temp_path = tempfile.mkstemp(suffix='.fasta', prefix=f'gap_{position}_')
            os.close(temp_fd)
            
            lines = content.strip().split('\n')
            if lines and lines[0].startswith('>'):
                original_header = lines[0][1:]
                new_header = f">gap{position}"
                
                if ' ' in original_header:
                    desc_parts = original_header.split(' ', 1)
                    if len(desc_parts) > 1:
                        new_header += f" {desc_parts[1]}"
                
                lines[0] = new_header
            
            with open(temp_path, 'w') as f:
                f.write('\n'.join(lines))
            
            return temp_path
            
        except Exception as e:
            self.log(f"    Failed to create temporary patch: {e}", "warning")
            return None
    
    def _apply_single_patch(self, patch_script: Path, chrom_file: Path, 
                          patch_file: Path, position: int, output_file: Path) -> bool:
        cmd = [
            "python", str(patch_script),
            "-r", str(chrom_file),
            "-p", str(patch_file),
            "--gap-position", str(position),
            "-o", str(output_file),
            "--flank-size", "10000",
            "--minimap2-mode", "asm5",
            "--min-score", "0.7",
            "--min-match-length", "100",
            "--min-mapq", "20",
            "--search-range", "500000"
        ]
        
        if self.quiet:
            cmd.append("--quiet")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode != 0:
                error_msg = result.stderr.lower()
                if "already covered" in error_msg or "already filled" in error_msg:
                    if chrom_file.exists():
                        shutil.copy2(chrom_file, output_file)
                    return True
                else:
                    return False
            
            return True
            
        except Exception as e:
            self.log(f"    Patch filling exception: {e}", "warning")
            return False
    
    def _apply_single_patch_with_temp_file(self, patch_script: Path, chrom_file: Path, 
                                         temp_patch: str, position: int, 
                                         output_file: Path) -> bool:
        cmd = [
            "python", str(patch_script),
            "-r", str(chrom_file),
            "-p", temp_patch,
            "--gap-position", str(position),
            "-o", str(output_file),
            "--flank-size", "10000",
            "--minimap2-mode", "asm5",
            "--min-score", "0.7",
            "--min-match-length", "100",
            "--min-mapq", "20",
            "--search-range", "500000"
        ]
        
        if self.quiet:
            cmd.append("--quiet")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode != 0:
                error_msg = result.stderr.lower()
                if "already covered" in error_msg or "already filled" in error_msg:
                    if chrom_file.exists():
                        shutil.copy2(chrom_file, output_file)
                    return True
                else:
                    return False
            
            return True
            
        except Exception as e:
            self.log(f"    Patch filling exception: {e}", "warning")
            return False
    
    def _find_script(self, script_name: str) -> Optional[Path]:
        possible_locations = [
            self.script_dir / script_name,
            self.script_dir.parent / script_name,
            Path(".") / script_name,
        ]
        
        for location in possible_locations:
            if location.exists():
                return location
        
        return None
    
    def run_complete_pipeline(self) -> bool:
        start_time = time.time()
        
        self.log("="*60, "info")
        self.log("Modular Enhanced GAP Complete Process Controller Starting", "info")
        self.log("="*60, "info")
        self.log(f"Query genome: {self.query_fasta.name}", "info")
        self.log(f"Reference genome: {self.reference_fasta.name}", "info")
        self.log(f"Threads: {self.threads}", "info")
        self.log(f"Repair mode: {self.repair_mode}", "info")
        self.log(f"Output directory: {self.output_dir}", "info")
        
        if self.coords_file:
            self.log(f"Using existing coords file: {self.coords_file.name}", "info")
        
        self.log(f"Execution mode: Modular call + post-fixer patching", "info")
        self.log("="*60, "info")
        
        try:
            self.log("\nStarting Step 1: Running gap_patcher", "info")
            if not self._run_patcher_module():
                self.log("Step 1 failed, process terminated", "error")
                return False
            
            self._collect_all_chromosome_patches()
            
            self.log("\nStarting Step 2: Processing chromosomes by category", "info")
            
            if self.no_patch_chromosomes:
                for chrom_name in self.no_patch_chromosomes:
                    self.log(f"\nProcessing no-patch chromosome: {chrom_name}", "info")
                    result = self._process_no_patch_chromosome(chrom_name)
                    if result["status"] != "success":
                        self.log(f"  Chromosome {chrom_name} processing failed: {result.get('error')}", "warning")
            else:
                self.log("No no-patch chromosomes to process", "info")
            
            if self.failed_completely_chromosomes:
                for chrom_name in self.failed_completely_chromosomes:
                    self.log(f"\nProcessing all-failed chromosome: {chrom_name}", "info")
                    result = self._process_failed_completely_chromosome(chrom_name)
                    if result["status"] != "success":
                        self.log(f"  Chromosome {chrom_name} processing failed: {result.get('error')}", "warning")
            else:
                self.log("No all-failed chromosomes to process", "info")
            
            if self.partially_successful_chromosomes:
                for chrom_name in self.partially_successful_chromosomes:
                    self.log(f"\nProcessing partially successful chromosome: {chrom_name}", "info")
                    result = self._process_partial_success_chromosome(chrom_name)
                    if result["status"] != "success":
                        self.log(f"  Chromosome {chrom_name} processing failed: {result.get('error')}", "warning")
            else:
                self.log("No partially successful chromosomes", "info")
            
            self.log("\nStarting Step 3: Creating final genome", "info")
            final_genome = self._create_final_genome_all()
            
            elapsed_time = time.time() - start_time
            
            if final_genome and final_genome.exists():
                self.log("\n" + "="*60, "success")
                self.log("✅ Modular Enhanced Process Completed!", "success")
                self.log(f"   Total time: {elapsed_time:.1f}s ({elapsed_time/60:.1f}min)", "info")
                self.log(f"   Final genome: {final_genome}", "info")
                self.log(f"   No-patch chromosomes: {len(self.no_patch_chromosomes)}", "info")
                self.log(f"   All-failed chromosomes: {len(self.failed_completely_chromosomes)}", "info")
                self.log(f"   Partial success chromosomes: {len(self.partially_successful_chromosomes)}", "info")
                self.log(f"   Log file: {self.log_file}", "info")
                self.log("="*60, "info")
                
                self._generate_final_report(elapsed_time, final_genome)
                
                return True
            else:
                self.log("\n❌ Process completed but no final genome generated", "error")
                return False
            
        except KeyboardInterrupt:
            self.log("\n❌ Process interrupted by user", "error")
            return False
        except Exception as e:
            self.log(f"\n❌ Process execution failed: {e}", "error")
            import traceback
            self.log(traceback.format_exc(), "error")
            return False
    
    def _create_final_genome_all(self) -> Optional[Path]:
        final_file = self.output_dir / "final_genome_module_enhanced.fasta"
        
        if not BIOPYTHON_AVAILABLE:
            self.log("BioPython not available, cannot create final genome", "error")
            return None
        
        all_records = []
        processed_chromosomes = set()
        
        for chrom_name, chrom_file in self.fixed_chromosomes.items():
            if chrom_file.exists() and self._is_valid_fasta(chrom_file):
                try:
                    records = list(SeqIO.parse(chrom_file, "fasta"))
                    for record in records:
                        record.id = chrom_name
                        strategy = self.chromosome_strategies.get(chrom_name, "unknown")
                        record.description = f"processed_{strategy}"
                        all_records.append(record)
                        processed_chromosomes.add(chrom_name)
                        self.log(f"  Adding processed chromosome: {chrom_name} ({strategy})", "info")
                except Exception as e:
                    self.log(f"  Failed to read chromosome {chrom_name}: {e}", "warning")
        
        try:
            for record in SeqIO.parse(self.query_fasta, "fasta"):
                chrom_name = record.id
                
                if chrom_name not in processed_chromosomes:
                    similar_exists = False
                    for processed_name in processed_chromosomes:
                        if chrom_name in processed_name or processed_name in chrom_name:
                            self.log(f"  Skipping possibly duplicate original chromosome: {chrom_name} (already have: {processed_name})", "info")
                            similar_exists = True
                            break
                    
                    if not similar_exists:
                        all_records.append(record)
                        self.log(f"  Adding original chromosome: {chrom_name}", "info")
        except Exception as e:
            self.log(f"  Failed to read original genome: {e}", "warning")
        
        try:
            if all_records:
                SeqIO.write(all_records, final_file, "fasta")
                
                records = list(SeqIO.parse(final_file, "fasta"))
                total_length = sum(len(rec.seq) for rec in records)
                
                self.log(f"\nFinal genome created: {final_file}", "success")
                self.log(f"  Contains {len(records)} chromosomes", "info")
                self.log(f"  Total length: {total_length:,} bp", "info")
                
                return final_file
            else:
                self.log(f"  No chromosome data to write", "error")
                return None
                
        except Exception as e:
            self.log(f"Failed to create final genome: {e}", "error")
            import traceback
            self.log(traceback.format_exc(), "error")
            return None
    
    def _generate_final_report(self, elapsed_time: float, final_genome: Path):
        report = {
            "pipeline": "GapCompleteControllerModuleEnhanced",
            "version": "2.0",
            "timestamp": datetime.now().isoformat(),
            "processing_time_seconds": round(elapsed_time, 2),
            "execution_mode": "pure_module_with_post_fixer_patching",
            "parameters": {
                "query_fasta": str(self.query_fasta),
                "reference_fasta": str(self.reference_fasta),
                "coords_file": str(self.coords_file) if self.coords_file else None,
                "threads": self.threads,
                "repair_mode": self.repair_mode,
                "gap_length": self.gap_length,
                "max_search_distance": self.max_search_distance,
                "search_step": self.search_step,
                "min_gap_size": self.min_gap_size,
                "output_dir": str(self.output_dir)
            },
            "results": {
                "final_genome": str(final_genome),
                "log_file": str(self.log_file),
                "chromosome_analysis": str(self.output_dir / "chromosome_analysis.json")
            },
            "strategy_summary": {
                "no_patch_chromosomes": {
                    "count": len(self.no_patch_chromosomes),
                    "list": self.no_patch_chromosomes,
                    "strategy": "full_fixer_and_repatch"
                },
                "failed_completely_chromosomes": {
                    "count": len(self.failed_completely_chromosomes),
                    "list": self.failed_completely_chromosomes,
                    "strategy": "full_fixer_and_reuse_patches"
                },
                "partially_successful_chromosomes": {
                    "count": len(self.partially_successful_chromosomes),
                    "list": self.partially_successful_chromosomes,
                    "strategy": "no_fixer_use_as_is"
                }
            },
            "chromosome_details": self.chromosome_info,
            "patch_counts": {k: len(v) for k, v in self.chromosome_patches.items()},
            "fixed_chromosomes": list(self.fixed_chromosomes.keys())
        }
        
        report_file = self.output_dir / "module_enhanced_complete_report.json"
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False, default=str)
        
        self.log(f"Complete report saved: {report_file}", "info")

def run_gap_complete(query_fasta, reference_fasta, **kwargs):
    config = {
        'query_fasta': str(query_fasta),
        'reference_fasta': str(reference_fasta),
        **kwargs
    }
    
    try:
        controller = GapCompleteControllerModuleEnhanced(config_dict=config)
        success = controller.run_complete_pipeline()
        
        if success:
            results = controller.get_results()
            results['status'] = 'success'
            return results
        else:
            return {
                'status': 'error',
                'error': 'Process execution failed',
                'output_dir': str(controller.output_dir),
                'log_file': str(controller.log_file)
            }
            
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }

def analyze_strategies(query_fasta, reference_fasta, **kwargs):
    config = {
        'query_fasta': str(query_fasta),
        'reference_fasta': str(reference_fasta),
        **kwargs
    }
    
    try:
        controller = GapCompleteControllerModuleEnhanced(config_dict=config)
        return controller.run_patcher_analysis()
        
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }

def fill_specific_chromosome(chrom_name, query_fasta, reference_fasta, strategy=None, **kwargs):
    config = {
        'query_fasta': str(query_fasta),
        'reference_fasta': str(reference_fasta),
        **kwargs
    }
    
    try:
        controller = GapCompleteControllerModuleEnhanced(config_dict=config)
        return controller.run_for_chromosome(chrom_name, strategy)
        
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e),
            'chromosome': chrom_name
        }

def create_controller(config_dict):
    return GapCompleteControllerModuleEnhanced(config_dict=config_dict)

def create_argument_parser():
    parser = argparse.ArgumentParser(
        description='Modular GAP Complete Process Controller - Enhanced Version',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  
  As command line tool:
    # Full process
    python gap_complete_controller_module_enhanced.py -q query.fasta -c ref.fasta
    
    # Strategy analysis only
    python gap_complete_controller_module_enhanced.py -q query.fasta -c ref.fasta --analyze-only
    
    # Process specific chromosome
    python gap_complete_controller_module_enhanced.py -q query.fasta -c ref.fasta --chromosome chr1
  
  As Python module:
    from gap_complete_controller_module_enhanced import run_gap_complete
    
    result = run_gap_complete(
        query_fasta="query.fasta",
        reference_fasta="ref.fasta",
        threads=16,
        output_dir="results"
    )
        """
    )
    
    parser.add_argument('-q', '--query-fasta', required=True,
                       help='Query genome FASTA file')
    parser.add_argument('-c', '--reference-fasta', required=True,
                       help='Reference genome FASTA file')
    
    parser.add_argument('-coords', '--coords-file',
                       help='Existing nucmer alignment coords file')
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='Total threads, default 8')
    parser.add_argument('-o', '--output-dir', default='gap_complete_module_enhanced_results',
                       help='Output directory, default gap_complete_module_enhanced_results')
    
    parser.add_argument('--repair-mode', choices=['conservative', 'aggressive'],
                       default='aggressive', help='Repair mode, default aggressive')
    
    parser.add_argument('--gap-length', type=int, default=100,
                       help='Final gap length (bp), default 100')
    parser.add_argument('--search', '--max-search-distance', type=int, 
                       default=500000, dest='max_search_distance',
                       help='Maximum search distance (bp), default 500000')
    parser.add_argument('--step', '--search-step', type=int, 
                       default=100000, dest='search_step',
                       help='Search step (bp), default 100000')
    
    parser.add_argument('--min-gap-size', type=int, default=100,
                       help='Minimum GAP size (bp), default 100')
    
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files')
    parser.add_argument('--quiet', action='store_true',
                       help='Quiet mode')
    
    parser.add_argument('--analyze-only', action='store_true',
                       help='Analyze strategies only, no filling')
    parser.add_argument('--chromosome', 
                       help='Specific chromosome to process')
    parser.add_argument('--strategy',
                       choices=['no_patch_fixer', 'full_fixer_reuse_patches', 'partial_success_no_fixer'],
                       help='Specify processing strategy (auto-selected if not specified)')
    
    return parser

def main():
    parser = create_argument_parser()
    args = parser.parse_args()
    
    if not os.path.exists(args.query_fasta):
        print(f"❌ Error: Query genome file does not exist: {args.query_fasta}")
        sys.exit(1)
    
    if not os.path.exists(args.reference_fasta):
        print(f"❌ Error: Reference genome file does not exist: {args.reference_fasta}")
        sys.exit(1)
    
    if args.coords_file and not os.path.exists(args.coords_file):
        print(f"⚠️  Warning: Coords file does not exist: {args.coords_file}")
        args.coords_file = None
    
    print("🚀 Starting Modular GAP Complete Process Controller...")
    
    try:
        controller = GapCompleteControllerModuleEnhanced(args)
        
        if args.analyze_only:
            print("Running strategy analysis mode...")
            result = controller.run_patcher_analysis()
            
        elif args.chromosome:
            print(f"Processing chromosome: {args.chromosome}")
            result = controller.run_for_chromosome(args.chromosome, args.strategy)
            
        else:
            success = controller.run_complete_workflow()
            result = {
                "status": "success" if success else "error",
                "success": success,
                "output_dir": str(controller.output_dir),
                "log_file": str(controller.log_file)
            }
        
        if isinstance(result, dict):
            if result.get('status') == 'success':
                print("\n🎉 Operation completed successfully!")
                if 'output_dir' in result:
                    print(f"   Output directory: {result['output_dir']}")
                if 'log_file' in result:
                    print(f"   Log file: {result['log_file']}")
                sys.exit(0)
            else:
                print("\n💥 Operation failed")
                if 'error' in result:
                    print(f"   Error: {result['error']}")
                sys.exit(1)
        else:
            if result:
                print("\n🎉 Process completed successfully!")
                sys.exit(0)
            else:
                print("\n💥 Process execution failed")
                sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n🛑 Process interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 Uncaught exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()