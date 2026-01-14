#!/usr/bin/env python3
"""
Genome GAP Patching Main Controller - Logger Error Fixed Version
Added adapter pattern support
Optimization: Perform nucmer analysis only on chromosomes with GAPs
New feature: Support using existing coords files directly, skipping nucmer alignment
"""

import os
import sys
import json
import time
import logging
import argparse
import subprocess
import concurrent.futures
import shutil
from typing import List, Dict, Tuple, Optional, Any

# ==================== Path Handling Utilities ====================

def get_script_dir():
    """Get absolute path of script directory"""
    return os.path.dirname(os.path.abspath(__file__))

def find_script(script_name):
    """Find script file path"""
    # 1. Try current directory
    current_dir = get_script_dir()
    script_path = os.path.join(current_dir, script_name)
    
    if os.path.exists(script_path):
        return script_path
    
    # 2. Try to find in PATH
    try:
        full_path = shutil.which(script_name)
        if full_path:
            return full_path
    except:
        pass
    
    # 3. Try relative paths
    for path in ['.', '..', '../..']:
        script_path = os.path.join(current_dir, path, script_name)
        if os.path.exists(script_path):
            return script_path
    
    raise FileNotFoundError(f"Script file not found: {script_name}")

# ==================== Main Controller Class ====================

class GenomeGapPatcherController:
    """Genome GAP Patching Main Controller"""
    
    def __init__(
        self,
        ref_fasta: str,
        query_fasta: str,
        threads: int = 8,
        output_dir: str = "results",
        min_gap_size: int = 100,
        keep_temp: bool = False,
        verbose: bool = True,
        coords_file: Optional[str] = None  # New: existing coords file path
    ):
        """
        Initialize controller
        
        Args:
            ref_fasta: Reference contigs file path
            query_fasta: Target genome file path
            threads: Number of threads
            output_dir: Output directory
            min_gap_size: Minimum GAP size
            keep_temp: Whether to keep temporary files
            verbose: Whether to show detailed information
            coords_file: Existing coords file path (if provided, skip nucmer alignment)
        """
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.query_fasta = os.path.abspath(query_fasta)
        self.threads = threads
        self.output_dir = os.path.abspath(output_dir)
        self.min_gap_size = min_gap_size
        self.keep_temp = keep_temp
        self.verbose = verbose
        self.coords_file = os.path.abspath(coords_file) if coords_file else None
        
        # Get script directory
        self.script_dir = get_script_dir()
        
        # Setup logging system first
        self._setup_logging()
        
        # Now can safely use logger
        self.logger.info(f"Script directory: {self.script_dir}")
        
        # Initialize status
        self.gap_positions = {}
        self.chromosome_dirs = {}
        self.status = {}
        self.results = {}
        
        # Find dependent script paths
        self.script_paths = self._find_script_paths()
        
        # Validate files
        self._validate_files()
    
    def _setup_logging(self):
        """Setup logging system"""
        self.logger = logging.getLogger('GenomeGapPatcher')
        self.logger.setLevel(logging.INFO if self.verbose else logging.WARNING)
        
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # Console handler
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        
        # File handler
        log_dir = os.path.join(self.output_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, "gap_patcher.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        
        self.logger.info("=" * 60)
        self.logger.info("Genome GAP Patching Controller Started")
        self.logger.info("=" * 60)
        
        # Record whether using existing coords file
        if self.coords_file:
            self.logger.info(f"Will use existing coords file, skipping nucmer alignment: {os.path.basename(self.coords_file)}")
    
    def _find_script_paths(self) -> Dict[str, str]:
        """Find paths of all dependent scripts"""
        scripts = [
            "find_gaps.py",
            "run_nucmer_parallel.py", 
            "extract_gap_patches.py",
            "patch_gap.py"
        ]
        
        script_paths = {}
        for script in scripts:
            try:
                script_paths[script] = find_script(script)
                self.logger.info(f"✓ Found script: {script} -> {os.path.basename(script_paths[script])}")
            except FileNotFoundError as e:
                self.logger.error(f"✗ Script not found: {script}")
                self.logger.error("Please ensure the following scripts are in the same directory:")
                for s in scripts:
                    self.logger.error(f"  - {s}")
                raise
        
        return script_paths
    
    def _validate_files(self):
        """Validate input files"""
        required_files = [
            (self.ref_fasta, "Reference contigs file"),
            (self.query_fasta, "Target genome file")
        ]
        
        for filepath, desc in required_files:
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"{desc} does not exist: {filepath}")
            self.logger.info(f"✓ {desc}: {os.path.basename(filepath)}")
        
        # Validate coords file (if provided)
        if self.coords_file:
            if not os.path.exists(self.coords_file):
                raise FileNotFoundError(f"Provided coords file does not exist: {self.coords_file}")
            self.logger.info(f"✓ Using existing coords file: {os.path.basename(self.coords_file)}")
    
    def run_workflow(self) -> Dict[str, Any]:
        """
        Run complete GAP patching workflow
        """
        start_time = time.time()
        
        try:
            # Step 1: Find all GAP positions
            self.logger.info("Step 1: Finding GAP positions in the genome")
            gap_results = self.find_all_gaps()
            
            # If no GAPs found, exit directly
            if not self.gap_positions or sum(len(gaps) for gaps in self.gap_positions.values()) == 0:
                self.logger.info("No GAPs found to patch, workflow completed")
                return self._create_report(start_time)
            
            # Step 2: Use existing coords file if provided; otherwise run nucmer alignment
            if self.coords_file and os.path.exists(self.coords_file):
                self.logger.info("Step 2: Using existing coords file, skipping nucmer alignment")
                nucmer_results = self.process_existing_coords_file()
            else:
                self.logger.info("Step 2: Running nucmer comparison for chromosomes with GAPs")
                nucmer_results = self.run_nucmer_comparison()
            
            # Step 3: Process chromosomes with GAPs
            self.logger.info("Step 3: Extracting patches and patching GAPs")
            patch_results = self.process_chromosomes_with_gaps()
            
            # Step 4: Merge patched chromosomes
            self.logger.info("Step 4: Merging all chromosomes")
            merge_result = self.merge_patched_chromosomes()
            
            # Generate final report
            self.logger.info("Step 5: Generating final report")
            final_report = self._create_report(start_time, {
                "gap_finding": gap_results,
                "nucmer_comparison": nucmer_results,
                "gap_patching": patch_results,
                "chromosome_merging": merge_result
            })
            
            # Clean up temporary files
            if not self.keep_temp:
                self.cleanup_temp_files()
            
            self.logger.info("=" * 60)
            self.logger.info("✅ GAP patching workflow completed!")
            self.logger.info("=" * 60)
            
            return final_report
            
        except Exception as e:
            self.logger.error(f"Workflow execution failed: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            raise
    
    def find_all_gaps(self) -> Dict[str, Any]:
        """
        Use find_gaps.py to find all GAP positions
        """
        self.logger.info("Running find_gaps.py for GAP analysis...")
        
        # Call find_gaps.py API
        try:
            # Add script directory to Python path
            sys.path.insert(0, self.script_dir)
            
            from find_gaps import GapFinderAPI
            
            finder = GapFinderAPI(verbose=self.verbose)
            result = finder.find_gaps(
                fasta_file=self.query_fasta,
                output_dir=os.path.join(self.output_dir, "gap_analysis"),
                output_prefix="gap_locations",
                min_gap_size=self.min_gap_size,
                parallel_threshold_mb=100,
                save_json=True,
                save_summary=True
            )
            
            # Extract GAP positions and group by chromosome
            if "gaps" in result:
                for gap in result["gaps"]:
                    chrom = gap["chrom"]
                    gap_pos = gap["start"] + gap["length"] // 2  # Use middle position
                    
                    if chrom not in self.gap_positions:
                        self.gap_positions[chrom] = []
                    
                    self.gap_positions[chrom].append(gap_pos)
            
            # Record statistics
            total_gaps = sum(len(gaps) for gaps in self.gap_positions.values())
            chromosomes_with_gaps = len(self.gap_positions)
            
            self.logger.info(f"Found {total_gaps} GAPs distributed across {chromosomes_with_gaps} chromosomes")
            
            for chrom, gaps in list(self.gap_positions.items())[:5]:
                self.logger.info(f"  {chrom}: {len(gaps)} GAPs")
            
            if chromosomes_with_gaps > 5:
                self.logger.info(f"  ... and {chromosomes_with_gaps - 5} more chromosomes")
            
            return result
            
        except ImportError as e:
            # If cannot import API, use command line call
            self.logger.warning(f"Cannot import GapFinderAPI: {e}")
            self.logger.warning("Using command line to call find_gaps.py")
            return self._run_find_gaps_commandline()
    
    def _run_find_gaps_commandline(self) -> Dict[str, Any]:
        """Run find_gaps.py via command line"""
        output_json = os.path.join(self.output_dir, "gap_analysis", "gap_locations.json")
        os.makedirs(os.path.dirname(output_json), exist_ok=True)
        
        cmd = [
            "python", self.script_paths["find_gaps.py"],
            self.query_fasta,
            "-o", output_json,
            "--gap-char", "N",
            "--min-gap-size", str(self.min_gap_size),
            "--output-dir", os.path.join(self.output_dir, "gap_analysis"),
            "--api-mode"
        ]
        
        if not self.verbose:
            cmd.append("--quiet")
        
        self.logger.info(f"Executing command: {' '.join(cmd[:6])}...")
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=3600
            )
            
            # Load results
            if os.path.exists(output_json):
                with open(output_json, 'r') as f:
                    result_data = json.load(f)
                
                # Extract GAP positions
                if "gaps" in result_data:
                    for gap in result_data["gaps"]:
                        chrom = gap["chrom"]
                        gap_pos = gap["start"] + gap["length"] // 2
                        
                        if chrom not in self.gap_positions:
                            self.gap_positions[chrom] = []
                        
                        self.gap_positions[chrom].append(gap_pos)
                
                return result_data
            else:
                raise FileNotFoundError(f"GAP analysis result file not found: {output_json}")
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"find_gaps.py execution failed: {e.stderr[:500]}")
            raise
    
    def process_existing_coords_file(self) -> Dict[str, Any]:
        """
        Process existing coords file, simulating nucmer comparison results
        Create directory for each chromosome with GAPs and copy coords file
        """
        self.logger.info(f"Using existing coords file: {self.coords_file}")
        
        if not os.path.exists(self.coords_file):
            raise FileNotFoundError(f"Coords file does not exist: {self.coords_file}")
        
        # Validate coords file format
        self._validate_coords_file(self.coords_file)
        
        chromosomes_with_gaps = list(self.gap_positions.keys())
        
        if not chromosomes_with_gaps:
            self.logger.warning("No chromosomes with GAPs found")
            return {"status": "skipped", "reason": "no chromosomes with gaps"}
        
        self.logger.info(f"Preparing coords files for {len(chromosomes_with_gaps)} chromosomes with GAPs")
        
        nucmer_results = {
            "status": "completed",
            "source": "existing_coords_file",
            "coords_file": self.coords_file,
            "chromosome_results": {}
        }
        
        # Create directory for each chromosome with GAPs
        for chrom in chromosomes_with_gaps:
            clean_chrom_name = self._sanitize_filename(chrom)
            chrom_dir = os.path.join(self.output_dir, f"chromosome_{clean_chrom_name}")
            os.makedirs(chrom_dir, exist_ok=True)
            
            # 1. Copy coords file to chromosome directory
            chrom_coords = os.path.join(chrom_dir, f"{clean_chrom_name}_merged.coords")
            shutil.copy2(self.coords_file, chrom_coords)
            
            # 2. Create chromosome FASTA file (if needed)
            # Assuming query genome is already single chromosome file, or user provides split files
            # Try to extract sequence for this chromosome from query genome
            chr_fasta = self._extract_chromosome_fasta(chrom, chrom_dir)
            
            # Record chromosome directory
            self.chromosome_dirs[chrom] = chrom_dir
            nucmer_results["chromosome_results"][chrom] = {
                "dir": chrom_dir,
                "coords_file": chrom_coords,
                "fasta_file": chr_fasta,
                "status": "existing_file_used"
            }
            
            self.logger.info(f"  Chromosome {chrom}: Directory prepared")
            if chr_fasta and os.path.exists(chr_fasta):
                self.logger.info(f"    FASTA file: {os.path.basename(chr_fasta)}")
        
        self.logger.info(f"Coords files set up for {len(chromosomes_with_gaps)} chromosomes with GAPs")
        return nucmer_results
    
    def _validate_coords_file(self, coords_file: str):
        """Validate coords file format"""
        self.logger.info(f"Validating coords file format: {coords_file}")
        
        try:
            with open(coords_file, 'r') as f:
                lines = f.readlines()
            
            # Check if it's a valid coords file
            valid_lines = 0
            for line in lines:
                line = line.strip()
                if line and not line.startswith('/') and not line.startswith('NUCMER') \
                   and not line.startswith('[') and not line.startswith('='):
                    parts = line.split()
                    if len(parts) >= 12:
                        valid_lines += 1
            
            if valid_lines == 0:
                raise ValueError(f"Invalid coords file format: {coords_file}")
            
            self.logger.info(f"  Coords file validation passed, contains {valid_lines} alignment records")
            
        except Exception as e:
            self.logger.error(f"Coords file validation failed: {e}")
            raise ValueError(f"Invalid coords file: {coords_file}")
    
    def _extract_chromosome_fasta(self, chrom_name: str, chrom_dir: str) -> Optional[str]:
        """
        Extract sequence for specific chromosome from query genome
        If query genome is multi-sequence file, try to extract specified chromosome
        """
        clean_chrom_name = self._sanitize_filename(chrom_name)
        output_file = os.path.join(chrom_dir, f"{clean_chrom_name}.fa")
        
        # If already exists, return directly
        if os.path.exists(output_file):
            return output_file
        
        try:
            # Try to extract chromosome sequence using samtools
            # First create index for query genome (if not exists)
            if not os.path.exists(f"{self.query_fasta}.fai"):
                self.logger.debug(f"Creating index for query genome: {self.query_fasta}")
                subprocess.run(["samtools", "faidx", self.query_fasta], 
                             check=True, capture_output=True)
            
            # Try to extract chromosome sequence
            cmd = ["samtools", "faidx", self.query_fasta, chrom_name]
            with open(output_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
            
            if result.returncode == 0 and os.path.getsize(output_file) > 0:
                self.logger.info(f"  Successfully extracted sequence for chromosome {chrom_name}")
                return output_file
            else:
                # If direct extraction fails, try partial matching with sequence IDs
                self.logger.debug(f"Direct extraction for chromosome {chrom_name} failed, trying alternative methods")
                
                # Try reading entire file and finding matching sequences
                try:
                    from Bio import SeqIO
                    sequences = list(SeqIO.parse(self.query_fasta, "fasta"))
                    
                    # Find matching sequences
                    for seq in sequences:
                        if chrom_name in seq.id or seq.id in chrom_name:
                            with open(output_file, 'w') as f:
                                SeqIO.write(seq, f, "fasta")
                            self.logger.info(f"  Using sequence {seq.id} for chromosome {chrom_name}")
                            return output_file
                except ImportError:
                    pass
                
                # If all fail, create empty file or copy entire file
                self.logger.warning(f"  Cannot extract sequence for chromosome {chrom_name}, will use entire query genome")
                shutil.copy2(self.query_fasta, output_file)
                return output_file
                
        except Exception as e:
            self.logger.warning(f"Failed to extract sequence for chromosome {chrom_name}: {e}")
            # If extraction fails, copy entire query genome
            shutil.copy2(self.query_fasta, output_file)
            return output_file
    
    def run_nucmer_comparison(self) -> Dict[str, Any]:
        """
        Run nucmer comparison using run_nucmer_parallel.py
        Only analyze chromosomes with GAPs
        """
        chromosomes_with_gaps = list(self.gap_positions.keys())
        
        if not chromosomes_with_gaps:
            self.logger.warning("No chromosomes with GAPs found, skipping nucmer analysis")
            return {"status": "skipped", "reason": "no chromosomes with gaps"}
        
        self.logger.info(f"Running nucmer comparison for {len(chromosomes_with_gaps)} chromosomes with GAPs")
        self.logger.info(f"Chromosomes: {', '.join(chromosomes_with_gaps)}")
        
        # Call run_nucmer_parallel.py API, passing chromosome list
        try:
            sys.path.insert(0, self.script_dir)
            from run_nucmer_parallel import NucmerAnalyzerAPI
            
            analyzer = NucmerAnalyzerAPI(verbose=self.verbose)
            
            # Key modification: add chromosomes parameter, only analyze chromosomes with GAPs
            result = analyzer.run_analysis(
                reference_fasta=self.ref_fasta,
                query_fasta=self.query_fasta,
                output_dir=self.output_dir,
                output_prefix="nucmer_analysis",
                threads=self.threads,
                overwrite=False,
                chromosomes=chromosomes_with_gaps  # Only analyze these chromosomes
            )
            
            # Record chromosome directory information
            if "chromosome_results" in result:
                for chrom_name, chrom_data in result["chromosome_results"].items():
                    if chrom_name in chromosomes_with_gaps:
                        self.chromosome_dirs[chrom_name] = chrom_data.get("dir", "")
                        self.logger.debug(f"Recorded chromosome directory: {chrom_name} -> {self.chromosome_dirs[chrom_name]}")
            
            # Check if all chromosomes with GAPs have results
            missing_chromosomes = []
            for chrom in chromosomes_with_gaps:
                if chrom not in self.chromosome_dirs:
                    missing_chromosomes.append(chrom)
            
            if missing_chromosomes:
                self.logger.warning(f"These chromosomes with GAPs did not generate nucmer results: {', '.join(missing_chromosomes)}")
                self.logger.warning("Trying to auto-discover chromosome directories...")
                self._discover_chromosome_dirs()
            
            return result
            
        except ImportError as e:
            self.logger.warning(f"Cannot import NucmerAnalyzerAPI: {e}")
            self.logger.warning("Using command line to call run_nucmer_parallel.py, only analyzing chromosomes with GAPs")
            return self._run_nucmer_commandline(chromosomes_with_gaps)
    
    def _run_nucmer_commandline(self, chromosomes_with_gaps: List[str]) -> Dict[str, Any]:
        """Run run_nucmer_parallel.py via command line, only analyze specified chromosomes"""
        # Create temporary configuration file specifying chromosomes
        config_dir = os.path.join(self.output_dir, "nucmer_config")
        os.makedirs(config_dir, exist_ok=True)
        
        config_file = os.path.join(config_dir, "chromosomes_to_process.txt")
        with open(config_file, 'w') as f:
            for chrom in chromosomes_with_gaps:
                f.write(f"{chrom}\n")
        
        cmd = [
            "python", self.script_paths["run_nucmer_parallel.py"],
            "-g", self.ref_fasta,
            "-q", self.query_fasta,
            "-o", self.output_dir,
            "-T", str(self.threads),
            "--chromosomes-file", config_file,  # Use chromosome list file
            "--api-mode"
        ]
        
        if not self.verbose:
            cmd.append("--quiet")
        
        self.logger.info(f"Executing command: {' '.join(cmd[:8])}...")
        self.logger.info(f"Only analyzing these chromosomes: {', '.join(chromosomes_with_gaps)}")
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=86400
            )
            
            # Load report file
            report_file = os.path.join(self.output_dir, "nucmer_analysis_api_report.json")
            if os.path.exists(report_file):
                with open(report_file, 'r') as f:
                    result_data = json.load(f)
                
                # Record chromosome directory information
                if "chromosome_results" in result_data:
                    for chrom_name, chrom_data in result_data["chromosome_results"].items():
                        if chrom_name in chromosomes_with_gaps:
                            self.chromosome_dirs[chrom_name] = chrom_data.get("dir", "")
                            self.logger.debug(f"Recorded chromosome directory: {chrom_name} -> {self.chromosome_dirs[chrom_name]}")
                
                return result_data
            
            # Try to find chromosome directories
            self._discover_chromosome_dirs()
            return {"status": "completed", "note": "Using discovered chromosome directories"}
                    
        except subprocess.CalledProcessError as e:
            self.logger.error(f"run_nucmer_parallel.py execution failed: {e.stderr[:500]}")
            raise
    
    def _discover_chromosome_dirs(self):
        """Discover chromosome directories"""
        self.logger.info("Auto-discovering chromosome directories...")
        
        if not os.path.exists(self.output_dir):
            raise FileNotFoundError(f"Output directory does not exist: {self.output_dir}")
        
        # Find all subdirectories, assuming they are chromosome directories
        chromosomes_with_gaps = list(self.gap_positions.keys())
        
        for item in os.listdir(self.output_dir):
            item_path = os.path.join(self.output_dir, item)
            if os.path.isdir(item_path):
                # Check if it's a chromosome directory (contains .fa file)
                for chrom_name in chromosomes_with_gaps:
                    # Clean chromosome name for matching
                    clean_chrom_name = self._sanitize_filename(chrom_name)
                    if clean_chrom_name in item or chrom_name in item:
                        chr_fasta = os.path.join(item_path, f"{clean_chrom_name}.fa")
                        if os.path.exists(chr_fasta):
                            self.chromosome_dirs[chrom_name] = item_path
                            self.logger.info(f"Discovered chromosome directory: {chrom_name} -> {item_path}")
                            break
        
        if not self.chromosome_dirs:
            self.logger.warning("No chromosome directories found, check if nucmer comparison succeeded")
        else:
            missing_chromosomes = [c for c in chromosomes_with_gaps if c not in self.chromosome_dirs]
            if missing_chromosomes:
                self.logger.warning(f"These chromosomes with GAPs have no directory: {', '.join(missing_chromosomes)}")
    
    def _sanitize_filename(self, filename: str) -> str:
        """Sanitize filename"""
        unsafe_chars = ['<', '>', ':', '"', '/', '\\', '|', '?', '*', ' ']
        for char in unsafe_chars:
            filename = filename.replace(char, '_')
        return filename
    
    def process_chromosomes_with_gaps(self) -> Dict[str, Dict[str, Any]]:
        """
        Process all chromosomes with GAPs
        """
        chromosomes_to_process = []
        
        # Prepare list of chromosomes to process
        for chrom_name, gaps in self.gap_positions.items():
            if chrom_name in self.chromosome_dirs:
                chr_dir = self.chromosome_dirs[chrom_name]
                if os.path.exists(chr_dir):
                    chromosomes_to_process.append((chrom_name, chr_dir, gaps))
                else:
                    self.logger.warning(f"Chromosome directory does not exist: {chr_dir}")
            else:
                self.logger.warning(f"Directory information missing for chromosome {chrom_name}")
        
        if not chromosomes_to_process:
            self.logger.error("No chromosome directories found to process")
            return {}
        
        self.logger.info(f"Preparing to process {len(chromosomes_to_process)} chromosomes with GAPs")
        
        # Process chromosomes in parallel
        results = {}
        
        # Calculate thread allocation
        max_workers = max(1, self.threads * 2 // 3)
        self.logger.info(f"Using {max_workers} threads for parallel chromosome processing")
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_chrom = {}
            for chrom_name, chr_dir, gaps in chromosomes_to_process:
                future = executor.submit(
                    self.process_single_chromosome,
                    chrom_name, chr_dir, gaps
                )
                future_to_chrom[future] = chrom_name
            
            for future in concurrent.futures.as_completed(future_to_chrom):
                chrom_name = future_to_chrom[future]
                try:
                    result = future.result(timeout=86400)
                    results[chrom_name] = result
                    self.logger.info(f"Chromosome {chrom_name} processing completed: {result.get('status', 'unknown')}")
                except Exception as e:
                    self.logger.error(f"Chromosome {chrom_name} processing failed: {e}")
                    results[chrom_name] = {"status": "failed", "error": str(e)}
        
        return results
    
    def process_single_chromosome(
        self, 
        chrom_name: str, 
        chr_dir: str, 
        gap_positions: List[int]
    ) -> Dict[str, Any]:
        """
        Process single chromosome
        """
        self.logger.info(f"Starting to process chromosome {chrom_name}, has {len(gap_positions)} GAPs")
        
        result = {
            "chromosome": chrom_name,
            "directory": chr_dir,
            "total_gaps": len(gap_positions),
            "successful_patches": 0,
            "failed_patches": 0,
            "processed_gaps": [],
            "final_file": None,
            "status": "processing"
        }
        
        try:
            # 1. Prepare file paths
            clean_chrom_name = self._sanitize_filename(chrom_name)
            chr_fasta = os.path.join(chr_dir, f"{clean_chrom_name}.fa")
            
            # Try other possible filenames
            if not os.path.exists(chr_fasta):
                possible_files = [
                    f"{chrom_name}.fa",
                    f"{chrom_name}.fasta",
                    f"{clean_chrom_name}.fasta",
                    "chromosome.fa",
                    "sequence.fa"
                ]
                
                for filename in possible_files:
                    test_path = os.path.join(chr_dir, filename)
                    if os.path.exists(test_path):
                        chr_fasta = test_path
                        self.logger.debug(f"Using alternative file: {filename}")
                        break
            
            # Find coords file
            coords_file = self._find_coords_file(chrom_name, chr_dir)
            
            if not os.path.exists(chr_fasta):
                raise FileNotFoundError(f"Chromosome FASTA file does not exist: {chr_fasta}")
            
            # 2. Create patches directory
            patches_dir = os.path.join(chr_dir, "gap_patches")
            os.makedirs(patches_dir, exist_ok=True)
            
            # 3. Process GAPs from end to beginning
            current_chr_file = chr_fasta
            gap_positions_sorted = sorted(gap_positions, reverse=True)
            
            for i, gap_pos in enumerate(gap_positions_sorted):
                self.logger.debug(f"  Processing GAP {i+1}/{len(gap_positions)}: position {gap_pos:,}")
                
                # Extract patch
                patch_file = self.extract_gap_patch(
                    chrom_name=chrom_name,
                    chr_dir=chr_dir,
                    patches_dir=patches_dir,
                    coords_file=coords_file,
                    gap_position=gap_pos,
                    gap_index=i
                )
                
                if not patch_file or not os.path.exists(patch_file):
                    result["failed_patches"] += 1
                    self.logger.warning(f"    GAP {gap_pos:,}: Patch extraction failed")
                    continue
                
                # Patch GAP
                patched_file = self.patch_single_gap(
                    chrom_name=chrom_name,
                    chr_dir=chr_dir,
                    current_chr_file=current_chr_file,
                    patch_file=patch_file,
                    gap_position=gap_pos,
                    gap_index=i
                )
                
                if patched_file and os.path.exists(patched_file):
                    current_chr_file = patched_file
                    result["successful_patches"] += 1
                    
                    gap_result = {
                        "position": gap_pos,
                        "patch_file": patch_file,
                        "patched_file": patched_file,
                        "status": "success"
                    }
                    result["processed_gaps"].append(gap_result)
                    
                    self.logger.debug(f"    GAP {gap_pos:,}: Patching successful")
                else:
                    result["failed_patches"] += 1
                    self.logger.warning(f"    GAP {gap_pos:,}: Patching failed")
            
            # 4. Record final file
            final_file = os.path.join(chr_dir, f"{clean_chrom_name}_patched.fasta")
            
            if result["successful_patches"] > 0 and os.path.exists(current_chr_file):
                shutil.copy2(current_chr_file, final_file)
                result["final_file"] = final_file
                result["status"] = "completed"
                
                orig_length = self._get_sequence_length(chr_fasta)
                final_length = self._get_sequence_length(final_file)
                result["length_change"] = final_length - orig_length
                
                self.logger.info(f"Chromosome {chrom_name} processing completed: "
                               f"Success {result['successful_patches']}/{result['total_gaps']}, "
                               f"Length change: {result.get('length_change', 0):+,} bp")
            else:
                result["final_file"] = chr_fasta
                result["status"] = "no_patches_applied"
                self.logger.info(f"Chromosome {chrom_name}: No successful patches applied, using original file")
            
            return result
            
        except Exception as e:
            self.logger.error(f"Chromosome {chrom_name} processing exception: {e}")
            result["status"] = "failed"
            result["error"] = str(e)
            return result
    
    def _find_coords_file(self, chrom_name: str, chr_dir: str) -> str:
        """Find coords file"""
        clean_chrom_name = self._sanitize_filename(chrom_name)
        
        # First look for coords file in chromosome directory
        coords_patterns = [
            f"{clean_chrom_name}_merged.coords",
            f"{chrom_name}_merged.coords",
            "*.coords",
            "merged.coords"
        ]
        
        for pattern in coords_patterns:
            import glob
            matches = glob.glob(os.path.join(chr_dir, pattern))
            if matches:
                return matches[0]
        
        # If not in chromosome directory, try using externally provided coords file
        if self.coords_file and os.path.exists(self.coords_file):
            self.logger.info(f"Using external coords file: {os.path.basename(self.coords_file)}")
            return self.coords_file
        
        # If all fail, try looking in parent directory
        parent_dir = os.path.dirname(chr_dir)
        for pattern in coords_patterns:
            import glob
            matches = glob.glob(os.path.join(parent_dir, pattern))
            if matches:
                return matches[0]
        
        raise FileNotFoundError(f"Cannot find coords file in directory {chr_dir}")
    
    def extract_gap_patch(
        self,
        chrom_name: str,
        chr_dir: str,
        patches_dir: str,
        coords_file: str,
        gap_position: int,
        gap_index: int
    ) -> Optional[str]:
        """
        Extract patch sequence for single GAP
        """
        patch_file = os.path.join(patches_dir, f"gap_{gap_position}_patch.fasta")
        
        # Check if valid patch file already exists
        if os.path.exists(patch_file):
            if self._is_valid_patch_file(patch_file):
                self.logger.debug(f"    Using existing patch file: {os.path.basename(patch_file)}")
                return patch_file
        
        # Use absolute script path
        extract_script = self.script_paths["extract_gap_patches.py"]
        
        # Build extraction command
        cmd = [
            "python", extract_script,
            "-c", self.ref_fasta,
            "-coords", coords_file,
            "--gap-position", str(gap_position),
            "-o", patch_file,
            "-f", "10000",
            "--quiet"
        ]
        
        try:
            self.logger.debug(f"    Patch extraction command: {' '.join(cmd[:6])}...")
            
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=600,
                cwd=chr_dir
            )
            
            if os.path.exists(patch_file) and self._is_valid_patch_file(patch_file):
                return patch_file
            else:
                self.logger.warning(f"    Patch file invalid or not generated: {patch_file}")
                return None
                
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr[:200] if e.stderr else str(e)[:200]
            self.logger.warning(f"    Patch extraction failed: {error_msg}")
            
            # Try to check output
            if os.path.exists(patch_file):
                with open(patch_file, 'r') as f:
                    content = f.read()
                    if "error" in content.lower() or "traceback" in content.lower():
                        self.logger.warning(f"    Patch file contains error message")
            
            return None
        except subprocess.TimeoutExpired:
            self.logger.warning(f"    Patch extraction timeout")
            return None
    
    def patch_single_gap(
        self,
        chrom_name: str,
        chr_dir: str,
        current_chr_file: str,
        patch_file: str,
        gap_position: int,
        gap_index: int
    ) -> Optional[str]:
        """
        Patch single GAP
        """
        output_file = os.path.join(chr_dir, f"temp_gap{gap_position}_iter{gap_index}.fasta")
        
        # Use absolute script path
        patch_script = self.script_paths["patch_gap.py"]
        
        # Build patching command
        cmd = [
            "python", patch_script,
            "-r", current_chr_file,
            "-p", patch_file,
            "--gap-position", str(gap_position),
            "-o", output_file,
            "--flank-size", "10000",
            "--minimap2-mode", "asm5",
            "--min-score", "0.7",
            "--min-match-length", "100",
            "--min-mapq", "20",
            "--search-range", "500000",
            "--quiet"
        ]
        
        try:
            self.logger.debug(f"    GAP patching command: {' '.join(cmd[:8])}...")
            
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=600,
                cwd=chr_dir
            )
            
            if os.path.exists(output_file) and self._is_valid_fasta(output_file):
                orig_length = self._get_sequence_length(current_chr_file)
                new_length = self._get_sequence_length(output_file)
                
                if abs(new_length - orig_length) > 100:
                    return output_file
                else:
                    self.logger.warning(f"    No significant length change after patching: {orig_length} -> {new_length}")
                    return None
            else:
                return None
                
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr if e.stderr else e.stdout
            error_str = str(error_msg).lower()
            
            # Check if GAP is already covered
            if "already covered" in error_str or "already filled" in error_str:
                self.logger.debug(f"    GAP {gap_position}: Already covered, skipping")
                return current_chr_file
            elif "no suitable patch" in error_str:
                self.logger.debug(f"    GAP {gap_position}: No suitable patch found")
                return None
            else:
                error_preview = error_msg[:200] if error_msg else str(e)[:200]
                self.logger.warning(f"    GAP patching failed: {error_preview}")
                return None
        except subprocess.TimeoutExpired:
            self.logger.warning(f"    GAP patching timeout")
            return None
    
    def _is_valid_patch_file(self, patch_file: str) -> bool:
        """Check if patch file is valid"""
        if not os.path.exists(patch_file):
            return False
        
        if os.path.getsize(patch_file) < 100:
            return False
        
        try:
            with open(patch_file, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    return False
                second_line = f.readline().strip()
                if not second_line or len(second_line) < 10:
                    return False
            return True
        except:
            return False
    
    def _is_valid_fasta(self, fasta_file: str) -> bool:
        """Check if FASTA file is valid"""
        if not os.path.exists(fasta_file):
            return False
        
        if os.path.getsize(fasta_file) < 100:
            return False
        
        try:
            with open(fasta_file, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>')
        except:
            return False
    
    def _get_sequence_length(self, fasta_file: str) -> int:
        """Get sequence length"""
        try:
            # Simple method
            with open(fasta_file, 'r') as f:
                length = 0
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('>'):
                        length += len(line)
                return length
        except:
            return 0
    
    def merge_patched_chromosomes(self) -> Dict[str, Any]:
        """
        Merge all patched chromosomes
        """
        self.logger.info("Merging all patched chromosomes...")
        
        # Final file placed in controller's parent directory
        base_dir = os.path.dirname(self.output_dir)
        merged_file = os.path.join(base_dir, "final_genome_patched.fasta")
        
        # Collect all patched chromosome files
        patched_files = []
        
        for chrom_name, chr_dir in self.chromosome_dirs.items():
            clean_chrom_name = self._sanitize_filename(chrom_name)
            patched_file = os.path.join(chr_dir, f"{clean_chrom_name}_patched.fasta")
            
            if os.path.exists(patched_file):
                patched_files.append(patched_file)
            else:
                # Try original file
                original_file = os.path.join(chr_dir, f"{clean_chrom_name}.fa")
                if os.path.exists(original_file):
                    patched_files.append(original_file)
                    self.logger.warning(f"Chromosome {chrom_name} using original file")
                else:
                    self.logger.warning(f"Chromosome {chrom_name} file does not exist")
        
        if not patched_files:
            self.logger.error("No chromosome files found for merging")
            return {"status": "failed", "error": "No chromosome files found"}
        
        # Simple FASTA file merging
        try:
            with open(merged_file, 'w') as out_f:
                for file in patched_files:
                    with open(file, 'r') as in_f:
                        out_f.write(in_f.read())
            
            # Count sequences
            seq_count = 0
            with open(merged_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        seq_count += 1
            
            self.logger.info(f"Merging completed: {merged_file} ({seq_count} sequences)")
            self.logger.info(f"File size: {os.path.getsize(merged_file) / (1024*1024):.2f} MB")
            
            return {
                "status": "success",
                "merged_file": merged_file,
                "chromosome_count": seq_count,
                "file_size_mb": os.path.getsize(merged_file) / (1024 * 1024)
            }
            
        except Exception as e:
            self.logger.error(f"Merging failed: {e}")
            return {"status": "failed", "error": str(e)}
    
    def cleanup_temp_files(self):
        """Clean up temporary files"""
        if self.keep_temp:
            self.logger.info("Keeping temporary files")
            return
        
        self.logger.info("Cleaning up temporary files...")
        
        temp_patterns = [
            "temp_*.fasta",
            "gap_*_patch.fasta",
            "*.error.log",
            "*.delta",
            "*.filtered.delta"
        ]
        
        deleted_count = 0
        for chr_dir in self.chromosome_dirs.values():
            if os.path.exists(chr_dir):
                import glob
                for pattern in temp_patterns:
                    for file in glob.glob(os.path.join(chr_dir, pattern)):
                        try:
                            os.remove(file)
                            deleted_count += 1
                        except:
                            pass
        
        self.logger.info(f"Cleaned up {deleted_count} temporary files")
    
    def _create_report(self, start_time: float, stage_results: Dict = None) -> Dict[str, Any]:
        """Create complete report"""
        elapsed_time = time.time() - start_time
        
        total_gaps = sum(len(gaps) for gaps in self.gap_positions.values())
        chromosomes_with_gaps = len(self.gap_positions)
        
        successful_patches = 0
        total_attempted = 0
        
        if stage_results and "gap_patching" in stage_results:
            for chrom_result in stage_results["gap_patching"].values():
                successful_patches += chrom_result.get("successful_patches", 0)
                total_attempted += chrom_result.get("total_gaps", 0)
        
        success_rate = (successful_patches / total_attempted * 100) if total_attempted > 0 else 0
        
        report = {
            "metadata": {
                "controller_version": "1.2",
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "processing_time_seconds": round(elapsed_time, 2),
                "input_files": {
                    "reference_contigs": os.path.basename(self.ref_fasta),
                    "target_genome": os.path.basename(self.query_fasta)
                },
                "parameters": {
                    "threads": self.threads,
                    "min_gap_size": self.min_gap_size,
                    "output_dir": self.output_dir,
                    "keep_temp": self.keep_temp,
                    "script_directory": self.script_dir,
                    "coords_file_provided": self.coords_file is not None
                }
            },
            "statistics": {
                "total_gaps_found": total_gaps,
                "chromosomes_with_gaps": chromosomes_with_gaps,
                "gaps_by_chromosome": {k: len(v) for k, v in self.gap_positions.items()},
                "total_gaps_patched": successful_patches,
                "total_gaps_attempted": total_attempted,
                "patch_success_rate": round(success_rate, 2),
                "chromosomes_processed": len(self.chromosome_dirs)
            },
            "stage_results": stage_results or {},
            "output_files": {
                "final_genome": os.path.join(os.path.dirname(self.output_dir), "final_genome_patched.fasta"),
                "log_file": os.path.join(self.output_dir, "logs", "gap_patcher.log"),
                "gap_analysis_dir": os.path.join(self.output_dir, "gap_analysis")
            },
            "status": "completed"
        }
        
        # Add relevant information if existing coords file was used
        if self.coords_file:
            report["metadata"]["coords_file"] = os.path.basename(self.coords_file)
            if stage_results and "nucmer_comparison" in stage_results:
                report["metadata"]["nucmer_skipped"] = True
                report["metadata"]["nucmer_source"] = "existing_file"
        
        # Save report
        base_dir = os.path.dirname(self.output_dir)
        report_file = os.path.join(base_dir, "gap_patching_report.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        # Print summary
        self._print_summary(report)
        
        return report
    
    def _print_summary(self, report: Dict[str, Any]):
        """Print processing summary"""
        stats = report["statistics"]
        
        self.logger.info("\n" + "="*60)
        self.logger.info("Processing Summary")
        self.logger.info("="*60)
        self.logger.info(f"Processing time: {report['metadata']['processing_time_seconds']:.1f} seconds")
        self.logger.info(f"Total GAPs found: {stats['total_gaps_found']:,}")
        self.logger.info(f"Chromosomes with GAPs: {stats['chromosomes_with_gaps']}")
        self.logger.info(f"Successfully patched GAPs: {stats['total_gaps_patched']:,}/{stats['total_gaps_attempted']:,}")
        self.logger.info(f"Patching success rate: {stats['patch_success_rate']:.1f}%")
        
        # Show whether nucmer alignment was skipped
        if report["metadata"].get("coords_file_provided", False):
            self.logger.info(f"nucmer alignment: Using existing coords file")
            if "coords_file" in report["metadata"]:
                self.logger.info(f"  Coords file: {report['metadata']['coords_file']}")
        
        self.logger.info(f"Final file: {report['output_files']['final_genome']}")
        self.logger.info(f"Detailed report: gap_patching_report.json")
        
        if "stage_results" in report and "gap_patching" in report["stage_results"]:
            self.logger.info("\nChromosome patching details:")
            for chrom_name, chr_result in report["stage_results"]["gap_patching"].items():
                if chr_result.get("status") == "completed":
                    success = chr_result.get("successful_patches", 0)
                    total = chr_result.get("total_gaps", 0)
                    change = chr_result.get("length_change", 0)
                    self.logger.info(f"  {chrom_name}: {success}/{total} successful, length change: {change:+,} bp")
                elif chr_result.get("status") == "failed":
                    self.logger.info(f"  {chrom_name}: Failed - {chr_result.get('error', 'Unknown error')}")
        
        self.logger.info("="*60)


# ==================== Adapter Class ====================

class GapPatcherAdapter:
    """
    Adapter class for gap_patcher_controller.py
    Provides unified controller interface
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        Initialize adapter
        
        Args:
            config: Configuration dictionary
        """
        self.config = config or {}
        self.controller = None
        self.result = None
        self.status = "initialized"
    
    def validate_config(self, config: Dict[str, Any]) -> Tuple[bool, str]:
        """
        Validate configuration parameters
        
        Args:
            config: Configuration dictionary
            
        Returns:
            (is_valid, error_message)
        """
        # Check required parameters
        required_params = ['query_fasta', 'reference_fasta']
        
        for param in required_params:
            if param not in config:
                return False, f"Missing required parameter: {param}"
            
            # Check if files exist
            if param in ['query_fasta', 'reference_fasta']:
                if not os.path.exists(config[param]):
                    return False, f"File does not exist: {config[param]}"
        
        # Check optional coords file
        if 'coords_file' in config and config['coords_file']:
            if not os.path.exists(config['coords_file']):
                return False, f"Coords file does not exist: {config['coords_file']}"
        
        # Check thread count
        if 'threads' in config and config['threads'] <= 0:
            return False, "Thread count must be greater than 0"
        
        return True, ""
    
    def run(self, config: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Run pipeline
        
        Args:
            config: Configuration dictionary, if None use initialized configuration
            
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
                "module": "gap_patcher",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
        
        try:
            self.status = "running"
            
            # Prepare parameters
            params = {
                'ref_fasta': self.config['reference_fasta'],
                'query_fasta': self.config['query_fasta'],
                'threads': self.config.get('threads', 8),
                'output_dir': self.config.get('output_dir', 'gap_patcher_results'),
                'min_gap_size': self.config.get('min_gap_size', 100),
                'keep_temp': self.config.get('keep_temp', False),
                'verbose': self.config.get('verbose', True),
                'coords_file': self.config.get('coords_file', None)  # New
            }
            
            # Create controller instance
            self.controller = GenomeGapPatcherController(**params)
            
            # Run workflow
            self.result = self.controller.run_workflow()
            
            # Standardize result format
            standardized_result = self._standardize_result(self.result)
            
            self.status = "completed"
            return standardized_result
            
        except Exception as e:
            self.status = "failed"
            error_result = {
                "status": "error",
                "error": str(e),
                "module": "gap_patcher",
                "timestamp": time.strftime("%Y%m%d_%H%M%S")
            }
            
            if hasattr(e, '__traceback__'):
                import traceback
                error_result["traceback"] = traceback.format_exc()
            
            return error_result
    
    def _standardize_result(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Standardize result format
        
        Args:
            result: Original result
            
        Returns:
            Standardized result
        """
        if result.get("status") != "completed":
            return {
                "status": "error",
                "error": "Pipeline execution not completed",
                "module": "gap_patcher",
                "raw_result": result
            }
        
        # Extract key information
        standardized = {
            "status": "success",
            "module": "gap_patcher",
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "pipeline_version": "1.2",
            "input_files": {
                "query_fasta": self.config.get('query_fasta', ''),
                "reference_fasta": self.config.get('reference_fasta', ''),
                "coords_file": self.config.get('coords_file', '')
            },
            "output_files": {
                "final_genome": result.get("output_files", {}).get("final_genome", ""),
                "log_file": result.get("output_files", {}).get("log_file", ""),
                "gap_analysis_dir": result.get("output_files", {}).get("gap_analysis_dir", ""),
                "report_file": "gap_patching_report.json"
            },
            "statistics": {
                "gap_finding": result.get("statistics", {}),
                "gap_patching": self._extract_patching_stats(result)
            },
            "processing_info": {
                "processing_time": result.get("metadata", {}).get("processing_time_seconds", 0),
                "threads_used": result.get("metadata", {}).get("parameters", {}).get("threads", 0),
                "timestamp": result.get("metadata", {}).get("timestamp", ""),
                "nucmer_skipped": result.get("metadata", {}).get("nucmer_skipped", False)
            }
        }
        
        # Add raw result (optional)
        if self.config.get('include_raw_result', False):
            standardized["raw_result"] = result
        
        return standardized
    
    def _extract_patching_stats(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Extract patching statistics"""
        stats = result.get("statistics", {})
        
        return {
            "total_gaps_found": stats.get("total_gaps_found", 0),
            "chromosomes_with_gaps": stats.get("chromosomes_with_gaps", 0),
            "total_gaps_patched": stats.get("total_gaps_patched", 0),
            "total_gaps_attempted": stats.get("total_gaps_attempted", 0),
            "patch_success_rate": stats.get("patch_success_rate", 0),
            "chromosomes_processed": stats.get("chromosomes_processed", 0)
        }
    
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
    
    def cleanup(self) -> bool:
        """
        Clean up temporary files
        
        Returns:
            Whether successful
        """
        try:
            if self.controller:
                self.controller.cleanup_temp_files()
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
            "reference_fasta": "",  # Required
            "coords_file": "",  # Optional, existing coords file
            "threads": 8,
            "output_dir": "gap_patcher_results",
            "min_gap_size": 100,
            "keep_temp": False,
            "verbose": True,
            "include_raw_result": False
        }
    
    def get_module_info(self) -> Dict[str, Any]:
        """
        Get module information
        
        Returns:
            Module information dictionary
        """
        return {
            "name": "GapPatcherAdapter",
            "version": "1.2",
            "description": "GAP patching controller adapter based on nucmer alignment, supports using existing coords files directly",
            "capabilities": [
                "GAP position detection",
                "nucmer synteny analysis (can be skipped using existing results)",
                "Parallel chromosome processing",
                "Patch extraction and patching",
                "Result merging and report generation"
            ],
            "dependencies": ["find_gaps", "run_nucmer_parallel", "extract_gap_patches", "patch_gap"]
        }


def get_controller_api(config: Dict[str, Any] = None):
    """
    Get controller API instance (factory function)
    
    Args:
        config: Configuration dictionary
        
    Returns:
        GapPatcherAdapter instance
    """
    return GapPatcherAdapter(config)


# ==================== Command Line Interface ====================

def main():
    parser = argparse.ArgumentParser(
        description='Genome GAP Patching Main Controller - Automated Workflow',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage Examples:
  # Basic usage (run complete workflow)
  python gap_patcher_controller.py -c ref_contigs.fasta -q target_genome.fasta -t 16
  
  # Use existing coords file, skip nucmer alignment
  python gap_patcher_controller.py -c ref.fa -q genome.fa --coords-file existing_nucmer.coords
  
  # Specify output directory, keep temporary files
  python gap_patcher_controller.py -c ref.fa -q genome.fa -o results --keep-temp
  
  # Minimum GAP size 500bp, quiet mode
  python gap_patcher_controller.py -c ref.fa -q genome.fa -t 8 --min-gap-size 500 --quiet
  
  # Adapter mode
  python gap_patcher_controller.py --adapter-mode --config config.json
  
  # Show module information
  python gap_patcher_controller.py --module-info
        """
    )
    
    # Required parameters
    parser.add_argument('-c', '--contigs', required=False, 
                       help='Reference contigs file path')
    parser.add_argument('-q', '--query', required=False,
                       help='Target genome file path to be patched')
    
    # New: coords file parameter (optional)
    parser.add_argument('--coords-file', type=str,
                       help='Existing nucmer alignment result coords file (if provided, skip nucmer alignment)')
    
    # Optional parameters
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='Total threads, default 8')
    parser.add_argument('-o', '--output-dir', default='results',
                       help='Output directory, default "results"')
    parser.add_argument('--min-gap-size', type=int, default=100,
                       help='Minimum GAP size (bp), default 100')
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files, default cleanup')
    parser.add_argument('--quiet', action='store_true',
                       help='Quiet mode, reduce output')
    
    # Adapter mode parameters
    parser.add_argument('--adapter-mode', action='store_true',
                       help='Use adapter mode')
    parser.add_argument('--config', type=str,
                       help='JSON configuration file path (adapter mode)')
    parser.add_argument('--module-info', action='store_true',
                       help='Show module information')
    
    args = parser.parse_args()
    
    # Show module information
    if args.module_info:
        adapter = GapPatcherAdapter()
        info = adapter.get_module_info()
        print(json.dumps(info, indent=2, ensure_ascii=False))
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
                print("Adapter mode requires --query and --contigs parameters")
                sys.exit(1)
                
            config = {
                'query_fasta': args.query,
                'reference_fasta': args.contigs,
                'coords_file': args.coords_file,  # New
                'threads': args.threads,
                'output_dir': args.output_dir,
                'min_gap_size': args.min_gap_size,
                'keep_temp': args.keep_temp,
                'verbose': not args.quiet
            }
        
        # Create adapter and run
        adapter = GapPatcherAdapter(config)
        result = adapter.run()
        
        # Output result
        print(json.dumps(result, indent=2, ensure_ascii=False))
        
        if result.get('status') == 'success':
            sys.exit(0)
        else:
            sys.exit(1)
    
    # Traditional command line mode
    else:
        # Check required parameters
        if not args.query or not args.contigs:
            print("Error: Must provide -q/--query and -c/--contigs parameters")
            sys.exit(1)
        
        try:
            controller = GenomeGapPatcherController(
                ref_fasta=args.contigs,
                query_fasta=args.query,
                threads=args.threads,
                output_dir=args.output_dir,
                min_gap_size=args.min_gap_size,
                keep_temp=args.keep_temp,
                verbose=not args.quiet,
                coords_file=args.coords_file  # New parameter
            )
            
            report = controller.run_workflow()
            
            final_file = report["output_files"]["final_genome"]
            if os.path.exists(final_file):
                print(f"\n✅ Patching completed!")
                print(f"   Final file: {final_file}")
                print(f"   Report file: gap_patching_report.json")
                
                if not args.quiet:
                    stats = report["statistics"]
                    print(f"\n📊 Statistics summary:")
                    print(f"   Total GAPs: {stats['total_gaps_found']:,}")
                    print(f"   Successfully patched: {stats['total_gaps_patched']:,} ({stats['patch_success_rate']:.1f}%)")
                    print(f"   Chromosomes: {stats['chromosomes_processed']}")
                    
                    # Show whether nucmer alignment was skipped
                    if args.coords_file:
                        print(f"   nucmer alignment: Using existing file {os.path.basename(args.coords_file)}")
                    else:
                        print(f"   nucmer alignment: Re-run")
                    
            else:
                print(f"\n⚠️  Processing completed but no final file generated")
                print(f"   Please check log: {args.output_dir}/logs/gap_patcher.log")
                
        except KeyboardInterrupt:
            print("\n❌ Processing interrupted by user")
            sys.exit(1)
        except Exception as e:
            print(f"\n❌ Error: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

if __name__ == "__main__":
    main()