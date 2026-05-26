#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import argparse
import re
import shutil
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
import json


@dataclass
class QualityConfig:
    """Quality assessment configuration"""
    reads_files: List[str] = field(default_factory=list)
    contigs_file: str = ""
    
    meryl_path: str = "meryl"
    merqury_path: str = "merqury.sh"
    craq_path: str = "craq"
    
    # External Meryl database path
    external_meryl_db: Optional[str] = None
    
    kmer_size: int = 21
    threads: int = 32
    memory_gb: int = 80
    
    min_qv: float = 30.0
    min_length: int = 100000
    
    output_prefix: str = "high_quality_contigs"
    keep_intermediate: bool = True
    force_rerun: bool = False
    
    working_dir: Optional[str] = None
    skip_meryl: bool = False
    skip_merqury: bool = False
    skip_craq: bool = False
    only_analyze: bool = False


class QualityPipeline:
    """Quality assessment pipeline class, supports modular calls"""
    
    def __init__(self, config: Optional[QualityConfig] = None):
        self.config = config or QualityConfig()
        self.results = {}
        self.output_dir = None
        
    def run(self) -> Dict[str, Any]:
        """Run complete quality assessment pipeline"""
        print("\n" + "="*80)
        print("CONTIG QUALITY PIPELINE - MODULE VERSION")
        print("="*80)
        
        if self.config.working_dir:
            os.chdir(self.config.working_dir)
            print(f"Working directory set to: {os.getcwd()}")
        
        if not self.validate_inputs():
            return {
                'success': False,
                'error': 'Input validation failed',
                'results': {}
            }
        
        self.output_dir = Path.cwd() / f"{self.config.output_prefix}_quality"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {self.output_dir}")
        
        try:
            # Now returns 5 values including meryl_db
            has_meryl, has_merqury, existing_qv_scores, existing_craq_scores, meryl_db = \
                self.check_existing_results()
            
            # Run meryl only if no existing database found
            if not self.config.skip_meryl and not self.config.only_analyze:
                if meryl_db is None:  # Only run if no existing database
                    meryl_db = self.run_meryl_count_if_needed(has_meryl)
            
            qv_scores = {}
            if not self.config.skip_merqury and not self.config.only_analyze:
                if meryl_db and Path(meryl_db).exists():
                    qv_scores = self.run_merqury_if_needed(meryl_db, has_merqury, existing_qv_scores)
            else:
                qv_scores = existing_qv_scores if existing_qv_scores else {}
            
            craq_scores = {}
            if not self.config.skip_craq and not self.config.only_analyze:
                craq_scores = self.run_craq_if_needed(existing_craq_scores)
            else:
                craq_scores = existing_craq_scores if existing_craq_scores else {}
            
            sequences = self.parse_fasta_sequences(Path(self.config.contigs_file))
            if not sequences:
                return {
                    'success': False,
                    'error': 'Failed to parse contigs file',
                    'results': {}
                }
            
            if not qv_scores and not craq_scores:
                print("\n⚠ WARNING: No quality scores found!")
                return {
                    'success': False,
                    'error': 'No quality scores found for filtering',
                    'results': {}
                }
            
            filtered_contigs = self.filter_and_score_contigs(
                sequences, qv_scores, craq_scores
            )
            
            output_file = self.output_dir / f"{self.config.output_prefix}.fasta"
            write_success = self.write_high_quality_contigs(filtered_contigs, output_file)
            
            self.generate_summary_report(
                sequences, filtered_contigs, qv_scores, craq_scores
            )
            
            final_output = Path.cwd() / f"{self.config.output_prefix}.fasta"
            final_output = final_output.resolve()
            
            if output_file.exists() and output_file.stat().st_size > 0:
                self.create_final_output_link(output_file, final_output)
            
            self.results = {
                'success': True,
                'output_dir': str(self.output_dir),
                'output_file': str(final_output),
                'filtered_contigs_count': len(filtered_contigs),
                'total_length': sum(c['length'] for c in filtered_contigs) if filtered_contigs else 0,
                'qv_scores_count': len(qv_scores),
                'craq_scores_count': len(craq_scores),
                'filtered_contigs': filtered_contigs[:10] if filtered_contigs else [],
                'meryl_db_used': str(meryl_db) if meryl_db else None
            }
            
            return self.results
            
        except Exception as e:
            import traceback
            print(f"\n✗ Pipeline failed with error: {e}")
            traceback.print_exc()
            return {
                'success': False,
                'error': str(e),
                'traceback': traceback.format_exc(),
                'results': self.results
            }
    
    def validate_inputs(self) -> bool:
        """Validate input files"""
        # Can skip reads files if external Meryl database is provided
        if not self.config.only_analyze and not self.config.reads_files and not self.config.external_meryl_db:
            print("Error: Reads files or external Meryl database required, or use --only-analyze mode")
            return False
        
        contigs_path = Path(self.config.contigs_file)
        if not contigs_path.exists():
            print(f"Error: Contigs file does not exist: {self.config.contigs_file}")
            return False
        
        # Validate external Meryl database
        if self.config.external_meryl_db:
            external_db = Path(self.config.external_meryl_db)
            if not external_db.exists():
                print(f"Warning: External Meryl database does not exist: {self.config.external_meryl_db}")
            else:
                print(f"External Meryl database: {external_db}")
        
        if self.config.reads_files:
            for rf in self.config.reads_files:
                if not os.path.exists(rf):
                    print(f"Warning: Reads file does not exist: {rf}")
        
        print(f"Contigs: {self.config.contigs_file}")
        print(f"Reads: {len(self.config.reads_files)} files")
        print(f"External Meryl DB: {self.config.external_meryl_db or 'None'}")
        print(f"Output prefix: {self.config.output_prefix}")
        return True
    
    def check_existing_results(self) -> Tuple[bool, bool, Dict, Dict, Optional[Path]]:
        """Check for existing result files"""
        print("\n" + "="*60)
        print("CHECKING FOR EXISTING RESULTS")
        print("="*60)
        
        has_meryl = False
        has_merqury = False
        qv_scores = {}
        craq_scores = {}
        meryl_db = None
        
        # Priority 1: Use external Meryl database if specified
        if self.config.external_meryl_db:
            external_db = Path(self.config.external_meryl_db)
            if external_db.exists() and not self.config.force_rerun:
                print(f"✓ Using external meryl database: {external_db}")
                meryl_db = external_db
                has_meryl = True
            elif not external_db.exists():
                print(f"✗ External meryl database does not exist: {external_db}")
        
        # Priority 2: Check locally generated database
        if not has_meryl:
            local_meryl_db = self.output_dir / f"reads_{self.config.kmer_size}mer.meryl"
            if local_meryl_db.exists() and not self.config.force_rerun:
                print(f"✓ Found existing local meryl database: {local_meryl_db}")
                meryl_db = local_meryl_db
                has_meryl = True
            else:
                print(f"✗ No local meryl database found or force rerun enabled")
        
        # Check Merqury QV file
        qv_file = self.output_dir / "merqury_qv.txt"
        if qv_file.exists() and not self.config.force_rerun:
            print(f"✓ Found existing merqury QV file: {qv_file}")
            has_merqury = True
            qv_scores = self.parse_merqury_qv_file(qv_file)
        else:
            for file in self.output_dir.glob("*.qv"):
                if file.exists() and not self.config.force_rerun:
                    print(f"✓ Found alternative QV file: {file}")
                    has_merqury = True
                    qv_scores = self.parse_merqury_qv_file(file)
                    shutil.copy2(file, qv_file)
                    break
        
        # Check CRAQ results
        craq_dir = self.output_dir / "craq_output"
        if craq_dir.exists() and not self.config.force_rerun:
            print(f"✓ Found existing CRAQ directory: {craq_dir}")
            craq_scores = self.find_and_parse_craq_report(craq_dir)
            if craq_scores:
                print(f"  Successfully parsed {len(craq_scores)} CRAQ scores")
            else:
                print(f"  Could not parse existing CRAQ scores")
        
        return has_meryl, has_merqury, qv_scores, craq_scores, meryl_db
    
    def run_meryl_count_if_needed(self, has_meryl: bool) -> Optional[Path]:
        """Run meryl if needed"""
        meryl_db = self.output_dir / f"reads_{self.config.kmer_size}mer.meryl"
        
        if has_meryl and not self.config.force_rerun:
            print(f"✓ Skipping meryl (using existing database)")
            return meryl_db.resolve()
        
        if not self.config.reads_files:
            print("✗ Cannot run meryl: no reads files provided")
            return None
        
        print(f"Running meryl count (k={self.config.kmer_size})...")
        
        cmd = [
            self.config.meryl_path,
            f"k={self.config.kmer_size}",
            "count",
            f"memory={self.config.memory_gb}G",
            f"threads={self.config.threads}",
            "output", str(meryl_db)
        ] + self.config.reads_files
        
        log_file = self.output_dir / "meryl_count.log"
        if self.run_command(cmd, "meryl count", log_file):
            return meryl_db.resolve()
        return None
    
    def run_merqury_if_needed(self, meryl_db: Path, has_merqury: bool, 
                             existing_qv_scores: Dict) -> Dict[str, float]:
        """Run merqury if needed"""
        if has_merqury and not self.config.force_rerun and existing_qv_scores:
            print(f"\n✓ Skipping merqury (using existing QV scores)")
            print(f"  Found {len(existing_qv_scores)} existing QV scores")
            return existing_qv_scores
        
        if meryl_db is None:
            print("✗ No meryl database available for merqury")
            return {}
        
        print("\n" + "="*60)
        print("RUNNING MERQURY")
        print("="*60)
        
        import time
        timestamp = int(time.time())
        temp_dir = self.output_dir / f"merqury_temp_{timestamp}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        meryl_link = temp_dir / "reads.meryl"
        contigs_link = temp_dir / "contigs.fasta"
        
        meryl_db = meryl_db.resolve()
        contigs_file = Path(self.config.contigs_file).resolve()
        
        if not meryl_db.exists():
            print(f"✗ Meryl database does not exist: {meryl_db}")
            return {}
        
        if not contigs_file.exists():
            print(f"✗ Contigs file does not exist: {contigs_file}")
            return {}
        
        print(f"Creating link: {meryl_db} -> {meryl_link}")
        print(f"Creating link: {contigs_file} -> {contigs_link}")
        
        meryl_success = self.create_symlink_with_fallback(meryl_db, meryl_link, "meryl database")
        contigs_success = self.create_symlink_with_fallback(contigs_file, contigs_link, "contigs file")
        
        if not meryl_success or not contigs_success:
            print("✗ Cannot create necessary links, merqury cannot run")
            return {}
        
        cmd = [
            self.config.merqury_path,
            "reads.meryl",
            "contigs.fasta",
            "merqury_out"
        ]
        
        log_file = self.output_dir / "merqury.log"
        
        print(f"\nRunning merqury in: {temp_dir}")
        print(f"Working directory: {temp_dir}")
        print(f"Current directory: {os.getcwd()}")
        
        original_cwd = os.getcwd()
        os.chdir(str(temp_dir))
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        os.chdir(original_cwd)
        
        with open(log_file, 'w') as f:
            f.write(f"Working directory: {temp_dir}\n")
            f.write(f"Return code: {result.returncode}\n")
            f.write(f"STDOUT:\n{result.stdout}\n")
            f.write(f"STDERR:\n{result.stderr}\n")
        
        if result.returncode == 0:
            print("✓ merqury completed successfully")
            
            qv_file = temp_dir / "merqury_out.contigs.qv"
            if qv_file.exists():
                print(f"Found QV file: {qv_file}")
                
                qv_file_final = self.output_dir / "merqury_qv.txt"
                shutil.copy2(qv_file, qv_file_final)
                print(f"Copied to: {qv_file_final}")
                
                qv_scores = self.parse_merqury_qv_file(qv_file)
                return qv_scores
            else:
                print("⚠ QV file not found in expected location")
                for file in temp_dir.rglob("*.qv"):
                    if file.exists():
                        print(f"Found alternative QV file: {file}")
                        qv_scores = self.parse_merqury_qv_file(file)
                        if qv_scores:
                            qv_file_final = self.output_dir / "merqury_qv.txt"
                            shutil.copy2(file, qv_file_final)
                            return qv_scores
        else:
            print(f"✗ merqury failed with return code {result.returncode}")
            if result.stderr:
                print(f"Error: {result.stderr[:500]}")
        
        if not self.config.keep_intermediate:
            try:
                shutil.rmtree(temp_dir, ignore_errors=True)
            except:
                pass
        
        return {}
    
    def run_craq_if_needed(self, existing_craq_scores: Dict) -> Dict[str, float]:
        """Run CRAQ if needed"""
        if existing_craq_scores and not self.config.force_rerun:
            print(f"\n✓ Skipping CRAQ (using existing scores)")
            print(f"  Found {len(existing_craq_scores)} existing CRAQ scores")
            return existing_craq_scores
        
        print("\n" + "="*60)
        print("RUNNING CRAQ")
        print("="*60)
        
        if '/' in self.config.output_prefix:
            base_prefix = Path(self.config.output_prefix).name
        else:
            base_prefix = self.config.output_prefix
        
        craq_dir = self.output_dir / "craq_output"
        
        if craq_dir.exists() and self.config.force_rerun:
            try:
                shutil.rmtree(craq_dir, ignore_errors=True)
            except Exception as e:
                print(f"✗ Cannot remove existing CRAQ directory: {e}")
        
        craq_dir.mkdir(parents=True, exist_ok=True)
        
        if not self.config.reads_files:
            print("✗ Cannot run CRAQ: no reads files")
            return {}
        
        reads_input = self.config.reads_files[0] if len(self.config.reads_files) == 1 else " ".join(self.config.reads_files)
        
        cmd = [
            self.config.craq_path,
            "-g", os.path.abspath(self.config.contigs_file),
            "-sms", reads_input,
            "-t", str(self.config.threads),
            "-o", str(craq_dir / base_prefix)
        ]
        
        log_file = self.output_dir / "craq.log"
        
        print(f"Running CRAQ command: {' '.join(cmd)}")
        print(f"Contigs file: {os.path.abspath(self.config.contigs_file)}")
        print(f"Reads input: {reads_input}")
        print(f"Output prefix: {base_prefix} (simplified from {self.config.output_prefix})")
        
        if self.run_command(cmd, "CRAQ", log_file):
            return self.find_and_parse_craq_report(craq_dir)
        
        return {}
    
    def filter_and_score_contigs(
        self,
        sequences: Dict[str, Dict],
        qv_scores: Dict[str, float],
        craq_scores: Dict[str, float]
    ) -> List[Dict]:
        """Filter contigs and calculate composite scores
        Modified: Keep contig if EITHER length OR QV meets threshold
        """
        filtered_contigs = []
        
        print(f"\n{'='*60}")
        print("FILTERING CONTIGS")
        print(f"Min QV: {self.config.min_qv}, Min length: {self.config.min_length:,} bp")
        print(f"Mode: Keep if EITHER length OR QV meets threshold")
        print("="*60)
        
        stats = {
            'total': len(sequences),
            'passed_length': 0,
            'passed_qv': 0,
            'final': 0
        }
        
        for contig_id, contig_data in sequences.items():
            length = contig_data['length']
            qv = qv_scores.get(contig_id, 0.0)
            
            # Check if criteria are met
            length_pass = length >= self.config.min_length
            qv_pass = qv >= self.config.min_qv
            
            # Only filter if BOTH criteria fail
            if not length_pass and not qv_pass:
                continue
            
            # Count passing criteria
            if length_pass:
                stats['passed_length'] += 1
            if qv_pass:
                stats['passed_qv'] += 1
            
            craq = craq_scores.get(contig_id, 0.0)
            total_score = qv + craq
            
            filtered_contigs.append({
                'id': contig_id,
                'header': contig_data['header'],
                'sequence': contig_data['sequence'],
                'length': length,
                'qv': qv,
                'craq': craq,
                'total_score': total_score,
                'passed_length': length_pass,  # For debugging
                'passed_qv': qv_pass            # For debugging
            })
            stats['final'] += 1
        
        filtered_contigs.sort(key=lambda x: x['total_score'], reverse=True)
        
        print(f"\nFiltering results:")
        print(f"  Total contigs: {stats['total']}")
        print(f"  Passed length filter: {stats['passed_length']}")
        print(f"  Passed QV filter: {stats['passed_qv']}")
        print(f"  Final retained (met at least one criterion): {stats['final']}")
        
        if filtered_contigs:
            print(f"  Best score: {filtered_contigs[0]['total_score']:.2f}")
            print(f"  Worst retained: {filtered_contigs[-1]['total_score']:.2f}")
            
            # Show breakdown of why contigs were retained
            only_length = sum(1 for c in filtered_contigs if c['passed_length'] and not c['passed_qv'])
            only_qv = sum(1 for c in filtered_contigs if not c['passed_length'] and c['passed_qv'])
            both = sum(1 for c in filtered_contigs if c['passed_length'] and c['passed_qv'])
            print(f"\nRetention breakdown:")
            print(f"  Met both criteria: {both}")
            print(f"  Met only length criterion: {only_length}")
            print(f"  Met only QV criterion: {only_qv}")
        
        return filtered_contigs
    
    def run_command(self, cmd: List[str], description: str, log_file: Optional[Path] = None) -> bool:
        """Run command"""
        print(f"\n{'='*60}")
        print(f"Running: {description}")
        print(f"Command: {' '.join(cmd)}")
        print("="*60)
        
        try:
            if log_file:
                with open(log_file, 'w') as log:
                    result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
            else:
                result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"✓ {description} completed successfully")
                return True
            else:
                print(f"✗ {description} failed with return code {result.returncode}")
                if result.stderr:
                    print(f"Error output: {result.stderr[:500]}")
                return False
        except Exception as e:
            print(f"✗ {description} failed with exception: {e}")
            return False
    
    def create_symlink_with_fallback(self, source: Path, target: Path, description: str = "file") -> bool:
        """Create symlink, fallback to copy"""
        try:
            source = source.resolve()
            if not source.exists():
                print(f"✗ {description} does not exist: {source}")
                return False
            
            if target.exists():
                try:
                    target.unlink()
                except Exception as e:
                    print(f"✗ Cannot remove existing target: {target}, error: {e}")
                    return False
            
            target.parent.mkdir(parents=True, exist_ok=True)
            
            try:
                os.symlink(str(source), str(target))
                print(f"✓ {description} symlink successful: {target} -> {source}")
                return True
            except Exception as symlink_error:
                print(f"✗ {description} symlink failed: {symlink_error}")
                
                try:
                    shutil.copy2(str(source), str(target))
                    print(f"✓ {description} copy successful: {target}")
                    return True
                except Exception as copy_error:
                    print(f"✗ {description} copy failed: {copy_error}")
                    return False
                    
        except Exception as e:
            print(f"✗ Error processing {description}: {e}")
            return False
    
    def parse_merqury_qv_file(self, qv_file: Path) -> Dict[str, float]:
        """Parse merqury QV file
        Modified: +inf is converted to 0.0 instead of 40.0
        """
        qv_scores = {}
        inf_count = 0
        error_count = 0
        
        try:
            with open(qv_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 4:
                        contig = parts[0]
                        qv_str = parts[3]
                        
                        # Handle +inf specially
                        if qv_str == "+inf":
                            qv = 0.0  # MODIFIED: +inf now becomes 0.0 instead of 40.0
                            inf_count += 1
                        else:
                            try:
                                qv = float(qv_str)
                            except ValueError:
                                qv = 0.0
                                error_count += 1
                        
                        qv_scores[contig] = qv
            
            print(f"Parsed QV scores for {len(qv_scores)} contigs")
            print(f"  Note: '+inf' values are now converted to 0.0 (was 40.0 in original)")
            
            if inf_count > 0:
                print(f"  {inf_count} contigs had '+inf' (converted to 0.0)")
            if error_count > 0:
                print(f"  {error_count} contigs had invalid values (converted to 0.0)")
            
            if qv_scores:
                print("Sample QV scores (first 5):")
                for i, (contig, qv) in enumerate(list(qv_scores.items())[:5]):
                    print(f"  {contig}: {qv:.2f}")
            
            return qv_scores
            
        except Exception as e:
            print(f"Error parsing QV file: {e}")
            return {}
    
    def find_and_parse_craq_report(self, craq_dir: Path) -> Dict[str, float]:
        """Find and parse CRAQ out_final.Report file - extracts scores from parentheses"""
        print(f"\nLooking for CRAQ report in: {craq_dir}")
        
        possible_paths = []
        
        # Priority: look for out_final.Report files
        for pattern in ["**/out_final.Report", "**/*.Report", "**/*final*report*"]:
            for file_path in craq_dir.glob(pattern):
                if file_path.exists() and file_path.stat().st_size > 0:
                    if file_path not in possible_paths:
                        possible_paths.append(file_path)
                        print(f"  Found report candidate: {file_path.relative_to(craq_dir)}")
        
        # Common paths
        common_paths = [
            craq_dir / "runAQI_out" / "out_final.Report",  # Primary path
            craq_dir / "out_final.Report",
            craq_dir / "runAQI_out" / "out_regional.Report",  # Alternative
            craq_dir / "out_regional.Report"
        ]
        
        for path in common_paths:
            if path.exists() and path.stat().st_size > 0 and path not in possible_paths:
                possible_paths.append(path)
                print(f"  Found common report: {path.relative_to(craq_dir)}")
        
        # Sort by priority: out_final.Report first
        possible_paths.sort(key=lambda x: 0 if "out_final.Report" in str(x) else 1)
        
        print(f"Checking {len(possible_paths)} possible report files...")
        
        for report_file in possible_paths:
            print(f"\nTrying to parse: {report_file.relative_to(craq_dir)}")
            scores = self.parse_craq_report_file(report_file)
            if scores:
                print(f"✓ Successfully parsed {len(scores)} scores from {report_file.name}")
                return scores
        
        print("⚠ No valid CRAQ report file found or all files empty!")
        
        # Debug: list directory contents
        print(f"\nDirectory contents of {craq_dir}:")
        for item in craq_dir.rglob("*"):
            if item.is_file():
                size = item.stat().st_size
                print(f"  {item.relative_to(craq_dir)} ({size} bytes)")
        
        return {}
    
    def parse_craq_report_file(self, report_file: Path) -> Dict[str, float]:
        """Parse CRAQ Report file - extracts score from parentheses in second last column"""
        # Import re inside the function to ensure it's available
        import re
        craq_scores = {}
        
        try:
            with open(report_file, 'r') as f:
                lines = f.readlines()
            
            print(f"Parsing CRAQ report: {report_file}")
            print(f"Total lines: {len(lines)}")
            
            # Find where data starts
            start_processing = False
            for line_num, line in enumerate(lines):
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Find header line
                if line.startswith('#Chr') or 'Covered.Rate' in line or 'Avg.CRE(R-AQI)' in line:
                    start_processing = True
                    print(f"Found header at line {line_num + 1}: {line}")
                    continue
                
                # Process data lines
                if start_processing:
                    # Skip comment lines
                    if line.startswith('#'):
                        continue
                    
                    # Split line by whitespace (tabs or spaces)
                    parts = re.split(r'\s+', line)
                    if len(parts) < 6:  # Need at least 6 columns
                        continue
                    
                    contig = parts[0]
                    
                    # Skip summary lines
                    if contig in ["Genome", "Total", "Summary"]:
                        continue
                    
                    # Get second last column (Avg.CRE(R-AQI) column)
                    # Format: "0.115956336230135(98.8471336632426)"
                    second_last_col = parts[-2].strip()
                    
                    # Extract score from inside parentheses
                    score = None
                    
                    # Method 1: Extract number inside parentheses
                    match = re.search(r'\(([^)]+)\)', second_last_col)
                    if match:
                        try:
                            score = float(match.group(1))
                            if len(craq_scores) < 3:  # Print first few examples
                                print(f"  Extracted from parentheses: {second_last_col} -> {score}")
                        except ValueError:
                            pass
                    
                    # Method 2: If method 1 fails, try to find any number in parentheses
                    if score is None:
                        # Look for pattern like (number)
                        paren_match = re.findall(r'\(([^)]+)\)', second_last_col)
                        if paren_match:
                            try:
                                score = float(paren_match[-1])  # Take the last parentheses content
                                print(f"  Extracted using fallback: {second_last_col} -> {score}")
                            except ValueError:
                                pass
                    
                    # Save score if successfully extracted
                    if score is not None:
                        craq_scores[contig] = score
                        if len(craq_scores) <= 5:  # Print first 5 examples
                            print(f"  Example: {contig} -> {score:.4f} (from '{second_last_col}')")
            
            if craq_scores:
                print(f"\n✓ Parsed CRAQ scores for {len(craq_scores)} contigs")
                
                # Statistics
                scores_list = list(craq_scores.values())
                min_score = min(scores_list)
                max_score = max(scores_list)
                avg_score = sum(scores_list) / len(scores_list)
                print(f"  Score range: {min_score:.4f} - {max_score:.4f}, average: {avg_score:.4f}")
                
                # Show sample scores
                print("\nSample CRAQ scores (from parentheses, rounded to 2 decimal places):")
                for contig, score in list(craq_scores.items())[:5]:
                    print(f"  {contig}: {score:.2f}")
                
                return craq_scores
            else:
                print("⚠ No valid CRAQ scores found in report")
                
                # Print first 10 lines for debugging
                print("\nFirst 10 lines of file for debugging:")
                for i, line in enumerate(lines[:10]):
                    print(f"  Line {i+1}: {line.strip()}")
                    
        except Exception as e:
            print(f"Error parsing CRAQ report {report_file}: {e}")
            import traceback
            traceback.print_exc()
        
        return {}
    
    def parse_fasta_sequences(self, fasta_file: Path) -> Dict[str, Dict]:
        """Parse FASTA sequences"""
        sequences = {}
        
        try:
            with open(fasta_file, 'r') as f:
                current_contig = None
                current_seq = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_contig:
                            seq_str = ''.join(current_seq)
                            sequences[current_contig] = {
                                'header': current_contig,
                                'sequence': seq_str,
                                'length': len(seq_str)
                            }
                        
                        current_contig = line[1:].split()[0] if ' ' in line[1:] else line[1:]
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                if current_contig:
                    seq_str = ''.join(current_seq)
                    sequences[current_contig] = {
                        'header': current_contig,
                        'sequence': seq_str,
                        'length': len(seq_str)
                    }
            
            print(f"Parsed {len(sequences)} contigs from {fasta_file}")
            return sequences
            
        except Exception as e:
            print(f"Error parsing FASTA file {fasta_file}: {e}")
            return {}
    
    def write_high_quality_contigs(
        self,
        contigs: List[Dict],
        output_file: Path
    ) -> bool:
        """Write high-quality contigs"""
        if not contigs:
            print("\n⚠ No contigs passed filtering!")
            try:
                output_file.touch()
                print(f"Created empty file: {output_file}")
            except Exception as e:
                print(f"Error creating empty file: {e}")
            return False
        
        print(f"\nWriting {len(contigs)} contigs to: {output_file}")
        
        total_length = sum(c['length'] for c in contigs)
        print(f"Total length: {total_length:,} bp")
        print(f"Average length: {total_length/len(contigs):,.0f} bp")
        
        try:
            with open(output_file, 'w') as f:
                for contig in contigs:
                    # Add flags to header indicating why contig was retained
                    criteria = []
                    if contig.get('passed_length', False):
                        criteria.append("LENGTH_OK")
                    if contig.get('passed_qv', False):
                        criteria.append(f"QV_OK({contig['qv']:.2f})")
                    criteria_str = ",".join(criteria) if criteria else "NONE"
                    
                    header = f">{contig['id']} total_score={contig['total_score']:.2f} QV={contig['qv']:.2f} CRAQ={contig['craq']:.2f} length={contig['length']} [{criteria_str}]"
                    f.write(f"{header}\n")
                    
                    sequence = contig['sequence']
                    for i in range(0, len(sequence), 80):
                        f.write(f"{sequence[i:i+80]}\n")
            
            print(f"✓ File written successfully ({output_file.stat().st_size} bytes)")
            return True
        except Exception as e:
            print(f"✗ Error writing to {output_file}: {e}")
            return False
    
    def generate_summary_report(
        self,
        sequences: Dict[str, Dict],
        filtered_contigs: List[Dict],
        qv_scores: Dict[str, float],
        craq_scores: Dict[str, float]
    ) -> None:
        """Generate summary report"""
        report_file = self.output_dir / "filtering_summary.txt"
        
        try:
            with open(report_file, 'w') as f:
                f.write("="*80 + "\n")
                f.write("CONTIG QUALITY FILTERING SUMMARY\n")
                f.write("="*80 + "\n\n")
                
                f.write("FILTERING MODE:\n")
                f.write("  Keep contig if EITHER length OR QV meets threshold\n\n")
                
                f.write("PARAMETERS:\n")
                f.write(f"  Minimum QV: {self.config.min_qv}\n")
                f.write(f"  Minimum length: {self.config.min_length:,} bp\n")
                if self.config.external_meryl_db:
                    f.write(f"  External Meryl DB: {self.config.external_meryl_db}\n")
                f.write("\n")
                
                f.write("INPUT STATISTICS:\n")
                f.write(f"  Total contigs: {len(sequences)}\n")
                
                if sequences:
                    lengths = [s['length'] for s in sequences.values()]
                    total_length = sum(lengths)
                    f.write(f"  Total length: {total_length:,} bp\n")
                    
                    n50 = self.calculate_n50(lengths)
                    f.write(f"  N50: {n50:,} bp\n")
                    
                    f.write(f"  Longest contig: {max(lengths):,} bp\n")
                    f.write(f"  Shortest contig: {min(lengths):,} bp\n")
                    f.write(f"  Average length: {total_length/len(sequences):,.0f} bp\n\n")
                
                f.write("QUALITY SCORE COVERAGE:\n")
                f.write(f"  Contigs with QV scores: {len(qv_scores)}\n")
                f.write(f"  Contigs with CRAQ scores: {len(craq_scores)}\n")
                
                if qv_scores:
                    avg_qv = sum(qv_scores.values()) / len(qv_scores)
                    f.write(f"  Average QV: {avg_qv:.2f}\n")
                    # Count +inf values (now 0.0)
                    inf_count = sum(1 for v in qv_scores.values() if v == 0.0)
                    if inf_count > 0:
                        f.write(f"  Note: {inf_count} contigs had '+inf' (converted to 0.0)\n")
                
                if craq_scores:
                    avg_craq = sum(craq_scores.values()) / len(craq_scores)
                    f.write(f"  Average CRAQ: {avg_craq:.2f}\n")
                f.write("\n")
                
                f.write("FILTERING RESULTS:\n")
                if filtered_contigs:
                    filtered_count = len(filtered_contigs)
                    filtered_length = sum(c['length'] for c in filtered_contigs)
                    
                    f.write(f"  Contigs retained: {filtered_count} ({filtered_count/len(sequences)*100:.1f}%)\n")
                    f.write(f"  Length retained: {filtered_length:,} bp ({filtered_length/total_length*100:.1f}%)\n\n")
                    
                    # Retention breakdown
                    only_length = sum(1 for c in filtered_contigs if c.get('passed_length', False) and not c.get('passed_qv', False))
                    only_qv = sum(1 for c in filtered_contigs if not c.get('passed_length', False) and c.get('passed_qv', False))
                    both = sum(1 for c in filtered_contigs if c.get('passed_length', False) and c.get('passed_qv', False))
                    
                    f.write("RETENTION BREAKDOWN:\n")
                    f.write(f"  Met both criteria: {both}\n")
                    f.write(f"  Met only length criterion: {only_length}\n")
                    f.write(f"  Met only QV criterion: {only_qv}\n\n")
                    
                    filtered_lengths = [c['length'] for c in filtered_contigs]
                    filtered_n50 = self.calculate_n50(filtered_lengths)
                    f.write(f"  Filtered N50: {filtered_n50:,} bp\n")
                    
                    avg_qv = sum(c['qv'] for c in filtered_contigs) / filtered_count
                    avg_craq = sum(c['craq'] for c in filtered_contigs) / filtered_count
                    avg_total = sum(c['total_score'] for c in filtered_contigs) / filtered_count
                    
                    f.write(f"  Average QV: {avg_qv:.2f}\n")
                    f.write(f"  Average CRAQ: {avg_craq:.2f}\n")
                    f.write(f"  Average total score: {avg_total:.2f}\n\n")
                    
                    f.write("TOP 10 HIGH-QUALITY CONTIGS:\n")
                    f.write("-"*80 + "\n")
                    f.write("Contig\tLength\tQV\tCRAQ\tTotal Score\tRetention Reason\n")
                    
                    for contig in filtered_contigs[:10]:
                        reason = []
                        if contig.get('passed_length', False):
                            reason.append("LENGTH")
                        if contig.get('passed_qv', False):
                            reason.append("QV")
                        reason_str = "+".join(reason) if reason else "UNKNOWN"
                        f.write(f"{contig['id']}\t{contig['length']:,}\t{contig['qv']:.2f}\t{contig['craq']:.2f}\t{contig['total_score']:.2f}\t{reason_str}\n")
                else:
                    f.write("  No contigs passed filtering!\n")
            
            print(f"✓ Summary report: {report_file}")
        except Exception as e:
            print(f"✗ Error generating summary report: {e}")
    
    def calculate_n50(self, lengths: List[int]) -> int:
        """Calculate N50"""
        if not lengths:
            return 0
        
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        half_length = total_length / 2
        
        cumulative = 0
        for length in sorted_lengths:
            cumulative += length
            if cumulative >= half_length:
                return length
        
        return sorted_lengths[-1]
    
    def create_final_output_link(self, source_file: Path, target_file: Path) -> bool:
        """Create final output link"""
        print(f"\nCreating final output: {target_file} -> {source_file}")
        
        try:
            source_file = source_file.resolve()
            if not source_file.exists():
                print(f"✗ Source file does not exist: {source_file}")
                return False
            
            if target_file.exists():
                try:
                    target_file.unlink()
                    print(f"✓ Removed existing target file: {target_file}")
                except Exception as e:
                    print(f"✗ Cannot remove existing target file: {e}")
                    return False
            
            target_file.parent.mkdir(parents=True, exist_ok=True)
            
            try:
                os.symlink(str(source_file), str(target_file))
                print(f"✓ Symlink created successfully: {target_file} -> {source_file}")
                return True
            except Exception as symlink_error:
                print(f"✗ Symlink creation failed: {symlink_error}")
                
                try:
                    shutil.copy2(str(source_file), str(target_file))
                    print(f"✓ File copy successful: {target_file}")
                    return True
                except Exception as copy_error:
                    print(f"✗ File copy failed: {copy_error}")
                    return False
                    
        except Exception as e:
            print(f"✗ Error creating final output: {e}")
            return False


def run_quality_pipeline(
    contigs_file: str,
    reads_files: Optional[List[str]] = None,
    output_prefix: str = "high_quality_contigs",
    min_qv: float = 30.0,
    min_length: int = 100000,
    threads: int = 32,
    memory_gb: int = 80,
    kmer_size: int = 21,
    force_rerun: bool = False,
    skip_meryl: bool = False,
    skip_merqury: bool = False,
    skip_craq: bool = False,
    only_analyze: bool = False,
    working_dir: Optional[str] = None,
    external_meryl_db: Optional[str] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Main function to run quality assessment pipeline
    
    Args:
        contigs_file: Path to assembled contigs file
        reads_files: List of sequencing reads files
        output_prefix: Output file prefix
        min_qv: Minimum QV value
        min_length: Minimum length (bp)
        threads: Number of threads
        memory_gb: Memory (GB)
        kmer_size: k-mer size
        force_rerun: Force rerun all steps
        skip_meryl: Skip meryl step
        skip_merqury: Skip merqury step
        skip_craq: Skip CRAQ step
        only_analyze: Use only existing results for analysis, don't run any tools
        working_dir: Working directory
        external_meryl_db: Path to external Meryl database
        **kwargs: Other parameters
        
    Returns:
        Dictionary containing results
    """
    config = QualityConfig(
        contigs_file=contigs_file,
        reads_files=reads_files or [],
        output_prefix=output_prefix,
        min_qv=min_qv,
        min_length=min_length,
        threads=threads,
        memory_gb=memory_gb,
        kmer_size=kmer_size,
        force_rerun=force_rerun,
        skip_meryl=skip_meryl,
        skip_merqury=skip_merqury,
        skip_craq=skip_craq,
        only_analyze=only_analyze,
        working_dir=working_dir,
        external_meryl_db=external_meryl_db
    )
    
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)
    
    pipeline = QualityPipeline(config)
    return pipeline.run()


def run_quality_pipeline_from_dict(config_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Run quality assessment pipeline from dictionary configuration"""
    config = QualityConfig(**config_dict)
    pipeline = QualityPipeline(config)
    return pipeline.run()


def run_quality_pipeline_from_json(json_file: str) -> Dict[str, Any]:
    """Run quality assessment pipeline from JSON file"""
    with open(json_file, 'r') as f:
        config_dict = json.load(f)
    return run_quality_pipeline_from_dict(config_dict)


def test_craq_parser():
    """Test CRAQ parser functionality - verifies extraction of scores from parentheses"""
    test_content = """#Chr	Covered.Rate	Low-confident.Rate	Avg.CRH	Avg.CSH	Avg.CRE(R-AQI)	Avg.CSE(S-AQI)
Genome	0.995782863	0.007482992	0.014962108	0.007481054	0.115956336230135(98.8471336632426)	0.0299242158013251(97.0519080772449)
Chr3_RagTag	0.999792065	0.002824466	0	0	0.0765047782206809(99.2378712595736)	0(100)
Chr1_RagTag	0.998388014	0.006531438	0.030735327	0	0.09220598201676(99.0821781159545)	0(100)
Chr5_RagTag	0.995915193	0.006757991	0	0	0.0999421967647978(99.0055556574048)	0(100)
Chr4_RagTag	0.989851745	0.01291176	0	0	0.136222639739568(98.6470099189802)	0.0454075465798559(95.5607947711069)
Chr2_RagTag	0.993043469	0.009873145	0.04466446	0.04466446	0.200990068142779(98.0101631762718)	0.133993378761853(87.4595855491604)"""
    
    # Create temporary file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.Report', delete=False) as f:
        f.write(test_content)
        temp_file = Path(f.name)
    
    try:
        # Test parsing
        pipeline = QualityPipeline()
        scores = pipeline.parse_craq_report_file(temp_file)
        
        print("\n" + "="*60)
        print("CRAQ PARSER TEST RESULTS")
        print("="*60)
        print(f"Parsed {len(scores)} contigs")
        
        # Expected values (scores inside parentheses)
        expected = {
            'Chr3_RagTag': 99.2378712595736,
            'Chr1_RagTag': 99.0821781159545,
            'Chr5_RagTag': 99.0055556574048,
            'Chr4_RagTag': 98.6470099189802,
            'Chr2_RagTag': 98.0101631762718
        }
        
        print("\nValidation (should match scores inside parentheses):")
        print("-" * 60)
        print(f"{'Contig':<15} {'Parsed':>12} {'Expected':>12} {'Diff':>12}")
        print("-" * 60)
        
        all_correct = True
        for contig, expected_score in expected.items():
            if contig in scores:
                diff = abs(scores[contig] - expected_score)
                status = "✓" if diff < 0.0001 else "✗"
                print(f"{status} {contig:<13} {scores[contig]:>12.4f} {expected_score:>12.4f} {diff:>12.6f}")
                if diff >= 0.0001:
                    all_correct = False
            else:
                print(f"✗ {contig:<13} {'NOT FOUND':>12} {expected_score:>12.4f} {'N/A':>12}")
                all_correct = False
        
        print("-" * 60)
        if all_correct:
            print("✓ All scores correctly parsed (values from parentheses)")
        else:
            print("✗ Some scores don't match expected values")
        
        # Print parsed scores rounded to 2 decimal places
        print("\nCRAQ Scores (from parentheses, rounded to 2 decimal places):")
        for contig, score in list(scores.items())[:5]:
            print(f"  {contig}: {score:.2f}")
        
    finally:
        # Clean up temporary file
        temp_file.unlink()


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Contig quality assessment and filtering - supports skipping completed steps and using external Meryl database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # First run (complete pipeline)
  python contigs_quality_smart.py -r reads.fastq.gz -c contigs.fasta -o filtered
  
  # Use external Meryl database
  python contigs_quality_smart.py --meryl-db /path/to/existing/reads.meryl -c contigs.fasta -o filtered
  
  # Force rerun all steps
  python contigs_quality_smart.py -r reads.fastq.gz -c contigs.fasta --force-rerun
  
  # Test CRAQ parser
  python contigs_quality_smart.py --test
  
  # As a module
  from contigs_quality_smart import run_quality_pipeline
  results = run_quality_pipeline(
      contigs_file="contigs.fasta",
      reads_files=["reads.fastq.gz"],
      external_meryl_db="/path/to/reads.meryl",
      output_prefix="filtered_contigs"
  )
        """
    )
    
    parser.add_argument('-r', '--reads', nargs='+', required=False,
                       help='Sequencing reads files (can be omitted if using external Meryl database)')
    parser.add_argument('-c', '--contigs', required=True,
                       help='Assembled contigs file')
    
    # External Meryl database parameter
    parser.add_argument('--meryl-db', dest='external_meryl_db',
                       help='Path to existing meryl database (skip meryl counting step)')
    
    parser.add_argument('--meryl-path', default='meryl',
                       help='Path to meryl')
    parser.add_argument('--merqury-path', default='merqury.sh',
                       help='Path to merqury.sh')
    parser.add_argument('--craq-path', default='craq',
                       help='Path to CRAQ')
    
    parser.add_argument('-k', '--kmer-size', type=int, default=21,
                       help='k-mer size')
    parser.add_argument('-t', '--threads', type=int, default=32,
                       help='Number of threads')
    parser.add_argument('-m', '--memory', type=int, default=80,
                       help='Memory (GB)')
    
    parser.add_argument('--min-qv', type=float, default=30.0,
                       help='Minimum QV value')
    parser.add_argument('--min-length', type=int, default=100000,
                       help='Minimum length (bp)')
    
    parser.add_argument('-o', '--output', default='high_quality_contigs',
                       help='Output prefix')
    
    parser.add_argument('--force-rerun', action='store_true',
                       help='Force rerun all steps')
    parser.add_argument('--only-analyze', action='store_true',
                       help='Use only existing results for analysis, don\'t run any tools')
    parser.add_argument('--skip-meryl', action='store_true',
                       help='Skip meryl step')
    parser.add_argument('--skip-merqury', action='store_true',
                       help='Skip merqury step')
    parser.add_argument('--skip-craq', action='store_true',
                       help='Skip CRAQ step')
    
    # Test option
    parser.add_argument('--test', action='store_true',
                       help='Run CRAQ parser test and exit')
    
    return parser.parse_args()


def main():
    """Command line main function"""
    args = parse_arguments()
    
    # Run test if requested
    if args.test:
        test_craq_parser()
        sys.exit(0)
    
    results = run_quality_pipeline(
        contigs_file=args.contigs,
        reads_files=args.reads,
        output_prefix=args.output,
        min_qv=args.min_qv,
        min_length=args.min_length,
        threads=args.threads,
        memory_gb=args.memory,
        kmer_size=args.kmer_size,
        force_rerun=args.force_rerun,
        skip_meryl=args.skip_meryl,
        skip_merqury=args.skip_merqury,
        skip_craq=args.skip_craq,
        only_analyze=args.only_analyze,
        meryl_path=args.meryl_path,
        merqury_path=args.merqury_path,
        craq_path=args.craq_path,
        external_meryl_db=args.external_meryl_db
    )
    
    print("\n" + "="*80)
    if results['success']:
        print("✓ PIPELINE COMPLETED SUCCESSFULLY")
        print(f"  Produced {results['filtered_contigs_count']} high-quality contigs")
        print(f"  Total length: {results['total_length']:,} bp")
        print(f"  Output: {results['output_file']}")
        if results.get('meryl_db_used'):
            print(f"  Meryl DB used: {results['meryl_db_used']}")
    else:
        print("✗ PIPELINE FAILED")
        print(f"  Error: {results.get('error', 'Unknown error')}")
    print("="*80)
    
    sys.exit(0 if results['success'] else 1)


if __name__ == "__main__":
    main()
