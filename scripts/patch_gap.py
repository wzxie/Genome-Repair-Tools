#!/usr/bin/env python3
"""
Genome Gap Patching Tool - Unified API Version based on minimap2
Provides consistent API design with extract_gap_patches.py
"""

import sys
import os
import argparse
import re
import subprocess
import tempfile
import json
import time
import logging
from typing import List, Optional, Dict, Any, Tuple
from dataclasses import dataclass, field
from collections import defaultdict
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

@dataclass
class MatchResult:
    """Match result"""
    position: int
    score: float
    match_length: int
    total_kmers: int
    matched_kmers: int
    flank_start_offset: int
    flank_end_offset: int
    genome_match_start: int = 0
    genome_match_end: int = 0
    match_type: str = "minimap2"
    alignment_length: int = 0
    matches_in_alignment: int = 0
    anchor_count: int = 0
    chain_length: int = 0
    metadata: Dict[str, Any] = field(default_factory=dict)
    cigar_string: str = ""
    
    def __str__(self):
        return (f"Position: {self.position:,} | "
                f"Score: {self.score:.3f} | "
                f"Length: {self.match_length:,} | "
                f"Type: {self.match_type} | "
                f"Flank position: {self.flank_start_offset:,}-{self.flank_end_offset:,}")

class GenomeGapPatcherAPI:
    """Unified API for genome gap patching"""
    
    def __init__(self, verbose: bool = True, log_file: Optional[str] = None):
        """
        Initialize GenomeGapPatcher API
        
        Args:
            verbose: Show detailed output
            log_file: Log file path (optional)
        """
        self.verbose = verbose
        self.log_file = log_file
        self.temp_files = []
        self._setup_logging()
    
    def _setup_logging(self):
        """Configure logging system"""
        self.logger = logging.getLogger('GenomeGapPatcher')
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
    
    def patch_gap(
        self,
        reference_fasta: str,
        patch_fasta: str,
        gap_position: int,
        output_fasta: str = "patched_genome.fasta",
        chromosome: Optional[str] = None,
        minimap2_mode: str = 'asm5',
        min_score: float = 0.7,
        min_match_length: int = 100,
        min_mapq: int = 20,
        search_range: int = 500000,
        flank_size: int = 10000,
        require_both_matches: bool = True,
        output_json: Optional[str] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Patch gap in genome
        
        Args:
            reference_fasta: Target chromosome/genome FASTA file
            patch_fasta: Patch sequence FASTA file
            gap_position: Gap start coordinate
            output_fasta: Output patched FASTA file
            chromosome: Specify chromosome name (optional, if None use first sequence)
            minimap2_mode: minimap2 alignment mode
            min_score: Minimum match score
            min_match_length: Minimum match length
            min_mapq: Minimum mapping quality
            search_range: Search range (bp)
            flank_size: Flank size (bp)
            require_both_matches: Require both left and right matches
            output_json: JSON report file path (optional)
            **kwargs: Other parameters
            
        Returns:
            Dictionary containing patching results
        """
        self.log(f"Starting genome gap patching", "info")
        self.log(f"Reference genome: {reference_fasta}", "info")
        self.log(f"Patch file: {patch_fasta}", "info")
        self.log(f"Gap position: {gap_position:,}", "info")
        
        start_time = time.time()
        
        self._check_dependencies()
        
        output_dir = os.path.dirname(output_fasta) or "."
        os.makedirs(output_dir, exist_ok=True)
        
        self.log("Reading chromosome sequence...", "info")
        chr_record, chr_seq = self._load_chromosome_sequence(reference_fasta, chromosome)
        
        self.log("Parsing patch sequence...", "info")
        selected_patch = self._load_patch_sequence(patch_fasta, gap_position)
        
        if not selected_patch:
            raise ValueError(f"No patch sequence matching gap position {gap_position} found")
        
        config = {
            'minimap2_mode': minimap2_mode,
            'min_match_score': min_score,
            'min_matching_length': min_match_length,
            'search_range': search_range,
            'flank_size': flank_size,
            'min_mapq': min_mapq,
            'require_both_matches': require_both_matches
        }
        config.update(kwargs)
        
        patcher = Minimap2GapPatcher(config, logger=self.logger)
        
        self.log(f"\n{'='*60}", "info")
        self.log(f"Starting to patch gap position: {gap_position:,}", "info")
        self.log(f"{'='*60}", "info")
        
        result = patcher.find_matching_regions(
            gap_position,
            selected_patch['sequence'],
            chr_seq
        )
        
        if not result['success']:
            error_msg = result.get('error', 'Unknown error')
            raise RuntimeError(f"Patching failed: {error_msg}")
        
        try:
            patched_seq, patch_details = patcher.apply_patch(
                chr_seq, gap_position, result['matches']
            )
            
            self._save_patched_sequence(chr_record.id, patched_seq, output_fasta)
            
            result_summary = self._create_result_summary(
                chr_seq=chr_seq,
                patched_seq=patched_seq,
                chr_record=chr_record,
                selected_patch=selected_patch,
                patch_details=patch_details,
                parameters=config,
                gap_position=gap_position,
                elapsed_time=time.time() - start_time
            )
            
            if output_json:
                self._save_results_to_json(result_summary, output_json)
                self.log(f"JSON report saved: {output_json}", "info")
            
            self._print_summary(result_summary)
            
            self._cleanup_temp_files()
            
            return result_summary
            
        except Exception as e:
            self.log(f"Failed to apply patch: {e}", "error")
            raise
    
    def _check_dependencies(self):
        """Check required dependencies"""
        try:
            result = subprocess.run(['minimap2', '--version'], 
                                  capture_output=True, text=True, check=False)
            if result.returncode == 0:
                version = result.stdout.strip()
                self.log(f"✓ minimap2 installed: {version}", "info")
            else:
                raise RuntimeError("minimap2 not properly installed")
        except FileNotFoundError:
            raise RuntimeError("minimap2 not installed. Install: conda install -c bioconda minimap2")
        
        if not BIOPYTHON_AVAILABLE:
            raise RuntimeError("BioPython not installed. Install: conda install -c bioconda biopython")
    
    def _load_chromosome_sequence(self, fasta_file: str, chromosome: Optional[str] = None) -> Tuple[Any, str]:
        """Load chromosome sequence"""
        try:
            records = list(SeqIO.parse(fasta_file, "fasta"))
            if len(records) == 0:
                raise ValueError("No sequences in FASTA file")
            
            if chromosome:
                target_record = None
                for record in records:
                    if chromosome in record.id or chromosome in record.description:
                        target_record = record
                        break
                
                if not target_record:
                    raise ValueError(f"Chromosome not found: {chromosome}")
            else:
                target_record = records[0]
                self.log(f"Using first sequence: {target_record.id}", "info")
            
            chr_seq = str(target_record.seq).upper()
            self.log(f"Chromosome ID: {target_record.id}", "info")
            self.log(f"Chromosome length: {len(chr_seq):,} bp", "info")
            
            return target_record, chr_seq
            
        except Exception as e:
            self.log(f"Failed to read chromosome sequence: {e}", "error")
            raise
    
    def _load_patch_sequence(self, patch_file: str, gap_position: int) -> Optional[Dict]:
        """Load patch sequence"""
        try:
            patch_records = list(SeqIO.parse(patch_file, "fasta"))
            
            for record in patch_records:
                header = record.description
                patch_info = self._parse_patch_header(header)
                
                if patch_info.get('gap_position') == gap_position:
                    selected_patch = {
                        'sequence': str(record.seq).upper(),
                        'header': header,
                        'info': patch_info,
                        'record': record
                    }
                    
                    self.log(f"Found matching patch: {header[:80]}...", "info")
                    self.log(f"  Patch length: {len(selected_patch['sequence']):,} bp", "info")
                    self.log(f"  Chromosome: {patch_info.get('chromosome', 'unknown')}", "info")
                    
                    return selected_patch
            
            return None
            
        except Exception as e:
            self.log(f"Failed to read patch sequence: {e}", "error")
            raise
    
    def _parse_patch_header(self, header: str) -> Dict[str, Any]:
        """Parse patch sequence header information"""
        info = {
            'gap_position': 0,
            'gap_length': 0,
            'chromosome': 'unknown',
            'source': 'unknown'
        }
        
        if 'gap' in header.lower():
            gap_match = re.search(r'gap[_-]?(\d+)', header, re.IGNORECASE)
            if gap_match:
                info['gap_position'] = int(gap_match.group(1))
            
            len_match = re.search(r'length[_-]?(\d+)', header, re.IGNORECASE)
            if len_match:
                info['gap_length'] = int(len_match.group(1))
            
            chr_match = re.search(r'chr(?:omosome)?[_-]?([\w\.]+)', header, re.IGNORECASE)
            if chr_match:
                info['chromosome'] = chr_match.group(1)
        
        return info
    
    def _save_patched_sequence(self, seq_id: str, sequence: str, output_file: str):
        """Save patched sequence"""
        try:
            with open(output_file, 'w') as f:
                f.write(f">{seq_id}_patched\n")
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")
            
            self.log(f"Patched sequence saved: {output_file}", "info")
            
        except Exception as e:
            self.log(f"Failed to save sequence: {e}", "error")
            raise
    
    def _create_result_summary(self, chr_seq: str, patched_seq: str, chr_record: Any, 
                             selected_patch: Dict, patch_details: Dict, 
                             parameters: Dict, gap_position: int, 
                             elapsed_time: float) -> Dict[str, Any]:
        """Create result summary"""
        
        result_summary = {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'parameters': {
                'reference_file': chr_record.id,
                'patch_file': selected_patch['header'],
                'gap_position': gap_position,
                'output_file': '',
                'minimap2_mode': parameters.get('minimap2_mode', 'asm5'),
                'min_score': parameters.get('min_match_score', 0.7),
                'min_match_length': parameters.get('min_matching_length', 100),
                'search_range': parameters.get('search_range', 500000),
                'flank_size': parameters.get('flank_size', 10000),
                'require_both_matches': parameters.get('require_both_matches', True)
            },
            'statistics': {
                'original_length': len(chr_seq),
                'patched_length': len(patched_seq),
                'net_change': len(patched_seq) - len(chr_seq),
                'removed_bases': patch_details.get('statistics', {}).get('removed_bases', 0),
                'inserted_bases': patch_details.get('statistics', {}).get('inserted_bases', 0),
                'processing_time_seconds': round(elapsed_time, 2)
            },
            'patch_details': patch_details,
            'patch_info': selected_patch['info'],
            'chromosome_info': {
                'id': chr_record.id,
                'description': chr_record.description,
                'original_length': len(chr_seq)
            }
        }
        
        return result_summary
    
    def _save_results_to_json(self, result_summary: Dict[str, Any], json_file: str):
        """Save result summary to JSON file"""
        try:
            with open(json_file, 'w') as f:
                json.dump(result_summary, f, indent=2, ensure_ascii=False)
        except Exception as e:
            self.log(f"Failed to save JSON report: {e}", "error")
            raise
    
    def _print_summary(self, result_summary: Dict[str, Any]):
        """Print analysis summary"""
        stats = result_summary['statistics']
        
        self.log("\n" + "="*60, "info")
        self.log("Genome gap patching completed!", "info")
        self.log("="*60, "info")
        
        self.log(f"Processing time: {stats['processing_time_seconds']:.2f} seconds", "info")
        self.log(f"Original chromosome length: {stats['original_length']:,} bp", "info")
        self.log(f"Patched length: {stats['patched_length']:,} bp", "info")
        self.log(f"Net length change: {stats['net_change']:+,} bp", "info")
        self.log(f"Removed genome sequence: {stats['removed_bases']:,} bp", "info")
        self.log(f"Inserted patch sequence: {stats['inserted_bases']:,} bp", "info")
        
        patch_details = result_summary['patch_details']
        match_quality = patch_details.get('match_quality', {})
        if match_quality:
            self.log(f"Left flank match score: {match_quality.get('left_score', 0):.3f} "
                  f"(MAPQ: {match_quality.get('left_mapq', 0)})", "info")
            self.log(f"Right flank match score: {match_quality.get('right_score', 0):.3f} "
                  f"(MAPQ: {match_quality.get('right_mapq', 0)})", "info")
        
        self.log(f"Patching details: {json_file}", "info") if result_summary['parameters']['output_file'] else None
    
    def _cleanup_temp_files(self):
        """Clean up temporary files"""
        for temp_file in self.temp_files:
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                    self.log(f"Cleaned temporary file: {temp_file}", "debug")
                except Exception as e:
                    self.log(f"Failed to clean temporary file {temp_file}: {e}", "warning")
        self.temp_files.clear()
    
    def batch_patch(
        self,
        reference_fasta: str,
        patch_fastas: List[str],
        gap_positions: List[int],
        output_dir: str = ".",
        batch_prefix: str = "batch",
        **kwargs
    ) -> Dict[str, Any]:
        """
        Batch patch multiple gaps
        
        Args:
            reference_fasta: Reference genome file
            patch_fastas: List of patch files
            gap_positions: List of gap positions
            output_dir: Output directory
            batch_prefix: Batch processing prefix
            **kwargs: Other parameters passed to patch_gap
            
        Returns:
            Batch processing summary
        """
        if len(patch_fastas) != len(gap_positions):
            raise ValueError("Patch file list and gap position list must have same length")
        
        self.log(f"Starting batch patching for {len(gap_positions)} gaps", "info")
        
        results = {}
        for i, (patch_fasta, gap_pos) in enumerate(zip(patch_fastas, gap_positions)):
            self.log(f"Processing gap {i+1}/{len(gap_positions)}: position {gap_pos:,}", "info")
            
            output_file = os.path.join(output_dir, f"{batch_prefix}_{i+1:03d}_patched.fasta")
            output_json = os.path.join(output_dir, f"{batch_prefix}_{i+1:03d}_report.json")
            
            try:
                result = self.patch_gap(
                    reference_fasta=reference_fasta,
                    patch_fasta=patch_fasta,
                    gap_position=gap_pos,
                    output_fasta=output_file,
                    output_json=output_json,
                    **kwargs
                )
                results[f"gap_{gap_pos}"] = {
                    "success": True,
                    "result": result,
                    "output_files": {
                        "patched_fasta": output_file,
                        "json_report": output_json
                    }
                }
            except Exception as e:
                results[f"gap_{gap_pos}"] = {
                    "success": False,
                    "error": str(e),
                    "output_files": {
                        "patched_fasta": output_file,
                        "json_report": output_json
                    }
                }
                self.log(f"Patching failed: {e}", "warning")
        
        summary = {
            "total": len(gap_positions),
            "successful": sum(1 for r in results.values() if r["success"]),
            "failed": sum(1 for r in results.values() if not r["success"]),
            "results": results,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        summary_file = os.path.join(output_dir, f"{batch_prefix}_summary.json")
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        self.log(f"Batch patching completed, successful: {summary['successful']}/{summary['total']}", "info")
        self.log(f"Summary report saved: {summary_file}", "info")
        
        return summary
    
    def cleanup(self):
        """Clean up all temporary files"""
        self._cleanup_temp_files()
    
    def __del__(self):
        """Destructor, automatically clean temporary files"""
        self.cleanup()

class Minimap2Matcher:
    """minimap2 aligner wrapper"""
    
    def __init__(self, mode: str = 'asm5', min_score: float = 0.7,
                 min_match_length: int = 100, min_mapq: int = 20,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize minimap2 aligner
        """
        self.mode = mode
        self.min_score = min_score
        self.min_match_length = min_match_length
        self.min_mapq = min_mapq
        self.logger = logger or logging.getLogger(__name__)
        
        self._check_minimap2()
        
        self.stats = {
            'minimap2_calls': 0,
            'successful_matches': 0,
            'failed_matches': 0,
            'total_time': 0
        }
    
    def _check_minimap2(self):
        """Check if minimap2 is installed"""
        try:
            result = subprocess.run(['minimap2', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                version = result.stdout.strip()
                self.logger.info(f"minimap2 installed: {version}")
            else:
                self.logger.error("minimap2 not properly installed")
                sys.exit(1)
        except FileNotFoundError:
            self.logger.error("minimap2 not installed")
            sys.exit(1)
    
    def _run_minimap2(self, query_seq: str, target_seq: str, 
                     extra_args: List[str] = None) -> str:
        """Run minimap2 and return result"""
        if extra_args is None:
            extra_args = []
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_file:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as target_file:
                query_file.write(f">query\n{query_seq}\n")
                target_file.write(f">target\n{target_seq}\n")
                query_file.flush()
                target_file.flush()
                
                cmd = ['minimap2', '-x', self.mode, '-c'] + extra_args + [
                    target_file.name, query_file.name
                ]
                
                start_time = time.time()
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, 
                                          check=True, timeout=300)
                    self.stats['minimap2_calls'] += 1
                    self.stats['total_time'] += time.time() - start_time
                    
                    return result.stdout
                except subprocess.CalledProcessError as e:
                    self.logger.error(f"minimap2 run failed: {e}")
                    return ""
                except subprocess.TimeoutExpired:
                    self.logger.error("minimap2 run timeout")
                    return ""
                finally:
                    try:
                        os.unlink(query_file.name)
                        os.unlink(target_file.name)
                    except:
                        pass
    
    def _parse_paf_line(self, line: str, target_start: int = 0) -> Optional[MatchResult]:
        """Parse single PAF format result line"""
        fields = line.strip().split('\t')
        
        if len(fields) < 12:
            return None
        
        query_start = int(fields[2])
        query_end = int(fields[3])
        strand = fields[4]
        target_match_start = int(fields[7])
        target_match_end = int(fields[8])
        matches = int(fields[9])
        alignment_length = int(fields[10])
        mapq = int(fields[11])
        
        if mapq < self.min_mapq:
            return None
        
        cigar_string = ""
        anchor_count = 0
        for i in range(12, len(fields)):
            if fields[i].startswith('cg:Z:'):
                cigar_string = fields[i][5:]
                anchor_count = self._estimate_anchors_from_cigar(cigar_string)
                break
        
        if alignment_length > 0:
            match_score = matches / alignment_length
        else:
            match_score = 0.0
        
        match_length = query_end - query_start
        genome_match_start = target_start + target_match_start
        genome_match_end = target_start + target_match_end
        
        result = MatchResult(
            position=genome_match_start,
            score=match_score,
            match_length=match_length,
            total_kmers=max(1, match_length // 15),
            matched_kmers=matches,
            flank_start_offset=query_start,
            flank_end_offset=query_end,
            genome_match_start=genome_match_start,
            genome_match_end=genome_match_end,
            match_type=f"minimap2_{self.mode}",
            alignment_length=alignment_length,
            matches_in_alignment=matches,
            anchor_count=anchor_count,
            chain_length=anchor_count,
            metadata={
                'mapq': mapq,
                'strand': strand,
                'cigar': cigar_string,
                'query_coverage': (query_end - query_start) / int(fields[1]) if int(fields[1]) > 0 else 0.0
            },
            cigar_string=cigar_string
        )
        
        return result
    
    def _estimate_anchors_from_cigar(self, cigar: str) -> int:
        """Estimate anchor count from CIGAR string"""
        if not cigar:
            return 0
        
        pattern = re.compile(r'(\d+)([MIDNSHPX=])')
        matches = list(pattern.finditer(cigar))
        
        anchors = 0
        for match in matches:
            length = int(match.group(1))
            op = match.group(2)
            if op == 'M' and length >= 15:
                anchors += 1
        
        return max(1, anchors)
    
    def find_best_match(self, query_seq: str, target_seq: str,
                       target_start: int = 0,
                       expected_position: Optional[int] = None,
                       use_secondary: bool = False) -> Optional[MatchResult]:
        """
        Find best match using minimap2
        """
        query_len = len(query_seq)
        target_len = len(target_seq)
        
        self.logger.info(f"Aligning {query_len:,}bp vs {target_len:,}bp (mode: {self.mode})")
        
        if query_len < 50 or target_len < 50:
            self.logger.warning("Sequence too short, skipping")
            return None
        
        extra_args = []
        
        if self.mode.startswith('asm'):
            extra_args.extend(['-r', '100', '-N', '5'])
        else:
            extra_args.extend(['-N', '3'])
        
        if expected_position is not None and target_len > 0:
            expected_rel_pos = expected_position - target_start
            if 0 <= expected_rel_pos < target_len:
                search_radius = min(200000, target_len // 2)
                search_start = max(0, expected_rel_pos - search_radius)
                search_end = min(target_len, expected_rel_pos + search_radius)
                
                if search_end - search_start > 1000:
                    self.logger.info(f"Searching near expected position: ±{search_radius:,}bp")
                    target_seq = target_seq[search_start:search_end]
                    target_start = target_start + search_start
        
        output = self._run_minimap2(query_seq, target_seq, extra_args)
        
        if not output:
            self.logger.error("Alignment failed")
            self.stats['failed_matches'] += 1
            return None
        
        lines = output.strip().split('\n')
        results = []
        
        for line in lines:
            if not line or line.startswith('#'):
                continue
            
            result = self._parse_paf_line(line, target_start)
            if result:
                results.append(result)
        
        if not results:
            self.logger.warning("No matches found")
            self.stats['failed_matches'] += 1
            return None
        
        def score_key(r: MatchResult) -> float:
            length_factor = min(1.0, r.match_length / 1000.0)
            return r.score * (0.7 + 0.3 * length_factor)
        
        results.sort(key=score_key, reverse=True)
        best_result = results[0]
        
        if (best_result.score >= self.min_score and 
            best_result.match_length >= self.min_match_length):
            
            self.logger.info(f"Match found: score={best_result.score:.3f}, "
                           f"length={best_result.match_length:,}bp, "
                           f"MAPQ={best_result.metadata['mapq']}")
            
            if best_result.metadata['strand'] == '-':
                self.logger.warning("Match on reverse strand")
            
            self.stats['successful_matches'] += 1
            return best_result
        else:
            self.logger.warning(f"Match quality insufficient: score={best_result.score:.3f} "
                              f"(requires>={self.min_score}), length={best_result.match_length:,}bp "
                              f"(requires>={self.min_match_length})")
            self.stats['failed_matches'] += 1
            return None

class Minimap2GapPatcher:
    """minimap2-based gap patcher"""
    
    def __init__(self, config: Optional[Dict] = None, logger: Optional[logging.Logger] = None):
        default_config = {
            'minimap2_mode': 'asm5',
            'min_match_score': 0.7,
            'min_matching_length': 100,
            'search_range': 500000,
            'flank_size': 10000,
            'min_mapq': 20,
            'require_both_matches': True
        }
        
        if config:
            default_config.update(config)
        
        self.config = default_config
        self.logger = logger or logging.getLogger(__name__)
        
        self.matcher = Minimap2Matcher(
            mode=self.config['minimap2_mode'],
            min_score=self.config['min_match_score'],
            min_match_length=self.config['min_matching_length'],
            min_mapq=self.config['min_mapq'],
            logger=self.logger
        )
    
    def extract_flanks(self, patch_seq: str) -> Tuple[str, str, str, Dict]:
        """Extract flank and core sequences"""
        patch_len = len(patch_seq)
        flank_size = self.config['flank_size']
        
        if patch_len <= 2 * flank_size:
            flank_size = patch_len // 4
            if flank_size < 100:
                flank_size = min(100, patch_len // 2)
                if flank_size < 20:
                    raise ValueError(f"Patch sequence too short: {patch_len}")
        
        left_flank = patch_seq[:flank_size]
        core_patch = patch_seq[flank_size:-flank_size] if patch_len > 2 * flank_size else ""
        right_flank = patch_seq[-flank_size:]
        
        flank_info = {
            'left_size': flank_size,
            'right_size': flank_size,
            'core_size': len(core_patch),
            'full_patch_length': patch_len
        }
        
        return left_flank, core_patch, right_flank, flank_info
    
    def find_matching_regions(self, gap_start: int, patch_seq: str, 
                            genome_seq: str) -> Dict:
        """Find matching regions for gap"""
        
        result = {
            'gap_start': gap_start,
            'success': False,
            'error': None,
            'matches': {}
        }
        
        try:
            left_flank, core_patch, right_flank, flank_info = self.extract_flanks(patch_seq)
            
            self.logger.info(f"Processing gap start position: {gap_start:,}")
            self.logger.info(f"Patch info: total={len(patch_seq):,}bp, flank={flank_info['left_size']:,}bp, "
                          f"core={flank_info['core_size']:,}bp")
            
            search_range = self.config['search_range']
            genome_len = len(genome_seq)
            
            left_context_start = max(0, gap_start - search_range)
            left_context_end = min(genome_len, gap_start + search_range)
            left_context = genome_seq[left_context_start:left_context_end]
            
            right_context_start = max(0, gap_start - search_range)
            right_context_end = min(genome_len, gap_start + search_range)
            right_context = genome_seq[right_context_start:right_context_end]
            
            self.logger.info(f"Search range: left[{left_context_start:,}-{left_context_end:,}], "
                          f"right[{right_context_start:,}-{right_context_end:,}]")
            
            self.logger.info("Finding left match...")
            left_match = self.matcher.find_best_match(
                left_flank, left_context, left_context_start,
                expected_position=gap_start
            )
            
            self.logger.info("Finding right match...")
            right_match = self.matcher.find_best_match(
                right_flank, right_context, right_context_start,
                expected_position=gap_start
            )
            
            if left_match and right_match:
                left_min_score = self.config['min_match_score']
                right_min_score = self.config['min_match_score']
                
                if left_match.metadata.get('mapq', 0) >= 60:
                    left_min_score = max(0.5, left_min_score * 0.8)
                if right_match.metadata.get('mapq', 0) >= 60:
                    right_min_score = max(0.5, right_min_score * 0.8)
                
                left_valid = (left_match.score >= left_min_score and 
                             left_match.match_length >= self.config['min_matching_length'])
                right_valid = (right_match.score >= right_min_score and 
                              right_match.match_length >= self.config['min_matching_length'])
                
                if left_valid and right_valid:
                    self.logger.info("Match validation passed!")
                    
                    result['success'] = True
                    result['matches'] = {
                        'left': {
                            'match': left_match,
                            'flank': left_flank,
                            'context_start': left_context_start,
                            'flank_size': flank_info['left_size'],
                            'full_patch': patch_seq
                        },
                        'right': {
                            'match': right_match,
                            'flank': right_flank,
                            'context_start': right_context_start,
                            'flank_size': flank_info['right_size'],
                            'full_patch': patch_seq
                        },
                        'core': core_patch,
                        'flank_info': flank_info
                    }
                    
                    self.logger.info(f"Left flank match: {left_match.match_length:,}bp, score={left_match.score:.3f}, "
                                  f"MAPQ={left_match.metadata.get('mapq', 'N/A')}")
                    self.logger.info(f"Right flank match: {right_match.match_length:,}bp, score={right_match.score:.3f}, "
                                  f"MAPQ={right_match.metadata.get('mapq', 'N/A')}")
                    
                    left_distance = abs(left_match.position - gap_start)
                    right_distance = abs(right_match.position - gap_start)
                    self.logger.info(f"Match distance from gap: left flank{left_distance:,}bp, right flank{right_distance:,}bp")
                    
                    if left_match.metadata.get('strand') == '-' or right_match.metadata.get('strand') == '-':
                        self.logger.warning("At least one match on reverse strand, may need adjustment")
                    
                else:
                    self.logger.error("Match quality insufficient")
                    result['error'] = "Match quality insufficient"
            else:
                missing = []
                if not left_match:
                    missing.append("left")
                if not right_match:
                    missing.append("right")
                self.logger.error(f"Incomplete match: missing {', '.join(missing)} match")
                result['error'] = f"Missing {', '.join(missing)} match"
                
        except Exception as e:
            result['error'] = str(e)
            self.logger.error(f"Processing failed: {e}")
            import traceback
            traceback.print_exc()
        
        return result
    
    def apply_patch(self, genome_seq: str, gap_start: int,
                  match_results: Dict) -> Tuple[str, Dict]:
        """Apply patch to genome"""
        
        left_match = match_results['left']['match']
        right_match = match_results['right']['match']
        left_context_start = match_results['left']['context_start']
        left_flank_size = match_results['left']['flank_size']
        right_flank_size = match_results['right']['flank_size']
        full_patch = match_results['left']['full_patch']
        flank_info = match_results['flank_info']
        
        self.logger.info("Calculating patch positions...")
        
        left_genome_start = left_match.genome_match_start
        left_genome_end = left_match.genome_match_end
        right_genome_start = right_match.genome_match_start
        right_genome_end = right_match.genome_match_end
        
        self.logger.info(f"Left genome match: {left_genome_start:,}-{left_genome_end:,} ({left_match.match_length:,}bp)")
        self.logger.info(f"Right genome match: {right_genome_start:,}-{right_genome_end:,} ({right_match.match_length:,}bp)")
        
        left_is_reverse = left_match.metadata.get('strand') == '-'
        right_is_reverse = right_match.metadata.get('strand') == '-'
        
        if left_is_reverse or right_is_reverse:
            self.logger.warning("Reverse strand match detected, needs special handling")
            if left_is_reverse and right_is_reverse:
                self.logger.info("Both matches on reverse strand, using reverse complement of patch")
                complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
                full_patch = ''.join(complement.get(base, base) for base in reversed(full_patch))
        
        left_cut = left_genome_end
        right_cut = right_genome_start
        
        if left_cut > right_cut:
            self.logger.warning(f"Cut points cross ({left_cut:,} > {right_cut:,}), adjusting")
            
            left_cigar = left_match.cigar_string
            right_cigar = right_match.cigar_string
            
            if left_cigar and right_cigar:
                self.logger.info("Using CIGAR information to adjust cut points")
                midpoint = (left_cut + right_cut) // 2
                left_cut = min(midpoint, gap_start)
                right_cut = max(midpoint, gap_start + 100)
            else:
                left_cut = min(left_cut, gap_start)
                right_cut = max(right_cut, gap_start + 100)
        
        genome_len = len(genome_seq)
        left_cut = max(0, min(left_cut, gap_start + self.config['search_range']))
        right_cut = max(gap_start - self.config['search_range'], min(right_cut, genome_len))
        
        left_patch_start = left_match.flank_start_offset
        left_patch_end = left_match.flank_end_offset
        
        right_patch_offset = len(full_patch) - right_flank_size
        right_patch_start = right_patch_offset + right_match.flank_start_offset
        right_patch_end = right_patch_offset + right_match.flank_end_offset
        
        self.logger.info(f"Left patch match: {left_patch_start:,}-{left_patch_end:,}")
        self.logger.info(f"Right patch match: {right_patch_start:,}-{right_patch_end:,}")
        
        patch_to_insert_start = left_patch_end
        patch_to_insert_end = right_patch_start
        
        if patch_to_insert_start > patch_to_insert_end:
            self.logger.info("Patch match regions overlap, merging regions")
            patch_to_insert_start = min(left_patch_start, right_patch_start)
            patch_to_insert_end = max(left_patch_end, right_patch_end)
        
        patch_to_insert_start = max(0, min(patch_to_insert_start, len(full_patch)))
        patch_to_insert_end = max(0, min(patch_to_insert_end, len(full_patch)))
        
        if patch_to_insert_start >= patch_to_insert_end:
            self.logger.warning("Patch insertion region invalid, using core region")
            patch_to_insert_start = left_flank_size
            patch_to_insert_end = len(full_patch) - right_flank_size
        
        patch_to_insert = full_patch[patch_to_insert_start:patch_to_insert_end]
        
        removed_length = max(0, right_cut - left_cut)
        inserted_length = len(patch_to_insert)
        
        self.logger.info("Patch operation:")
        self.logger.info(f"Delete genome region: {left_cut:,}-{right_cut:,}")
        self.logger.info(f"Delete length: {removed_length:,}bp")
        self.logger.info(f"Insert patch region: {patch_to_insert_start:,}-{patch_to_insert_end:,}")
        self.logger.info(f"Insert length: {inserted_length:,}bp")
        self.logger.info(f"Match direction: left flank {'reverse' if left_is_reverse else 'forward'}, "
                      f"right flank {'reverse' if right_is_reverse else 'forward'}")
        
        left_part = genome_seq[:left_cut]
        right_part = genome_seq[right_cut:]
        patched_seq = left_part + patch_to_insert + right_part
        
        net_change = inserted_length - removed_length
        
        self.logger.info(f"Net change: {net_change:+,}bp")
        
        patch_info = {
            'gap_position': gap_start,
            'deleted_genome_region': f"{left_cut:,}-{right_cut:,}",
            'inserted_patch_region': f"{patch_to_insert_start:,}-{patch_to_insert_end:,}",
            'patch_insert_length': inserted_length,
            'merge_positions': {
                'left_cut': left_cut,
                'right_cut': right_cut,
                'removed_length': removed_length,
                'patch_to_insert_length': inserted_length,
                'left_match_type': left_match.match_type,
                'right_match_type': right_match.match_type,
                'left_match_score': left_match.score,
                'right_match_score': right_match.score,
                'left_mapq': left_match.metadata.get('mapq', 0),
                'right_mapq': right_match.metadata.get('mapq', 0),
                'left_strand': left_match.metadata.get('strand', '+'),
                'right_strand': right_match.metadata.get('strand', '+')
            },
            'match_quality': {
                'left_score': left_match.score,
                'right_score': right_match.score,
                'combined_score': (left_match.score + right_match.score) / 2,
                'left_mapq': left_match.metadata.get('mapq', 0),
                'right_mapq': right_match.metadata.get('mapq', 0)
            },
            'statistics': {
                'removed_bases': removed_length,
                'inserted_bases': inserted_length,
                'net_change': net_change
            },
            'minimap2_stats': self.matcher.stats
        }
        
        return patched_seq, patch_info

def patch_genome_gap(
    reference_fasta: str,
    patch_fasta: str,
    gap_position: int,
    output_fasta: str = "patched_genome.fasta",
    verbose: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Simplified interface for patching genome gap
    
    Example:
        >>> from patch_gap import patch_genome_gap
        >>> result = patch_genome_gap(
        ...     reference_fasta="chromosome.fasta",
        ...     patch_fasta="patches.fasta",
        ...     gap_position=11533087,
        ...     output_fasta="patched.fasta",
        ...     minimap2_mode="asm10",
        ...     flank_size=5000,
        ...     verbose=False
        ... )
    """
    patcher = GenomeGapPatcherAPI(verbose=verbose)
    return patcher.patch_gap(
        reference_fasta=reference_fasta,
        patch_fasta=patch_fasta,
        gap_position=gap_position,
        output_fasta=output_fasta,
        **kwargs
    )

def main():
    """Command line interface main function"""
    parser = argparse.ArgumentParser(
        description='minimap2-based genome gap patching tool (Unified API Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Basic usage
  python patch_gap.py -r chromosome.fasta -p patches.fasta --gap-position 11533087
  
  # Specify chromosome and multiple parameters
  python patch_gap.py -r genome.fasta -p patches.fasta --gap-position 11533087 --chromosome Chr1 --minimap2-mode asm10 --flank-size 5000
  
  # Python module call
  from patch_gap import GenomeGapPatcherAPI
  patcher = GenomeGapPatcherAPI(verbose=True)
  result = patcher.patch_gap(
      reference_fasta="chromosome.fasta",
      patch_fasta="patches.fasta",
      gap_position=11533087,
      output_fasta="patched.fasta",
      minimap2_mode="asm5",
      flank_size=10000
  )
        """
    )
    
    parser.add_argument("-r", "--reference", dest="reference_fasta", required=True,
                       help="Target chromosome/genome FASTA file")
    parser.add_argument("-p", "--patch", dest="patch_fasta", required=True,
                       help="Patch sequence FASTA file")
    parser.add_argument("--gap-position", type=int, required=True,
                       help="Gap position coordinate")
    
    parser.add_argument("-o", "--output", default="patched_genome.fasta",
                       help="Output FASTA file, default: patched_genome.fasta")
    parser.add_argument("--chromosome", help="Specify chromosome name (if reference has multiple sequences)")
    parser.add_argument("--json", default="patch_report.json",
                       help="Result JSON file, default: patch_report.json")
    parser.add_argument("--output-dir", default=".", help="Output directory, default current directory")
    
    parser.add_argument("--minimap2-mode", default='asm5', 
                       choices=['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb', 'splice'],
                       help='minimap2 alignment mode (default: asm5)')
    parser.add_argument("--min-score", type=float, default=0.7, 
                       help='Minimum match score (default: 0.7)')
    parser.add_argument("--min-mapq", type=int, default=20, 
                       help='Minimum mapping quality (default: 20)')
    parser.add_argument("--min-match-length", type=int, default=100, 
                       help='Minimum match length (default: 100)')
    
    parser.add_argument("--search-range", type=int, default=500000, 
                       help='Search range (default: 500,000 bp)')
    parser.add_argument("--flank-size", type=int, default=10000, 
                       help='Flank size (default: 10,000 bp)')
    parser.add_argument("--require-both-matches", action="store_true", default=True,
                       help='Require both left and right matches (default)')
    parser.add_argument("--allow-single-match", action="store_false", dest="require_both_matches",
                       help='Allow single match')
    
    parser.add_argument("--verbose", action="store_true", help="Show detailed output")
    parser.add_argument("--quiet", action="store_true", help="Quiet mode")
    parser.add_argument("--log-file", help="Log file path")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.reference_fasta):
        print(f"Error: Reference genome file {args.reference_fasta} not found", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.patch_fasta):
        print(f"Error: Patch file {args.patch_fasta} not found", file=sys.stderr)
        sys.exit(1)
    
    verbose = args.verbose and not args.quiet
    
    try:
        output_file = os.path.join(args.output_dir, args.output)
        output_json = os.path.join(args.output_dir, args.json)
        
        patcher = GenomeGapPatcherAPI(verbose=verbose, log_file=args.log_file)
        
        result = patcher.patch_gap(
            reference_fasta=args.reference_fasta,
            patch_fasta=args.patch_fasta,
            gap_position=args.gap_position,
            output_fasta=output_file,
            chromosome=args.chromosome,
            minimap2_mode=args.minimap2_mode,
            min_score=args.min_score,
            min_match_length=args.min_match_length,
            min_mapq=args.min_mapq,
            search_range=args.search_range,
            flank_size=args.flank_size,
            require_both_matches=args.require_both_matches,
            output_json=output_json
        )
        
        print(f"\nPatching completed!")
        print(f"Patched sequence saved to: {output_file}")
        print(f"Detailed report saved to: {output_json}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()