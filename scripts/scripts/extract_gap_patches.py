#!/usr/bin/env python3
"""
Gap Patch Sequence Extractor - Unified API Version
Extract gap patch sequences from coords files, can be called as submodule
"""

import sys
import os
import subprocess
import argparse
import json
import time
import logging
from typing import List, Dict, Tuple, Optional, Any, Union
from collections import defaultdict
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

class GapPatchesExtractorAPI:
    """Unified API for gap patch extraction"""
    
    def __init__(self, verbose: bool = True, log_file: Optional[str] = None):
        """
        Initialize GapPatchesExtractor API
        
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
        self.logger = logging.getLogger('GapPatchesExtractor')
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
    
    def extract_patches(
        self,
        reference_fasta: str,
        coords_file: str,
        gap_positions: List[int],
        output_file: str = "gap_patches.fasta",
        flank_size: int = 10000,
        require_both_anchors: bool = True,
        min_anchor_distance: int = 1000,
        max_anchor_distance: int = 1000000,
        output_json: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Extract gap patch sequences from reference genome
        
        Args:
            reference_fasta: Reference genome file path
            coords_file: Nucmer alignment result coords file
            gap_positions: List of gap positions (based on query coordinates)
            output_file: Output FASTA file path
            flank_size: Flanking sequence length (bp)
            require_both_anchors: Require both left and right anchors
            min_anchor_distance: Minimum anchor distance (bp)
            max_anchor_distance: Maximum anchor distance (bp)
            output_json: JSON report file path (optional)
            
        Returns:
            Dictionary containing extraction results
        """
        self.log(f"Starting gap patch sequence extraction", "info")
        self.log(f"Reference genome: {reference_fasta}", "info")
        self.log(f"Alignment file: {coords_file}", "info")
        self.log(f"Processing {len(gap_positions)} gap positions", "info")
        
        start_time = time.time()
        
        self._check_dependencies()
        
        self._ensure_reference_index(reference_fasta)
        
        output_dir = os.path.dirname(output_file) or "."
        os.makedirs(output_dir, exist_ok=True)
        
        self.log("Parsing alignment file...", "info")
        alignments = self._parse_coords_file(coords_file)
        
        if not alignments:
            raise ValueError("No alignment information parsed")
        
        qry_contigs = set(align[5] for align in alignments)
        if len(qry_contigs) > 1:
            self.log(f"Warning: Multiple query contigs found: {qry_contigs}", "warning")
        
        qry_contig = list(qry_contigs)[0]
        self.log(f"Query contig: {qry_contig}", "info")
        
        ref_contig_groups = self._group_alignments_by_ref_contig(alignments)
        self.log(f"Alignment involves {len(ref_contig_groups)} reference contigs", "info")
        
        result_summary = {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'parameters': {
                'reference_file': reference_fasta,
                'coords_file': coords_file,
                'flank_size': flank_size,
                'output_file': output_file,
                'require_both_anchors': require_both_anchors,
                'min_anchor_distance': min_anchor_distance,
                'max_anchor_distance': max_anchor_distance,
                'gap_positions': gap_positions
            },
            'statistics': {
                'total_gaps': len(gap_positions),
                'successful_gaps': 0,
                'total_patches': 0,
                'patches_per_gap': defaultdict(int)
            },
            'gaps': []
        }
        
        if os.path.exists(output_file):
            os.remove(output_file)
        
        patch_counter = 1
        
        for gap_pos in gap_positions:
            self.log(f"Processing gap position: {gap_pos:,}", "info")
            
            gap_result = self._process_single_gap(
                reference_fasta=reference_fasta,
                alignments=alignments,
                ref_contig_groups=ref_contig_groups,
                qry_contig=qry_contig,
                gap_pos=gap_pos,
                flank_size=flank_size,
                require_both_anchors=require_both_anchors,
                min_anchor_distance=min_anchor_distance,
                max_anchor_distance=max_anchor_distance,
                output_file=output_file,
                start_patch_counter=patch_counter
            )
            
            gap_result['gap_position'] = gap_pos
            result_summary['gaps'].append(gap_result)
            
            successful_patches = gap_result['successful_patches']
            if successful_patches > 0:
                result_summary['statistics']['successful_gaps'] += 1
                result_summary['statistics']['total_patches'] += successful_patches
                result_summary['statistics']['patches_per_gap'][successful_patches] += 1
            
            patch_counter += gap_result['total_patches_attempted']
        
        elapsed_time = time.time() - start_time
        result_summary['statistics']['processing_time_seconds'] = round(elapsed_time, 2)
        
        result_summary['api_info'] = {
            "execution_time": elapsed_time,
            "input_files": {
                "reference_fasta": reference_fasta,
                "coords_file": coords_file,
                "gap_positions": gap_positions
            },
            "parameters": result_summary['parameters'],
            "output_files": {
                "patches_fasta": output_file,
                "json_report": output_json
            },
            "timestamp": result_summary['timestamp']
        }
        
        if output_json:
            self._save_results_to_json(result_summary, output_json)
            self.log(f"JSON report saved: {output_json}", "info")
        
        self._print_summary(result_summary)
        
        self._cleanup_temp_files()
        
        return result_summary
    
    def _check_dependencies(self):
        """Check required dependencies"""
        try:
            subprocess.run(["samtools", "--version"], 
                         stdout=subprocess.DEVNULL, 
                         stderr=subprocess.DEVNULL,
                         check=False)
        except FileNotFoundError:
            raise RuntimeError("samtools not installed or not in PATH. Install samtools: conda install -c bioconda samtools")
        
        if not BIOPYTHON_AVAILABLE:
            raise RuntimeError("BioPython not installed. Install: conda install -c bioconda biopython")
    
    def _ensure_reference_index(self, reference_fasta: str):
        """Ensure reference genome has index"""
        if not os.path.exists(reference_fasta):
            raise FileNotFoundError(f"Reference genome file not found: {reference_fasta}")
        
        index_file = f"{reference_fasta}.fai"
        if not os.path.exists(index_file):
            self.log(f"Creating index for reference genome: {reference_fasta}", "info")
            try:
                subprocess.run(["samtools", "faidx", reference_fasta], 
                             check=True, 
                             stderr=subprocess.PIPE)
                self.log(f"Index created successfully: {index_file}", "info")
            except subprocess.CalledProcessError as e:
                error_msg = e.stderr.decode() if e.stderr else str(e)
                raise RuntimeError(f"Failed to create index for reference genome: {error_msg}")
    
    def _parse_coords_file(self, coords_file: str) -> List[Tuple]:
        """Parse coords file and extract alignment information"""
        alignments = []
        
        try:
            with open(coords_file, 'r') as f:
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
                            qry_start = int(parts[3])
                            qry_end = int(parts[4])
                            ref_contig = parts[-2]
                            qry_contig = parts[-1]
                            
                            alignments.append((ref_start, ref_end, qry_start, qry_end, ref_contig, qry_contig))
                        except (ValueError, IndexError):
                            continue
        
        except Exception as e:
            self.log(f"Failed to parse coords file: {e}", "error")
            raise
        
        return alignments
    
    def _group_alignments_by_ref_contig(self, alignments: List[Tuple]) -> Dict:
        """Group syntenic regions by reference contig"""
        ref_contig_groups = {}
        for align in alignments:
            *_, ref_contig, qry_contig = align
            if ref_contig not in ref_contig_groups:
                ref_contig_groups[ref_contig] = []
            ref_contig_groups[ref_contig].append(align)
        
        return ref_contig_groups
    
    def _find_surrounding_alignments(self, alignments: List[Tuple], gap_pos: int, 
                                   require_both: bool = True) -> Tuple[Optional[Tuple], Optional[Tuple]]:
        """Find closest left and right syntenic blocks around gap"""
        left_candidates = []
        right_candidates = []
        
        for align in alignments:
            _, _, qry_start, qry_end, _, _ = align
            q_min = min(qry_start, qry_end)
            q_max = max(qry_start, qry_end)
            
            if q_max <= gap_pos:
                left_candidates.append((q_max, align))
            elif q_min > gap_pos:
                right_candidates.append((q_min, align))
        
        before_align = max(left_candidates, key=lambda x: x[0])[1] if left_candidates else None
        after_align = min(right_candidates, key=lambda x: x[0])[1] if right_candidates else None
        
        if require_both and (before_align is None or after_align is None):
            return None, None
        
        return before_align, after_align
    
    def _get_patch_region_for_ref_contig(self, before_align: Tuple, after_align: Tuple, 
                                       flank_size: int) -> Tuple[Optional[int], Optional[int], int, Optional[str]]:
        """Build patch region with flank for single reference contig"""
        r1s, r1e, q1s, q1e, ref_ctg1, _ = before_align
        r2s, r2e, q2s, q2e, ref_ctg2, _ = after_align
        
        if ref_ctg1 != ref_ctg2:
            return None, None, 0, None
        
        left_ref_min, left_ref_max = min(r1s, r1e), max(r1s, r1e)
        right_ref_min, right_ref_max = min(r2s, r2e), max(r2s, r2e)
        
        if left_ref_max < right_ref_min:
            direction = 1
            gap_adj_left = left_ref_max
            gap_adj_right = right_ref_min
        elif right_ref_max < left_ref_min:
            direction = -1
            gap_adj_left = right_ref_max
            gap_adj_right = left_ref_min
        else:
            return None, None, 0, None
        
        patch_start = max(1, gap_adj_left - flank_size + 1)
        patch_end = gap_adj_right + flank_size - 1
        
        return patch_start, patch_end, direction, ref_ctg1
    
    def _extract_patch_sequence(self, ref_fasta: str, contig: str, start: int, end: int, 
                              direction: int, flank_size: int, patch_index: int, 
                              gap_pos: int, output_file: str) -> Tuple[bool, Optional[int]]:
        """Extract patch sequence and write to specified output file"""
        temp_file = f"temp_patch_{patch_index}_{int(time.time())}.fa"
        self.temp_files.append(temp_file)
        
        try:
            extract_cmd = ["samtools", "faidx", ref_fasta, f"{contig}:{start}-{end}"]
            
            with open(temp_file, 'w') as outfile:
                result = subprocess.run(
                    extract_cmd,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    check=True
                )
            
            if os.path.getsize(temp_file) == 0:
                self.log(f"  Patch{patch_index}: Extracted sequence is empty", "warning")
                return False, None
            
            with open(temp_file, 'r') as f:
                lines = f.read().strip().split('\n')
            
            if not lines or not lines[0].startswith('>'):
                self.log(f"  Patch{patch_index}: Invalid FASTA format", "warning")
                return False, None
            
            sequence = ''.join(lines[1:]).replace('\n', '').replace(' ', '')
            
            if direction == -1:
                seq_obj = Seq(sequence)
                rev_comp_seq = str(seq_obj.reverse_complement())
                new_header = f">gap{gap_pos}_patch{patch_index}_contig{contig}_{start}_{end}_revcomp_flank{flank_size}"
                formatted_seq = '\n'.join([rev_comp_seq[i:i+80] for i in range(0, len(rev_comp_seq), 80)])
            else:
                new_header = f">gap{gap_pos}_patch{patch_index}_contig{contig}_{start}_{end}_flank{flank_size}"
                formatted_seq = '\n'.join([sequence[i:i+80] for i in range(0, len(sequence), 80)])
            
            final_content = f"{new_header}\n{formatted_seq}\n"
            
            with open(output_file, 'a') as f:
                f.write(final_content)
            
            self.log(f"  Patch{patch_index}: Successfully extracted {len(sequence):,} bp sequence", "info")
            
            return True, len(sequence)
            
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            self.log(f"  Patch{patch_index}: Failed to extract sequence: {error_msg}", "warning")
            return False, None
        except Exception as e:
            self.log(f"  Patch{patch_index}: Unknown error: {e}", "warning")
            return False, None
    
    def _process_single_gap(self, reference_fasta: str, alignments: List[Tuple], 
                          ref_contig_groups: Dict, qry_contig: str, gap_pos: int, 
                          flank_size: int, require_both_anchors: bool,
                          min_anchor_distance: int, max_anchor_distance: int,
                          output_file: str, start_patch_counter: int) -> Dict[str, Any]:
        """Process single gap position"""
        gap_result = {
            'query_contig': qry_contig,
            'gap_position': gap_pos,
            'valid_reference_contigs': 0,
            'successful_patches': 0,
            'total_patches_attempted': 0,
            'patches': []
        }
        
        valid_pairs = []
        
        for ref_contig, aligns in ref_contig_groups.items():
            qry_aligns = [align for align in aligns if align[5] == qry_contig]
            if len(qry_aligns) < 2:
                continue
            
            before_align, after_align = self._find_surrounding_alignments(
                qry_aligns, gap_pos, require_both_anchors
            )
            
            if before_align and after_align:
                _, _, q1s, q1e, _, _ = before_align
                _, _, q2s, q2e, _, _ = after_align
                
                q1_max = max(q1s, q1e)
                q2_min = min(q2s, q2e)
                
                anchor_distance = q2_min - q1_max
                
                if (min_anchor_distance <= anchor_distance <= max_anchor_distance):
                    valid_pairs.append((ref_contig, before_align, after_align, anchor_distance))
                else:
                    self.log(f"  Reference contig {ref_contig}: Anchor distance {anchor_distance:,} bp out of range", "debug")
        
        gap_result['valid_reference_contigs'] = len(valid_pairs)
        self.log(f"  Found {len(valid_pairs)} valid reference contigs", "info")
        
        if not valid_pairs:
            return gap_result
        
        success_count = 0
        patch_counter = start_patch_counter
        
        for i, (ref_contig, before_align, after_align, anchor_distance) in enumerate(valid_pairs):
            patch_num = patch_counter + i
            
            patch_start, patch_end, direction, ref_contig_name = self._get_patch_region_for_ref_contig(
                before_align, after_align, flank_size
            )
            
            if patch_start is None:
                continue
            
            success, seq_length = self._extract_patch_sequence(
                reference_fasta, ref_contig_name, patch_start, patch_end, 
                direction, flank_size, patch_num, gap_pos, output_file
            )
            
            if success:
                success_count += 1
                
                patch_info = {
                    'patch_id': patch_num,
                    'reference_contig': ref_contig_name,
                    'patch_coordinates': f"{patch_start}-{patch_end}",
                    'sequence_length': seq_length,
                    'direction': 'reverse_complement' if direction == -1 else 'forward',
                    'flank_size': flank_size,
                    'anchor_distance': anchor_distance,
                    'anchor_left': {
                        'query': f"{before_align[2]}-{before_align[3]}",
                        'reference': f"{before_align[0]}-{before_align[1]}",
                        'contig': before_align[4]
                    },
                    'anchor_right': {
                        'query': f"{after_align[2]}-{after_align[3]}",
                        'reference': f"{after_align[0]}-{after_align[1]}",
                        'contig': after_align[4]
                    }
                }
                gap_result['patches'].append(patch_info)
        
        gap_result['successful_patches'] = success_count
        gap_result['total_patches_attempted'] = len(valid_pairs)
        
        self.log(f"  Successfully extracted {success_count}/{len(valid_pairs)} patch sequences", "info")
        
        return gap_result
    
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
        self.log("Gap patch extraction completed!", "info")
        self.log("="*60, "info")
        
        self.log(f"Total processing time: {stats['processing_time_seconds']:.2f} seconds", "info")
        self.log(f"Successfully processed: {stats['successful_gaps']}/{stats['total_gaps']} gaps", "info")
        self.log(f"Total patches: {stats['total_patches']}", "info")
        
        if stats['total_patches'] > 0:
            self.log(f"Patch distribution:", "info")
            for count, gaps in sorted(stats['patches_per_gap'].items()):
                self.log(f"  {count} patches/gap: {gaps} gaps", "info")
        
        output_file = result_summary['parameters']['output_file']
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            self.log(f"Output file: {output_file} ({file_size:,} bytes)", "info")
    
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
    
    def batch_extract(
        self,
        reference_fastas: List[str],
        coords_files: List[str],
        gap_positions_list: List[List[int]],
        output_dir: str = ".",
        batch_prefix: str = "batch",
        **kwargs
    ) -> Dict[str, Any]:
        """
        Batch extract multiple gap patch sets
        
        Args:
            reference_fastas: List of reference genome files
            coords_files: List of alignment result files
            gap_positions_list: Corresponding gap position lists
            output_dir: Output directory
            batch_prefix: Batch processing prefix
            **kwargs: Other parameters passed to extract_patches
            
        Returns:
            Batch processing summary
        """
        if len(reference_fastas) != len(coords_files) or len(reference_fastas) != len(gap_positions_list):
            raise ValueError("All input lists must have same length")
        
        self.log(f"Starting batch extraction for {len(reference_fastas)} gap patch sets", "info")
        
        results = {}
        for i, (ref_fasta, coords_file, gap_positions) in enumerate(zip(reference_fastas, coords_files, gap_positions_list)):
            self.log(f"Processing extraction task {i+1}/{len(reference_fastas)}: {os.path.basename(coords_file)}", "info")
            
            output_file = os.path.join(output_dir, f"{batch_prefix}_{i+1:03d}_patches.fasta")
            output_json = os.path.join(output_dir, f"{batch_prefix}_{i+1:03d}_report.json")
            
            try:
                result = self.extract_patches(
                    reference_fasta=ref_fasta,
                    coords_file=coords_file,
                    gap_positions=gap_positions,
                    output_file=output_file,
                    output_json=output_json,
                    **kwargs
                )
                results[f"{os.path.basename(coords_file)}"] = {
                    "success": True,
                    "result": result,
                    "output_files": {
                        "patches_fasta": output_file,
                        "json_report": output_json
                    }
                }
            except Exception as e:
                results[f"{os.path.basename(coords_file)}"] = {
                    "success": False,
                    "error": str(e),
                    "output_files": {
                        "patches_fasta": output_file,
                        "json_report": output_json
                    }
                }
                self.log(f"Extraction failed: {e}", "warning")
        
        summary = {
            "total": len(reference_fastas),
            "successful": sum(1 for r in results.values() if r["success"]),
            "failed": sum(1 for r in results.values() if not r["success"]),
            "results": results,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        summary_file = os.path.join(output_dir, f"{batch_prefix}_summary.json")
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        self.log(f"Batch extraction completed, successful: {summary['successful']}/{summary['total']}", "info")
        self.log(f"Summary report saved: {summary_file}", "info")
        
        return summary
    
    def extract_from_analysis_result(
        self,
        nucmer_analysis_result: Dict[str, Any],
        gap_positions: List[int],
        chromosome: Optional[str] = None,
        output_file: str = "extracted_patches.fasta",
        flank_size: int = 10000,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Extract gap patches from nucmer analysis results
        
        Args:
            nucmer_analysis_result: Nucmer analysis result dictionary
            gap_positions: List of gap positions
            chromosome: Specify chromosome (optional, None=first chromosome)
            output_file: Output file path
            flank_size: Flanking sequence length
            **kwargs: Other parameters
            
        Returns:
            Extraction result dictionary
        """
        reference_fasta = nucmer_analysis_result.get("parameters", {}).get("ref_fasta")
        if not reference_fasta or not os.path.exists(reference_fasta):
            raise ValueError(f"Reference genome file not found or not specified: {reference_fasta}")
        
        chromosome_results = nucmer_analysis_result.get("chromosome_results", {})
        
        if chromosome:
            if chromosome not in chromosome_results:
                raise ValueError(f"Chromosome {chromosome} not found in analysis results")
            chr_data = chromosome_results[chromosome]
        else:
            if not chromosome_results:
                raise ValueError("No chromosome data in analysis results")
            chromosome = list(chromosome_results.keys())[0]
            chr_data = chromosome_results[chromosome]
        
        chr_dir = chr_data.get("dir", "")
        coords_file = None
        
        for file in os.listdir(chr_dir):
            if file.endswith('_merged.coords'):
                coords_file = os.path.join(chr_dir, file)
                break
        
        if not coords_file or not os.path.exists(coords_file):
            raise ValueError(f"Coords file not found: {chr_dir}")
        
        self.log(f"Extracting gap patches from chromosome {chromosome}", "info")
        self.log(f"Using coords file: {coords_file}", "info")
        
        return self.extract_patches(
            reference_fasta=reference_fasta,
            coords_file=coords_file,
            gap_positions=gap_positions,
            output_file=output_file,
            flank_size=flank_size,
            **kwargs
        )
    
    def cleanup(self):
        """Clean up all temporary files"""
        self._cleanup_temp_files()
    
    def __del__(self):
        """Destructor, automatically clean temporary files"""
        self.cleanup()

def _parse_coord_file_original(coord_file: str) -> List[Tuple]:
    """Parse coords file and extract alignment information (original implementation)"""
    alignments = []
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
                    qry_start = int(parts[3])
                    qry_end = int(parts[4])
                    ref_contig = parts[-2]
                    qry_contig = parts[-1]
                    alignments.append((ref_start, ref_end, qry_start, qry_end, ref_contig, qry_contig))
                except:
                    continue
    return alignments

def _group_alignments_by_ref_contig_original(alignments: List[Tuple]) -> Dict:
    """Group syntenic regions by reference contig (original implementation)"""
    ref_contig_groups = {}
    for align in alignments:
        *_, ref_contig, qry_contig = align
        if ref_contig not in ref_contig_groups:
            ref_contig_groups[ref_contig] = []
        ref_contig_groups[ref_contig].append(align)
    return ref_contig_groups

def _find_surrounding_alignments_original(alignments: List[Tuple], gap_pos: int):
    """Find closest left and right syntenic blocks around gap (original implementation)"""
    left_candidates = []
    right_candidates = []

    for align in alignments:
        _, _, qry_start, qry_end, _, _ = align
        q_min = min(qry_start, qry_end)
        q_max = max(qry_start, qry_end)

        if q_max <= gap_pos:
            left_candidates.append((q_max, align))
        elif q_min > gap_pos:
            right_candidates.append((q_min, align))

    before_align = None
    after_align = None

    if left_candidates:
        before_align = max(left_candidates, key=lambda x: x[0])[1]

    if right_candidates:
        after_align = min(right_candidates, key=lambda x: x[0])[1]

    return before_align, after_align

def extract_gap_patches(
    reference_fasta: str,
    coords_file: str,
    gap_positions: List[int],
    output_file: str = "gap_patches.fasta",
    flank_size: int = 10000,
    verbose: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Simplified interface for extracting gap patches
    
    Example:
        >>> from extract_gap_patches import extract_gap_patches
        >>> result = extract_gap_patches(
        ...     reference_fasta="reference.fasta",
        ...     coords_file="alignment.coords",
        ...     gap_positions=[1000000, 2000000],
        ...     output_file="patches.fasta",
        ...     flank_size=5000,
        ...     verbose=False
        ... )
    """
    extractor = GapPatchesExtractorAPI(verbose=verbose)
    return extractor.extract_patches(
        reference_fasta=reference_fasta,
        coords_file=coords_file,
        gap_positions=gap_positions,
        output_file=output_file,
        flank_size=flank_size,
        **kwargs
    )

def main():
    """Command line interface main function"""
    parser = argparse.ArgumentParser(
        description='Extract gap patch sequences from coords files (Unified API Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Basic usage, original mode
  python extract_gap_patches.py -c reference.fasta -coords alignment.coords --gap-position 1804444 12272394
  
  # API mode, more parameters
  python extract_gap_patches.py -c reference.fasta -coords alignment.coords --gap-position 1804444 12272394 --api-mode --output-dir results --flank-size 20000
  
  # Python module call
  from extract_gap_patches import GapPatchesExtractorAPI
  extractor = GapPatchesExtractorAPI(verbose=True)
  result = extractor.extract_patches(
      reference_fasta="reference.fasta",
      coords_file="alignment.coords",
      gap_positions=[1804444, 12272394],
      output_file="patches.fasta",
      flank_size=10000
  )
        """
    )
    
    parser.add_argument("-c", "--contigs", dest="reference_fasta", required=True,
                       help="Reference contigs file")
    parser.add_argument("-coords", required=True, help="Coords file (nucmer alignment results)")
    parser.add_argument("--gap-position", dest="gap_positions", type=int, nargs='+', required=True,
                       help="Gap position coordinates (multiple, space separated)")
    
    parser.add_argument("-f", "--flank-size", type=int, default=10000,
                       help="Flanking extension distance on both sides (bp), default 10000")
    parser.add_argument("-o", "--output", default="gap_patches.fasta",
                       help="Output fasta file, default: gap_patches.fasta")
    parser.add_argument("-j", "--json", default="extraction_results.json",
                       help="Result JSON file, default: extraction_results.json")
    parser.add_argument("--output-dir", default=".", help="Output directory, default current directory")
    parser.add_argument("--require-both-anchors", action="store_true", default=True,
                       help="Require both left and right anchors (default)")
    parser.add_argument("--allow-single-anchor", action="store_false", dest="require_both_anchors",
                       help="Allow single anchor")
    parser.add_argument("--min-anchor-distance", type=int, default=1000,
                       help="Minimum anchor distance (bp), default 1000")
    parser.add_argument("--max-anchor-distance", type=int, default=1000000,
                       help="Maximum anchor distance (bp), default 1000000")
    parser.add_argument("--api-mode", action="store_true", help="Use API mode")
    parser.add_argument("--verbose", action="store_true", help="Show detailed output")
    parser.add_argument("--quiet", action="store_true", help="Quiet mode")
    parser.add_argument("--log-file", help="Log file path")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.reference_fasta):
        print(f"Error: Reference contigs file {args.reference_fasta} not found", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.coords):
        print(f"Error: Coords file {args.coords} not found", file=sys.stderr)
        sys.exit(1)
    
    verbose = args.verbose and not args.quiet
    
    try:
        if args.api_mode:
            extractor = GapPatchesExtractorAPI(verbose=verbose, log_file=args.log_file)
            
            output_file = os.path.join(args.output_dir, args.output)
            output_json = os.path.join(args.output_dir, args.json)
            
            result = extractor.extract_patches(
                reference_fasta=args.reference_fasta,
                coords_file=args.coords,
                gap_positions=args.gap_positions,
                output_file=output_file,
                flank_size=args.flank_size,
                require_both_anchors=args.require_both_anchors,
                min_anchor_distance=args.min_anchor_distance,
                max_anchor_distance=args.max_anchor_distance,
                output_json=output_json
            )
            
            if verbose:
                print(f"\nDetailed report saved to: {output_json}")
                print(f"Patch sequences saved to: {output_file}")
            
        else:
            from Bio import SeqIO
            from Bio.Seq import Seq
            
            if args.flank_size < 0:
                print("Error: flank_size must be non-negative integer", file=sys.stderr)
                sys.exit(1)
            
            if not os.path.exists(args.reference_fasta + ".fai"):
                try:
                    subprocess.run(["samtools", "faidx", args.reference_fasta], 
                                 check=True, stderr=subprocess.DEVNULL)
                except Exception as e:
                    print(f"Error: Failed to create index for reference genome: {e}", file=sys.stderr)
                    sys.exit(1)
            
            print(f"Processing {len(args.gap_positions)} gap positions...")
            print(f"Reference file: {args.reference_fasta}")
            print(f"Coords file: {args.coords}")
            print(f"Flanking length: {args.flank_size} bp")
            print(f"Output file: {args.output}")
            
            extractor = GapPatchesExtractorAPI(verbose=verbose)
            result = extractor.extract_patches(
                reference_fasta=args.reference_fasta,
                coords_file=args.coords,
                gap_positions=args.gap_positions,
                output_file=args.output,
                flank_size=args.flank_size,
                output_json=args.json
            )
            
            stats = result['statistics']
            print(f"\nProcessing completed!")
            print(f"Processing time: {stats['processing_time_seconds']:.2f} seconds")
            print(f"Successfully processed: {stats['successful_gaps']}/{stats['total_gaps']} gaps")
            print(f"Total patches: {stats['total_patches']}")
            print(f"Output file: {args.output} ({os.path.getsize(args.output) if os.path.exists(args.output) else 0} bytes)")
            print(f"Result details: {args.json}")
        
        stats = result['statistics']
        if stats['successful_gaps'] == 0:
            print(f"Warning: No patches successfully extracted", file=sys.stderr)
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()