#!/usr/bin/env python3
"""
Genome Gap Analysis and Repair Tool based on Synteny - Unified API Version (Enhanced Repair Mode)
Fixed version: Fixed return value unpacking error
"""

import sys
import os
import json
import logging
import time
from collections import defaultdict, deque
from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple, Optional, Set, Any, Union
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

### Unified API Interface ###

class GapAnalyzerAPI:
    """Unified API interface for Gap analysis"""
    
    def __init__(self, verbose: bool = True, log_file: Optional[str] = None):
        """
        Initialize GapAnalyzer API
        
        Args:
            verbose: Whether to show detailed output
            log_file: Log file path (optional)
        """
        self.verbose = verbose
        self.log_file = log_file
        self._setup_logging()
    
    def _setup_logging(self):
        """Setup logging system"""
        self.logger = logging.getLogger('GapAnalyzer')
        self.logger.setLevel(logging.INFO if self.verbose else logging.WARNING)
        
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO if self.verbose else logging.WARNING)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        # File handler (if specified)
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
    
    def log_repair_decision(self, gap_pos: int, repair_mode: str, 
                           decision: bool, details: str = ""):
        """Log repair decision"""
        decision_text = "Repair" if decision else "No repair"
        message = f"Gap {gap_pos:,}: {repair_mode} mode -> {decision_text}"
        if details:
            message += f" ({details})"
        self.log(message, "info")
    
    def analyze_gaps(
        self,
        coords_file: str,
        gap_positions: List[int],
        query_fasta: Optional[str] = None,
        output_dir: str = ".",
        output_prefix: str = "gap_analysis",
        max_search_distance: int = 500000,
        search_step: int = 100000,
        final_gap_length: int = 100,
        min_confidence: str = "low",
        repair_mode: str = "conservative"  # "conservative" or "aggressive"
    ) -> Dict[str, Any]:
        """
        Analyze gap positions and generate repair suggestions
        
        Args:
            coords_file: nucmer alignment result coords file
            gap_positions: List of gap positions (based on query coordinates)
            query_fasta: Query genome FASTA file (optional, for repair)
            output_dir: Output directory
            output_prefix: Output file prefix
            max_search_distance: Maximum search distance (bp)
            search_step: Search step size (bp)
            final_gap_length: Final unified gap length (bp)
            min_confidence: Minimum confidence threshold ("low", "medium", "high")
            repair_mode: Repair mode
            
        Returns:
            Dictionary containing analysis results
        """
        self.log(f"Starting gap analysis, repair mode: {repair_mode}", "info")
        self.log(f"Alignment file: {coords_file}", "info")
        self.log(f"Analyzing {len(gap_positions)} gap positions", "info")
        
        start_time = time.time()
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Create analyzer instance
        analyzer = _GapAnalyzer()
        
        # Modify step size parameter in find_surrounding_alignments function
        original_method = analyzer.find_surrounding_alignments
        def new_find_surrounding_alignments(query_contig, gap_pos, max_distance):
            return original_method(query_contig, gap_pos, max_distance, search_step)
        analyzer.find_surrounding_alignments = new_find_surrounding_alignments
        
        # Parse coords file
        self.log("Parsing alignment file and building synteny blocks...", "info")
        alignments = analyzer.parse_coords_file(coords_file)
        
        if not alignments:
            raise ValueError("No alignment information parsed")
        
        # Analyze all gaps
        self.log(f"Analyzing all gaps (including synteny analysis)...", "info")
        gap_analyses = analyzer.analyze_all_gaps(gap_positions, max_search_distance)
        
        if not gap_analyses:
            raise ValueError("No gaps successfully analyzed")
        
        # Apply repairs
        repair_log = []
        repair_stats = {}
        filled_fasta = None
        
        if query_fasta and os.path.exists(query_fasta):
            self.log(f"Applying repairs (mode: {repair_mode})...", "info")
            
            # Generate output file paths
            output_fasta = os.path.join(output_dir, f"{output_prefix}_repaired.fasta")
            
            try:
                filled_fasta, repair_log, repair_stats = analyzer.apply_unified_repairs(
                    query_fasta, gap_analyses, output_fasta, final_gap_length, repair_mode
                )
                
                self.log(f"Repaired genome saved: {filled_fasta}", "info")
                self.log(f"Performed {len(repair_log)} replacements", "info")
                
                # Record repair statistics
                if repair_stats:
                    self.log(f"Repair statistics: {repair_stats}", "info")
                    
            except Exception as e:
                self.log(f"Error during repair process: {e}", "error")
                raise
        
        # Save report
        report_file = os.path.join(output_dir, f"{output_prefix}_report.json")
        report = analyzer.save_report(gap_analyses, repair_log, report_file, query_fasta, repair_mode)
        
        # Add repair statistics to report
        if repair_stats:
            report.setdefault('repair_stats', {}).update(repair_stats)
        
        # Add API-specific information
        report["api_info"] = {
            "execution_time": time.time() - start_time,
            "input_files": {
                "coords_file": coords_file,
                "query_fasta": query_fasta,
                "gap_positions": gap_positions
            },
            "parameters": {
                "max_search_distance": max_search_distance,
                "search_step": search_step,
                "final_gap_length": final_gap_length,
                "min_confidence": min_confidence,
                "repair_mode": repair_mode
            },
            "output_files": {
                "report": report_file,
                "filled_fasta": filled_fasta
            }
        }
        
        # Print summary
        self._print_summary(report)
        
        return report
    
    def _print_summary(self, report: Dict[str, Any]):
        """Print analysis summary"""
        summary = report.get("summary", {})
        
        self.log("\n" + "="*60, "info")
        self.log("Gap analysis completed!", "info")
        self.log("="*60, "info")
        
        self.log(f"Successfully analyzed: {summary.get('total_gaps_analyzed', 0)} gaps", "info")
        self.log(f"Gaps in error regions: {summary.get('gaps_in_error_region', 0)}", "info")
        self.log(f"Gaps not in error regions: {summary.get('gaps_not_in_error_region', 0)}", "info")
        
        repair_mode = report.get("api_info", {}).get("parameters", {}).get("repair_mode", "conservative")
        self.log(f"Repair mode: {repair_mode}", "info")
        
        if "gaps_with_large_distance_anchor" in summary:
            self.log(f"Large distance anchor warnings: {summary['gaps_with_large_distance_anchor']}", "info")
        
        if "repair_summary" in summary:
            repair_summary = summary["repair_summary"]
            self.log(f"\nRepair decision statistics ({repair_mode} mode):", "info")
            self.log(f"  Recommended repairs: {repair_summary.get('recommended_repairs', 0)}", "info")
            self.log(f"  Skipped (conservative conditions not met): {repair_summary.get('skipped_conservative', 0)}", "info")
            if repair_mode == "aggressive":
                self.log(f"  Additional repairs in aggressive mode: {repair_summary.get('aggressive_repairs', 0)}", "info")
                self.log(f"  Skipped in aggressive mode (insufficient anchors): {repair_summary.get('skipped_aggressive', 0)}", "info")
        
        self.log("\nError type distribution:", "info")
        if "error_type_counts" in summary:
            for err_type, count in summary["error_type_counts"].items():
                self.log(f"  {err_type}: {count}", "info")
        
        if "repair_details" in report and report["repair_details"]:
            total_replaced = sum(r.get('original_length', 0) for r in report["repair_details"])
            self.log(f"\nTotal replacement length: {total_replaced:,} bp", "info")
            self.log(f"All gaps finally unified to: {report['api_info']['parameters']['final_gap_length']} bp", "info")
    
    def batch_analyze(
        self,
        coords_files: List[str],
        gap_positions_list: List[List[int]],
        query_fastas: Optional[List[str]] = None,
        output_dir: str = ".",
        batch_prefix: str = "batch",
        **kwargs
    ) -> Dict[str, Any]:
        """
        Batch analyze multiple alignment results
        
        Args:
            coords_files: List of alignment result files
            gap_positions_list: Corresponding gap position lists
            query_fastas: List of query genome files (optional)
            output_dir: Output directory
            batch_prefix: Batch processing prefix
            **kwargs: Other parameters passed to analyze_gaps
            
        Returns:
            Batch processing result summary
        """
        if len(coords_files) != len(gap_positions_list):
            raise ValueError("coords_files and gap_positions_list must have same length")
        
        if query_fastas and len(query_fastas) != len(coords_files):
            raise ValueError("query_fastas length must match coords_files")
        
        self.log(f"Starting batch analysis of {len(coords_files)} alignment results", "info")
        
        results = {}
        for i, (coords_file, gap_positions) in enumerate(zip(coords_files, gap_positions_list)):
            self.log(f"Processing {i+1}/{len(coords_files)}: {os.path.basename(coords_file)}", "info")
            
            query_fasta = query_fastas[i] if query_fastas else None
            output_prefix = f"{batch_prefix}_{i+1:03d}"
            
            try:
                result = self.analyze_gaps(
                    coords_file=coords_file,
                    gap_positions=gap_positions,
                    query_fasta=query_fasta,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    **kwargs
                )
                results[os.path.basename(coords_file)] = {
                    "success": True,
                    "result": result,
                    "output_prefix": output_prefix
                }
            except Exception as e:
                results[os.path.basename(coords_file)] = {
                    "success": False,
                    "error": str(e),
                    "output_prefix": output_prefix
                }
                self.log(f"Analysis failed: {e}", "warning")
        
        # Generate summary report
        summary = {
            "total": len(coords_files),
            "successful": sum(1 for r in results.values() if r["success"]),
            "failed": sum(1 for r in results.values() if not r["success"]),
            "results": results,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        # Save summary report
        summary_file = os.path.join(output_dir, f"{batch_prefix}_summary.json")
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        self.log(f"Batch analysis completed, successful: {summary['successful']}/{summary['total']}", "info")
        self.log(f"Summary report saved: {summary_file}", "info")
        
        return summary

### Internal implementation classes (preserving original functionality) ###

@dataclass
class _Alignment:
    """Store alignment information"""
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    ref_contig: str
    query_contig: str
    identity: float
    
    @property
    def ref_min(self):
        return min(self.ref_start, self.ref_end)
    
    @property
    def ref_max(self):
        return max(self.ref_start, self.ref_end)
    
    @property
    def query_min(self):
        return min(self.query_start, self.query_end)
    
    @property
    def query_max(self):
        return max(self.query_start, self.query_end)
    
    @property
    def direction(self):
        """1 for forward, -1 for reverse"""
        return 1 if self.query_start < self.query_end else -1
    
    @property
    def length(self):
        return abs(self.query_end - self.query_start) + 1
    
    def __str__(self):
        dir_str = "→" if self.direction == 1 else "←"
        return f"{self.query_contig}:{self.query_min:,}-{self.query_max:,}{dir_str}{self.ref_contig}:{self.ref_min:,}-{self.ref_max:,}"

@dataclass
class _SyntenyBlock:
    """Synteny block information"""
    ref_contig: str
    start_pos: int
    end_pos: int
    direction: int
    alignments: List[_Alignment]
    
    @property
    def length(self):
        return self.end_pos - self.start_pos + 1
    
    @property
    def query_contigs(self) -> Set[str]:
        """Included query contigs"""
        return set(align.query_contig for align in self.alignments)
    
    @property
    def avg_identity(self) -> float:
        """Average alignment identity"""
        if not self.alignments:
            return 0.0
        return sum(align.identity for align in self.alignments) / len(self.alignments)

@dataclass
class _GapAnalysis:
    """Gap analysis result"""
    gap_pos: int
    query_contig: str
    supporting_ref_contigs: Dict[str, Tuple[Optional[_Alignment], Optional[_Alignment]]]
    error_type: str
    replace_start: int
    replace_end: int
    confidence: str
    evidence_details: Dict[str, List[str]]
    gap_in_error_region: bool = True
    synteny_analysis: Optional[Dict] = None
    large_distance_anchor: bool = False
    total_anchor_length: int = 0  # New: total anchor length
    repair_reason: Optional[str] = None  # New: repair reason
    
    @property
    def replace_length(self):
        return self.replace_end - self.replace_start + 1 if self.replace_start > 0 else 0

class _GapAnalyzer:
    """Internal Gap analyzer (enhanced repair mode)"""
    
    def __init__(self):
        self.coords_file: Optional[str] = None
        self.alignments: List[_Alignment] = []
        self.query_contig_groups = defaultdict(list)
        self.ref_contig_groups = defaultdict(list)
        self.synteny_blocks: Dict[str, List[_SyntenyBlock]] = defaultdict(list)
    
    def parse_coords_file(self, coords_file: str):
        """Parse coords file"""
        self.coords_file = coords_file
        alignments = []
        
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
                        identity = float(parts[6]) if len(parts) > 6 else 100.0
                        
                        alignment = _Alignment(
                            ref_start=ref_start,
                            ref_end=ref_end,
                            query_start=qry_start,
                            query_end=qry_end,
                            ref_contig=ref_contig,
                            query_contig=qry_contig,
                            identity=identity
                        )
                        alignments.append(alignment)
                    except:
                        continue
        
        self.alignments = alignments
        
        # Group by query contig and ref contig
        for align in alignments:
            self.query_contig_groups[align.query_contig].append(align)
            self.ref_contig_groups[align.ref_contig].append(align)
        
        # Sort
        for contig in self.query_contig_groups:
            self.query_contig_groups[contig].sort(key=lambda x: x.query_min)
        for contig in self.ref_contig_groups:
            self.ref_contig_groups[contig].sort(key=lambda x: x.ref_min)
        
        # Build synteny blocks
        self._build_synteny_blocks()
        
        print(f"Successfully parsed {len(alignments)} alignment regions")
        print(f"Involving {len(self.ref_contig_groups)} reference contigs")
        print(f"Built {sum(len(blocks) for blocks in self.synteny_blocks.values())} synteny blocks")
        return alignments
    
    def _build_synteny_blocks(self):
        """Build synteny blocks on reference contigs"""
        for ref_contig, aligns in self.ref_contig_groups.items():
            if not aligns:
                continue
            
            # Sort by reference coordinates
            aligns_sorted = sorted(aligns, key=lambda x: x.ref_min)
            
            blocks = []
            current_block = None
            
            for align in aligns_sorted:
                if current_block is None:
                    # Start new block
                    current_block = _SyntenyBlock(
                        ref_contig=ref_contig,
                        start_pos=align.ref_min,
                        end_pos=align.ref_max,
                        direction=align.direction,
                        alignments=[align]
                    )
                else:
                    # Check if can merge into current block
                    gap_size = align.ref_min - current_block.end_pos
                    if (align.direction == current_block.direction and 
                        gap_size < 100000 and 
                        gap_size > -10000):
                        # Merge into current block
                        current_block.end_pos = max(current_block.end_pos, align.ref_max)
                        current_block.alignments.append(align)
                    else:
                        # Save current block, start new block
                        blocks.append(current_block)
                        current_block = _SyntenyBlock(
                            ref_contig=ref_contig,
                            start_pos=align.ref_min,
                            end_pos=align.ref_max,
                            direction=align.direction,
                            alignments=[align]
                        )
            
            if current_block:
                blocks.append(current_block)
            
            self.synteny_blocks[ref_contig] = blocks
    
    def find_surrounding_alignments(self, query_contig: str, gap_pos: int, 
                                  max_distance: int = 500000, step: int = 100000):
        """Find alignments around gap, dynamically expanding search range"""
        contig_aligns = self.query_contig_groups.get(query_contig, [])
        
        for window in range(step, max_distance + step, step):
            left_candidates = []
            right_candidates = []
            crossing_candidates = []
            
            for align in contig_aligns:
                # Check if crosses gap
                if align.query_min <= gap_pos <= align.query_max:
                    crossing_candidates.append(align)
                # Check if within window
                elif abs(align.query_min - gap_pos) <= window or abs(align.query_max - gap_pos) <= window:
                    if align.query_max < gap_pos:
                        left_candidates.append(align)
                    elif align.query_min > gap_pos:
                        right_candidates.append(align)
            
            # If found alignments on both sides or crossing gap, return results
            if crossing_candidates or (left_candidates and right_candidates):
                return left_candidates, right_candidates, crossing_candidates, window
        
        # If not found, return empty results
        return [], [], [], max_distance
    
    def analyze_synteny_for_ref_contig(self, ref_contig: str, gap_pos: int, 
                                     query_contig: str) -> Dict:
        """
        Analyze synteny on specific reference contig
        """
        result = {
            'has_synteny': False,
            'blocks_near_gap': [],
            'anchor_distance': float('inf'),
            'assessment': 'no_data',
            'synteny_blocks': []
        }
        
        # Get synteny blocks for this reference contig
        blocks = self.synteny_blocks.get(ref_contig, [])
        if not blocks:
            result['assessment'] = 'no_synteny_blocks'
            return result
        
        result['synteny_blocks'] = [{
            'start': block.start_pos,
            'end': block.end_pos,
            'direction': block.direction,
            'length': block.length,
            'avg_identity': block.avg_identity,
            'query_contigs': list(block.query_contigs)
        } for block in blocks]
        
        # Find reference coordinates corresponding to gap position
        query_aligns = [a for a in self.alignments 
                       if a.query_contig == query_contig and a.ref_contig == ref_contig]
        
        if not query_aligns:
            result['assessment'] = 'no_alignment_for_query'
            return result
        
        # Try to infer gap position on reference genome
        ref_positions = []
        for align in query_aligns:
            if align.direction == 1:
                query_to_ref_ratio = (align.ref_max - align.ref_min) / (align.query_max - align.query_min)
                if align.query_min <= gap_pos <= align.query_max:
                    offset = gap_pos - align.query_min
                    ref_pos = align.ref_min + offset * query_to_ref_ratio
                    ref_positions.append(int(ref_pos))
            else:
                query_to_ref_ratio = (align.ref_max - align.ref_min) / (align.query_max - align.query_min)
                if align.query_min <= gap_pos <= align.query_max:
                    offset = gap_pos - align.query_min
                    ref_pos = align.ref_max - offset * query_to_ref_ratio
                    ref_positions.append(int(ref_pos))
        
        if not ref_positions:
            # If cannot directly map, use nearest alignments for inference
            left_aligns = [a for a in query_aligns if a.query_max < gap_pos]
            right_aligns = [a for a in query_aligns if a.query_min > gap_pos]
            
            if left_aligns and right_aligns:
                left_align = max(left_aligns, key=lambda x: x.query_max)
                right_align = min(right_aligns, key=lambda x: x.query_min)
                
                query_dist_left = gap_pos - left_align.query_max
                query_dist_right = right_align.query_min - gap_pos
                
                result['anchor_distance'] = min(query_dist_left, query_dist_right)
                
                for block in blocks:
                    block_align_refs = [a.ref_min for a in block.alignments]
                    if (left_align.ref_min in block_align_refs and 
                        right_align.ref_min in block_align_refs):
                        result['has_synteny'] = True
                        result['assessment'] = 'synteny_present'
                        result['blocks_near_gap'].append({
                            'block': block,
                            'distance_to_gap': min(query_dist_left, query_dist_right)
                        })
                        
                        if result['anchor_distance'] > 500000:
                            result['assessment'] = 'large_distance_with_synteny'
                        break
                else:
                    result['assessment'] = 'no_synteny_between_anchors'
            else:
                result['assessment'] = 'insufficient_anchors'
        else:
            avg_ref_pos = sum(ref_positions) / len(ref_positions)
            result['has_synteny'] = True
            
            for block in blocks:
                if block.start_pos <= avg_ref_pos <= block.end_pos:
                    result['blocks_near_gap'].append({
                        'block': block,
                        'distance_to_gap': 0,
                        'ref_position': avg_ref_pos
                    })
                    result['assessment'] = 'gap_in_synteny_block'
                    break
            else:
                result['assessment'] = 'gap_outside_synteny_blocks'
        
        return result
    
    def check_gap_in_error_region(self, gap_pos: int, supporting_refs: Dict, 
                                error_type: str, crossing_aligns: List[_Alignment]) -> bool:
        """Check if gap is in error region"""
        # Case 1: Has alignments crossing gap
        if crossing_aligns:
            return True
        
        # Case 2: Not enough supporting alignments
        if not supporting_refs:
            return False
        
        # Analyze each error type
        if error_type == "type1":
            return True
            
        elif error_type == "type2":
            for left_align, right_align in supporting_refs.values():
                if left_align.direction != right_align.direction:
                    conflict_region_start = left_align.query_max
                    conflict_region_end = right_align.query_min
                    
                    if conflict_region_start <= gap_pos <= conflict_region_end:
                        return True
                    else:
                        return False
                        
        elif error_type == "type3":
            for left_align, right_align in supporting_refs.values():
                if left_align.ref_contig != right_align.ref_contig:
                    junction_region_start = left_align.query_max
                    junction_region_end = right_align.query_min
                    
                    distance_to_junction = min(abs(gap_pos - junction_region_start),
                                              abs(gap_pos - junction_region_end))
                    
                    if distance_to_junction <= 5000:
                        return True
                    else:
                        return False
                        
        elif error_type == "type4":
            for left_align, right_align in supporting_refs.values():
                if self.check_overlap(left_align, right_align):
                    if (left_align.query_min <= gap_pos <= left_align.query_max or
                        right_align.query_min <= gap_pos <= right_align.query_max):
                        return True
                    else:
                        return False
        
        return False
    
    def analyze_gap(self, gap_pos: int, max_search_distance: int = 500000) -> Optional[_GapAnalysis]:
        """Analyze single gap"""
        if not self.query_contig_groups:
            return None
        
        query_contig = list(self.query_contig_groups.keys())[0]
        
        # Find surrounding alignments
        left_aligns, right_aligns, crossing_aligns, window_used = self.find_surrounding_alignments(
            query_contig, gap_pos, max_search_distance
        )
        
        # Check if there are alignments crossing gap
        if crossing_aligns:
            supporting_refs = defaultdict(lambda: {"left": [], "right": []})
            for align in crossing_aligns:
                supporting_refs[align.ref_contig]["left"].append(align)
                supporting_refs[align.ref_contig]["right"].append(align)
        elif left_aligns and right_aligns:
            supporting_refs = defaultdict(lambda: {"left": [], "right": []})
            for align in left_aligns:
                supporting_refs[align.ref_contig]["left"].append(align)
            for align in right_aligns:
                supporting_refs[align.ref_contig]["right"].append(align)
        else:
            # Perform synteny analysis
            synteny_results = {}
            ref_contigs_to_analyze = list(self.ref_contig_groups.keys())
            
            for ref_contig in ref_contigs_to_analyze:
                synteny_result = self.analyze_synteny_for_ref_contig(ref_contig, gap_pos, query_contig)
                synteny_results[ref_contig] = synteny_result
            
            # Analyze synteny results
            valid_ref_contigs = []
            large_distance_refs = []
            no_synteny_refs = []
            
            for ref_contig, result in synteny_results.items():
                if result['assessment'] in ['synteny_present', 'gap_in_synteny_block']:
                    valid_ref_contigs.append(ref_contig)
                    if result.get('anchor_distance', 0) > 500000:
                        large_distance_refs.append(ref_contig)
                elif result['assessment'] == 'large_distance_with_synteny':
                    large_distance_refs.append(ref_contig)
                elif result['assessment'] in ['no_synteny_blocks', 'no_alignment_for_query', 
                                            'no_synteny_between_anchors']:
                    no_synteny_refs.append(ref_contig)
            
            # Build supporting alignments
            supporting_refs = defaultdict(lambda: {"left": [], "right": []})
            
            for ref_contig in valid_ref_contigs + large_distance_refs:
                ref_aligns = [a for a in self.alignments 
                            if a.ref_contig == ref_contig and a.query_contig == query_contig]
                
                if not ref_aligns:
                    continue
                
                left_for_ref = [a for a in ref_aligns if a.query_max < gap_pos]
                right_for_ref = [a for a in ref_aligns if a.query_min > gap_pos]
                
                if left_for_ref:
                    closest_left = max(left_for_ref, key=lambda x: x.query_max)
                    supporting_refs[ref_contig]["left"].append(closest_left)
                
                if right_for_ref:
                    closest_right = min(right_for_ref, key=lambda x: x.query_min)
                    supporting_refs[ref_contig]["right"].append(closest_right)
            
            if not supporting_refs:
                return None
        
        # Select closest left/right alignments for each reference contig
        final_supporting = {}
        evidence_details = defaultdict(list)
        large_distance_detected = False
        total_anchor_length = 0  # Calculate total anchor length
        
        for ref_contig, aligns in supporting_refs.items():
            left_best = max(aligns["left"], key=lambda x: x.query_max) if aligns["left"] else None
            right_best = min(aligns["right"], key=lambda x: x.query_min) if aligns["right"] else None
            
            if left_best and right_best:
                left_distance = gap_pos - left_best.query_max
                right_distance = right_best.query_min - gap_pos
                min_distance = min(left_distance, right_distance)
                
                final_supporting[ref_contig] = (left_best, right_best)
                evidence_details[ref_contig].append(f"Left: {left_best} (distance: {left_distance:,} bp)")
                evidence_details[ref_contig].append(f"Right: {right_best} (distance: {right_distance:,} bp)")
                
                # Calculate anchor length for this reference contig
                anchor_length = left_best.length + right_best.length
                total_anchor_length += anchor_length
                evidence_details[ref_contig].append(f"Anchor length: {anchor_length:,} bp")
                
                if min_distance > 500000:
                    large_distance_detected = True
                    evidence_details[ref_contig].append(f"Warning: Anchor distance too large ({min_distance:,} bp)")
            elif left_best:
                distance = gap_pos - left_best.query_max
                evidence_details[ref_contig].append(f"Left only: {left_best} (distance: {distance:,} bp)")
            elif right_best:
                distance = right_best.query_min - gap_pos
                evidence_details[ref_contig].append(f"Right only: {right_best} (distance: {distance:,} bp)")
        
        # Determine error type
        error_type, confidence = self.determine_error_type(final_supporting, crossing_aligns)
        
        if not error_type:
            return None
        
        # Check if gap is in error region
        gap_in_error_region = self.check_gap_in_error_region(
            gap_pos, final_supporting, error_type, crossing_aligns
        )
        
        # If no left/right anchors, perform more detailed synteny analysis
        synteny_analysis = None
        if not left_aligns or not right_aligns:
            synteny_analysis = {}
            for ref_contig in list(self.ref_contig_groups.keys())[:5]:
                result = self.analyze_synteny_for_ref_contig(ref_contig, gap_pos, query_contig)
                synteny_analysis[ref_contig] = {
                    'assessment': result['assessment'],
                    'has_synteny': result['has_synteny'],
                    'anchor_distance': result.get('anchor_distance', 'N/A')
                }
        
        if not gap_in_error_region:
            analysis = _GapAnalysis(
                gap_pos=gap_pos,
                query_contig=query_contig,
                supporting_ref_contigs=final_supporting,
                error_type=error_type,
                replace_start=0,
                replace_end=0,
                confidence=confidence,
                evidence_details=evidence_details,
                gap_in_error_region=False,
                synteny_analysis=synteny_analysis,
                large_distance_anchor=large_distance_detected,
                total_anchor_length=total_anchor_length
            )
            return analysis
        
        # Determine replacement region
        replace_start, replace_end = self.determine_replace_region(
            gap_pos, final_supporting, error_type, crossing_aligns
        )
        
        analysis = _GapAnalysis(
            gap_pos=gap_pos,
            query_contig=query_contig,
            supporting_ref_contigs=final_supporting,
            error_type=error_type,
            replace_start=replace_start,
            replace_end=replace_end,
            confidence=confidence,
            evidence_details=evidence_details,
            gap_in_error_region=True,
            synteny_analysis=synteny_analysis,
            large_distance_anchor=large_distance_detected,
            total_anchor_length=total_anchor_length
        )
        
        return analysis
    
    def _should_repair_in_aggressive_mode(self, analysis: _GapAnalysis) -> bool:
        """
        Determine if should repair in aggressive mode
        Condition: If total anchor region length > region length to replace, can replace
        """
        if not analysis.gap_in_error_region:
            return False
        
        if not analysis.supporting_ref_contigs:
            return False
        
        # Use pre-calculated total anchor length
        total_anchor_length = analysis.total_anchor_length
        
        # Calculate region length to replace
        replace_length = analysis.replace_length
        
        if replace_length <= 0:
            return False
        
        # Aggressive mode judgment: total anchor length > replacement region length
        if total_anchor_length > replace_length:
            ratio = total_anchor_length / replace_length
            print(f"Aggressive mode: Total anchor length({total_anchor_length:,}bp) > Replacement region({replace_length:,}bp), ratio: {ratio:.2f}, recommended repair")
            return True
        else:
            ratio = total_anchor_length / replace_length if replace_length > 0 else 0
            print(f"Aggressive mode: Total anchor length({total_anchor_length:,}bp) ≤ Replacement region({replace_length:,}bp), ratio: {ratio:.2f}, no repair")
            return False
    
    def _should_repair_in_conservative_mode(self, analysis: _GapAnalysis) -> bool:
        """
        Determine if should repair in conservative mode
        Condition: gap in error region and no large distance anchor warnings
        """
        # Basic condition: gap must be in error region
        if not analysis.gap_in_error_region:
            return False
        
        # Conservative mode: no repair if large distance anchor warning
        if analysis.large_distance_anchor:
            print(f"Conservative mode: Large distance anchor warning at {analysis.gap_pos:,}, no repair")
            return False
        
        # Need certain confidence level
        if analysis.confidence == "low":
            print(f"Conservative mode: Confidence too low ({analysis.confidence}), no repair")
            return False
        
        return True
    
    def filter_analyses_by_repair_mode(self, gap_analyses: List[_GapAnalysis],
                                     repair_mode: str = "conservative") -> List[_GapAnalysis]:
        """
        Filter analysis results based on repair mode
        """
        filtered_analyses = []
        
        for analysis in gap_analyses:
            should_repair = False
            repair_reason = ""
            
            if repair_mode == "conservative":
                if self._should_repair_in_conservative_mode(analysis):
                    should_repair = True
                    repair_reason = "Conservative mode conditions met"
                else:
                    repair_reason = "Conservative mode conditions not met"
            
            elif repair_mode == "aggressive":
                # Aggressive mode: first check if conservative mode conditions met
                conservative_ok = self._should_repair_in_conservative_mode(analysis)
                
                if conservative_ok:
                    # Conservative mode conditions met, repair directly
                    should_repair = True
                    repair_reason = "Conservative mode conditions met"
                else:
                    # Conservative mode not met, use aggressive mode rules
                    aggressive_ok = self._should_repair_in_aggressive_mode(analysis)
                    if aggressive_ok:
                        should_repair = True
                        repair_reason = "Aggressive mode sufficient anchors"
                    else:
                        repair_reason = "Aggressive mode insufficient anchors"
            
            if should_repair:
                # Mark repair reason
                analysis.repair_reason = repair_reason
                filtered_analyses.append(analysis)
                print(f"Gap {analysis.gap_pos:,}: Recommended repair ({repair_reason})")
            else:
                print(f"Gap {analysis.gap_pos:,}: Not recommended for repair ({repair_reason})")
        
        return filtered_analyses
    
    def determine_error_type(self, supporting_refs: Dict, crossing_aligns: List[_Alignment]) -> Tuple[str, str]:
        """Determine error type and confidence"""
        if crossing_aligns:
            return "type1", "high"
        
        if not supporting_refs:
            return None, "low"
        
        type_counts = defaultdict(int)
        for ref_contig, (left_align, right_align) in supporting_refs.items():
            if left_align.ref_contig != right_align.ref_contig:
                type_counts["type3"] += 1
            elif left_align.direction != right_align.direction:
                type_counts["type2"] += 1
            elif self.check_overlap(left_align, right_align):
                type_counts["type4"] += 1
            else:
                type_counts["type1"] += 1
        
        if not type_counts:
            return None, "low"
        
        consensus_type = max(type_counts.items(), key=lambda x: x[1])[0]
        
        total_refs = len(supporting_refs)
        if total_refs >= 3 and type_counts[consensus_type] >= 2:
            confidence = "high"
        elif total_refs >= 2:
            confidence = "medium"
        else:
            confidence = "low"
        
        return consensus_type, confidence
    
    def check_overlap(self, left_align: _Alignment, right_align: _Alignment, threshold: float = 0.1) -> bool:
        """Check if two alignments overlap on reference genome"""
        if left_align.ref_contig != right_align.ref_contig:
            return False
        
        overlap_start = max(left_align.ref_min, right_align.ref_min)
        overlap_end = min(left_align.ref_max, right_align.ref_max)
        
        if overlap_start <= overlap_end:
            overlap_size = overlap_end - overlap_start + 1
            min_length = min(left_align.length, right_align.length)
            return overlap_size / min_length > threshold
        
        return False
    
    def determine_replace_region(self, gap_pos: int, supporting_refs: Dict, 
                               error_type: str, crossing_aligns: List[_Alignment]) -> Tuple[int, int]:
        """Determine replacement region"""
        replace_start = gap_pos
        replace_end = gap_pos
        
        if error_type == "type1":
            if crossing_aligns:
                min_pos = min(a.query_min for a in crossing_aligns)
                max_pos = max(a.query_max for a in crossing_aligns)
                replace_start = min_pos
                replace_end = max_pos
            elif supporting_refs:
                min_left = float('inf')
                max_right = -float('inf')
                
                for left_align, right_align in supporting_refs.values():
                    if left_align.query_max < gap_pos:
                        min_left = min(min_left, left_align.query_max)
                    if right_align.query_min > gap_pos:
                        max_right = max(max_right, right_align.query_min)
                
                if min_left != float('inf') and max_right != -float('inf'):
                    replace_start = min_left + 1
                    replace_end = max_right - 1
        
        elif error_type == "type2":
            for left_align, right_align in supporting_refs.values():
                if left_align.direction != right_align.direction:
                    if left_align.direction == -1:
                        replace_start = min(replace_start, left_align.query_min)
                        replace_end = max(replace_end, left_align.query_max)
                    elif right_align.direction == -1:
                        replace_start = min(replace_start, right_align.query_min)
                        replace_end = max(replace_end, right_align.query_max)
        
        elif error_type == "type3":
            for left_align, right_align in supporting_refs.values():
                if left_align.ref_contig != right_align.ref_contig:
                    left_len = left_align.length
                    right_len = right_align.length
                    
                    if left_len < right_len:
                        replace_start = min(replace_start, left_align.query_min)
                        replace_end = max(replace_end, left_align.query_max)
                    else:
                        replace_start = min(replace_start, right_align.query_min)
                        replace_end = max(replace_end, right_align.query_max)
        
        elif error_type == "type4":
            for left_align, right_align in supporting_refs.values():
                if self.check_overlap(left_align, right_align):
                    left_len = left_align.length
                    right_len = right_align.length
                    
                    if left_len < right_len:
                        replace_start = min(replace_start, left_align.query_min)
                        replace_end = max(replace_end, left_align.query_max)
                    else:
                        replace_start = min(replace_start, right_align.query_min)
                        replace_end = max(replace_end, right_align.query_max)
        
        replace_start = min(replace_start, gap_pos)
        replace_end = max(replace_end, gap_pos)
        
        return replace_start, replace_end
    
    def analyze_all_gaps(self, gap_positions: List[int], max_search_distance: int = 500000) -> List[_GapAnalysis]:
        """Analyze all gap positions"""
        print("=" * 60)
        print(f"Starting analysis of {len(gap_positions)} gap positions")
        print(f"Maximum search distance: {max_search_distance:,} bp")
        print("=" * 60)
        
        all_analyses = []
        gaps_in_error = 0
        gaps_not_in_error = 0
        gaps_with_large_distance = 0
        
        for gap_pos in gap_positions:
            analysis = self.analyze_gap(gap_pos, max_search_distance)
            if analysis:
                all_analyses.append(analysis)
                if analysis.gap_in_error_region:
                    gaps_in_error += 1
                    if analysis.large_distance_anchor:
                        gaps_with_large_distance += 1
                else:
                    gaps_not_in_error += 1
        
        print(f"\nSuccessfully analyzed {len(all_analyses)}/{len(gap_positions)} gaps")
        print(f"  - {gaps_in_error} gaps in error regions")
        print(f"    * {gaps_with_large_distance} with large distance anchor warnings")
        print(f"  - {gaps_not_in_error} gaps not in error regions")
        return all_analyses
    
    def apply_unified_repairs(self, query_fasta: str, gap_analyses: List[_GapAnalysis], 
                            output_fasta: str, final_gap_length: int = 100,
                            repair_mode: str = "conservative") -> Tuple[str, List[Dict], Dict]:
        """Apply all repairs uniformly, filtered by repair_mode"""
        print(f"\n" + "=" * 60)
        print(f"Applying repairs (mode: {repair_mode})")
        print("=" * 60)
        
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython not installed, cannot perform sequence repair")
        
        # Filter analysis results based on repair mode
        print(f"Filtering gap analysis results based on {repair_mode} mode...")
        analyses_to_repair = self.filter_analyses_by_repair_mode(gap_analyses, repair_mode)
        
        print(f"Total analyses: {len(gap_analyses)}")
        print(f"Recommended repairs: {len(analyses_to_repair)}")
        
        # Statistics
        repair_stats = {
            'recommended_repairs': len(analyses_to_repair),
            'skipped_conservative': 0,
            'aggressive_repairs': 0,
            'skipped_aggressive': 0
        }
        
        if repair_mode == "conservative":
            large_distance_count = sum(1 for a in gap_analyses if a.large_distance_anchor and a.gap_in_error_region)
            repair_stats['skipped_conservative'] = large_distance_count
            print(f"Conservative mode skipped (large distance anchors): {large_distance_count}")
        elif repair_mode == "aggressive":
            # Count additional repairs in aggressive mode
            aggressive_repairs = []
            for analysis in analyses_to_repair:
                if not self._should_repair_in_conservative_mode(analysis):
                    aggressive_repairs.append(analysis)
            
            repair_stats['aggressive_repairs'] = len(aggressive_repairs)
            
            # Count skipped in aggressive mode (insufficient anchors)
            aggressive_skipped = 0
            for analysis in gap_analyses:
                if (analysis.gap_in_error_region and 
                    not self._should_repair_in_conservative_mode(analysis) and
                    not self._should_repair_in_aggressive_mode(analysis)):
                    aggressive_skipped += 1
            
            repair_stats['skipped_aggressive'] = aggressive_skipped
            print(f"Aggressive mode additional repairs: {len(aggressive_repairs)}")
            print(f"Aggressive mode skipped (insufficient anchors): {aggressive_skipped}")
        
        records = list(SeqIO.parse(query_fasta, "fasta"))
        if len(records) != 1:
            print(f"Warning: Expected single contig, but found {len(records)}")
        
        record = records[0]
        seq = str(record.seq)
        seq_length = len(seq)
        
        # Collect all regions to replace (only using recommended repair analyses)
        replace_regions = []
        for analysis in analyses_to_repair:
            if (analysis.query_contig == record.id and
                analysis.replace_start > 0 and analysis.replace_end > 0):
                
                start = max(1, analysis.replace_start)
                end = min(seq_length, analysis.replace_end)
                if start <= end:
                    replace_regions.append({
                        'start': start,
                        'end': end,
                        'gap_pos': analysis.gap_pos,
                        'error_type': analysis.error_type,
                        'large_distance': analysis.large_distance_anchor,
                        'confidence': analysis.confidence,
                        'total_anchor_length': analysis.total_anchor_length,
                        'replace_length': analysis.replace_length,
                        'repair_reason': analysis.repair_reason
                    })
        
        if not replace_regions:
            print("No regions need repair")
            return query_fasta, [], repair_stats
        
        print(f"Found {len(replace_regions)} regions to replace")
        
        # Sort by start position
        replace_regions.sort(key=lambda x: x['start'])
        
        # Merge overlapping replacement regions
        merged_regions = []
        if replace_regions:
            current = replace_regions[0].copy()
            
            for region in replace_regions[1:]:
                if region['start'] <= current['end'] + 1:
                    # Merge regions
                    current['end'] = max(current['end'], region['end'])
                    current['error_type'] = f"{current['error_type']}+{region['error_type']}"
                    current['large_distance'] = current['large_distance'] or region['large_distance']
                    current['confidence'] = min(current['confidence'], region['confidence'], 
                                              key=lambda x: ['low', 'medium', 'high'].index(x))
                    current['total_anchor_length'] += region['total_anchor_length']
                else:
                    merged_regions.append(current)
                    current = region.copy()
            
            merged_regions.append(current)
        
        print(f"Merged into {len(merged_regions)} regions")
        
        # Apply replacements from end to beginning
        seq_list = list(seq)
        repair_log = []
        
        for region in reversed(merged_regions):
            start, end = region['start'], region['end']
            original_length = end - start + 1
            replacement = 'N' * original_length
            
            repair_log.append({
                'gap_position': region['gap_pos'],
                'replace_region': f"{start:,}-{end:,}",
                'original_length': original_length,
                'error_type': region['error_type'],
                'confidence': region['confidence'],
                'large_distance_anchor': region['large_distance'],
                'total_anchor_length': region['total_anchor_length'],
                'replace_length': region['replace_length'],
                'repair_mode': repair_mode,
                'repair_reason': region['repair_reason'],
                'anchor_to_replace_ratio': region['total_anchor_length'] / region['replace_length'] if region['replace_length'] > 0 else 0
            })
            
            seq_list[start-1:end] = list(replacement)
            print(f"Replace region: {start:,}-{end:,} ({original_length:,}bp) - {region['error_type']} - {region['repair_reason']}")
        
        repaired_seq = ''.join(seq_list)
        
        # Adjust all consecutive N regions to unified length
        final_seq = self.normalize_gaps(repaired_seq, final_gap_length)
        
        # Save result
        new_record = SeqRecord(Seq(final_seq), id=record.id,
                              description=f"repaired_{repair_mode}_{len(repair_log)}_gaps")
        SeqIO.write([new_record], output_fasta, "fasta")
        
        # Print repair statistics
        total_replaced = sum(r['original_length'] for r in repair_log)
        total_anchor = sum(r.get('total_anchor_length', 0) for r in repair_log)
        print(f"\nRepair statistics:")
        print(f"  - Total replacement regions: {len(repair_log)}")
        print(f"  - Total replacement length: {total_replaced:,} bp")
        print(f"  - Total anchor length: {total_anchor:,} bp")
        print(f"  - Anchor/replacement ratio: {total_anchor/total_replaced:.2f}" if total_replaced > 0 else "  - Anchor/replacement ratio: N/A")
        print(f"  - Repair mode: {repair_mode}")
        print(f"  - Final gap length: {final_gap_length} bp")
        
        return output_fasta, repair_log, repair_stats
    
    def normalize_gaps(self, sequence: str, gap_length: int = 100) -> str:
        """Adjust all consecutive N regions to specified length"""
        result = []
        i = 0
        n = len(sequence)
        
        while i < n:
            if sequence[i] == 'N':
                start = i
                while i < n and sequence[i] == 'N':
                    i += 1
                result.append('N' * gap_length)
            else:
                result.append(sequence[i])
                i += 1
        
        return ''.join(result)
    
    def save_report(self, gap_analyses: List[_GapAnalysis], repair_log: List[Dict], 
                   output_json: str, query_fasta: str = None, repair_mode: str = "conservative") -> Dict[str, Any]:
        """Save analysis report to JSON file"""
        # Collect repair statistics
        repair_stats = {
            'recommended_repairs': len(repair_log),
            'skipped_conservative': 0,
            'aggressive_repairs': 0,
            'skipped_aggressive': 0
        }
        
        # Calculate repair statistics
        if repair_mode == "conservative":
            large_distance_count = sum(1 for a in gap_analyses if a.large_distance_anchor and a.gap_in_error_region)
            repair_stats['skipped_conservative'] = large_distance_count
        elif repair_mode == "aggressive":
            aggressive_repairs = 0
            aggressive_skipped = 0
            
            for analysis in gap_analyses:
                if analysis.gap_in_error_region:
                    conservative_ok = self._should_repair_in_conservative_mode(analysis)
                    if not conservative_ok:
                        aggressive_ok = self._should_repair_in_aggressive_mode(analysis)
                        if aggressive_ok:
                            aggressive_repairs += 1
                        else:
                            aggressive_skipped += 1
            
            repair_stats['aggressive_repairs'] = aggressive_repairs
            repair_stats['skipped_aggressive'] = aggressive_skipped
        
        report = {
            'summary': {
                'total_gaps_analyzed': len(gap_analyses),
                'gaps_in_error_region': sum(1 for a in gap_analyses if a.gap_in_error_region),
                'gaps_not_in_error_region': sum(1 for a in gap_analyses if not a.gap_in_error_region),
                'gaps_with_large_distance_anchor': sum(1 for a in gap_analyses if a.large_distance_anchor),
                'gaps_repaired': len(repair_log),
                'total_bp_replaced': sum(r.get('original_length', 0) for r in repair_log),
                'total_anchor_bp': sum(r.get('total_anchor_length', 0) for r in repair_log),
                'final_gap_length': 100,
                'repair_mode': repair_mode,
                'repair_summary': repair_stats,
                'error_type_counts': defaultdict(int),
                'error_type_counts_in_error_region': defaultdict(int)
            },
            'gap_analyses': [],
            'repair_details': repair_log,
            'synteny_summary': defaultdict(int)
        }
        
        # Add sequence information
        if query_fasta and os.path.exists(query_fasta):
            records = list(SeqIO.parse(query_fasta, "fasta"))
            if records:
                report['summary']['query_contig'] = records[0].id
                report['summary']['query_length'] = len(records[0].seq)
        
        for analysis in gap_analyses:
            gap_report = {
                'gap_position': analysis.gap_pos,
                'error_type': analysis.error_type,
                'confidence': analysis.confidence,
                'gap_in_error_region': analysis.gap_in_error_region,
                'large_distance_anchor': analysis.large_distance_anchor,
                'total_anchor_length': analysis.total_anchor_length,
                'replace_region': f"{analysis.replace_start}-{analysis.replace_end}" if analysis.replace_start > 0 else "None",
                'replace_length': analysis.replace_length if analysis.replace_start > 0 else 0,
                'supporting_ref_contigs': len(analysis.supporting_ref_contigs),
                'repair_recommended': analysis.repair_reason is not None,
                'repair_reason': analysis.repair_reason,
                'evidence_details': analysis.evidence_details
            }
            
            if analysis.synteny_analysis:
                gap_report['synteny_analysis'] = analysis.synteny_analysis
                for ref_contig, result in analysis.synteny_analysis.items():
                    if isinstance(result, dict):
                        report['synteny_summary'][result.get('assessment', 'unknown')] += 1
            
            report['gap_analyses'].append(gap_report)
            report['summary']['error_type_counts'][analysis.error_type] += 1
            
            if analysis.gap_in_error_region:
                report['summary']['error_type_counts_in_error_region'][analysis.error_type] += 1
        
        with open(output_json, 'w') as f:
            json.dump(report, f, indent=2, default=lambda x: dict(x) if hasattr(x, '__dict__') else str(x))
        
        return report

### Command Line Interface ###
def main():
    """Command line interface main function"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Genome Gap Analysis and Repair Tool based on Synteny (Enhanced: Enhanced Repair Mode)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage Examples:
  # Basic analysis, using default search distance(500kb)
  python gap_analyzer.py -coords scaffold_5_vs_ref.coords --gaps 1804444 12272394 14567005 15036440
  
  # Conservative mode analysis and repair
  python gap_analyzer.py -coords scaffold_5_vs_ref.coords --gaps 1804444 12272394 --query scaffold_5.fasta --repair-mode conservative
  
  # Aggressive mode analysis and repair
  python gap_analyzer.py -coords scaffold_5_vs_ref.coords --gaps 1804444 12272394 --query scaffold_5.fasta --repair-mode aggressive --verbose
  
  # Python module call
  from gap_analyzer import GapAnalyzerAPI
  analyzer = GapAnalyzerAPI(verbose=True)
  result = analyzer.analyze_gaps(
      coords_file="scaffold_5_vs_ref.coords",
      gap_positions=[1804444, 12272394],
      query_fasta="scaffold_5.fasta",
      repair_mode="aggressive"
  )
        """
    )
    
    # Required parameters
    parser.add_argument('-coords', required=True, help='Alignment coordinates file (coords format)')
    parser.add_argument('--gaps', nargs='+', type=int, required=True, help='List of gap positions')
    
    # Optional parameters
    parser.add_argument('--query', help='Query genome FASTA file (for repair)')
    parser.add_argument('--output-dir', default='.', help='Output directory, default current directory')
    parser.add_argument('--prefix', default='gap_analysis', help='Output file prefix, default gap_analysis')
    parser.add_argument('--json', default='gap_analysis_report.json', help='JSON report file')
    parser.add_argument('--gap-length', type=int, default=100, help='Final gap length (default 100bp)')
    parser.add_argument('--search', type=int, default=500000, help='Maximum search distance (bp), default 500kb')
    parser.add_argument('--step', type=int, default=100000, help='Search step size (bp), default 100kb')
    parser.add_argument('--min-confidence', choices=['low', 'medium', 'high'], default='low', 
                       help='Minimum confidence threshold')
    parser.add_argument('--repair-mode', choices=['conservative', 'aggressive'], default='conservative',
                       help='Repair mode: conservative(no repair with large distance anchor warnings) | aggressive(repair when anchors long enough)')
    parser.add_argument('--api-mode', action='store_true', help='Use API mode')
    parser.add_argument('--verbose', action='store_true', help='Show detailed output')
    parser.add_argument('--log-file', help='Log file path')
    
    args = parser.parse_args()
    
    # Check files
    if not os.path.exists(args.coords):
        print(f"Error: coords file does not exist: {args.coords}", file=sys.stderr)
        sys.exit(1)
    
    if args.query and not os.path.exists(args.query):
        print(f"Error: query file does not exist: {args.query}", file=sys.stderr)
        sys.exit(1)
    
    try:
        if args.api_mode:
            # Use API mode
            analyzer = GapAnalyzerAPI(verbose=args.verbose, log_file=args.log_file)
            
            result = analyzer.analyze_gaps(
                coords_file=args.coords,
                gap_positions=args.gaps,
                query_fasta=args.query,
                output_dir=args.output_dir,
                output_prefix=args.prefix,
                max_search_distance=args.search,
                search_step=args.step,
                final_gap_length=args.gap_length,
                min_confidence=args.min_confidence,
                repair_mode=args.repair_mode
            )
            
            print(f"\nDetailed report saved to: {result['api_info']['output_files']['report']}")
            if result['api_info']['output_files']['filled_fasta']:
                print(f"Repaired genome: {result['api_info']['output_files']['filled_fasta']}")
            
        else:
            # Use original mode
            analyzer = _GapAnalyzer()
            
            # Modify step size parameter in find_surrounding_alignments function
            original_method = analyzer.find_surrounding_alignments
            def new_find_surrounding_alignments(query_contig, gap_pos, max_distance):
                return original_method(query_contig, gap_pos, max_distance, args.step)
            analyzer.find_surrounding_alignments = new_find_surrounding_alignments
            
            # Parse coords file
            print("Parsing alignment file and building synteny blocks...")
            alignments = analyzer.parse_coords_file(args.coords)
            
            if not alignments:
                print("Error: No alignment information parsed", file=sys.stderr)
                sys.exit(1)
            
            # Analyze all gaps
            print("Analyzing all gaps (including synteny analysis)...")
            gap_analyses = analyzer.analyze_all_gaps(args.gaps, args.search)
            
            if not gap_analyses:
                print("Error: No gaps successfully analyzed", file=sys.stderr)
                sys.exit(1)
            
            # Apply repairs
            repair_log = []
            repair_stats = {}
            if args.query:
                print(f"Applying repairs (mode: {args.repair_mode})...")
                output_fasta = os.path.join(args.output_dir, f"{args.prefix}_repaired.fasta")
                output_fasta, repair_log, repair_stats = analyzer.apply_unified_repairs(
                    args.query, gap_analyses, output_fasta, args.gap_length, args.repair_mode
                )
                print(f"Repaired genome saved: {output_fasta}")
                print(f"Performed {len(repair_log)} replacements ({args.repair_mode} mode)")
                if repair_stats:
                    print(f"Repair statistics: {repair_stats}")
            
            # Save report
            report_file = os.path.join(args.output_dir, args.json)
            report = analyzer.save_report(gap_analyses, repair_log, report_file, args.query, args.repair_mode)
            print(f"Report saved: {report_file}")
            
            # Print summary
            print("\n" + "="*60)
            print("Analysis Summary")
            print("="*60)
            print(f"Successfully analyzed: {len(gap_analyses)}/{len(args.gaps)} gaps")
            print(f"Gaps in error regions: {report['summary']['gaps_in_error_region']}")
            print(f"Gaps not in error regions: {report['summary']['gaps_not_in_error_region']}")
            print(f"Repair mode: {args.repair_mode}")
            
            if 'gaps_with_large_distance_anchor' in report['summary']:
                print(f"Large distance anchor warnings: {report['summary']['gaps_with_large_distance_anchor']}")
            
            if repair_log:
                total_replaced = sum(r['original_length'] for r in repair_log)
                total_anchor = sum(r.get('total_anchor_length', 0) for r in repair_log)
                print(f"\nTotal replacement length: {total_replaced:,} bp")
                print(f"Total anchor length: {total_anchor:,} bp")
                print(f"Anchor/replacement ratio: {total_anchor/total_replaced:.2f}" if total_replaced > 0 else "Anchor/replacement ratio: N/A")
                print(f"All gaps finally unified to: {args.gap_length} bp")
            
            # Show repair statistics
            if 'repair_summary' in report['summary']:
                repair_summary = report['summary']['repair_summary']
                print(f"\nRepair decision statistics ({args.repair_mode} mode):")
                print(f"  Recommended repairs: {repair_summary.get('recommended_repairs', 0)}")
                if args.repair_mode == "conservative":
                    print(f"  Skipped (large distance anchors): {repair_summary.get('skipped_conservative', 0)}")
                elif args.repair_mode == "aggressive":
                    print(f"  Additional repairs in aggressive mode: {repair_summary.get('aggressive_repairs', 0)}")
                    print(f"  Skipped in aggressive mode (insufficient anchors): {repair_summary.get('skipped_aggressive', 0)}")
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()