#!/usr/bin/env python3
"""
Genome Gap Analysis and Repair Tool based on Synteny - Enhanced Version
Modified: Type5 now based on REFERENCE overlap with EXTRA MARGIN for gap filling
"""

import sys
import os
import json
import logging
import time
from collections import defaultdict, deque
from dataclasses import dataclass, asdict, field
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
        repair_mode: str = "conservative",
        type5_replace_side: str = "right",
        type5_extra_margin: int = 100
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
            type5_replace_side: For Type5 errors, which side to replace ("left" or "right")
            type5_extra_margin: Extra bp to delete for Type5 to facilitate gap filling (default 100)
            
        Returns:
            Dictionary containing analysis results
        """
        self.log(f"Starting gap analysis, repair mode: {repair_mode}", "info")
        self.log(f"Alignment file: {coords_file}", "info")
        self.log(f"Analyzing {len(gap_positions)} gap positions", "info")
        self.log(f"Type5 replace side: {type5_replace_side}", "info")
        self.log(f"Type5 extra margin: {type5_extra_margin} bp", "info")
        
        start_time = time.time()
        
        os.makedirs(output_dir, exist_ok=True)
        
        analyzer = _GapAnalyzer(type5_replace_side=type5_replace_side, type5_extra_margin=type5_extra_margin)
        
        original_method = analyzer.find_surrounding_alignments
        def new_find_surrounding_alignments(query_contig, gap_pos, max_distance):
            return original_method(query_contig, gap_pos, max_distance, search_step)
        analyzer.find_surrounding_alignments = new_find_surrounding_alignments
        
        self.log("Parsing alignment file and building synteny blocks...", "info")
        alignments = analyzer.parse_coords_file(coords_file)
        
        if not alignments:
            raise ValueError("No alignment information parsed")
        
        self.log(f"Analyzing all gaps (including synteny analysis)...", "info")
        gap_analyses = analyzer.analyze_all_gaps(gap_positions, max_search_distance)
        
        if not gap_analyses:
            raise ValueError("No gaps successfully analyzed")
        
        repair_log = []
        repair_stats = {}
        filled_fasta = None
        
        if query_fasta and os.path.exists(query_fasta):
            self.log(f"Applying repairs (mode: {repair_mode})...", "info")
            
            output_fasta = os.path.join(output_dir, f"{output_prefix}_repaired.fasta")
            
            try:
                filled_fasta, repair_log, repair_stats = analyzer.apply_unified_repairs(
                    query_fasta, gap_analyses, output_fasta, final_gap_length, repair_mode
                )
                
                self.log(f"Repaired genome saved: {filled_fasta}", "info")
                self.log(f"Performed {len(repair_log)} replacements", "info")
                
                if repair_stats:
                    self.log(f"Repair statistics: {repair_stats}", "info")
                    
            except Exception as e:
                self.log(f"Error during repair process: {e}", "error")
                raise
        
        report_file = os.path.join(output_dir, f"{output_prefix}_report.json")
        report = analyzer.save_report(gap_analyses, repair_log, report_file, query_fasta, repair_mode)
        
        if repair_stats:
            report.setdefault('repair_stats', {}).update(repair_stats)
        
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
                "repair_mode": repair_mode,
                "type5_replace_side": type5_replace_side,
                "type5_extra_margin": type5_extra_margin
            },
            "output_files": {
                "report": report_file,
                "filled_fasta": filled_fasta
            }
        }
        
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
        type5_replace_side = report.get("api_info", {}).get("parameters", {}).get("type5_replace_side", "right")
        type5_extra_margin = report.get("api_info", {}).get("parameters", {}).get("type5_extra_margin", 100)
        self.log(f"Repair mode: {repair_mode}", "info")
        self.log(f"Type5 replace side: {type5_replace_side}", "info")
        self.log(f"Type5 extra margin: {type5_extra_margin} bp", "info")
        
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
        
        self.log("\nError type distribution (main type):", "info")
        if "error_type_counts" in summary:
            for err_type, count in summary["error_type_counts"].items():
                self.log(f"  {err_type}: {count}", "info")
        
        self.log("\nError subtype distribution:", "info")
        if "error_subtype_counts" in summary:
            for subtype, count in summary["error_subtype_counts"].items():
                self.log(f"  {subtype}: {count}", "info")
        
        self.log("\nFeature distribution:", "info")
        if "feature_counts" in summary:
            for feature, count in list(summary["feature_counts"].items())[:10]:
                self.log(f"  {feature}: {count}", "info")
        
        if "type5_details" in summary and summary["type5_details"]:
            self.log("\nType5 (Reference Overlap with extra margin) details:", "info")
            for detail in summary["type5_details"][:5]:
                self.log(f"  Gap {detail['gap_pos']}: ref overlap {detail['overlap_length']}bp, "
                        f"extra margin {detail.get('extra_margin', 100)}bp, "
                        f"total deleted {detail.get('total_deleted', 0)}bp, "
                        f"replaced {detail['replace_side']} side", "info")
        
        if "repair_details" in report and report["repair_details"]:
            total_deleted = sum(r.get('deleted_length', 0) for r in report["repair_details"] if r.get('action') == 'delete_with_margin')
            total_replaced = sum(r.get('original_length', 0) for r in report["repair_details"] if r.get('action') != 'delete_with_margin')
            self.log(f"\nTotal deleted (Type5 with margin): {total_deleted:,} bp", "info")
            self.log(f"Total replaced (other types): {total_replaced:,} bp", "info")
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
        
        summary = {
            "total": len(coords_files),
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

### Internal implementation classes ###

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
        dir_str = "fwd" if self.direction == 1 else "rev"
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
    """Enhanced Gap analysis result with detailed classification"""
    gap_pos: int
    query_contig: str
    supporting_ref_contigs: Dict[str, Tuple[Optional[_Alignment], Optional[_Alignment]]]
    error_type: str
    error_subtype: str
    error_features: List[str]
    confidence: str
    replace_start: int
    replace_end: int
    gap_in_error_region: bool = True
    synteny_analysis: Optional[Dict] = None
    large_distance_anchor: bool = False
    total_anchor_length: int = 0
    repair_reason: Optional[str] = None
    
    query_overlap_length: int = 0
    query_overlap_region: Tuple[int, int] = (0, 0)
    ref_overlap_length: int = 0
    ref_overlap_region: Tuple[int, int] = (0, 0)
    ref_overlap_in_query: Tuple[int, int] = (0, 0)
    direction_conflict: bool = False
    ref_contig_conflict: bool = False
    multiple_ref_support: int = 0
    anchor_distances: Dict[str, int] = field(default_factory=dict)
    synteny_block_info: Dict[str, Any] = field(default_factory=dict)
    type_confidence_score: float = 0.0
    type5_replace_side: str = "right"
    type5_extra_margin: int = 100
    
    @property
    def replace_length(self):
        return self.replace_end - self.replace_start + 1 if self.replace_start > 0 else 0

class _GapAnalyzer:
    """Enhanced Gap analyzer with comprehensive error type classification"""
    
    def __init__(self, type5_replace_side: str = "right", type5_extra_margin: int = 100):
        self.coords_file: Optional[str] = None
        self.alignments: List[_Alignment] = []
        self.query_contig_groups = defaultdict(list)
        self.ref_contig_groups = defaultdict(list)
        self.synteny_blocks: Dict[str, List[_SyntenyBlock]] = defaultdict(list)
        self.type5_replace_side = type5_replace_side
        self.type5_extra_margin = type5_extra_margin
        self.current_type5_metadata = {}
        self.seq_length = 0
        self.type5_details = []
    
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
                    except Exception as e:
                        continue
        
        self.alignments = alignments
        
        for align in alignments:
            self.query_contig_groups[align.query_contig].append(align)
            self.ref_contig_groups[align.ref_contig].append(align)
        
        for contig in self.query_contig_groups:
            self.query_contig_groups[contig].sort(key=lambda x: x.query_min)
        for contig in self.ref_contig_groups:
            self.ref_contig_groups[contig].sort(key=lambda x: x.ref_min)
        
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
            
            aligns_sorted = sorted(aligns, key=lambda x: x.ref_min)
            
            blocks = []
            current_block = None
            
            for align in aligns_sorted:
                if current_block is None:
                    current_block = _SyntenyBlock(
                        ref_contig=ref_contig,
                        start_pos=align.ref_min,
                        end_pos=align.ref_max,
                        direction=align.direction,
                        alignments=[align]
                    )
                else:
                    gap_size = align.ref_min - current_block.end_pos
                    if (align.direction == current_block.direction and 
                        gap_size < 100000 and 
                        gap_size > -10000):
                        current_block.end_pos = max(current_block.end_pos, align.ref_max)
                        current_block.alignments.append(align)
                    else:
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
                if align.query_min <= gap_pos <= align.query_max:
                    crossing_candidates.append(align)
                elif abs(align.query_min - gap_pos) <= window or abs(align.query_max - gap_pos) <= window:
                    if align.query_max < gap_pos:
                        left_candidates.append(align)
                    elif align.query_min > gap_pos:
                        right_candidates.append(align)
            
            if crossing_candidates or (left_candidates and right_candidates):
                return left_candidates, right_candidates, crossing_candidates, window
        
        return [], [], [], max_distance
    
    def analyze_synteny_for_ref_contig(self, ref_contig: str, gap_pos: int, 
                                     query_contig: str) -> Dict:
        """Analyze synteny on specific reference contig"""
        result = {
            'has_synteny': False,
            'blocks_near_gap': [],
            'anchor_distance': float('inf'),
            'assessment': 'no_data',
            'synteny_blocks': []
        }
        
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
        
        query_aligns = [a for a in self.alignments 
                       if a.query_contig == query_contig and a.ref_contig == ref_contig]
        
        if not query_aligns:
            result['assessment'] = 'no_alignment_for_query'
            return result
        
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
    
    def analyze_gap_features(self, gap_pos: int, left_align: _Alignment, 
                           right_align: _Alignment) -> Dict[str, Any]:
        """
        Detailed analysis of gap features
        """
        features = {
            'query_overlap': False,
            'query_overlap_length': 0,
            'query_overlap_region': (0, 0),
            'query_overlap_ratio_left': 0,
            'query_overlap_ratio_right': 0,
            'overlap_type': 'none',
            'ref_overlap': False,
            'ref_overlap_length': 0,
            'ref_overlap_region': (0, 0),
            'ref_overlap_ratio': 0.0,
            'ref_overlap_in_query': (0, 0),
            'direction_match': True,
            'ref_contig_match': True,
            'left_distance': 0,
            'right_distance': 0,
            'left_in_synteny': False,
            'right_in_synteny': False,
            'synteny_block_consistent': False,
            'left_anchor_quality': 0,
            'right_anchor_quality': 0,
            'total_gap_size': 0,
            'ref_contig': left_align.ref_contig if left_align else None,
            'left_anchor_length': 0,
            'right_anchor_length': 0
        }
        
        if not left_align or not right_align:
            return features
        
        q_overlap_start = max(left_align.query_min, right_align.query_min)
        q_overlap_end = min(left_align.query_max, right_align.query_max)
        if q_overlap_start <= q_overlap_end:
            features['query_overlap'] = True
            features['query_overlap_length'] = q_overlap_end - q_overlap_start + 1
            features['query_overlap_region'] = (q_overlap_start, q_overlap_end)
            
            left_len = left_align.length
            right_len = right_align.length
            features['query_overlap_ratio_left'] = features['query_overlap_length'] / left_len if left_len > 0 else 0
            features['query_overlap_ratio_right'] = features['query_overlap_length'] / right_len if right_len > 0 else 0
            
            if left_align.query_max < right_align.query_min:
                features['overlap_type'] = 'adjacent'
            elif left_align.query_max >= right_align.query_min:
                features['overlap_type'] = 'overlapping'
        
        if left_align.ref_contig == right_align.ref_contig:
            r_overlap_start = max(left_align.ref_min, right_align.ref_min)
            r_overlap_end = min(left_align.ref_max, right_align.ref_max)
            if r_overlap_start <= r_overlap_end:
                features['ref_overlap'] = True
                features['ref_overlap_length'] = r_overlap_end - r_overlap_start + 1
                features['ref_overlap_region'] = (r_overlap_start, r_overlap_end)
                
                min_ref_len = min(left_align.length, right_align.length)
                features['ref_overlap_ratio'] = features['ref_overlap_length'] / min_ref_len if min_ref_len > 0 else 0
                
                if left_align.direction == 1:
                    left_query_of_overlap_start = left_align.query_min + (r_overlap_start - left_align.ref_min)
                    left_query_of_overlap_end = left_align.query_min + (r_overlap_end - left_align.ref_min)
                else:
                    left_query_of_overlap_start = left_align.query_max - (r_overlap_end - left_align.ref_min)
                    left_query_of_overlap_end = left_align.query_max - (r_overlap_start - left_align.ref_min)
                
                if right_align.direction == 1:
                    right_query_of_overlap_start = right_align.query_min + (r_overlap_start - right_align.ref_min)
                    right_query_of_overlap_end = right_align.query_min + (r_overlap_end - right_align.ref_min)
                else:
                    right_query_of_overlap_start = right_align.query_max - (r_overlap_end - right_align.ref_min)
                    right_query_of_overlap_end = right_align.query_max - (r_overlap_start - right_align.ref_min)
                
                features['ref_overlap_in_query'] = (
                    min(left_query_of_overlap_start, left_query_of_overlap_end),
                    max(right_query_of_overlap_start, right_query_of_overlap_end)
                )
        
        features['direction_match'] = (left_align.direction == right_align.direction)
        
        features['ref_contig_match'] = (left_align.ref_contig == right_align.ref_contig)
        
        features['left_distance'] = gap_pos - left_align.query_max
        features['right_distance'] = right_align.query_min - gap_pos
        features['total_gap_size'] = features['left_distance'] + features['right_distance']
        
        features['left_anchor_quality'] = left_align.length * (left_align.identity / 100)
        features['right_anchor_quality'] = right_align.length * (right_align.identity / 100)
        features['left_anchor_length'] = left_align.length
        features['right_anchor_length'] = right_align.length
        
        if left_align.ref_contig in self.synteny_blocks:
            for block in self.synteny_blocks[left_align.ref_contig]:
                if block.start_pos <= left_align.ref_min <= block.end_pos:
                    features['left_in_synteny'] = True
                    features['left_synteny_block'] = {
                        'start': block.start_pos,
                        'end': block.end_pos,
                        'length': block.length,
                        'avg_identity': block.avg_identity
                    }
        
        if right_align.ref_contig in self.synteny_blocks:
            for block in self.synteny_blocks[right_align.ref_contig]:
                if block.start_pos <= right_align.ref_min <= block.end_pos:
                    features['right_in_synteny'] = True
                    features['right_synteny_block'] = {
                        'start': block.start_pos,
                        'end': block.end_pos,
                        'length': block.length,
                        'avg_identity': block.avg_identity
                    }
        
        if (features['left_in_synteny'] and features['right_in_synteny'] and
            features.get('left_synteny_block') == features.get('right_synteny_block')):
            features['synteny_block_consistent'] = True
        
        return features
    
    def check_query_overlap(self, left_align: _Alignment, right_align: _Alignment) -> Tuple[bool, int, Tuple[int, int]]:
        """Check if two alignments overlap on query sequence"""
        if left_align.query_contig != right_align.query_contig:
            return False, 0, (0, 0)
        
        overlap_start = max(left_align.query_min, right_align.query_min)
        overlap_end = min(left_align.query_max, right_align.query_max)
        
        if overlap_start <= overlap_end:
            overlap_length = overlap_end - overlap_start + 1
            return True, overlap_length, (overlap_start, overlap_end)
        
        return False, 0, (0, 0)
    
    def check_ref_overlap(self, left_align: _Alignment, right_align: _Alignment, threshold: float = 0.1) -> Tuple[bool, int, Tuple[int, int]]:
        """Check if two alignments overlap on reference genome"""
        if left_align.ref_contig != right_align.ref_contig:
            return False, 0, (0, 0)
        
        overlap_start = max(left_align.ref_min, right_align.ref_min)
        overlap_end = min(left_align.ref_max, right_align.ref_max)
        
        if overlap_start <= overlap_end:
            overlap_size = overlap_end - overlap_start + 1
            min_length = min(left_align.length, right_align.length)
            if overlap_size / min_length > threshold:
                return True, overlap_size, (overlap_start, overlap_end)
        
        return False, 0, (0, 0)
    
    def determine_error_type_detailed(self, supporting_refs: Dict, 
                                     crossing_aligns: List[_Alignment],
                                     gap_pos: int) -> Tuple[str, str, List[str], float, Dict]:
        """
        Detailed error type analysis with MODIFIED Type5 (now based on reference overlap)
        Returns: (main_type, subtype, feature_list, confidence_score, details)
        """
        if crossing_aligns:
            return self._analyze_crossing_gap(crossing_aligns, gap_pos)
        
        if not supporting_refs:
            return "unknown", "no_evidence", ["insufficient_data"], 0.0, {}
        
        type5_detected = False
        type5_details = {}
        type5_overlap_length = 0
        
        for ref_contig, (left_align, right_align) in supporting_refs.items():
            if left_align and right_align:
                if left_align.ref_contig == right_align.ref_contig:
                    ref_overlap_start = max(left_align.ref_min, right_align.ref_min)
                    ref_overlap_end = min(left_align.ref_max, right_align.ref_max)
                    
                    if ref_overlap_start <= ref_overlap_end:
                        has_overlap = True
                        overlap_len = ref_overlap_end - ref_overlap_start + 1
                        overlap_region = (ref_overlap_start, ref_overlap_end)
                        
                        if has_overlap and overlap_len > 0:
                            type5_detected = True
                            type5_overlap_length = max(type5_overlap_length, overlap_len)
                            type5_details = {
                                'overlap_length': overlap_len,
                                'overlap_region': overlap_region,
                                'ref_contig': ref_contig,
                                'type': 'reference_overlap'
                            }
        
        if type5_detected:
            features = [f"ref_overlap_{type5_overlap_length}"]
            
            if type5_overlap_length < 10000:
                subtype = "small_ref_overlap"
            elif type5_overlap_length < 50000:
                subtype = "medium_ref_overlap"
            else:
                subtype = "large_ref_overlap"
            
            confidence_score = min(0.9 + (type5_overlap_length / 1000000), 0.99)
            
            return "type5", subtype, features, confidence_score, type5_details
        
        all_features = []
        all_feature_lists = []
        
        for ref_contig, (left_align, right_align) in supporting_refs.items():
            if left_align and right_align:
                features = self.analyze_gap_features(gap_pos, left_align, right_align)
                all_features.append(features)
                
                feature_list = []
                if features['query_overlap']:
                    feature_list.append(f"query_overlap_{features['query_overlap_length']}")
                if features['ref_overlap']:
                    feature_list.append(f"ref_overlap_{features['ref_overlap_length']}")
                if not features['direction_match']:
                    feature_list.append("direction_conflict")
                if not features['ref_contig_match']:
                    feature_list.append("ref_contig_conflict")
                if features['left_distance'] > 100000:
                    feature_list.append(f"large_left_gap_{features['left_distance']}")
                if features['right_distance'] > 100000:
                    feature_list.append(f"large_right_gap_{features['right_distance']}")
                if features['synteny_block_consistent']:
                    feature_list.append("synteny_consistent")
                else:
                    feature_list.append("synteny_inconsistent")
                
                if features['total_gap_size'] < 1000:
                    feature_list.append("small_gap")
                elif features['total_gap_size'] < 50000:
                    feature_list.append("medium_gap")
                else:
                    feature_list.append("large_gap")
                
                all_feature_lists.append(feature_list)
        
        main_type, type_score = self._determine_main_type(all_features, all_feature_lists)
        subtype = self._determine_subtype(all_features, all_feature_lists, main_type)
        
        all_features_combined = set()
        for fl in all_feature_lists:
            all_features_combined.update(fl)
        
        confidence_score = self._calculate_confidence_score(all_features, len(supporting_refs))
        
        details = {
            'num_supporting_refs': len(supporting_refs),
            'feature_summary': list(all_features_combined),
            'anchor_qualities': [f['left_anchor_quality'] + f['right_anchor_quality'] 
                                for f in all_features],
            'distances': [(f['left_distance'], f['right_distance']) for f in all_features],
            'type_score': type_score,
            'confidence_score': confidence_score
        }
        
        return main_type, subtype, list(all_features_combined), confidence_score, details
    
    def _analyze_crossing_gap(self, crossing_aligns: List[_Alignment], 
                            gap_pos: int) -> Tuple[str, str, List[str], float, Dict]:
        """Analyze gap that is within an alignment region"""
        features = ["crossing_alignment"]
        
        if len(crossing_aligns) > 1:
            features.append("multiple_crossings")
            
            ref_contigs = set(a.ref_contig for a in crossing_aligns)
            if len(ref_contigs) > 1:
                features.append("cross_contig_crossing")
        
        for align in crossing_aligns:
            if align.query_min <= gap_pos <= align.query_max:
                pos_in_align = gap_pos - align.query_min
                align_length = align.length
                if pos_in_align < 100:
                    features.append("gap_at_start")
                elif pos_in_align > align_length - 100:
                    features.append("gap_at_end")
                else:
                    features.append("gap_in_middle")
                break
        
        return "type1", "crossing_alignment", features, 0.9, {"crossing_aligns": len(crossing_aligns)}
    
    def _determine_main_type(self, all_features: List[Dict], 
                           all_feature_lists: List[List[str]]) -> Tuple[str, float]:
        """
        Determine main error type based on features
        MODIFIED: Type5 now based on reference overlap
        Returns: (main_type, confidence_score)
        """
        type_scores = defaultdict(float)
        
        for i, features in enumerate(all_features):
            feature_list = all_feature_lists[i]
            
            if features['ref_overlap']:
                overlap_ratio = features['ref_overlap_length'] / min(
                    features['left_anchor_length'], features['right_anchor_length']
                )
                type_scores['type5'] += 2.5 + overlap_ratio * 2.5
                
                if not features['direction_match']:
                    type_scores['type5'] += 1.0
            
            if features['ref_overlap'] and not features.get('type5_override', False):
                overlap_ratio = features['ref_overlap_length'] / min(
                    features['left_anchor_length'], features['right_anchor_length']
                )
                type_scores['type4'] += 1.0 + overlap_ratio * 1.0
            
            if not features['direction_match']:
                type_scores['type2'] += 2.0
            
            if not features['ref_contig_match']:
                type_scores['type3'] += 2.0
            
            conflict_count = 0
            if not features['direction_match']:
                conflict_count += 1
            if not features['ref_contig_match']:
                conflict_count += 1
            if features['query_overlap'] and features['ref_overlap']:
                conflict_count += 1
            if conflict_count >= 2:
                type_scores['type6'] += 2.5
            
            if (not features['query_overlap'] and not features['ref_overlap'] and
                features['direction_match'] and features['ref_contig_match']):
                gap_size = features['total_gap_size']
                if gap_size < 1000:
                    type_scores['type1'] += 1.0
                elif gap_size < 50000:
                    type_scores['type1'] += 0.7
                else:
                    type_scores['type1'] += 0.3
        
        if not type_scores:
            return "type1", 0.3
        
        max_score = max(type_scores.values())
        normalized_scores = {t: s/max_score for t, s in type_scores.items()}
        
        main_type = max(normalized_scores.items(), key=lambda x: x[1])[0]
        confidence_score = normalized_scores[main_type]
        
        return main_type, confidence_score
    
    def _determine_subtype(self, all_features: List[Dict], 
                         all_feature_lists: List[List[str]], 
                         main_type: str) -> str:
        """Determine subtype based on main type and features"""
        
        if main_type == "type1":
            avg_gap_size = sum(f['total_gap_size'] for f in all_features) / len(all_features)
            if avg_gap_size < 1000:
                return "small_gap"
            elif avg_gap_size < 50000:
                return "medium_gap"
            else:
                return "large_gap"
        
        elif main_type == "type2":
            if any(not f['direction_match'] for f in all_features):
                return "direction_conflict"
            else:
                return "orientation_issue"
        
        elif main_type == "type3":
            ref_contigs = set()
            for f in all_features:
                if f.get('ref_contig'):
                    ref_contigs.add(f['ref_contig'])
            if len(ref_contigs) > 2:
                return "multi_contig_translocation"
            else:
                return "simple_translocation"
        
        elif main_type == "type4":
            avg_overlap = sum(f['ref_overlap_length'] for f in all_features) / len(all_features)
            if avg_overlap < 10000:
                return "small_ref_overlap"
            elif avg_overlap < 100000:
                return "medium_ref_overlap"
            else:
                return "large_ref_overlap"
        
        elif main_type == "type5":
            avg_overlap = sum(f['ref_overlap_length'] for f in all_features) / len(all_features)
            if avg_overlap < 10000:
                return "small_ref_overlap"
            elif avg_overlap < 50000:
                return "medium_ref_overlap"
            else:
                return "large_ref_overlap"
        
        elif main_type == "type6":
            features_present = set()
            for fl in all_feature_lists:
                for f in fl:
                    if 'conflict' in f or 'overlap' in f:
                        features_present.add(f.split('_')[0] if '_' in f else f)
            return "+".join(sorted(list(features_present))[:3])
        
        return "unspecified"
    
    def _calculate_confidence_score(self, all_features: List[Dict], num_supporting_refs: int) -> float:
        """Calculate confidence score (0-1)"""
        if not all_features:
            return 0.0
        
        if num_supporting_refs >= 3:
            base_confidence = 0.8
        elif num_supporting_refs == 2:
            base_confidence = 0.6
        else:
            base_confidence = 0.4
        
        feature_consistency = 1.0
        if len(all_features) > 1:
            consistent_count = 0
            for i in range(len(all_features)-1):
                if (all_features[i]['direction_match'] == all_features[i+1]['direction_match'] and
                    all_features[i]['ref_contig_match'] == all_features[i+1]['ref_contig_match']):
                    consistent_count += 1
            feature_consistency = 0.5 + 0.5 * (consistent_count / (len(all_features) - 1))
        
        avg_quality = sum(f['left_anchor_quality'] + f['right_anchor_quality'] 
                         for f in all_features) / (len(all_features) * 2)
        quality_factor = min(avg_quality / 1000000, 1.0)
        
        avg_gap_size = sum(f['total_gap_size'] for f in all_features) / len(all_features)
        if avg_gap_size < 1000:
            gap_factor = 1.0
        elif avg_gap_size < 10000:
            gap_factor = 0.9
        elif avg_gap_size < 100000:
            gap_factor = 0.7
        else:
            gap_factor = 0.5
        
        synteny_factor = 1.0
        synteny_consistent = sum(1 for f in all_features if f['synteny_block_consistent'])
        if synteny_consistent > 0:
            synteny_factor = 0.8 + 0.2 * (synteny_consistent / len(all_features))
        
        final_confidence = (base_confidence * 0.4 + 
                          feature_consistency * 0.2 + 
                          quality_factor * 0.2 + 
                          gap_factor * 0.1 + 
                          synteny_factor * 0.1)
        
        return min(final_confidence, 1.0)
    
    def _score_to_confidence_level(self, score: float) -> str:
        """Convert confidence score to level"""
        if score >= 0.7:
            return "high"
        elif score >= 0.4:
            return "medium"
        else:
            return "low"
    
    def check_gap_in_error_region(self, gap_pos: int, supporting_refs: Dict, 
                                error_type: str, crossing_aligns: List[_Alignment]) -> bool:
        """Check if gap is in error region"""
        if crossing_aligns:
            return True
        
        if not supporting_refs:
            return False
        
        if error_type in ["type1", "type2", "type3", "type4", "type5", "type6"]:
            return True
        
        return False
    
    def determine_replace_region(self, gap_pos: int, supporting_refs: Dict, 
                               error_type: str, error_subtype: str,
                               crossing_aligns: List[_Alignment]) -> Tuple[int, int]:
        """
        Determine replacement region based on error type
        MODIFIED: Type5 now based on reference overlap with EXTRA MARGIN for gap filling
        """
        
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
                    if left_align and left_align.query_max < gap_pos:
                        min_left = min(min_left, left_align.query_max)
                    if right_align and right_align.query_min > gap_pos:
                        max_right = max(max_right, right_align.query_min)
                
                if min_left != float('inf') and max_right != -float('inf'):
                    replace_start = min_left + 1
                    replace_end = max_right - 1
        
        elif error_type == "type2":
            for left_align, right_align in supporting_refs.values():
                if left_align and right_align and left_align.direction != right_align.direction:
                    if left_align.direction == -1:
                        replace_start = min(replace_start, left_align.query_min)
                        replace_end = max(replace_end, left_align.query_max)
                    elif right_align.direction == -1:
                        replace_start = min(replace_start, right_align.query_min)
                        replace_end = max(replace_end, right_align.query_max)
        
        elif error_type == "type3":
            for left_align, right_align in supporting_refs.values():
                if left_align and right_align and left_align.ref_contig != right_align.ref_contig:
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
                if left_align and right_align and self.check_ref_overlap(left_align, right_align)[0]:
                    left_len = left_align.length
                    right_len = right_align.length
                    
                    if left_len < right_len:
                        replace_start = min(replace_start, left_align.query_min)
                        replace_end = max(replace_end, left_align.query_max)
                    else:
                        replace_start = min(replace_start, right_align.query_min)
                        replace_end = max(replace_end, right_align.query_max)
        
        elif error_type == "type5":
            for ref_contig, (left_align, right_align) in supporting_refs.items():
                if left_align and right_align and left_align.ref_contig == right_align.ref_contig:
                    ref_overlap_start = max(left_align.ref_min, right_align.ref_min)
                    ref_overlap_end = min(left_align.ref_max, right_align.ref_max)
                    
                    if ref_overlap_start <= ref_overlap_end:
                        ref_overlap_length = ref_overlap_end - ref_overlap_start + 1
                        
                        if self.type5_replace_side == "right":
                            if right_align.direction == 1:
                                ref_len = right_align.ref_max - right_align.ref_min
                                query_len = right_align.query_max - right_align.query_min
                                ratio = query_len / ref_len if ref_len > 0 else 1
                                
                                overlap_offset_in_right = ref_overlap_start - right_align.ref_min
                                query_overlap_start = right_align.query_min + int(overlap_offset_in_right * ratio)
                                query_overlap_end = right_align.query_min + int((ref_overlap_end - right_align.ref_min) * ratio)
                                
                                replace_start = query_overlap_start
                                replace_end = query_overlap_end + self.type5_extra_margin
                                
                                print(f"Type5 right side: overlap {query_overlap_start}-{query_overlap_end} ({ref_overlap_length}bp) + margin {self.type5_extra_margin}bp = {replace_start}-{replace_end}")
                                
                            else:
                                ref_len = right_align.ref_max - right_align.ref_min
                                query_len = right_align.query_max - right_align.query_min
                                ratio = query_len / ref_len if ref_len > 0 else 1
                                
                                overlap_offset_in_right = ref_overlap_start - right_align.ref_min
                                query_overlap_end = right_align.query_max - int(overlap_offset_in_right * ratio)
                                query_overlap_start = right_align.query_max - int((ref_overlap_end - right_align.ref_min) * ratio)
                                
                                replace_start = query_overlap_start
                                replace_end = query_overlap_end + self.type5_extra_margin
                                
                                print(f"Type5 right side (reverse): overlap {query_overlap_start}-{query_overlap_end} ({ref_overlap_length}bp) + margin {self.type5_extra_margin}bp = {replace_start}-{replace_end}")
                        
                        else:
                            if left_align.direction == 1:
                                ref_len = left_align.ref_max - left_align.ref_min
                                query_len = left_align.query_max - left_align.query_min
                                ratio = query_len / ref_len if ref_len > 0 else 1
                                
                                overlap_offset_in_left = ref_overlap_start - left_align.ref_min
                                query_overlap_start = left_align.query_min + int(overlap_offset_in_left * ratio)
                                query_overlap_end = left_align.query_min + int((ref_overlap_end - left_align.ref_min) * ratio)
                                
                                replace_start = query_overlap_start - self.type5_extra_margin
                                replace_end = query_overlap_end
                                
                                print(f"Type5 left side: overlap {query_overlap_start}-{query_overlap_end} ({ref_overlap_length}bp) + margin {self.type5_extra_margin}bp = {replace_start}-{replace_end}")
                                
                            else:
                                ref_len = left_align.ref_max - left_align.ref_min
                                query_len = left_align.query_max - left_align.query_min
                                ratio = query_len / ref_len if ref_len > 0 else 1
                                
                                overlap_offset_in_left = ref_overlap_start - left_align.ref_min
                                query_overlap_end = left_align.query_max - int(overlap_offset_in_left * ratio)
                                query_overlap_start = left_align.query_max - int((ref_overlap_end - left_align.ref_min) * ratio)
                                
                                replace_start = query_overlap_start - self.type5_extra_margin
                                replace_end = query_overlap_end
                                
                                print(f"Type5 left side (reverse): overlap {query_overlap_start}-{query_overlap_end} ({ref_overlap_length}bp) + margin {self.type5_extra_margin}bp = {replace_start}-{replace_end}")
                    break
        
        elif error_type == "type6":
            min_left = float('inf')
            max_right = -float('inf')
            for left_align, right_align in supporting_refs.values():
                if left_align and left_align.query_max < gap_pos:
                    min_left = min(min_left, left_align.query_max)
                if right_align and right_align.query_min > gap_pos:
                    max_right = max(max_right, right_align.query_min)
            if min_left != float('inf') and max_right != -float('inf'):
                replace_start = min_left + 1
                replace_end = max_right - 1
        
        replace_start = min(replace_start, gap_pos)
        replace_end = max(replace_end, gap_pos)
        
        replace_start = max(1, replace_start)
        
        return replace_start, replace_end
    
    def analyze_gap(self, gap_pos: int, max_search_distance: int = 500000) -> Optional[_GapAnalysis]:
        """Analyze single gap with enhanced classification"""
        if not self.query_contig_groups:
            return None
        
        query_contig = list(self.query_contig_groups.keys())[0]
        
        left_aligns, right_aligns, crossing_aligns, window_used = self.find_surrounding_alignments(
            query_contig, gap_pos, max_search_distance
        )
        
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
            return self._analyze_gap_with_synteny(gap_pos, query_contig, max_search_distance)
        
        final_supporting = {}
        evidence_details = defaultdict(list)
        large_distance_detected = False
        total_anchor_length = 0
        query_overlap_length = 0
        query_overlap_region = (0, 0)
        ref_overlap_length = 0
        ref_overlap_region = (0, 0)
        ref_overlap_in_query = (0, 0)
        
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
                
                has_q_overlap, q_overlap_len, q_overlap_region = self.check_query_overlap(left_best, right_best)
                if has_q_overlap:
                    query_overlap_length = q_overlap_len
                    query_overlap_region = q_overlap_region
                    evidence_details[ref_contig].append(f"Query overlap detected: {q_overlap_len:,} bp")
                
                if left_best.ref_contig == right_best.ref_contig:
                    ref_overlap_start = max(left_best.ref_min, right_best.ref_min)
                    ref_overlap_end = min(left_best.ref_max, right_best.ref_max)
                    if ref_overlap_start <= ref_overlap_end:
                        ref_overlap_length = ref_overlap_end - ref_overlap_start + 1
                        ref_overlap_region = (ref_overlap_start, ref_overlap_end)
                        evidence_details[ref_contig].append(f"Reference overlap detected: {ref_overlap_length:,} bp")
                        
                        features = self.analyze_gap_features(gap_pos, left_best, right_best)
                        ref_overlap_in_query = features.get('ref_overlap_in_query', (0, 0))
                
                anchor_length = left_best.length + right_best.length
                total_anchor_length += anchor_length
                evidence_details[ref_contig].append(f"Anchor length: {anchor_length:,} bp")
                
                if min_distance > 500000:
                    large_distance_detected = True
                    evidence_details[ref_contig].append(f"Warning: Large anchor distance ({min_distance:,} bp)")
        
        if not final_supporting:
            return None
        
        main_type, subtype, features, confidence_score, details = self.determine_error_type_detailed(
            final_supporting, crossing_aligns, gap_pos
        )
        
        confidence = self._score_to_confidence_level(confidence_score)
        
        gap_in_error_region = self.check_gap_in_error_region(
            gap_pos, final_supporting, main_type, crossing_aligns
        )
        
        if not gap_in_error_region:
            analysis = _GapAnalysis(
                gap_pos=gap_pos,
                query_contig=query_contig,
                supporting_ref_contigs=final_supporting,
                error_type=main_type,
                error_subtype=subtype,
                error_features=features,
                confidence=confidence,
                replace_start=0,
                replace_end=0,
                gap_in_error_region=False,
                large_distance_anchor=large_distance_detected,
                total_anchor_length=total_anchor_length,
                query_overlap_length=query_overlap_length,
                query_overlap_region=query_overlap_region,
                ref_overlap_length=ref_overlap_length,
                ref_overlap_region=ref_overlap_region,
                ref_overlap_in_query=ref_overlap_in_query,
                type_confidence_score=confidence_score,
                type5_replace_side=self.type5_replace_side,
                type5_extra_margin=self.type5_extra_margin,
                anchor_distances={'left': left_distance, 'right': right_distance}
            )
            return analysis
        
        replace_start, replace_end = self.determine_replace_region(
            gap_pos, final_supporting, main_type, subtype, crossing_aligns
        )
        
        analysis = _GapAnalysis(
            gap_pos=gap_pos,
            query_contig=query_contig,
            supporting_ref_contigs=final_supporting,
            error_type=main_type,
            error_subtype=subtype,
            error_features=features,
            confidence=confidence,
            replace_start=replace_start,
            replace_end=replace_end,
            gap_in_error_region=True,
            large_distance_anchor=large_distance_detected,
            total_anchor_length=total_anchor_length,
            query_overlap_length=query_overlap_length,
            query_overlap_region=query_overlap_region,
            ref_overlap_length=ref_overlap_length,
            ref_overlap_region=ref_overlap_region,
            ref_overlap_in_query=ref_overlap_in_query,
            type_confidence_score=confidence_score,
            type5_replace_side=self.type5_replace_side,
            type5_extra_margin=self.type5_extra_margin,
            anchor_distances={'left': left_distance, 'right': right_distance}
        )
        
        return analysis
    
    def _analyze_gap_with_synteny(self, gap_pos: int, query_contig: str, 
                                 max_search_distance: int) -> Optional[_GapAnalysis]:
        """Analyze gap using synteny information when direct alignments are insufficient"""
        synteny_results = {}
        for ref_contig in list(self.ref_contig_groups.keys())[:5]:
            result = self.analyze_synteny_for_ref_contig(ref_contig, gap_pos, query_contig)
            synteny_results[ref_contig] = result
        
        supporting_refs = defaultdict(lambda: {"left": [], "right": []})
        
        for ref_contig, result in synteny_results.items():
            if result['assessment'] in ['synteny_present', 'gap_in_synteny_block']:
                ref_aligns = [a for a in self.alignments 
                            if a.ref_contig == ref_contig and a.query_contig == query_contig]
                
                left_for_ref = [a for a in ref_aligns if a.query_max < gap_pos]
                right_for_ref = [a for a in ref_aligns if a.query_min > gap_pos]
                
                if left_for_ref:
                    supporting_refs[ref_contig]["left"].append(max(left_for_ref, key=lambda x: x.query_max))
                if right_for_ref:
                    supporting_refs[ref_contig]["right"].append(min(right_for_ref, key=lambda x: x.query_min))
        
        if not supporting_refs:
            return None
        
        return self.analyze_gap(gap_pos, max_search_distance)
    
    def _should_repair_in_aggressive_mode(self, analysis: _GapAnalysis) -> bool:
        """Determine if should repair in aggressive mode"""
        if not analysis.gap_in_error_region:
            return False
        
        if not analysis.supporting_ref_contigs:
            return False
        
        total_anchor_length = analysis.total_anchor_length
        replace_length = analysis.replace_length
        
        if replace_length <= 0:
            return False
        
        if total_anchor_length > replace_length:
            ratio = total_anchor_length / replace_length
            print(f"Aggressive mode: Anchor length({total_anchor_length:,}bp) > Replacement({replace_length:,}bp), ratio: {ratio:.2f}")
            return True
        else:
            ratio = total_anchor_length / replace_length if replace_length > 0 else 0
            print(f"Aggressive mode: Anchor length({total_anchor_length:,}bp) <= Replacement({replace_length:,}bp), ratio: {ratio:.2f}")
            return False
    
    def _should_repair_in_conservative_mode(self, analysis: _GapAnalysis) -> bool:
        """Determine if should repair in conservative mode"""
        if not analysis.gap_in_error_region:
            return False
        
        if analysis.large_distance_anchor:
            print(f"Conservative mode: Large distance anchor at {analysis.gap_pos:,}, skipping")
            return False
        
        if analysis.confidence == "low" and analysis.type_confidence_score < 0.5:
            print(f"Conservative mode: Low confidence ({analysis.confidence}, score={analysis.type_confidence_score:.2f})")
            return False
        
        return True
    
    def filter_analyses_by_repair_mode(self, gap_analyses: List[_GapAnalysis],
                                     repair_mode: str = "conservative") -> List[_GapAnalysis]:
        """Filter analysis results based on repair mode"""
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
                conservative_ok = self._should_repair_in_conservative_mode(analysis)
                
                if conservative_ok:
                    should_repair = True
                    repair_reason = "Conservative mode conditions met"
                else:
                    aggressive_ok = self._should_repair_in_aggressive_mode(analysis)
                    if aggressive_ok:
                        should_repair = True
                        repair_reason = "Aggressive mode sufficient anchors"
                    else:
                        repair_reason = "Aggressive mode insufficient anchors"
            
            if should_repair:
                analysis.repair_reason = repair_reason
                filtered_analyses.append(analysis)
                print(f"Gap {analysis.gap_pos:,}: Repair ({repair_reason}) [{analysis.error_type}/{analysis.error_subtype}]")
            else:
                print(f"Gap {analysis.gap_pos:,}: Skip ({repair_reason}) [{analysis.error_type}/{analysis.error_subtype}]")
        
        return filtered_analyses
    
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
        type_counts = defaultdict(int)
        
        for gap_pos in gap_positions:
            analysis = self.analyze_gap(gap_pos, max_search_distance)
            if analysis:
                all_analyses.append(analysis)
                if analysis.gap_in_error_region:
                    gaps_in_error += 1
                    if analysis.large_distance_anchor:
                        gaps_with_large_distance += 1
                    type_counts[analysis.error_type] += 1
                else:
                    gaps_not_in_error += 1
        
        print(f"\nSuccessfully analyzed {len(all_analyses)}/{len(gap_positions)} gaps")
        print(f"  - {gaps_in_error} gaps in error regions")
        print(f"    * {gaps_with_large_distance} with large distance anchor warnings")
        print("\nError type distribution:")
        for err_type, count in type_counts.items():
            print(f"    {err_type}: {count}")
        print(f"  - {gaps_not_in_error} gaps not in error regions")
        
        return all_analyses
    
    def apply_unified_repairs(self, query_fasta: str, gap_analyses: List[_GapAnalysis], 
                            output_fasta: str, final_gap_length: int = 100,
                            repair_mode: str = "conservative") -> Tuple[str, List[Dict], Dict]:
        """Apply all repairs uniformly"""
        print(f"\n" + "=" * 60)
        print(f"Applying repairs (mode: {repair_mode})")
        print("=" * 60)
        
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython not installed, cannot perform sequence repair")
        
        records = list(SeqIO.parse(query_fasta, "fasta"))
        if len(records) != 1:
            print(f"Warning: Expected single contig, but found {len(records)}")
        
        record = records[0]
        self.seq_length = len(record.seq)
        print(f"Sequence length: {self.seq_length:,} bp")
        
        print(f"Filtering gap analysis results based on {repair_mode} mode...")
        analyses_to_repair = self.filter_analyses_by_repair_mode(gap_analyses, repair_mode)
        
        print(f"Total analyses: {len(gap_analyses)}")
        print(f"Recommended repairs: {len(analyses_to_repair)}")
        
        repair_stats = {
            'recommended_repairs': len(analyses_to_repair),
            'skipped_conservative': 0,
            'aggressive_repairs': 0,
            'skipped_aggressive': 0,
            'type5_repairs': 0
        }
        
        if repair_mode == "conservative":
            large_distance_count = sum(1 for a in gap_analyses if a.large_distance_anchor and a.gap_in_error_region)
            repair_stats['skipped_conservative'] = large_distance_count
            print(f"Conservative mode skipped (large distance anchors): {large_distance_count}")
        
        seq = str(record.seq)
        seq_length = len(seq)
        
        replace_regions = []
        for analysis in analyses_to_repair:
            if (analysis.query_contig == record.id and
                analysis.replace_start > 0 and analysis.replace_end > 0):
                
                start = max(1, analysis.replace_start)
                end = min(seq_length, analysis.replace_end)
                
                if analysis.error_type == "type5":
                    repair_stats['type5_repairs'] += 1
                    
                if start <= end:
                    region = {
                        'start': start,
                        'end': end,
                        'gap_pos': analysis.gap_pos,
                        'error_type': analysis.error_type,
                        'error_subtype': analysis.error_subtype,
                        'error_features': analysis.error_features,
                        'large_distance': analysis.large_distance_anchor,
                        'confidence': analysis.confidence,
                        'confidence_score': analysis.type_confidence_score,
                        'total_anchor_length': analysis.total_anchor_length,
                        'replace_length': analysis.replace_length,
                        'repair_reason': analysis.repair_reason,
                    }
                    
                    if analysis.error_type == "type5":
                        region['ref_overlap_length'] = analysis.ref_overlap_length
                        region['ref_overlap_region'] = analysis.ref_overlap_region
                        region['extra_margin'] = analysis.type5_extra_margin
                        region['replace_side'] = analysis.type5_replace_side
                    
                    replace_regions.append(region)
        
        if not replace_regions:
            print("No regions need repair")
            return query_fasta, [], repair_stats
        
        print(f"Found {len(replace_regions)} regions to repair")
        
        replace_regions.sort(key=lambda x: x['start'])
        
        merged_regions = []
        if replace_regions:
            current = replace_regions[0].copy()
            
            for region in replace_regions[1:]:
                if region['start'] <= current['end'] + 1:
                    current['end'] = max(current['end'], region['end'])
                    current['error_type'] = f"{current['error_type']}+{region['error_type']}"
                    current['error_subtype'] = f"{current['error_subtype']}+{region['error_subtype']}"
                    if 'ref_overlap_length' in current and 'ref_overlap_length' in region:
                        current['ref_overlap_length'] += region['ref_overlap_length']
                else:
                    merged_regions.append(current)
                    current = region.copy()
            
            merged_regions.append(current)
        
        print(f"Merged into {len(merged_regions)} regions")
        
        seq_list = list(seq)
        repair_log = []
        total_deleted = 0
        total_replaced = 0
        
        for region in reversed(merged_regions):
            start, end = region['start'], region['end']
            original_length = end - start + 1
            
            if region['error_type'] == 'type5' or 'type5' in str(region['error_type']):
                ref_overlap = region.get('ref_overlap_length', 0)
                margin = region.get('extra_margin', 100)
                
                print(f"\nType5 deletion:")
                print(f"  - Region: {start:,}-{end:,} ({original_length:,} bp)")
                print(f"  - Reference overlap: {ref_overlap:,} bp")
                print(f"  - Extra margin: {margin} bp (for gap filling)")
                print(f"  - Total deleted: {original_length:,} bp")
                
                del seq_list[start-1:end]
                total_deleted += original_length
                
                self.type5_details.append({
                    'gap_pos': region['gap_pos'],
                    'overlap_length': ref_overlap,
                    'extra_margin': margin,
                    'total_deleted': original_length,
                    'replace_side': region.get('replace_side', 'right'),
                    'region': f"{start}-{end}"
                })
                
                repair_entry = {
                    'gap_position': region['gap_pos'],
                    'delete_region': f"{start:,}-{end:,}",
                    'deleted_length': original_length,
                    'ref_overlap_length': ref_overlap,
                    'extra_margin': margin,
                    'error_type': region['error_type'],
                    'error_subtype': region['error_subtype'],
                    'action': 'delete_with_margin',
                    'repair_reason': region['repair_reason']
                }
            else:
                replacement = 'N' * original_length
                seq_list[start-1:end] = list(replacement)
                total_replaced += original_length
                
                print(f"Replace: {start:,}-{end:,} ({original_length:,}bp) - {region['error_type']}")
                
                repair_entry = {
                    'gap_position': region['gap_pos'],
                    'replace_region': f"{start:,}-{end:,}",
                    'original_length': original_length,
                    'error_type': region['error_type'],
                    'error_subtype': region['error_subtype'],
                    'action': 'replace',
                    'repair_reason': region['repair_reason']
                }
            
            repair_log.append(repair_entry)
        
        repaired_seq = ''.join(seq_list)
        
        final_seq = self.normalize_gaps(repaired_seq, final_gap_length)
        
        length_change = len(final_seq) - self.seq_length
        print(f"\nLength change summary:")
        print(f"  - Original length: {self.seq_length:,} bp")
        print(f"  - Final length: {len(final_seq):,} bp")
        print(f"  - Net change: {length_change:+,} bp")
        print(f"  - Total deleted (Type5 with margin): {total_deleted:,} bp")
        print(f"  - Total replaced (other types): {total_replaced:,} bp")
        
        new_record = SeqRecord(Seq(final_seq), id=record.id,
                              description=f"repaired_{repair_mode}_{len(repair_log)}_gaps")
        SeqIO.write([new_record], output_fasta, "fasta")
        
        repair_stats.update({
            'total_deleted_with_margin': total_deleted,
            'total_replaced': total_replaced,
            'net_length_change': length_change
        })
        
        print(f"\nRepair statistics:")
        print(f"  - Total repair regions: {len(repair_log)}")
        print(f"  - Type5 repairs (with margin): {repair_stats['type5_repairs']}")
        print(f"  - Total deleted with margin: {total_deleted:,} bp")
        print(f"  - Total replaced: {total_replaced:,} bp")
        print(f"  - Repair mode: {repair_mode}")
        print(f"  - Final gap length: {final_gap_length} bp")
        
        return output_fasta, repair_log, repair_stats
    
    def normalize_gaps(self, sequence: str, gap_length: int = 100) -> str:
        """Adjust all consecutive N regions to specified length"""
        result = []
        i = 0
        n = len(sequence)
        gap_count = 0
        total_ns_original = 0
        total_ns_final = 0
        
        while i < n:
            if sequence[i] == 'N':
                start = i
                while i < n and sequence[i] == 'N':
                    i += 1
                original_gap_length = i - start
                gap_count += 1
                total_ns_original += original_gap_length
                
                result.append('N' * gap_length)
                total_ns_final += gap_length
            else:
                result.append(sequence[i])
                i += 1
        
        final_seq = ''.join(result)
        
        if gap_count > 0:
            print(f"\nGap normalization:")
            print(f"  - Found {gap_count} gap regions")
            print(f"  - Original total Ns: {total_ns_original:,} bp")
            print(f"  - Final total Ns: {total_ns_final:,} bp")
            print(f"  - Ns reduced by: {total_ns_original - total_ns_final:,} bp")
        
        return final_seq
    
    def save_report(self, gap_analyses: List[_GapAnalysis], repair_log: List[Dict], 
                   output_json: str, query_fasta: str = None, repair_mode: str = "conservative") -> Dict[str, Any]:
        """Save analysis report to JSON file"""
        
        error_type_counts = defaultdict(int)
        error_subtype_counts = defaultdict(int)
        feature_counts = defaultdict(int)
        confidence_counts = defaultdict(int)
        
        for analysis in gap_analyses:
            error_type_counts[analysis.error_type] += 1
            error_subtype_counts[analysis.error_subtype] += 1
            confidence_counts[analysis.confidence] += 1
            
            for feature in analysis.error_features:
                feature_counts[feature] += 1
        
        repair_stats = {
            'recommended_repairs': len(repair_log),
            'skipped_conservative': 0,
            'aggressive_repairs': 0,
            'skipped_aggressive': 0,
            'type5_repairs': sum(1 for a in gap_analyses if a.error_type == "type5" and a.gap_in_error_region)
        }
        
        total_deleted = sum(r.get('deleted_length', 0) for r in repair_log if r.get('action') == 'delete_with_margin')
        total_replaced = sum(r.get('original_length', 0) for r in repair_log if r.get('action') == 'replace')
        
        report = {
            'summary': {
                'total_gaps_analyzed': len(gap_analyses),
                'gaps_in_error_region': sum(1 for a in gap_analyses if a.gap_in_error_region),
                'gaps_not_in_error_region': sum(1 for a in gap_analyses if not a.gap_in_error_region),
                'gaps_with_large_distance_anchor': sum(1 for a in gap_analyses if a.large_distance_anchor),
                'gaps_repaired': len(repair_log),
                'total_deleted_with_margin': total_deleted,
                'total_replaced': total_replaced,
                'final_gap_length': 100,
                'repair_mode': repair_mode,
                'repair_summary': repair_stats,
                'error_type_counts': dict(error_type_counts),
                'error_subtype_counts': dict(error_subtype_counts),
                'feature_counts': dict(sorted(feature_counts.items(), key=lambda x: x[1], reverse=True)),
                'confidence_counts': dict(confidence_counts),
                'type5_details': self.type5_details if hasattr(self, 'type5_details') else []
            },
            'gap_analyses': [],
            'repair_details': repair_log
        }
        
        if query_fasta and os.path.exists(query_fasta):
            records = list(SeqIO.parse(query_fasta, "fasta"))
            if records:
                report['summary']['query_contig'] = records[0].id
                report['summary']['query_length'] = len(records[0].seq)
        
        for analysis in gap_analyses:
            gap_report = {
                'gap_position': analysis.gap_pos,
                'error_type': analysis.error_type,
                'error_subtype': analysis.error_subtype,
                'error_features': analysis.error_features,
                'confidence': analysis.confidence,
                'confidence_score': analysis.type_confidence_score,
                'gap_in_error_region': analysis.gap_in_error_region,
                'large_distance_anchor': analysis.large_distance_anchor,
                'total_anchor_length': analysis.total_anchor_length,
                'replace_region': f"{analysis.replace_start}-{analysis.replace_end}" if analysis.replace_start > 0 else "None",
                'replace_length': analysis.replace_length if analysis.replace_start > 0 else 0,
                'supporting_ref_contigs': len(analysis.supporting_ref_contigs),
                'repair_recommended': analysis.repair_reason is not None,
                'repair_reason': analysis.repair_reason,
                'anchor_distances': analysis.anchor_distances
            }
            
            if analysis.error_type == "type5":
                gap_report['ref_overlap_length'] = analysis.ref_overlap_length
                gap_report['ref_overlap_region'] = f"{analysis.ref_overlap_region[0]}-{analysis.ref_overlap_region[1]}"
                gap_report['ref_overlap_in_query'] = f"{analysis.ref_overlap_in_query[0]}-{analysis.ref_overlap_in_query[1]}"
                gap_report['type5_replace_side'] = analysis.type5_replace_side
                gap_report['type5_extra_margin'] = analysis.type5_extra_margin
            
            report['gap_analyses'].append(gap_report)
        
        with open(output_json, 'w') as f:
            json.dump(report, f, indent=2, default=lambda x: dict(x) if hasattr(x, '__dict__') else str(x))
        
        return report

### Command Line Interface ###
def main():
    """Command line interface main function"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Genome Gap Analysis and Repair Tool - Enhanced Version with Type1-Type6 Classification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Error Types (MODIFIED - Type5 now based on REFERENCE overlap with EXTRA MARGIN):
  Type1: Simple gap - normal gap between alignments
  Type2: Direction conflict - alignments on same reference with opposite directions
  Type3: Reference conflict - alignments on different reference contigs
  Type4: Reference overlap - alignments overlap on reference genome (original)
  Type5: Reference overlap with margin - alignments overlap on reference genome, delete one side + extra margin
  Type6: Complex - multiple conflicting features

The extra margin in Type5 (default 100bp) creates space for gap filling software.

Usage Examples:
  # Basic analysis
  python gap_analyzer.py -coords align.coords --gaps 1804444 12272394
  
  # Analysis and repair with Type5 handling (based on reference overlap)
  python gap_analyzer.py -coords align.coords --gaps 7053263 --query genome.fasta --type5-side right --type5-margin 100
  
  # API mode
  python gap_analyzer.py -coords align.coords --gaps 1804444 --query genome.fasta --api-mode --verbose
        """
    )
    
    parser.add_argument('-coords', required=True, help='Alignment coordinates file (coords format)')
    parser.add_argument('--gaps', nargs='+', type=int, required=True, help='List of gap positions')
    
    parser.add_argument('--query', help='Query genome FASTA file (for repair)')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--prefix', default='gap_analysis', help='Output file prefix')
    parser.add_argument('--json', default='gap_analysis_report.json', help='JSON report file')
    parser.add_argument('--gap-length', type=int, default=100, help='Final gap length (default 100bp)')
    parser.add_argument('--search', type=int, default=500000, help='Maximum search distance (bp)')
    parser.add_argument('--step', type=int, default=100000, help='Search step size (bp)')
    parser.add_argument('--min-confidence', choices=['low', 'medium', 'high'], default='low')
    parser.add_argument('--repair-mode', choices=['conservative', 'aggressive'], default='conservative')
    parser.add_argument('--type5-side', choices=['left', 'right'], default='right',
                       help='For Type5 errors (based on reference overlap), which side to delete')
    parser.add_argument('--type5-margin', type=int, default=100,
                       help='Extra bp to delete for Type5 to facilitate gap filling (default 100)')
    parser.add_argument('--api-mode', action='store_true', help='Use API mode')
    parser.add_argument('--verbose', action='store_true', help='Show detailed output')
    parser.add_argument('--log-file', help='Log file path')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.coords):
        print(f"Error: coords file does not exist: {args.coords}", file=sys.stderr)
        sys.exit(1)
    
    if args.query and not os.path.exists(args.query):
        print(f"Error: query file does not exist: {args.query}", file=sys.stderr)
        sys.exit(1)
    
    try:
        if args.api_mode:
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
                repair_mode=args.repair_mode,
                type5_replace_side=args.type5_side,
                type5_extra_margin=args.type5_margin
            )
            
            print(f"\nDetailed report saved to: {result['api_info']['output_files']['report']}")
            if result['api_info']['output_files']['filled_fasta']:
                print(f"Repaired genome: {result['api_info']['output_files']['filled_fasta']}")
            
        else:
            analyzer = _GapAnalyzer(type5_replace_side=args.type5_side, type5_extra_margin=args.type5_margin)
            
            original_method = analyzer.find_surrounding_alignments
            def new_find_surrounding_alignments(query_contig, gap_pos, max_distance):
                return original_method(query_contig, gap_pos, max_distance, args.step)
            analyzer.find_surrounding_alignments = new_find_surrounding_alignments
            
            print("Parsing alignment file and building synteny blocks...")
            alignments = analyzer.parse_coords_file(args.coords)
            
            if not alignments:
                print("Error: No alignment information parsed", file=sys.stderr)
                sys.exit(1)
            
            print("Analyzing all gaps with enhanced classification...")
            gap_analyses = analyzer.analyze_all_gaps(args.gaps, args.search)
            
            if not gap_analyses:
                print("Error: No gaps successfully analyzed", file=sys.stderr)
                sys.exit(1)
            
            repair_log = []
            repair_stats = {}
            if args.query:
                print(f"Applying repairs (mode: {args.repair_mode})...")
                output_fasta = os.path.join(args.output_dir, f"{args.prefix}_repaired.fasta")
                output_fasta, repair_log, repair_stats = analyzer.apply_unified_repairs(
                    args.query, gap_analyses, output_fasta, args.gap_length, args.repair_mode
                )
                print(f"Repaired genome saved: {output_fasta}")
                print(f"Performed {len(repair_log)} repairs ({args.repair_mode} mode)")
            
            report_file = os.path.join(args.output_dir, args.json)
            report = analyzer.save_report(gap_analyses, repair_log, report_file, args.query, args.repair_mode)
            print(f"Report saved: {report_file}")
            
            print("\n" + "="*60)
            print("Analysis Summary")
            print("="*60)
            print(f"Successfully analyzed: {len(gap_analyses)}/{len(args.gaps)} gaps")
            print(f"Gaps in error regions: {report['summary']['gaps_in_error_region']}")
            print(f"Gaps not in error regions: {report['summary']['gaps_not_in_error_region']}")
            print(f"Repair mode: {args.repair_mode}")
            print(f"Type5 extra margin: {args.type5_margin} bp")
            
            print("\nError type distribution:")
            for err_type, count in report['summary']['error_type_counts'].items():
                print(f"  {err_type}: {count}")
            
            if 'error_subtype_counts' in report['summary']:
                print("\nError subtype distribution:")
                for subtype, count in list(report['summary']['error_subtype_counts'].items())[:5]:
                    print(f"  {subtype}: {count}")
            
            if 'type5_details' in report['summary'] and report['summary']['type5_details']:
                print("\nType5 details (based on reference overlap with margin):")
                for detail in report['summary']['type5_details'][:5]:
                    print(f"  Gap {detail['gap_pos']}: ref overlap {detail['overlap_length']}bp, "
                          f"margin {detail.get('extra_margin', 100)}bp, "
                          f"total deleted {detail.get('total_deleted', 0)}bp, "
                          f"replaced {detail['replace_side']} side")
            
            if repair_log:
                total_deleted = sum(r.get('deleted_length', 0) for r in repair_log if r.get('action') == 'delete_with_margin')
                total_replaced = sum(r.get('original_length', 0) for r in repair_log if r.get('action') == 'replace')
                print(f"\nTotal deleted (Type5 with margin): {total_deleted:,} bp")
                print(f"Total replaced (other types): {total_replaced:,} bp")
                print(f"Final gap length: {args.gap_length} bp")
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()