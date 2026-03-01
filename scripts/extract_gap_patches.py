#!/usr/bin/env python3
"""
Gap Patch Sequence Extractor - 修复边界版本
Extract gap patch sequences from coords files, 修复了边界判断问题
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
    """Unified API for gap patch extraction - 修复边界版本"""
    
    def __init__(self, verbose: bool = True, log_file: Optional[str] = None, debug: bool = False):
        """
        Initialize GapPatchesExtractor API
        
        Args:
            verbose: Show detailed output
            log_file: Log file path (optional)
            debug: Enable debug mode with extra output
        """
        self.verbose = verbose
        self.log_file = log_file
        self.debug = debug
        self.temp_files = []
        self._setup_logging()
    
    def _setup_logging(self):
        """Configure logging system"""
        self.logger = logging.getLogger('GapPatchesExtractor')
        log_level = logging.DEBUG if self.debug else (logging.INFO if self.verbose else logging.WARNING)
        self.logger.setLevel(log_level)
        
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(log_level)
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
        output_json: Optional[str] = None,
        boundary_tolerance: int = 10,  # 新增参数：边界容差
        allow_adjacent_anchors: bool = True  # 新增参数：允许相邻锚点
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
            boundary_tolerance: Tolerance for boundary matching (bp)
            allow_adjacent_anchors: Allow anchors that are adjacent to gap
            
        Returns:
            Dictionary containing extraction results
        """
        self.log(f"Starting gap patch sequence extraction (修复边界版本)", "info")
        self.log(f"边界容差: {boundary_tolerance} bp", "info")
        self.log(f"允许相邻锚点: {allow_adjacent_anchors}", "info")
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
        
        # 显示前几个比对用于调试
        if self.debug:
            self.log("前5个比对记录:", "debug")
            for i, align in enumerate(alignments[:5]):
                ref_s, ref_e, qry_s, qry_e, ref_c, qry_c = align
                self.log(f"  {i+1}: ref={ref_c}:{min(ref_s,ref_e)}-{max(ref_s,ref_e)}, qry={qry_c}:{min(qry_s,qry_e)}-{max(qry_s,qry_e)}", "debug")
        
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
                'boundary_tolerance': boundary_tolerance,
                'allow_adjacent_anchors': allow_adjacent_anchors,
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
        
        # 对每个gap位置进行排序处理，便于调试
        sorted_gaps = sorted(gap_positions)
        self.log(f"排序后的gap位置: {sorted_gaps}", "debug")
        
        for gap_pos in sorted_gaps:
            self.log(f"\n{'='*60}", "info")
            self.log(f"Processing gap position: {gap_pos:,}", "info")
            self.log(f"{'='*60}", "info")
            
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
                boundary_tolerance=boundary_tolerance,
                allow_adjacent_anchors=allow_adjacent_anchors,
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
        line_count = 0
        
        try:
            with open(coords_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    line_count += 1
                    
                    # 跳过注释行和空行
                    if not line or line.startswith('/') or line.startswith('NUCMER') or \
                       line.startswith('[') or line.startswith('=') or line.startswith('#'):
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
                            
                            # 确保坐标是正向的
                            ref_min, ref_max = min(ref_start, ref_end), max(ref_start, ref_end)
                            qry_min, qry_max = min(qry_start, qry_end), max(qry_start, qry_end)
                            
                            alignments.append((ref_min, ref_max, qry_min, qry_max, ref_contig, qry_contig))
                            
                            if self.debug and len(alignments) <= 10:
                                self.log(f"Parsed alignment: ref={ref_contig}:{ref_min}-{ref_max}, qry={qry_contig}:{qry_min}-{qry_max}", "debug")
                                
                        except (ValueError, IndexError) as e:
                            if self.debug:
                                self.log(f"Failed to parse line {line_count}: {line}", "debug")
                            continue
        
        except Exception as e:
            self.log(f"Failed to parse coords file: {e}", "error")
            raise
        
        self.log(f"Parsed {len(alignments)} alignments from {line_count} lines", "info")
        return alignments
    
    def _group_alignments_by_ref_contig(self, alignments: List[Tuple]) -> Dict:
        """Group syntenic regions by reference contig"""
        ref_contig_groups = {}
        for align in alignments:
            *_, ref_contig, qry_contig = align
            if ref_contig not in ref_contig_groups:
                ref_contig_groups[ref_contig] = []
            ref_contig_groups[ref_contig].append(align)
        
        # 按query坐标排序每个组内的比对
        for ref_contig in ref_contig_groups:
            ref_contig_groups[ref_contig].sort(key=lambda x: x[2])  # 按qry_min排序
        
        return ref_contig_groups
    
    def _find_surrounding_alignments(self, alignments: List[Tuple], gap_pos: int, 
                                   require_both: bool = True, boundary_tolerance: int = 10,
                                   allow_adjacent: bool = True) -> Tuple[Optional[Tuple], Optional[Tuple]]:
        """
        修复版本的锚点查找函数
        增加边界容差，允许相邻锚点
        """
        left_candidates = []
        right_candidates = []
        
        self.log(f"   查找gap位置 {gap_pos} 周围的锚点 (容差: {boundary_tolerance} bp)", "debug")
        
        for align in alignments:
            _, _, qry_min, qry_max, ref_contig, qry_contig = align
            
            # 调试输出
            if self.debug:
                self.log(f"    检查比对: {ref_contig}:{qry_min}-{qry_max}", "debug")
                self.log(f"      qry_max={qry_max}, gap_pos={gap_pos}, 差值={qry_max - gap_pos}", "debug")
                self.log(f"      qry_min={qry_min}, gap_pos={gap_pos}, 差值={qry_min - gap_pos}", "debug")
            
            # 左侧锚点：比对结束位置在gap之前（或允许容差范围内）
            if qry_max <= gap_pos + boundary_tolerance:
                distance = gap_pos - qry_max
                left_candidates.append((qry_max, distance, align))
                if self.debug:
                    self.log(f"      ✓ 左侧锚点候选: 距离gap {distance} bp", "debug")
            
            # 右侧锚点：比对开始位置在gap之后（或允许容差范围内）
            if qry_min >= gap_pos - boundary_tolerance:
                distance = qry_min - gap_pos
                right_candidates.append((qry_min, distance, align))
                if self.debug:
                    self.log(f"      ✓ 右侧锚点候选: 距离gap {distance} bp", "debug")
        
        # 选择最接近的锚点
        before_align = None
        after_align = None
        before_distance = None
        after_distance = None
        
        if left_candidates:
            # 选择最接近gap的左侧锚点（最大的qry_max）
            left_candidates.sort(key=lambda x: x[0], reverse=True)
            before_align = left_candidates[0][2]
            before_distance = left_candidates[0][1]
            self.log(f"    选择左侧锚点: 距离gap {before_distance} bp", "debug")
        
        if right_candidates:
            # 选择最接近gap的右侧锚点（最小的qry_min）
            right_candidates.sort(key=lambda x: x[0])
            after_align = right_candidates[0][2]
            after_distance = right_candidates[0][1]
            self.log(f"    选择右侧锚点: 距离gap {after_distance} bp", "debug")
        
        # 检查锚点有效性
        if require_both and (before_align is None or after_align is None):
            self.log(f"    未找到完整的双侧锚点", "debug")
            return None, None
        
        # 检查锚点是否相邻或重叠
        if before_align and after_align:
            _, _, before_qry_min, before_qry_max, _, _ = before_align
            _, _, after_qry_min, after_qry_max, _, _ = after_align
            
            # 计算锚点间距离
            gap_between = after_qry_min - before_qry_max
            self.log(f"    锚点间距离: {gap_between} bp", "debug")
            
            if gap_between < 0:
                self.log(f"    警告: 锚点重叠 {abs(gap_between)} bp", "debug")
            elif gap_between == 0:
                self.log(f"    锚点相邻", "debug")
        
        return before_align, after_align
    
    def _get_patch_region_for_ref_contig(self, before_align: Tuple, after_align: Tuple, 
                                       flank_size: int) -> Tuple[Optional[int], Optional[int], int, Optional[str]]:
        """Build patch region with flank for single reference contig"""
        r1s, r1e, q1s, q1e, ref_ctg1, _ = before_align
        r2s, r2e, q2s, q2e, ref_ctg2, _ = after_align
        
        if ref_ctg1 != ref_ctg2:
            self.log(f"    左右锚点在不同contig上: {ref_ctg1} vs {ref_ctg2}", "debug")
            return None, None, 0, None
        
        left_ref_min, left_ref_max = min(r1s, r1e), max(r1s, r1e)
        right_ref_min, right_ref_max = min(r2s, r2e), max(r2s, r2e)
        
        self.log(f"    左侧锚点参考坐标: {left_ref_min}-{left_ref_max}", "debug")
        self.log(f"    右侧锚点参考坐标: {right_ref_min}-{right_ref_max}", "debug")
        
        if left_ref_max < right_ref_min:
            direction = 1
            gap_adj_left = left_ref_max
            gap_adj_right = right_ref_min
            self.log(f"    正向排列", "debug")
        elif right_ref_max < left_ref_min:
            direction = -1
            gap_adj_left = right_ref_max
            gap_adj_right = left_ref_min
            self.log(f"    反向排列", "debug")
        else:
            self.log(f"    锚点重叠或包含关系", "debug")
            return None, None, 0, None
        
        # 计算补丁区域
        patch_start = max(1, gap_adj_left - flank_size + 1)
        patch_end = gap_adj_right + flank_size - 1
        
        self.log(f"    补丁区域: {patch_start}-{patch_end} (长度: {patch_end - patch_start + 1} bp)", "debug")
        
        return patch_start, patch_end, direction, ref_ctg1
    
    def _extract_patch_sequence(self, ref_fasta: str, contig: str, start: int, end: int, 
                              direction: int, flank_size: int, patch_index: int, 
                              gap_pos: int, output_file: str) -> Tuple[bool, Optional[int]]:
        """Extract patch sequence and write to specified output file"""
        temp_file = f"temp_patch_{patch_index}_{int(time.time())}.fa"
        self.temp_files.append(temp_file)
        
        try:
            extract_cmd = ["samtools", "faidx", ref_fasta, f"{contig}:{start}-{end}"]
            self.log(f"    提取命令: {' '.join(extract_cmd)}", "debug")
            
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
                self.log(f"    应用反向互补", "debug")
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
                          boundary_tolerance: int, allow_adjacent_anchors: bool,
                          output_file: str, start_patch_counter: int) -> Dict[str, Any]:
        """Process single gap position - 修复版本"""
        gap_result = {
            'query_contig': qry_contig,
            'gap_position': gap_pos,
            'valid_reference_contigs': 0,
            'successful_patches': 0,
            'total_patches_attempted': 0,
            'patches': []
        }
        
        valid_pairs = []
        
        self.log(f"  检查 {len(ref_contig_groups)} 个参考contig", "info")
        
        for ref_contig, aligns in ref_contig_groups.items():
            self.log(f"  处理参考contig: {ref_contig}", "debug")
            
            # 筛选出当前query contig的比对
            qry_aligns = [align for align in aligns if align[5] == qry_contig]
            if len(qry_aligns) < 2:
                self.log(f"    contig {ref_contig}: 只有 {len(qry_aligns)} 个比对，跳过", "debug")
                continue
            
            # 按query坐标排序
            qry_aligns.sort(key=lambda x: x[2])  # 按qry_min排序
            
            if self.debug:
                self.log(f"    contig {ref_contig} 有 {len(qry_aligns)} 个比对:", "debug")
                for i, align in enumerate(qry_aligns):
                    _, _, q_min, q_max, _, _ = align
                    self.log(f"      {i+1}: {q_min}-{q_max}", "debug")
            
            before_align, after_align = self._find_surrounding_alignments(
                qry_aligns, gap_pos, require_both_anchors, 
                boundary_tolerance, allow_adjacent_anchors
            )
            
            if before_align and after_align:
                _, _, q1s, q1e, _, _ = before_align
                _, _, q2s, q2e, _, _ = after_align
                
                q1_max = max(q1s, q1e)
                q2_min = min(q2s, q2e)
                
                anchor_distance = q2_min - q1_max
                
                self.log(f"    contig {ref_contig}: 锚点距离 = {anchor_distance:,} bp", "debug")
                
                # 检查锚点距离是否在允许范围内
                distance_valid = True
                if anchor_distance < min_anchor_distance:
                    self.log(f"      距离小于最小值 {min_anchor_distance}", "debug")
                    distance_valid = False
                if anchor_distance > max_anchor_distance:
                    self.log(f"      距离大于最大值 {max_anchor_distance}", "debug")
                    distance_valid = False
                
                if distance_valid:
                    valid_pairs.append((ref_contig, before_align, after_align, anchor_distance))
                    self.log(f"    ✓ contig {ref_contig} 有效", "info")
                else:
                    self.log(f"    ✗ contig {ref_contig} 距离不满足要求", "debug")
            else:
                self.log(f"    contig {ref_contig}: 未找到完整锚点对", "debug")
        
        gap_result['valid_reference_contigs'] = len(valid_pairs)
        self.log(f"  找到 {len(valid_pairs)} 个有效参考contig", "info")
        
        if not valid_pairs:
            self.log(f"  没有找到有效的锚点对", "warning")
            return gap_result
        
        success_count = 0
        patch_counter = start_patch_counter
        
        for i, (ref_contig, before_align, after_align, anchor_distance) in enumerate(valid_pairs):
            patch_num = patch_counter + i
            self.log(f"  处理补丁 {patch_num} (contig {ref_contig})", "info")
            
            patch_start, patch_end, direction, ref_contig_name = self._get_patch_region_for_ref_contig(
                before_align, after_align, flank_size
            )
            
            if patch_start is None:
                self.log(f"    无法确定补丁区域", "warning")
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
                        'query': f"{min(before_align[2], before_align[3])}-{max(before_align[2], before_align[3])}",
                        'reference': f"{min(before_align[0], before_align[1])}-{max(before_align[0], before_align[1])}",
                        'contig': before_align[4]
                    },
                    'anchor_right': {
                        'query': f"{min(after_align[2], after_align[3])}-{max(after_align[2], after_align[3])}",
                        'reference': f"{min(after_align[0], after_align[1])}-{max(after_align[0], after_align[1])}",
                        'contig': after_align[4]
                    }
                }
                gap_result['patches'].append(patch_info)
        
        gap_result['successful_patches'] = success_count
        gap_result['total_patches_attempted'] = len(valid_pairs)
        
        self.log(f"  成功提取 {success_count}/{len(valid_pairs)} 个补丁序列", "info")
        
        return gap_result
    
    def _save_results_to_json(self, result_summary: Dict[str, Any], json_file: str):
        """Save result summary to JSON file"""
        try:
            # 转换defaultdict为普通dict
            result_summary['statistics']['patches_per_gap'] = dict(result_summary['statistics']['patches_per_gap'])
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
    
    def cleanup(self):
        """Clean up all temporary files"""
        self._cleanup_temp_files()
    
    def __del__(self):
        """Destructor, automatically clean temporary files"""
        self.cleanup()

def extract_gap_patches(
    reference_fasta: str,
    coords_file: str,
    gap_positions: List[int],
    output_file: str = "gap_patches.fasta",
    flank_size: int = 10000,
    verbose: bool = True,
    debug: bool = False,
    boundary_tolerance: int = 10,
    allow_adjacent_anchors: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Simplified interface for extracting gap patches (修复版本)
    
    Example:
        >>> from extract_gap_patches import extract_gap_patches
        >>> result = extract_gap_patches(
        ...     reference_fasta="reference.fasta",
        ...     coords_file="alignment.coords",
        ...     gap_positions=[3744012, 3744112],
        ...     output_file="patches.fasta",
        ...     flank_size=5000,
        ...     boundary_tolerance=10,  # 允许10bp的边界容差
        ...     allow_adjacent_anchors=True,  # 允许相邻锚点
        ...     debug=True  # 开启调试输出
        ... )
    """
    extractor = GapPatchesExtractorAPI(verbose=verbose, debug=debug)
    return extractor.extract_patches(
        reference_fasta=reference_fasta,
        coords_file=coords_file,
        gap_positions=gap_positions,
        output_file=output_file,
        flank_size=flank_size,
        boundary_tolerance=boundary_tolerance,
        allow_adjacent_anchors=allow_adjacent_anchors,
        **kwargs
    )

def main():
    """Command line interface main function"""
    parser = argparse.ArgumentParser(
        description='Extract gap patch sequences from coords files (修复边界版本)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # 基本用法，处理你的特定gap
  python extract_gap_patches.py -c reference.fasta -coords alignment.coords \\
    --gap-position 3744012 3744112 --boundary-tolerance 10 --allow-adjacent-anchors \\
    --min-anchor-distance 0 --debug
  
  # 处理所有gap
  python extract_gap_patches.py -c reference.fasta -coords alignment.coords \\
    --gap-position 95508 246186 751477 1580197 2080759 2550820 3215153 3744012 4488592 4734877 \\
    --boundary-tolerance 10 --allow-adjacent-anchors --output-dir results --verbose
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
    
    # 新增参数
    parser.add_argument("--boundary-tolerance", type=int, default=10,
                       help="Tolerance for boundary matching (bp), default 10")
    parser.add_argument("--allow-adjacent-anchors", action="store_true", default=True,
                       help="Allow anchors that are adjacent to gap (default)")
    parser.add_argument("--disallow-adjacent-anchors", action="store_false", dest="allow_adjacent_anchors",
                       help="Disallow adjacent anchors")
    
    parser.add_argument("--require-both-anchors", action="store_true", default=True,
                       help="Require both left and right anchors (default)")
    parser.add_argument("--allow-single-anchor", action="store_false", dest="require_both_anchors",
                       help="Allow single anchor")
    parser.add_argument("--min-anchor-distance", type=int, default=0,  # 修改默认值为0
                       help="Minimum anchor distance (bp), default 0")
    parser.add_argument("--max-anchor-distance", type=int, default=1000000,
                       help="Maximum anchor distance (bp), default 1000000")
    
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
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
        extractor = GapPatchesExtractorAPI(verbose=verbose, log_file=args.log_file, debug=args.debug)
        
        output_file = os.path.join(args.output_dir, args.output)
        output_json = os.path.join(args.output_dir, args.json)
        
        # 特别处理你的gap
        if args.debug:
            print("\n" + "="*60)
            print("处理特定gap: 3744012-3744112")
            print("="*60)
            
            # 先解析并显示比对信息
            alignments = extractor._parse_coords_file(args.coords)
            print(f"\n找到 {len(alignments)} 个比对记录")
            
            # 显示相关比对
            target_gap = 3744012
            print(f"\n与gap {target_gap} 相关的比对:")
            for align in alignments:
                _, _, q_min, q_max, ref_c, qry_c = align
                if abs(q_max - target_gap) < 1000 or abs(q_min - target_gap) < 1000:
                    print(f"  比对: {ref_c}:{q_min}-{q_max}")
                    print(f"    q_max={q_max}, 差值={q_max-target_gap}")
                    print(f"    q_min={q_min}, 差值={q_min-target_gap}")
        
        result = extractor.extract_patches(
            reference_fasta=args.reference_fasta,
            coords_file=args.coords,
            gap_positions=args.gap_positions,
            output_file=output_file,
            flank_size=args.flank_size,
            require_both_anchors=args.require_both_anchors,
            min_anchor_distance=args.min_anchor_distance,
            max_anchor_distance=args.max_anchor_distance,
            output_json=output_json,
            boundary_tolerance=args.boundary_tolerance,
            allow_adjacent_anchors=args.allow_adjacent_anchors
        )
        
        if verbose:
            print(f"\nDetailed report saved to: {output_json}")
            print(f"Patch sequences saved to: {output_file}")
        
        stats = result['statistics']
        if stats['successful_gaps'] == 0:
            print(f"\n警告: 没有成功提取任何补丁", file=sys.stderr)
            print(f"建议: 尝试增加 --boundary-tolerance 或设置 --min-anchor-distance 0", file=sys.stderr)
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
