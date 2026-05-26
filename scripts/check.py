#!/usr/bin/env python3
"""
fill_quality_check.py - Gap filling quality assessment with allele-specific analysis

Evaluate the quality of filled gaps in a genome assembly using HiFi reads.
Features: allele-specific depth, boundary consistency, coverage uniformity.
"""

import os
import sys
import subprocess
import argparse
import re
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, Any
from Bio import SeqIO
from collections import defaultdict
import threading
import numpy as np


class SimpleGapMapperWithQuality:
    """Gap mapper with CRAQ-style quality assessment"""
    
    def __init__(self, verbose=True, threads_per_alignment=2, keep_intermediate=False,
                 hifi_reads=None, min_depth=5, min_coverage=0.5, flank_size=5000,
                 min_zero_stretch_bp=1, boundary_good_ratio=0.25, boundary_min_spanning=5,
                 min_alignment_coverage=0.7, low_coverage_threshold=3,
                 max_pdf_length=100000, pdf_step_auto=True):
        self.verbose = verbose
        self.threads_per_alignment = threads_per_alignment
        self.keep_intermediate = keep_intermediate
        self.hifi_reads = hifi_reads
        self.min_depth = min_depth
        self.min_coverage = min_coverage
        self.flank_size = flank_size
        self.min_zero_stretch_bp = min_zero_stretch_bp
        self.boundary_good_ratio = boundary_good_ratio
        self.boundary_min_spanning = boundary_min_spanning
        self.min_alignment_coverage = min_alignment_coverage
        self.low_coverage_threshold = low_coverage_threshold
        self.max_pdf_length = max_pdf_length
        self.pdf_step_auto = pdf_step_auto
        self.lock = threading.Lock()
        self.completed = 0
        self.total = 0
        
        self.bam_file = None
        
        self.unmapped_count = 0
        self.mapped_with_n_count = 0
        self.mapped_clean_count = 0
        self.high_quality_count = 0
        self.good_quality_count = 0
        self.fair_quality_count = 0
        self.poor_quality_count = 0
        
        self.unmapped_gaps = []
        self.mapped_with_n_gaps = []
    
    def log(self, msg, level="info"):
        if self.verbose:
            prefix = {"info": "[INFO] ", "success": "[OK] ", "warning": "[WARNING] ", 
                      "error": "[ERROR] ", "debug": "[DEBUG] ", "quality": "[QUALITY] "}.get(level, "      ")
            print(f"{prefix}{msg}")
    
    def extract_chromosome(self, fasta_file, chrom, output_file):
        """Extract single chromosome from FASTA"""
        cmd = f"samtools faidx {fasta_file} {chrom} > {output_file} 2>/dev/null"
        result = subprocess.run(cmd, shell=True, capture_output=True)
        if result.returncode == 0 and os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            return True
        
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    header_name = header.split()[0]
                    if header_name == chrom:
                        cmd = f'samtools faidx {fasta_file} "{header}" > {output_file} 2>/dev/null'
                        result = subprocess.run(cmd, shell=True, capture_output=True)
                        return result.returncode == 0 and os.path.exists(output_file) and os.path.getsize(output_file) > 0
        return False
    
    def run_nucmer_alignment(self, before_fasta, after_fasta, output_prefix):
        """Run nucmer alignment"""
        try:
            cmd = f"nucmer -c 1000 -l 100 --batch=500000000 -t {self.threads_per_alignment} --prefix={output_prefix} {before_fasta} {after_fasta} > {output_prefix}.log 2>&1"
            subprocess.run(cmd, shell=True, capture_output=True)
            
            delta_file = f"{output_prefix}.delta"
            if not os.path.exists(delta_file) or os.path.getsize(delta_file) == 0:
                return [], "No alignments found"
            
            cmd = f"delta-filter -i -r {delta_file} > {output_prefix}.filter.delta"
            subprocess.run(cmd, shell=True, capture_output=True)
            
            cmd = f"show-coords -THrcl {output_prefix}.filter.delta > {output_prefix}.coords"
            subprocess.run(cmd, shell=True, capture_output=True)
            
            alignments = self.parse_coords(f"{output_prefix}.coords")
            
            if not self.keep_intermediate:
                for f in [f"{output_prefix}.delta", f"{output_prefix}.filter.delta", f"{output_prefix}.log"]:
                    if os.path.exists(f):
                        os.remove(f)
            
            return alignments, None
        except Exception as e:
            return [], str(e)
    
    def parse_coords(self, coords_file):
        """Parse nucmer coords file"""
        alignments = []
        if not os.path.exists(coords_file):
            return alignments
        
        with open(coords_file, 'r') as f:
            for line in f:
                if line.startswith('=') or line.startswith('-') or line.startswith('['):
                    continue
                parts = line.strip().split()
                if len(parts) >= 13:
                    try:
                        alignments.append({
                            'before_start': int(parts[0]),
                            'before_end': int(parts[1]),
                            'after_start': int(parts[2]),
                            'after_end': int(parts[3]),
                            'before_len': int(parts[4]),
                            'after_len': int(parts[5]),
                            'identity': float(parts[6]),
                            'before_chrom': parts[11],
                            'after_chrom': parts[12]
                        })
                    except:
                        continue
        return alignments
    
    def find_left_right_alignments(self, alignments, before_contig, gap_start, gap_end, max_distance=500000):
        """Find alignments flanking the gap"""
        contig_aligns = [a for a in alignments if a['before_chrom'] == before_contig]
        
        if not contig_aligns:
            return None, None, [], None, None
        
        all_starts = [min(a['before_start'], a['before_end']) for a in contig_aligns]
        all_ends = [max(a['before_start'], a['before_end']) for a in contig_aligns]
        chrom_min = min(all_starts)
        chrom_max = max(all_ends)
        
        left_candidates = []
        right_candidates = []
        containing_candidates = []
        
        overlap_tolerance = 10
        
        for align in contig_aligns:
            b_start = min(align['before_start'], align['before_end'])
            b_end = max(align['before_start'], align['before_end'])
            
            if b_start < gap_start and gap_end < b_end:
                containing_candidates.append(align)
            elif b_end <= gap_start + overlap_tolerance:
                dist = max(0, gap_start - b_end)
                if dist <= max_distance:
                    left_candidates.append((align, dist, b_end))
            elif b_start >= gap_end - overlap_tolerance:
                dist = max(0, b_start - gap_end)
                if dist <= max_distance:
                    right_candidates.append((align, dist, b_start))
        
        left_candidates.sort(key=lambda x: x[1])
        right_candidates.sort(key=lambda x: x[1])
        
        left_best = left_candidates[0][0] if left_candidates else None
        right_best = right_candidates[0][0] if right_candidates else None
        
        if not left_best and gap_start < chrom_min + 1000:
            left_best = {
                'before_start': chrom_min - 1,
                'before_end': chrom_min - 1,
                'after_start': 0,
                'after_end': 0,
                'identity': 100.0,
                'after_chrom': before_contig,
                'before_chrom': before_contig
            }
        
        if not right_best and gap_end > chrom_max - 1000:
            right_best = {
                'before_start': chrom_max + 1,
                'before_end': chrom_max + 1,
                'after_start': chrom_max + 1,
                'after_end': chrom_max + 1,
                'identity': 100.0,
                'after_chrom': before_contig,
                'before_chrom': before_contig
            }
        
        return left_best, right_best, containing_candidates, None, None
    
    def map_gap_between_alignments(self, gap, left_aln, right_aln):
        """Map gap located between two alignments"""
        is_virtual_left = (left_aln.get('before_start') == left_aln.get('before_end') and 
                          left_aln.get('after_start') == 0)
        is_virtual_right = (right_aln.get('before_start') == right_aln.get('before_end') and 
                           right_aln.get('after_start') == right_aln.get('before_start'))
        
        if not is_virtual_left:
            if left_aln['after_start'] < left_aln['after_end']:
                left_after_end = left_aln['after_end']
            else:
                left_after_end = left_aln['after_start']
        else:
            left_after_end = 0
        
        if not is_virtual_right:
            if right_aln['after_start'] < right_aln['after_end']:
                right_after_start = right_aln['after_start']
            else:
                right_after_start = right_aln['after_end']
        else:
            right_after_start = left_after_end + 1
        
        if not is_virtual_left:
            b_left_end = max(left_aln['before_start'], left_aln['before_end'])
            if b_left_end > gap['start']:
                overlap = b_left_end - gap['start']
                ratio = left_aln['after_len'] / left_aln['before_len'] if left_aln['before_len'] > 0 else 1.0
                if left_aln['after_start'] < left_aln['after_end']:
                    left_after_end = max(left_after_end - int(overlap * ratio), 0)
                else:
                    left_after_end = left_after_end + int(overlap * ratio)
        
        if not is_virtual_right:
            b_right_start = min(right_aln['before_start'], right_aln['before_end'])
            if b_right_start < gap['end']:
                overlap = gap['end'] - b_right_start
                ratio = right_aln['after_len'] / right_aln['before_len'] if right_aln['before_len'] > 0 else 1.0
                if right_aln['after_start'] < right_aln['after_end']:
                    right_after_start = right_after_start + int(overlap * ratio)
                else:
                    right_after_start = max(right_after_start - int(overlap * ratio), 0)
        
        if left_after_end > right_after_start:
            left_after_end, right_after_start = right_after_start, left_after_end
        
        filled_length = right_after_start - left_after_end
        
        return {
            'chrom': left_aln['after_chrom'],
            'start': left_after_end,
            'end': right_after_start,
            'length': filled_length,
            'identity': (left_aln.get('identity', 100.0) + right_aln.get('identity', 100.0)) / 2,
            'method': 'between_alignments'
        }
    
    def map_gap_inside_alignment(self, gap, containing_aln):
        """Map gap located inside an alignment"""
        b_start = min(containing_aln['before_start'], containing_aln['before_end'])
        b_end = max(containing_aln['before_start'], containing_aln['before_end'])
        a_start = min(containing_aln['after_start'], containing_aln['after_end'])
        a_end = max(containing_aln['after_start'], containing_aln['after_end'])
        
        ratio = (a_end - a_start) / (b_end - b_start) if b_end > b_start else 1.0
        
        gap_start_offset = gap['start'] - b_start
        gap_end_offset = gap['end'] - b_start
        
        if containing_aln['before_start'] < containing_aln['before_end']:
            filled_start = a_start + int(gap_start_offset * ratio)
            filled_end = a_start + int(gap_end_offset * ratio)
        else:
            filled_start = a_end - int(gap_end_offset * ratio)
            filled_end = a_end - int(gap_start_offset * ratio)
        
        if filled_start > filled_end:
            filled_start, filled_end = filled_end, filled_start
        
        return {
            'chrom': containing_aln['after_chrom'],
            'start': filled_start,
            'end': filled_end,
            'length': filled_end - filled_start,
            'identity': containing_aln['identity'],
            'method': 'inside_alignment'
        }
    
    # ============ Allele-Specific Depth Calculation Module ============
    
    def calculate_allele_specific_depth(self, bam_file, chrom, start, end, ref_fasta=None):
        """
        Calculate allele-specific coverage depth
        
        Returns for each position:
        - total_depth: Raw coverage (all reads)
        - effective_depth: Coverage supporting the major allele
        - major_base: The most frequent base at this position
        - major_freq: Frequency of the major allele
        - genotype: HOMOZYGOUS (≥80%), HETEROZYGOUS (20-80%), LOW_FREQ (<20%)
        - base_counts: Counts for each base (A, C, G, T)
        """
        region = f"{chrom}:{start}-{end}"
        
        # Get reference sequence if provided
        ref_seq = None
        if ref_fasta:
            cmd = f"samtools faidx {ref_fasta} {region} 2>/dev/null | grep -v '^>' | tr -d '\\n'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            ref_seq = result.stdout.upper() if result.stdout else None
        
        # Use samtools mpileup for base-level counts
        cmd = f"samtools mpileup -A -Q 0 -r {region} {bam_file} 2>/dev/null"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        position_data = {}
        
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            
            try:
                pos_chrom = parts[0]
                pos = int(parts[1])
                ref_base = parts[2].upper()
                depth = int(parts[3])
                pileup = parts[4] if len(parts) > 4 else ""
            except (ValueError, IndexError):
                continue
            
            # Parse pileup string to count each base
            base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'del': 0}
            
            i = 0
            while i < len(pileup):
                c = pileup[i]
                
                # Handle indel markers
                if c == '+' or c == '-':
                    # Read indel length
                    j = i + 1
                    indel_len = 0
                    while j < len(pileup) and pileup[j].isdigit():
                        indel_len = indel_len * 10 + int(pileup[j])
                        j += 1
                    # Skip the indel sequence
                    i = j + indel_len
                    if c == '-':
                        base_counts['del'] += 1
                    continue
                
                # Handle base calls
                if c == '.' or c == ',':
                    # Match to reference base
                    base_counts[ref_base] = base_counts.get(ref_base, 0) + 1
                elif c.upper() in base_counts:
                    base_counts[c.upper()] += 1
                elif c == '^':  # Start of read marker, skip next character
                    i += 1
                elif c == '$':  # End of read marker
                    pass
                
                i += 1
            
            # Find major allele
            if depth > 0:
                major_base = max(base_counts, key=lambda k: base_counts[k] if k != 'del' else 0)
                major_depth = base_counts[major_base] if major_base != 'del' else 0
                major_freq = major_depth / depth if depth > 0 else 0
            else:
                major_base = 'N'
                major_depth = 0
                major_freq = 0
            
            # Determine genotype
            if major_freq >= 0.8:
                genotype = "HOMOZYGOUS"
            elif major_freq >= 0.2:
                genotype = "HETEROZYGOUS"
            elif depth > 0:
                genotype = "LOW_FREQ"
            else:
                genotype = "NO_COVERAGE"
            
            position_data[pos] = {
                'total_depth': depth,
                'effective_depth': major_depth,
                'ref_base': ref_seq[pos - start] if ref_seq and (pos - start) < len(ref_seq) else ref_base,
                'major_base': major_base,
                'major_freq': major_freq,
                'base_counts': base_counts,
                'genotype': genotype
            }
        
        return position_data
    
    def get_effective_depth_stats(self, bam_file, chrom, start, end, ref_fasta=None):
        """Get effective depth statistics for a region"""
        pos_data = self.calculate_allele_specific_depth(bam_file, chrom, start, end, ref_fasta)
        
        if not pos_data:
            return {
                'mean_effective_depth': 0, 
                'mean_total_depth': 0,
                'homozygous_ratio': 0,
                'heterozygous_ratio': 0,
                'coverage_efficiency': 0,
                'total_bases_analyzed': 0
            }
        
        depths_total = [d['total_depth'] for d in pos_data.values()]
        depths_effective = [d['effective_depth'] for d in pos_data.values()]
        genotypes = [d['genotype'] for d in pos_data.values()]
        
        mean_total = np.mean(depths_total) if depths_total else 0
        mean_effective = np.mean(depths_effective) if depths_effective else 0
        
        return {
            'mean_total_depth': mean_total,
            'mean_effective_depth': mean_effective,
            'std_effective_depth': np.std(depths_effective) if depths_effective else 0,
            'homozygous_ratio': genotypes.count('HOMOZYGOUS') / len(genotypes) if genotypes else 0,
            'heterozygous_ratio': genotypes.count('HETEROZYGOUS') / len(genotypes) if genotypes else 0,
            'low_freq_ratio': genotypes.count('LOW_FREQ') / len(genotypes) if genotypes else 0,
            'coverage_efficiency': (mean_effective / mean_total * 100) if mean_total > 0 else 0,
            'total_bases_analyzed': len(pos_data)
        }
    
    # ============ Quality Assessment Module ============
    
    def index_genome(self, genome_fasta):
        """Index genome"""
        if not os.path.exists(f"{genome_fasta}.fai"):
            self.log("Indexing genome...", "info")
            cmd = f"samtools faidx {genome_fasta}"
            result = subprocess.run(cmd, shell=True, capture_output=True)
            return result.returncode == 0
        return True
    
    def parse_cigar(self, cigar):
        """Parse CIGAR string to get aligned length"""
        align_len = 0
        matches = re.finditer(r'(\d+)([MIDNSHP=X])', cigar)
        for match in matches:
            count = int(match.group(1))
            op = match.group(2)
            if op in 'MX=':
                align_len += count
        return align_len
    
    def get_read_end_from_cigar(self, read_start, cigar):
        """Calculate read end position from CIGAR string"""
        end = read_start
        cigar_parts = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        for length, op in cigar_parts:
            length = int(length)
            if op in 'MDN=X':
                end += length
        return end
    
    def align_hifi_reads(self, genome_fasta, output_bam):
        """
        Align HiFi reads to filled genome using pipeline (no intermediate SAM file)
        Only keep primary alignments (filter out secondary and supplementary)
        """
        if not self.hifi_reads:
            return False
        
        self.log("Aligning HiFi reads to filled genome...", "info")
        self.log(f"Filtering: primary alignments only (no coverage filtering)", "info")
        
        if os.path.exists(output_bam) and os.path.exists(f"{output_bam}.bai"):
            self.log("Using existing BAM file", "info")
            return True
        
        # Use pipeline: minimap2 -> samtools sort -> BAM (no intermediate SAM file)
        cmd = f"""minimap2 -ax map-hifi \
            -Y \
            --secondary=no \
            -t {self.threads_per_alignment} \
            {genome_fasta} \
            {self.hifi_reads} 2> {output_bam}.log | \
            samtools sort -@ {self.threads_per_alignment} -m 4G -o {output_bam} - 2>> {output_bam}.log"""
        
        self.log(f"Running minimap2 with {self.threads_per_alignment} threads (primary alignments only)", "debug")
        self.log("Using pipeline: minimap2 | samtools sort (no intermediate SAM file)", "debug")
        result = subprocess.run(cmd, shell=True, capture_output=True)
        
        if result.returncode == 0 and os.path.exists(output_bam) and os.path.getsize(output_bam) > 0:
            subprocess.run(f"samtools index {output_bam}", shell=True, capture_output=True)
            
            total_cmd = f"samtools view -c {output_bam} 2>/dev/null"
            total_result = subprocess.run(total_cmd, shell=True, capture_output=True, text=True)
            total_reads = int(total_result.stdout.strip()) if total_result.stdout.strip() else 0
            
            self.log(f"Alignment completed: {total_reads} primary alignments kept", "success")
            return True
        
        self.log(f"Minimap2 alignment failed. Check {output_bam}.log for details", "error")
        # Print last few lines of log for debugging
        if os.path.exists(f"{output_bam}.log"):
            with open(f"{output_bam}.log", 'r') as f:
                lines = f.readlines()[-10:]
                self.log("Last 10 lines of log:", "debug")
                for line in lines:
                    self.log(f"  {line.strip()}", "debug")
        return False
    
    def expand_region_for_quality(self, start, end):
        """Expand region: add flanking bases on each side"""
        new_start = max(0, start - self.flank_size)
        new_end = end + self.flank_size
        return new_start, new_end
    
    def check_n_content(self, genome_fasta, chrom, start, end):
        """Check N content in region"""
        region = f"{chrom}:{start}-{end}"
        cmd = f"samtools faidx {genome_fasta} {region} 2>/dev/null | grep -v '^>' | tr -d '\\n'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.stdout:
            sequence = result.stdout.upper()
            n_count = sequence.count('N')
            total_bases = len(sequence)
            n_ratio = n_count / total_bases if total_bases > 0 else 1.0
            return n_count, total_bases, n_ratio
        return 0, 0, 1.0
    
    def detect_zero_coverage_stretches(self, bam_file, chrom, start, end):
        """Detect consecutive zero-coverage stretches at base-pair resolution"""
        region_length = end - start
        
        cmd = f"samtools depth -r {chrom}:{start}-{end} {bam_file} 2>/dev/null | awk '{{print $2}}'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        covered_positions = set()
        if result.stdout:
            covered_positions = set(int(x) for x in result.stdout.strip().split())
        
        zero_stretches = []
        current_stretch_start = None
        current_stretch_length = 0
        zero_positions = []
        
        for pos in range(start, end):
            if pos not in covered_positions:
                zero_positions.append(pos)
                if current_stretch_start is None:
                    current_stretch_start = pos
                    current_stretch_length = 1
                else:
                    current_stretch_length += 1
            else:
                if current_stretch_start is not None:
                    zero_stretches.append({
                        'start': current_stretch_start,
                        'end': pos,
                        'length': current_stretch_length
                    })
                    current_stretch_start = None
                    current_stretch_length = 0
        
        if current_stretch_start is not None:
            zero_stretches.append({
                'start': current_stretch_start,
                'end': end,
                'length': current_stretch_length
            })
        
        total_zero_bases = sum(s['length'] for s in zero_stretches)
        zero_ratio = total_zero_bases / region_length if region_length > 0 else 1.0
        
        flagged_stretches = [s for s in zero_stretches if s['length'] >= self.min_zero_stretch_bp]
        worst_stretch = max(flagged_stretches, key=lambda x: x['length']) if flagged_stretches else None
        
        return {
            'has_zero_stretches': len(flagged_stretches) > 0,
            'all_zero_stretches': zero_stretches,
            'consecutive_stretches_count': len(flagged_stretches),
            'total_zero_bases': total_zero_bases,
            'zero_ratio': zero_ratio,
            'worst_stretch': worst_stretch,
            'zero_positions_count': len(zero_positions),
            'zero_positions': zero_positions
        }
    
    def calculate_low_coverage_density_score(self, bam_file, chrom, start, end, effective_depth_data=None):
        """
        Calculate low coverage density score using effective depth if available
        Penalizes both continuous and fragmented low coverage regions
        Low coverage threshold is now 3x (modified from 5x)
        """
        region_length = end - start
        if region_length <= 0:
            return 0, {}
        
        # Use effective depth if provided, otherwise use raw depth
        if effective_depth_data:
            depth_dict = {pos: effective_depth_data[pos]['effective_depth'] 
                         for pos in effective_depth_data}
        else:
            cmd = f"samtools depth -r {chrom}:{start}-{end} {bam_file} 2>/dev/null | awk '{{print $2, $3}}'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            depth_dict = {}
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            depth_dict[int(parts[0])] = float(parts[1])
                        except (ValueError, IndexError):
                            continue
        
        low_stretches = []
        current_stretch = None
        
        for pos in range(start, end):
            depth = depth_dict.get(pos, 0)
            if depth < self.low_coverage_threshold:  # Now default is 3
                if current_stretch is None:
                    current_stretch = {'start': pos, 'length': 1, 'depths': [depth], 'min_depth': depth}
                else:
                    current_stretch['length'] += 1
                    current_stretch['depths'].append(depth)
                    current_stretch['min_depth'] = min(current_stretch['min_depth'], depth)
            else:
                if current_stretch is not None:
                    low_stretches.append(current_stretch)
                    current_stretch = None
        
        if current_stretch is not None:
            low_stretches.append(current_stretch)
        
        if not low_stretches:
            return 100, {'low_stretches_count': 0, 'total_low_bases': 0, 'low_density': 0, 
                         'fragmentation_penalty': 0, 'depth_penalty': 1.0, 'very_low_ratio': 0}
        
        total_low_bases = sum(s['length'] for s in low_stretches)
        low_density = total_low_bases / region_length
        stretch_count = len(low_stretches)
        
        # Stronger fragmentation penalty
        fragmentation_penalty = min(0.5, stretch_count / 100)
        
        # Stricter base score based on low density
        if low_density < 0.01:
            base_score = 100
        elif low_density < 0.03:
            base_score = 80
        elif low_density < 0.05:
            base_score = 60
        elif low_density < 0.10:
            base_score = 35
        elif low_density < 0.20:
            base_score = 15
        elif low_density < 0.35:
            base_score = 5
        else:
            base_score = 0
        
        # Calculate very low depth penalty (<1x)
        very_low_bases = sum(1 for s in low_stretches for d in s['depths'] if d < 1)
        very_low_ratio = very_low_bases / total_low_bases if total_low_bases > 0 else 0
        
        if very_low_ratio > 0.5:
            depth_penalty = 0.6
        elif very_low_ratio > 0.3:
            depth_penalty = 0.75
        elif very_low_ratio > 0.1:
            depth_penalty = 0.85
        else:
            depth_penalty = 0.95
        
        density_score = base_score * (1 - fragmentation_penalty) * depth_penalty
        
        if stretch_count > 10 and self.verbose:
            self.log(f"      Fragmentation: {stretch_count} low-coverage stretches (density: {low_density:.2%})", "debug")
        
        return max(0, min(100, density_score)), {
            'low_stretches_count': stretch_count,
            'total_low_bases': total_low_bases,
            'low_density': low_density,
            'fragmentation_penalty': fragmentation_penalty,
            'depth_penalty': depth_penalty,
            'very_low_ratio': very_low_ratio,
            'avg_low_depth': np.mean([d for s in low_stretches for d in s['depths']]) if low_stretches else 0
        }
    
    def calculate_region_quality_index(self, bam_file, chrom, start, end, ref_fasta=None):
        """
        Calculate region quality index (R-AQI) using effective depth
        
        Weight distribution:
        - Depth adequacy: 20%
        - Coverage uniformity: 15%
        - Coverage completeness: 10%
        - Low coverage density: 25%
        - Mapping quality: 15%
        - Read support: 15%
        """
        if start >= end:
            return 0
        
        region_length = end - start
        
        # Calculate allele-specific depth
        pos_data = self.calculate_allele_specific_depth(bam_file, chrom, start, end, ref_fasta)
        
        if not pos_data:
            return 0
        
        depths_total = [d['total_depth'] for d in pos_data.values()]
        depths_effective = [d['effective_depth'] for d in pos_data.values()]
        
        # Use effective depth for quality metrics
        mean_depth = np.mean(depths_effective) if depths_effective else 0
        std_depth = np.std(depths_effective) if depths_effective else 0
        cv_depth = std_depth / mean_depth if mean_depth > 0 else 1.0
        
        # 1. Coverage uniformity score (15% weight)
        uniformity_score = max(0, 1 - cv_depth) * 100
        
        # 2. Depth adequacy score (20% weight) - using effective depth
        if mean_depth >= 15:
            depth_score = 100
        elif mean_depth >= 10:
            depth_score = 80
        elif mean_depth >= 5:
            depth_score = 60
        elif mean_depth >= 3:
            depth_score = 40
        else:
            depth_score = 20 * (mean_depth / 3) if mean_depth > 0 else 0
        
        # Apply heterozygosity penalty
        het_ratio = sum(1 for d in pos_data.values() if d['genotype'] == 'HETEROZYGOUS') / len(pos_data)
        if het_ratio > 0.3:
            depth_score = depth_score * 0.7
        elif het_ratio > 0.1:
            depth_score = depth_score * 0.85
        
        # 3. Coverage completeness (10% weight)
        covered_positions = sum(1 for d in pos_data.values() if d['effective_depth'] > 0)
        covered_ratio = covered_positions / region_length if region_length > 0 else 0
        
        # 4. Low coverage density score (25% weight)
        low_cov_score, low_cov_metrics = self.calculate_low_coverage_density_score(
            bam_file, chrom, start, end, pos_data
        )
        
        # 5. Mapping quality score (15% weight)
        cmd = f"samtools view {bam_file} {chrom}:{start}-{end} 2>/dev/null | awk '{{print $5}}'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        mapq_scores = [float(x) for x in result.stdout.strip().split() if x]
        
        if mapq_scores:
            mean_mapq = np.mean(mapq_scores)
            mapq_score = min(100, (mean_mapq / 60) * 100)
        else:
            mapq_score = 0
        
        # 6. Read support score (15% weight) - using effective depth data
        if depths_effective and depths_total:
            avg_effective = np.mean(depths_effective)
            avg_total = np.mean(depths_total)
            identity_score = (avg_effective / avg_total * 100) if avg_total > 0 else 0
        else:
            identity_score = 0
        
        # Calculate R-AQI
        r_aqi = (
            depth_score * 0.20 +
            uniformity_score * 0.15 +
            covered_ratio * 100 * 0.10 +
            low_cov_score * 0.25 +
            mapq_score * 0.15 +
            identity_score * 0.15
        )
        
        return r_aqi
    
    def get_coverage_variability(self, bam_file, chrom, start, end):
        """Calculate coverage variability"""
        variability_window = 1000
        window_starts = range(start, end, variability_window)
        window_depths = []
        
        for w_start in window_starts:
            w_end = min(w_start + variability_window, end)
            cmd = f"samtools depth -r {chrom}:{w_start}-{w_end} {bam_file} 2>/dev/null | \
                   awk '{{sum+=$3; count++}} END {{if(count>0) print sum/count; else print 0}}'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            depth = float(result.stdout.strip()) if result.stdout.strip() else 0
            window_depths.append(depth)
        
        if not window_depths:
            return {'cv': 1.0, 'has_outliers': False, 'zero_windows_count': 0}
        
        window_depths = np.array(window_depths)
        mean_depth = np.mean(window_depths)
        std_depth = np.std(window_depths)
        cv = std_depth / mean_depth if mean_depth > 0 else 1.0
        zero_windows = [i for i, d in enumerate(window_depths) if d == 0]
        
        return {'cv': cv, 'has_outliers': len(zero_windows) > 0, 'zero_windows_count': len(zero_windows)}
    
    def get_identity_stats(self, bam_file, chrom, start, end):
        """Calculate identity statistics for reads in region"""
        cmd = f"samtools view {bam_file} {chrom}:{start}-{end} 2>/dev/null | \
               grep -o 'NM:i:[0-9]*' | awk -F':' '{{print $3}}'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        nm_values = [int(x) for x in result.stdout.strip().split() if x]
        
        if not nm_values:
            return {'mean_identity': 0, 'min_identity': 0, 'std_identity': 0,
                    'has_low_identity': False, 'low_identity_ratio': 0, 'read_count': 0,
                    'filtered_reads': 0}
        
        cmd_len = f"samtools view {bam_file} {chrom}:{start}-{end} 2>/dev/null | \
                   awk '{{print length($10)}}'"
        result_len = subprocess.run(cmd_len, shell=True, capture_output=True, text=True)
        read_lengths = [int(x) for x in result_len.stdout.strip().split() if x]
        
        identities = []
        filtered_count = 0
        min_read_length = 50
        
        for nm, length in zip(nm_values, read_lengths):
            if length >= min_read_length:
                identity = 100 - (nm / length * 100)
                if identity < 0:
                    if self.verbose and filtered_count < 5:
                        self.log(f"Warning: Negative identity ({identity:.1f}) for read: NM={nm}, length={length}", "debug")
                    identity = max(0, identity)
                identities.append(identity)
            else:
                filtered_count += 1
                if self.verbose and filtered_count <= 5:
                    self.log(f"Filtering short read (length={length}bp, NM={nm}) - below minimum {min_read_length}bp", "debug")
        
        if not identities:
            return {'mean_identity': 0, 'min_identity': 0, 'std_identity': 0,
                    'has_low_identity': False, 'low_identity_ratio': 0, 'read_count': 0,
                    'filtered_reads': filtered_count}
        
        identities = np.array(identities)
        return {
            'mean_identity': np.mean(identities),
            'min_identity': np.min(identities),
            'std_identity': np.std(identities),
            'has_low_identity': np.mean(identities) < 95,
            'low_identity_ratio': sum(1 for i in identities if i < 95) / len(identities),
            'read_count': len(identities),
            'filtered_reads': filtered_count
        }
    
    def calculate_boundary_consistency(self, bam_file, chrom, boundary_pos, gap_start, gap_end, window=5000):
        """Calculate boundary consistency"""
        boundary_region = f"{chrom}:{max(0, boundary_pos - window)}-{boundary_pos + window}"
        
        cmd = f"samtools view {bam_file} {boundary_region} 2>/dev/null | awk '{{print $4, length($10), $5}}'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        spanning_reads = 0
        total_reads = 0
        spanning_reads_high_quality = 0
        
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        read_start = int(parts[0])
                        read_len = int(parts[1])
                        mapq = int(parts[2])
                        read_end = read_start + read_len
                        if read_len >= 50:
                            total_reads += 1
                            if read_start < boundary_pos < read_end:
                                spanning_reads += 1
                                if mapq >= 30:
                                    spanning_reads_high_quality += 1
                    except (ValueError, IndexError):
                        continue
        
        if total_reads == 0:
            return 0
        
        ratio = spanning_reads / total_reads
        high_quality_ratio = spanning_reads_high_quality / total_reads if total_reads > 0 else 0
        
        if ratio >= self.boundary_good_ratio:
            ratio_score = 85 + min(15, (ratio - self.boundary_good_ratio) / (1 - self.boundary_good_ratio) * 15)
        elif ratio >= 0.15:
            ratio_score = 65 + (ratio - 0.15) / (self.boundary_good_ratio - 0.15) * 20
        elif ratio >= 0.08:
            ratio_score = 45 + (ratio - 0.08) / 0.07 * 20
        elif ratio >= 0.03:
            ratio_score = 25 + (ratio - 0.03) / 0.05 * 20
        else:
            ratio_score = ratio / 0.03 * 25
        
        quality_boost = min(10, high_quality_ratio * 20)
        ratio_score = min(100, ratio_score + quality_boost)
        
        if spanning_reads >= self.boundary_min_spanning:
            count_score = 100
        elif spanning_reads >= 5:
            count_score = 70 + (spanning_reads - 5) / (self.boundary_min_spanning - 5) * 30
        elif spanning_reads >= 2:
            count_score = 40 + (spanning_reads - 2) / 3 * 30
        elif spanning_reads == 1:
            count_score = 25
        else:
            count_score = 0
        
        consistency = ratio_score * 0.6 + count_score * 0.4
        
        region_size = gap_end - gap_start
        if region_size < 5000 and spanning_reads_high_quality >= 1:
            consistency = max(consistency, 60)
        elif region_size < 5000 and spanning_reads >= 1:
            consistency = max(consistency, 45)
        
        return min(100, consistency)
    
    def calculate_mapping_aqi(self, filled_region, ref_fasta=None):
        """Calculate comprehensive mapping quality index (M-AQI) using effective depth"""
        chrom = filled_region['filled_chrom']
        start = filled_region['filled_start']
        end = filled_region['filled_end']
        bam_file = self.bam_file
        
        if not bam_file or not os.path.exists(bam_file):
            return {'m_aqi': 50 if filled_region.get('n_ratio', 1) < 0.05 else 0,
                    'internal_quality': 50, 'left_consistency': 50, 'right_consistency': 50}
        
        internal_quality = self.calculate_region_quality_index(bam_file, chrom, start, end, ref_fasta)
        left_consistency = self.calculate_boundary_consistency(bam_file, chrom, start, start, end)
        right_consistency = self.calculate_boundary_consistency(bam_file, chrom, end, start, end)
        
        m_aqi = internal_quality * 0.50 + left_consistency * 0.25 + right_consistency * 0.25
        
        return {'m_aqi': m_aqi, 'internal_quality': internal_quality,
                'left_consistency': left_consistency, 'right_consistency': right_consistency}
    
    # ============ Simplified Coverage Plot Module with Downsampling ============
    
    def get_downsample_step(self, region_length):
        """Determine downsampling step based on region length"""
        if region_length > 500000:      # >500kb
            return 500, 50
        elif region_length > 200000:    # 200-500kb
            return 200, 60
        elif region_length > 100000:    # 100-200kb
            return 100, 72
        elif region_length > 50000:     # 50-100kb
            return 50, 86
        elif region_length > 20000:     # 20-50kb
            return 20, 100
        elif region_length > 10000:     # 10-20kb
            return 10, 120
        else:
            return 1, 150  # No downsampling for short regions
    
    def get_sampled_depth_data(self, bam_file, chrom, start, end, step):
        """Get depth data with downsampling"""
        positions = []
        depths = []
        
        # Use samtools depth to get coverage
        cmd = f"samtools depth -r {chrom}:{start}-{end} {bam_file} 2>/dev/null"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        # Build depth dictionary
        depth_dict = {}
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        pos = int(parts[1])
                        depth_dict[pos] = float(parts[2])
                    except:
                        pass
        
        # Sample at step intervals
        for pos in range(start, end + 1, step):
            positions.append(pos)
            depths.append(depth_dict.get(pos, 0))
        
        return positions, depths
    
    def generate_coverage_pdf_plot(self, bam_file, chrom, start, end, gap_id, output_dir, ref_fasta=None):
        """
        Generate single PDF with downsampling for large regions
        """
        region_length = end - start
        
        if region_length <= 0:
            return None
        
        # Skip extremely large regions
        if region_length > self.max_pdf_length:
            self.log(f"    Skipping PDF for {gap_id}: region too large ({region_length:,} bp > {self.max_pdf_length:,} bp)", "warning")
            return None
        
        # Determine downsampling step
        step, dpi = self.get_downsample_step(region_length)
        
        if step > 1:
            self.log(f"    Downsampling {gap_id}: {region_length:,} bp → ~{region_length//step} points (step={step})", "debug")
        
        with self.lock:
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                import matplotlib.patches as patches
            except ImportError:
                self.log(f"matplotlib not available. Install with: pip install matplotlib", "warning")
                return None
            
            plot_file = os.path.join(output_dir, f"{gap_id}_coverage.pdf")
            
            # Get sampled depth data
            positions, depths = self.get_sampled_depth_data(bam_file, chrom, start, end, step)
            
            if not positions or max(depths) == 0:
                self.log(f"    No depth data for {gap_id}", "debug")
                return None
            
            # Get read alignment data (limited to 30 reads for clarity)
            cmd_reads = f"samtools view {bam_file} {chrom}:{start}-{end} 2>/dev/null | \
                         awk -v OFS='\\t' '{{print $1, $4, $6, length($10), $5, $2}}' | head -30"
            result = subprocess.run(cmd_reads, shell=True, capture_output=True, text=True)
            
            reads_data = []
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        try:
                            read_start = int(parts[1])
                            cigar = parts[2]
                            mapq = int(parts[4])
                            flag = int(parts[5]) if len(parts) > 5 else 0
                            
                            strand = '-' if (flag & 16) else '+'
                            read_end = self.get_read_end_from_cigar(read_start, cigar)
                            
                            if read_end > start and read_start < end:
                                reads_data.append({
                                    'start': max(read_start, start),
                                    'end': min(read_end, end),
                                    'mapq': mapq,
                                    'strand': strand
                                })
                        except (ValueError, IndexError):
                            continue
            
            # Sort reads by start position
            reads_data.sort(key=lambda x: x['start'])
            
            max_depth = max(depths) if depths else 10
            
            # Create figure with 2 subplots for simplicity (faster)
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
            
            # Panel 1: Coverage depth
            ax1.fill_between(positions, depths, 0, alpha=0.3, color='blue')
            ax1.plot(positions, depths, 'b-', linewidth=0.5, alpha=0.7)
            ax1.set_ylabel('Coverage Depth (X)', fontsize=10, fontweight='bold')
            title = f'{gap_id}: {chrom}:{start:,}-{end:,} ({region_length:,} bp)'
            if step > 1:
                title += f' [downsampled: {step}bp step]'
            ax1.set_title(title, fontsize=9)
            ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            ax1.axhline(y=self.low_coverage_threshold, color='orange', linestyle='--', 
                       linewidth=0.8, alpha=0.7, label=f'Threshold ({self.low_coverage_threshold}x)')
            ax1.legend(loc='upper right', fontsize=8)
            
            # Panel 2: Individual reads (simplified)
            if reads_data:
                y_pos = 0
                y_step = 0.6
                
                for read in reads_data:
                    read_start = read['start']
                    read_end = read['end']
                    mapq = read['mapq']
                    strand = read['strand']
                    
                    # Color by MAPQ
                    if mapq >= 60:
                        color = 'darkblue' if strand == '+' else 'darkred'
                    elif mapq >= 30:
                        color = 'steelblue' if strand == '+' else 'salmon'
                    elif mapq >= 10:
                        color = 'lightblue' if strand == '+' else 'lightcoral'
                    else:
                        color = 'lightgray'
                    
                    rect = patches.Rectangle(
                        (read_start, y_pos), read_end - read_start, 0.4,
                        linewidth=0.5, edgecolor=color, facecolor=color, alpha=0.6
                    )
                    ax2.add_patch(rect)
                    y_pos += y_step
                
                ax2.set_ylim(-0.5, y_pos + 0.5)
                ax2.set_ylabel('Individual Reads', fontsize=10, fontweight='bold')
                ax2.set_title(f'Read Alignments ({len(reads_data)} reads, color by MAPQ)', fontsize=9)
            else:
                ax2.text(0.5, 0.5, 'No reads in this region', transform=ax2.transAxes,
                        ha='center', va='center', fontsize=12, style='italic')
                ax2.set_ylabel('Individual Reads', fontsize=10, fontweight='bold')
            
            ax2.set_xlim(start, end)
            ax2.set_xlabel(f'Position on {chrom} (bp)', fontsize=10)
            ax2.tick_params(axis='x', rotation=45, labelsize=8)
            ax2.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='x')
            
            plt.tight_layout()
            
            os.makedirs(output_dir, exist_ok=True)
            plt.savefig(plot_file, dpi=dpi, bbox_inches='tight', format='pdf')
            plt.close()
            
            self.log(f"    Coverage plot saved: {gap_id}_coverage.pdf ({len(positions)} points, step={step})", "debug")
            return plot_file
    
    # ============ Quality Assessment Main Function ============
    
    def assess_filled_region(self, filled_region, genome_fasta, bam_file, ref_fasta=None):
        """Assess quality of a single filled region using allele-specific depth"""
        if not bam_file or not os.path.exists(bam_file):
            return None
        
        chrom = filled_region['filled_chrom']
        start = filled_region['filled_start']
        end = filled_region['filled_end']
        
        if start == 'N/A' or end == 'N/A' or start >= end:
            return None
        
        original_length = end - start
        
        expanded_start, expanded_end = self.expand_region_for_quality(start, end)
        region_length = expanded_end - expanded_start
        
        n_count, total_bases, n_ratio = self.check_n_content(genome_fasta, chrom, start, end)
        
        N_THRESHOLD = 0.05
        
        if n_ratio > N_THRESHOLD:
            return {
                'filled_length': original_length, 'n_count': n_count, 'n_ratio': n_ratio,
                'total_score': 0, 'quality_grade': 'FAILED', 'classification': 'MAPPED_WITH_N',
                'is_clean_fill': False, 'is_high_quality': False,
                'expanded_start': expanded_start, 'expanded_end': expanded_end, 'region_length': region_length,
                'coverage_rate': 0, 'mean_depth': 0, 'mean_mapq': 0, 'mean_identity': 0,
                'm_aqi': 0, 'internal_quality': 0, 'left_consistency': 0, 'right_consistency': 0,
                'cv_depth': 1.0, 'has_coverage_gaps': True, 'has_zero_stretches': True,
                'zero_stretches_count': 1, 'zero_ratio_in_filled': 1.0, 'worst_zero_stretch_length': original_length,
                'zero_positions_count': 0, 'zero_positions': [], 'filtered_reads': 0,
                'low_cov_stretches_count': 0, 'low_cov_density': 1.0, 'low_cov_score': 0
            }
        
        zero_result = self.detect_zero_coverage_stretches(bam_file, chrom, start, end)
        maqi_result = self.calculate_mapping_aqi(filled_region, ref_fasta)
        variability = self.get_coverage_variability(bam_file, chrom, start, end)
        identity_stats = self.get_identity_stats(bam_file, chrom, start, end)
        
        # Calculate effective depth statistics
        eff_depth_stats = self.get_effective_depth_stats(bam_file, chrom, start, end, ref_fasta)
        
        # Use effective depth for coverage metrics
        mean_depth = eff_depth_stats.get('mean_effective_depth', 0) if eff_depth_stats else 0
        coverage_rate = eff_depth_stats.get('coverage_efficiency', 0) / 100 if eff_depth_stats else 0
        covered_bases = eff_depth_stats.get('total_bases_analyzed', 0) if eff_depth_stats else 0
        
        # Calculate low coverage score with effective depth
        pos_data = self.calculate_allele_specific_depth(bam_file, chrom, start, end, ref_fasta)
        low_cov_score, low_cov_metrics = self.calculate_low_coverage_density_score(
            bam_file, chrom, start, end, pos_data
        )
        
        cmd = f"samtools view {bam_file} {chrom}:{expanded_start}-{expanded_end} 2>/dev/null | awk '{{sum+=$5; count++}} END {{if(count>0) print sum/count; else print 0}}'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        mean_mapq = float(result.stdout.strip()) if result.stdout.strip() else 0
        
        m_aqi = maqi_result['m_aqi']
        
        # Additional penalty for low coverage density
        if low_cov_score < 50:
            penalty = 0.6 + (low_cov_score / 100) * 0.4
            m_aqi = m_aqi * penalty
        
        if zero_result['has_zero_stretches']:
            worst_stretch = zero_result['worst_stretch']
            if worst_stretch:
                stretch_len = worst_stretch['length']
                
                if stretch_len >= 1000:
                    penalty = 0.1
                    self.log(f"      Severe: {stretch_len}bp zero coverage stretch", "warning")
                elif stretch_len >= 100:
                    penalty = 0.4
                    self.log(f"      Significant zero coverage: {stretch_len}bp", "warning")
                elif stretch_len >= 10:
                    penalty = 0.7
                else:
                    penalty = 0.85
                
                m_aqi = m_aqi * penalty
        
        # Heterozygosity penalty
        het_ratio = eff_depth_stats.get('heterozygous_ratio', 0) if eff_depth_stats else 0
        if het_ratio > 0.3:
            het_penalty = 0.7
            self.log(f"      High heterozygosity: {het_ratio:.1%} -> penalty {het_penalty:.2f}", "debug")
        elif het_ratio > 0.1:
            het_penalty = 0.85
        else:
            het_penalty = 1.0
        
        m_aqi = m_aqi * het_penalty
        
        # Revised quality grading
        if m_aqi >= 90:
            grade = "EXCELLENT"
            is_high_quality = True
        elif m_aqi >= 75:
            grade = "GOOD"
            is_high_quality = True
        elif m_aqi >= 60:
            grade = "FAIR"
            is_high_quality = False
        elif m_aqi >= 40:
            grade = "POOR"
            is_high_quality = False
        else:
            grade = "FAILED"
            is_high_quality = False
        
        return {
            'filled_length': original_length,
            'expanded_start': expanded_start, 'expanded_end': expanded_end, 'region_length': region_length,
            'coverage_rate': coverage_rate, 'mean_depth': mean_depth, 
            'covered_bases': covered_bases,
            'mean_mapq': mean_mapq, 'mean_identity': identity_stats['mean_identity'],
            'n_count': n_count, 'n_ratio': n_ratio, 'total_score': m_aqi, 'quality_grade': grade,
            'classification': 'MAPPED_CLEAN', 'is_clean_fill': True, 'is_high_quality': is_high_quality,
            'm_aqi': m_aqi, 'internal_quality': maqi_result['internal_quality'],
            'left_consistency': maqi_result['left_consistency'], 'right_consistency': maqi_result['right_consistency'],
            'cv_depth': variability['cv'], 'has_coverage_gaps': variability['zero_windows_count'] > 0,
            'low_identity_ratio': identity_stats['low_identity_ratio'],
            'has_zero_stretches': zero_result['has_zero_stretches'],
            'zero_stretches_count': zero_result['consecutive_stretches_count'],
            'total_zero_bases': zero_result['total_zero_bases'],
            'zero_ratio_in_filled': zero_result['zero_ratio'],
            'worst_zero_stretch_length': zero_result['worst_stretch']['length'] if zero_result['worst_stretch'] else 0,
            'zero_positions_count': zero_result['zero_positions_count'],
            'zero_positions': zero_result['zero_positions'][:100],
            'filtered_reads': identity_stats.get('filtered_reads', 0),
            'low_cov_stretches_count': low_cov_metrics.get('low_stretches_count', 0),
            'low_cov_density': low_cov_metrics.get('low_density', 0),
            'low_cov_score': low_cov_score,
            'heterozygous_ratio': het_ratio,
            'coverage_efficiency': eff_depth_stats.get('coverage_efficiency', 0) if eff_depth_stats else 0
        }
    
    def process_chromosome(self, chrom, gaps, before_fasta, after_fasta, temp_dir, output_dir,
                          plot_coverage=False):
        """Process a single chromosome"""
        self.log(f"  Processing {chrom} ({len(gaps)} gaps)...")
        
        before_chrom_fasta = os.path.join(temp_dir, f"{chrom}_before.fa")
        after_chrom_fasta = os.path.join(temp_dir, f"{chrom}_after.fa")
        
        if not self.extract_chromosome(before_fasta, chrom, before_chrom_fasta):
            return chrom, [], f"Before chromosome not found"
        
        if not self.extract_chromosome(after_fasta, chrom, after_chrom_fasta):
            return chrom, [], f"After chromosome not found"
        
        prefix = os.path.join(temp_dir, f"{chrom}_align")
        alignments, error = self.run_nucmer_alignment(before_chrom_fasta, after_chrom_fasta, prefix)
        
        if error:
            self.log(f"    Nucmer failed: {error}", "warning")
            return chrom, [], error
        
        self.log(f"    Found {len(alignments)} alignments")
        
        mapped_gaps = []
        for gap in gaps:
            gap_start = gap['start']
            gap_end = gap['end']
            
            left_aln, right_aln, containing, _, _ = self.find_left_right_alignments(
                alignments, chrom, gap_start, gap_end, max_distance=500000
            )
            
            if containing:
                containing_aln = containing[0]
                mapping = self.map_gap_inside_alignment(gap, containing_aln)
                status = 'MAPPED'
            elif left_aln and right_aln:
                mapping = self.map_gap_between_alignments(gap, left_aln, right_aln)
                status = 'MAPPED' if mapping['length'] > 0 else 'MAPPING_FAILED_ZERO_LENGTH'
            else:
                mapping = None
                status = 'NOT_MAPPED_NO_FLANKING'
            
            if mapping and status == 'MAPPED' and mapping['length'] > 0:
                quality = None
                if self.bam_file:
                    quality = self.assess_filled_region(
                        {'filled_chrom': mapping['chrom'],
                         'filled_start': mapping['start'],
                         'filled_end': mapping['end'],
                         'original_length': gap['length']},
                        after_fasta,
                        self.bam_file,
                        ref_fasta=after_fasta
                    )
                
                result = {
                    'gap_id': gap['gap_id'],
                    'original_chrom': gap['chrom'],
                    'original_start': gap['start'],
                    'original_end': gap['end'],
                    'original_length': gap['length'],
                    'filled_chrom': mapping['chrom'],
                    'filled_start': mapping['start'],
                    'filled_end': mapping['end'],
                    'filled_length': mapping['length'],
                    'identity': mapping['identity'],
                    'method': mapping['method'],
                    'status': status
                }
                
                if quality:
                    result.update({
                        'expanded_start': quality['expanded_start'],
                        'expanded_end': quality['expanded_end'],
                        'region_length': quality['region_length'],
                        'coverage_rate': quality['coverage_rate'],
                        'mean_depth': quality['mean_depth'],
                        'mean_mapq': quality['mean_mapq'],
                        'mean_identity': quality['mean_identity'],
                        'n_count': quality['n_count'],
                        'n_ratio': quality['n_ratio'],
                        'quality_score': quality['total_score'],
                        'quality_grade': quality['quality_grade'],
                        'classification': quality['classification'],
                        'is_clean_fill': quality['is_clean_fill'],
                        'is_high_quality': quality['is_high_quality'],
                        'm_aqi': quality.get('m_aqi', 0),
                        'internal_quality': quality.get('internal_quality', 0),
                        'left_consistency': quality.get('left_consistency', 0),
                        'right_consistency': quality.get('right_consistency', 0),
                        'cv_depth': quality.get('cv_depth', 1.0),
                        'has_coverage_gaps': quality.get('has_coverage_gaps', False),
                        'low_identity_ratio': quality.get('low_identity_ratio', 0),
                        'has_zero_stretches': quality.get('has_zero_stretches', False),
                        'zero_stretches_count': quality.get('zero_stretches_count', 0),
                        'total_zero_bases': quality.get('total_zero_bases', 0),
                        'zero_ratio_in_filled': quality.get('zero_ratio_in_filled', 0),
                        'worst_zero_stretch_length': quality.get('worst_zero_stretch_length', 0),
                        'zero_positions_count': quality.get('zero_positions_count', 0),
                        'filtered_reads': quality.get('filtered_reads', 0),
                        'low_cov_stretches_count': quality.get('low_cov_stretches_count', 0),
                        'low_cov_density': quality.get('low_cov_density', 0),
                        'low_cov_score': quality.get('low_cov_score', 0),
                        'heterozygous_ratio': quality.get('heterozygous_ratio', 0),
                        'coverage_efficiency': quality.get('coverage_efficiency', 0)
                    })
                    
                    # Generate coverage plot (single PDF with 4 panels)
                    if plot_coverage and self.bam_file:
                        if mapping['chrom'] != 'N/A' and mapping['start'] != 'N/A' and mapping['end'] != 'N/A':
                            plot_dir = os.path.join(output_dir, "coverage_plots")
                            os.makedirs(plot_dir, exist_ok=True)
                            
                            flank = self.flank_size
                            plot_start = max(0, mapping['start'] - flank)
                            plot_end = mapping['end'] + flank
                            
                            self.generate_coverage_pdf_plot(
                                self.bam_file,
                                mapping['chrom'],
                                plot_start,
                                plot_end,
                                gap['gap_id'],
                                plot_dir,
                                ref_fasta=after_fasta
                            )
                    
                    with self.lock:
                        if quality['classification'] == 'MAPPED_WITH_N':
                            self.mapped_with_n_count += 1
                            self.mapped_with_n_gaps.append(gap['gap_id'])
                            self.log(f"      X {gap['gap_id']}: MAPPED_WITH_N - N ratio: {quality['n_ratio']:.1%} [FAILED]", "warning")
                        elif quality['classification'] == 'MAPPED_CLEAN':
                            self.mapped_clean_count += 1
                            if quality['quality_grade'] == 'EXCELLENT':
                                self.high_quality_count += 1
                                self.log(f"      ✓ {gap['gap_id']}: EXCELLENT - M-AQI: {quality['m_aqi']:.1f} (eff_depth: {quality['mean_depth']:.1f}x)", "success")
                            elif quality['quality_grade'] == 'GOOD':
                                self.good_quality_count += 1
                                self.log(f"      ✓ {gap['gap_id']}: GOOD - M-AQI: {quality['m_aqi']:.1f} (eff_depth: {quality['mean_depth']:.1f}x)", "info")
                            elif quality['quality_grade'] == 'FAIR':
                                self.fair_quality_count += 1
                                self.log(f"      ⚠ {gap['gap_id']}: FAIR - M-AQI: {quality['m_aqi']:.1f} (het: {quality.get('heterozygous_ratio', 0):.1%})", "warning")
                            elif quality['quality_grade'] == 'POOR':
                                self.poor_quality_count += 1
                                self.log(f"      ⚠ {gap['gap_id']}: POOR - M-AQI: {quality['m_aqi']:.1f} (het: {quality.get('heterozygous_ratio', 0):.1%})", "warning")
                            else:
                                self.log(f"      ? {gap['gap_id']}: {quality['quality_grade']} - M-AQI: {quality['m_aqi']:.1f}", "warning")
                else:
                    result.update({
                        'classification': 'MAPPED_NO_QUALITY',
                        'is_clean_fill': None,
                        'is_high_quality': None
                    })
                    with self.lock:
                        self.mapped_clean_count += 1
                    self.log(f"      ✓ {gap['gap_id']}: MAPPED (no quality data)", "success")
                
                mapped_gaps.append(result)
            else:
                with self.lock:
                    self.unmapped_count += 1
                    self.unmapped_gaps.append(gap['gap_id'])
                
                result = {
                    'gap_id': gap['gap_id'],
                    'original_chrom': gap['chrom'],
                    'original_start': gap['start'],
                    'original_end': gap['end'],
                    'original_length': gap['length'],
                    'filled_chrom': 'N/A',
                    'filled_start': 'N/A',
                    'filled_end': 'N/A',
                    'filled_length': 0,
                    'identity': 0,
                    'method': 'N/A',
                    'status': status,
                    'classification': 'UNMAPPED',
                    'is_clean_fill': False,
                    'is_high_quality': False
                }
                mapped_gaps.append(result)
                self.log(f"      ? {gap['gap_id']}: UNMAPPED - {status}", "warning")
        
        if mapped_gaps:
            chrom_output = os.path.join(output_dir, f"{chrom}_mapping_with_quality.tsv")
            with open(chrom_output, 'w') as f:
                f.write("gap_id\toriginal_chrom\toriginal_start\toriginal_end\toriginal_length\t"
                       "filled_chrom\tfilled_start\tfilled_end\tfilled_length\tidentity\tmethod\tstatus\tclassification\t"
                       "n_count\tn_ratio\tcoverage_rate\tmean_depth\tmean_mapq\tmean_identity\t"
                       "quality_score\tquality_grade\tm_aqi\tinternal_quality\tleft_consistency\t"
                       "right_consistency\tcv_depth\thas_coverage_gaps\tlow_identity_ratio\t"
                       "has_zero_stretches\tzero_stretches_count\ttotal_zero_bases\tzero_ratio_in_filled\t"
                       "worst_zero_stretch_length\tzero_positions_count\tfiltered_reads\t"
                       "low_cov_stretches_count\tlow_cov_density\tlow_cov_score\theterozygous_ratio\tcoverage_efficiency\n")
                for r in mapped_gaps:
                    f.write(f"{r['gap_id']}\t{r['original_chrom']}\t{r['original_start']}\t{r['original_end']}\t{r['original_length']}\t")
                    f.write(f"{r['filled_chrom']}\t{r['filled_start']}\t{r['filled_end']}\t{r['filled_length']}\t")
                    f.write(f"{r['identity']:.1f}\t{r['method']}\t{r['status']}\t{r.get('classification', 'N/A')}\t")
                    f.write(f"{r.get('n_count', 'N/A')}\t{r.get('n_ratio', 'N/A')}\t")
                    f.write(f"{r.get('coverage_rate', 'N/A')}\t{r.get('mean_depth', 'N/A')}\t")
                    f.write(f"{r.get('mean_mapq', 'N/A')}\t{r.get('mean_identity', 'N/A')}\t")
                    f.write(f"{r.get('quality_score', 'N/A')}\t{r.get('quality_grade', 'N/A')}\t")
                    f.write(f"{r.get('m_aqi', 'N/A')}\t{r.get('internal_quality', 'N/A')}\t")
                    f.write(f"{r.get('left_consistency', 'N/A')}\t{r.get('right_consistency', 'N/A')}\t")
                    f.write(f"{r.get('cv_depth', 'N/A')}\t{r.get('has_coverage_gaps', 'N/A')}\t")
                    f.write(f"{r.get('low_identity_ratio', 'N/A')}\t")
                    f.write(f"{r.get('has_zero_stretches', 'N/A')}\t{r.get('zero_stretches_count', 'N/A')}\t")
                    f.write(f"{r.get('total_zero_bases', 'N/A')}\t{r.get('zero_ratio_in_filled', 'N/A')}\t")
                    f.write(f"{r.get('worst_zero_stretch_length', 'N/A')}\t")
                    f.write(f"{r.get('zero_positions_count', 'N/A')}\t")
                    f.write(f"{r.get('filtered_reads', 'N/A')}\t")
                    f.write(f"{r.get('low_cov_stretches_count', 'N/A')}\t")
                    f.write(f"{r.get('low_cov_density', 'N/A')}\t")
                    f.write(f"{r.get('low_cov_score', 'N/A')}\t")
                    f.write(f"{r.get('heterozygous_ratio', 'N/A')}\t")
                    f.write(f"{r.get('coverage_efficiency', 'N/A')}\n")
        
        with self.lock:
            self.completed += 1
            self.log(f"  Completed {self.completed}/{self.total} chromosomes", "info")
        
        return chrom, mapped_gaps, None
    
    def find_gaps(self, fasta_file, min_gap_size=1):
        """Find gap regions in FASTA file"""
        gaps = []
        gap_id = 1
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id
            sequence = str(record.seq).upper()
            
            gap_pattern = re.compile('N+')
            for match in gap_pattern.finditer(sequence):
                start = match.start()
                end = match.end()
                length = end - start
                if length >= min_gap_size:
                    gaps.append({
                        'chrom': seq_id,
                        'start': start,
                        'end': end,
                        'length': length,
                        'gap_id': f"gap_{gap_id}"
                    })
                    gap_id += 1
        return gaps
    
    def print_summary(self, results):
        """Print summary information"""
        print("\n" + "="*70)
        print("GAP MAPPING SUMMARY")
        print("="*70)
        
        total_gaps = len(results)
        
        # 修复除零错误
        if total_gaps == 0:
            print("No gaps were successfully mapped!")
            print("Possible issues:")
            print("  1. No gaps found in before genome (check --before file)")
            print("  2. Chromosome names don't match between assemblies")
            print("  3. Minimap2 alignment failed (check log file)")
            print("  4. Nucmer alignment failed (check temp_nucmer/*.log)")
            return
        
        mapped_count = sum(1 for r in results if r.get('status') == 'MAPPED')
        
        print(f"Total gaps: {total_gaps}")
        print(f"  UNMAPPED: {self.unmapped_count} ({self.unmapped_count*100/total_gaps:.1f}%)")
        print(f"  MAPPED: {mapped_count} ({mapped_count*100/total_gaps:.1f}%)")
        
        if self.hifi_reads:
            print(f"\n" + "="*70)
            print("QUALITY ASSESSMENT RESULTS (Allele-specific coverage analysis)")
            print("="*70)
            print(f"  MAPPED_WITH_N (has Ns, fill failed): {self.mapped_with_n_count}")
            print(f"  MAPPED_CLEAN (no Ns, fill successful): {self.mapped_clean_count}")
            print(f"    ├─ EXCELLENT (M-AQI ≥ 90): {self.high_quality_count}")
            print(f"    ├─ GOOD (75-89): {self.good_quality_count}")
            print(f"    ├─ FAIR (60-74): {self.fair_quality_count}")
            print(f"    └─ POOR/FAILED (<60): {self.poor_quality_count}")
            
            # Calculate heterozygosity stats
            het_ratios = [r.get('heterozygous_ratio', 0) for r in results if r.get('heterozygous_ratio') is not None]
            if het_ratios:
                avg_het = np.mean(het_ratios)
                print(f"\n  Average heterozygosity across filled regions: {avg_het:.1%}")
            
            zero_stretch_gaps = [r for r in results if r.get('has_zero_stretches', False)]
            if zero_stretch_gaps:
                print(f"\n  WARNING: Gaps with zero coverage stretches (≥{self.min_zero_stretch_bp}bp): {len(zero_stretch_gaps)}")
                total_zero_bases = sum(r.get('total_zero_bases', 0) for r in zero_stretch_gaps)
                print(f"     Total zero coverage bases: {total_zero_bases:,} bp")
            
            print(f"\n" + "="*70)
            print("PERFORMANCE METRICS")
            print("="*70)
            
            fill_rate = self.mapped_clean_count / total_gaps if total_gaps > 0 else 0
            print(f"Fill Rate (clean fills / total gaps):")
            print(f"  {self.mapped_clean_count} / {total_gaps} = {fill_rate:.1%}")
            
            if self.mapped_clean_count > 0:
                high_quality_rate = self.high_quality_count / self.mapped_clean_count
                good_quality_rate = self.good_quality_count / self.mapped_clean_count
                fair_quality_rate = self.fair_quality_count / self.mapped_clean_count
                poor_quality_rate = self.poor_quality_count / self.mapped_clean_count
                print(f"\nQuality Distribution (among clean fills):")
                print(f"  EXCELLENT (≥90): {self.high_quality_count:>8} ({high_quality_rate:>10.1%})")
                print(f"  GOOD (75-89):    {self.good_quality_count:>8} ({good_quality_rate:>10.1%})")
                print(f"  FAIR (60-74):    {self.fair_quality_count:>8} ({fair_quality_rate:>10.1%})")
                print(f"  POOR/FAILED (<60): {self.poor_quality_count:>8} ({poor_quality_rate:>10.1%})")
            else:
                print(f"\nQuality Rate: No clean fills to evaluate")
            
            print(f"\n" + "-"*70)
            print("Summary Table:")
            print(f"  +------------------------+----------+-------------+")
            print(f"  | Metric                  | Value    | Percentage  |")
            print(f"  +------------------------+----------+-------------+")
            print(f"  | Fill Rate (clean/total) | {self.mapped_clean_count:>8} | {fill_rate:>10.1%} |")
            print(f"  | EXCELLENT (≥90)         | {self.high_quality_count:>8} | {high_quality_rate:>10.1%} |")
            print(f"  | GOOD (75-89)            | {self.good_quality_count:>8} | {good_quality_rate:>10.1%} |")
            print(f"  | FAIR (60-74)            | {self.fair_quality_count:>8} | {fair_quality_rate:>10.1%} |")
            print(f"  | POOR/FAILED (<60)       | {self.poor_quality_count:>8} | {poor_quality_rate:>10.1%} |")
            print(f"  +------------------------+----------+-------------+")
            
            if self.unmapped_gaps:
                print(f"\nWARNING: UNMAPPED gaps ({len(self.unmapped_gaps)}): {', '.join(self.unmapped_gaps[:10])}")
                if len(self.unmapped_gaps) > 10:
                    print(f"   ... and {len(self.unmapped_gaps) - 10} more")
            
            if self.mapped_with_n_gaps:
                print(f"\nWARNING: MAPPED_WITH_N gaps ({len(self.mapped_with_n_gaps)}): {', '.join(self.mapped_with_n_gaps[:10])}")
                if len(self.mapped_with_n_gaps) > 10:
                    print(f"   ... and {len(self.mapped_with_n_gaps) - 10} more")
            
            print(f"\n" + "="*70)
            print("ALLELE-SPECIFIC COVERAGE INTERPRETATION")
            print("="*70)
            print("This analysis uses EFFECTIVE DEPTH (support for the major allele):")
            print("  - Raw depth: Total number of reads covering a position")
            print("  - Effective depth: Reads supporting the most frequent base at that position")
            print("  - Coverage efficiency: Effective depth / Raw depth × 100%")
            print("\nM-AQI (Mapping Assembly Quality Index) components:")
            print("  - Internal region quality (R-AQI) - 50% weight")
            print("  - Left boundary consistency - 25% weight")
            print("  - Right boundary consistency - 25% weight")
            print("\nR-AQI components (using effective depth):")
            print("  - Depth adequacy: 20% (based on effective depth)")
            print("  - Coverage uniformity: 15%")
            print("  - Coverage completeness: 10%")
            print("  - Low coverage density: 25% (using effective depth)")
            print("  - Mapping quality (MAPQ): 15%")
            print("  - Read support (identity): 15%")
            print(f"\nGenotype interpretation (per base):")
            print(f"  - HOMOZYGOUS: Major allele frequency ≥80%")
            print(f"  - HETEROZYGOUS: Major allele frequency 20-80%")
            print(f"  - LOW FREQ: Major allele frequency <20% (possible variant or error)")
            print("\nScore interpretation:")
            print("  EXCELLENT (90-100): High quality fill, low heterozygosity, high efficiency")
            print("  GOOD (75-89): Good quality fill")
            print("  FAIR (60-74): Acceptable, may have moderate heterozygosity")
            print("  POOR (40-59): Low quality, high heterozygosity or low efficiency")
            print("  FAILED (0-39): Fill failed or very poor quality")
            print(f"\nMODIFICATIONS:")
            print(f"  - Pipeline alignment: minimap2 | samtools sort (no intermediate SAM file)")
            print(f"  - Downsampling: automatic step sizing for large regions (max {self.max_pdf_length:,} bp)")
            print(f"  - Low coverage threshold: {self.low_coverage_threshold}x (was 5x)")
    
    def run(self, before_fasta, after_fasta, output_dir, max_workers=8, min_gap_size=1, plot_coverage=False):
        """Main run method"""
        print("="*70)
        print("fill_quality_check.py - Gap Filling Quality Assessment")
        print("Core logic: Map gaps from before genome to after genome")
        print("Using effective depth (major allele support) for quality assessment")
        print(f"Zero coverage detection: threshold={self.min_zero_stretch_bp}bp for flagging")
        print(f"Low coverage detection: threshold={self.low_coverage_threshold}x (effective depth)")
        print(f"Boundary consistency: good ratio={self.boundary_good_ratio:.0%}, min spanning={self.boundary_min_spanning}")
        print(f"Alignment: primary alignments only (--secondary=no), pipeline mode (no intermediate SAM)")
        print(f"PDF downsampling: automatic (max region {self.max_pdf_length:,} bp)")
        if plot_coverage:
            print("PDF coverage plots: ENABLED (with automatic downsampling for large regions)")
        print("="*70)
        print(f"Before genome (with gaps): {before_fasta}")
        print(f"After genome (filled): {after_fasta}")
        print(f"Output: {output_dir}")
        if self.hifi_reads:
            print(f"HiFi reads: {self.hifi_reads}")
            print(f"Flank size for quality assessment: {self.flank_size} bp on each side")
        print("-"*70)
        
        self.log("Finding gaps in before genome...")
        gaps = self.find_gaps(before_fasta, min_gap_size)
        
        if not gaps:
            self.log("No gaps found in before genome!", "warning")
            self.log("Please check that your --before file contains N bases", "warning")
            return []
        
        gaps_by_chrom = {}
        for gap in gaps:
            chrom = gap['chrom']
            if chrom not in gaps_by_chrom:
                gaps_by_chrom[chrom] = []
            gaps_by_chrom[chrom].append(gap)
        
        self.total = len(gaps_by_chrom)
        total_gaps = len(gaps)
        self.log(f"Found {total_gaps} gaps in {self.total} chromosomes")
        
        os.makedirs(output_dir, exist_ok=True)
        temp_dir = os.path.join(output_dir, "temp_nucmer")
        os.makedirs(temp_dir, exist_ok=True)
        
        if self.hifi_reads:
            self.log("Preparing for quality assessment...", "info")
            self.index_genome(after_fasta)
            self.bam_file = os.path.join(output_dir, "aligned_reads.primary.bam")
            if not self.align_hifi_reads(after_fasta, self.bam_file):
                self.log("HiFi alignment failed, continuing without quality assessment", "warning")
                self.bam_file = None
            else:
                self.log(f"Using filtered BAM: {self.bam_file}", "info")
        
        all_results = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for chrom, gaps_list in gaps_by_chrom.items():
                future = executor.submit(
                    self.process_chromosome,
                    chrom, gaps_list, before_fasta, after_fasta, temp_dir, output_dir,
                    plot_coverage
                )
                futures[future] = chrom
            
            for future in as_completed(futures):
                chrom, mapped_gaps, error = future.result()
                if mapped_gaps:
                    all_results.extend(mapped_gaps)
        
        # 只有当有结果时才写入合并文件
        if all_results:
            merged_file = os.path.join(output_dir, "all_gaps_mapping_with_quality.tsv")
            with open(merged_file, 'w') as f:
                f.write("gap_id\toriginal_chrom\toriginal_start\toriginal_end\toriginal_length\t"
                       "filled_chrom\tfilled_start\tfilled_end\tfilled_length\tidentity\tmethod\tstatus\tclassification\t"
                       "n_count\tn_ratio\tcoverage_rate\tmean_depth\tmean_mapq\tmean_identity\t"
                       "quality_score\tquality_grade\tm_aqi\tinternal_quality\tleft_consistency\t"
                       "right_consistency\tcv_depth\thas_coverage_gaps\tlow_identity_ratio\t"
                       "has_zero_stretches\tzero_stretches_count\ttotal_zero_bases\tzero_ratio_in_filled\t"
                       "worst_zero_stretch_length\tzero_positions_count\tfiltered_reads\t"
                       "low_cov_stretches_count\tlow_cov_density\tlow_cov_score\theterozygous_ratio\tcoverage_efficiency\n")
                for r in all_results:
                    f.write(f"{r['gap_id']}\t{r['original_chrom']}\t{r['original_start']}\t{r['original_end']}\t{r['original_length']}\t")
                    f.write(f"{r['filled_chrom']}\t{r['filled_start']}\t{r['filled_end']}\t{r['filled_length']}\t")
                    f.write(f"{r['identity']:.1f}\t{r['method']}\t{r['status']}\t{r.get('classification', 'N/A')}\t")
                    f.write(f"{r.get('n_count', 'N/A')}\t{r.get('n_ratio', 'N/A')}\t")
                    f.write(f"{r.get('coverage_rate', 'N/A')}\t{r.get('mean_depth', 'N/A')}\t")
                    f.write(f"{r.get('mean_mapq', 'N/A')}\t{r.get('mean_identity', 'N/A')}\t")
                    f.write(f"{r.get('quality_score', 'N/A')}\t{r.get('quality_grade', 'N/A')}\t")
                    f.write(f"{r.get('m_aqi', 'N/A')}\t{r.get('internal_quality', 'N/A')}\t")
                    f.write(f"{r.get('left_consistency', 'N/A')}\t{r.get('right_consistency', 'N/A')}\t")
                    f.write(f"{r.get('cv_depth', 'N/A')}\t{r.get('has_coverage_gaps', 'N/A')}\t")
                    f.write(f"{r.get('low_identity_ratio', 'N/A')}\t")
                    f.write(f"{r.get('has_zero_stretches', 'N/A')}\t{r.get('zero_stretches_count', 'N/A')}\t")
                    f.write(f"{r.get('total_zero_bases', 'N/A')}\t{r.get('zero_ratio_in_filled', 'N/A')}\t")
                    f.write(f"{r.get('worst_zero_stretch_length', 'N/A')}\t")
                    f.write(f"{r.get('zero_positions_count', 'N/A')}\t")
                    f.write(f"{r.get('filtered_reads', 'N/A')}\t")
                    f.write(f"{r.get('low_cov_stretches_count', 'N/A')}\t")
                    f.write(f"{r.get('low_cov_density', 'N/A')}\t")
                    f.write(f"{r.get('low_cov_score', 'N/A')}\t")
                    f.write(f"{r.get('heterozygous_ratio', 'N/A')}\t")
                    f.write(f"{r.get('coverage_efficiency', 'N/A')}\n")
        else:
            self.log("No gaps were successfully mapped!", "warning")
        
        self.print_summary(all_results)
        
        if self.hifi_reads and all_results:
            stats_file = os.path.join(output_dir, "assessment_stats.json")
            stats = {
                'total_gaps': total_gaps,
                'unmapped_count': self.unmapped_count,
                'mapped_with_n_count': self.mapped_with_n_count,
                'mapped_clean_count': self.mapped_clean_count,
                'excellent_count': self.high_quality_count,
                'good_count': self.good_quality_count,
                'fair_count': self.fair_quality_count,
                'poor_count': self.poor_quality_count,
                'fill_rate': self.mapped_clean_count / total_gaps if total_gaps > 0 else 0,
                'flank_size': self.flank_size,
                'min_zero_stretch_bp': self.min_zero_stretch_bp,
                'low_coverage_threshold': self.low_coverage_threshold,
                'min_alignment_coverage': self.min_alignment_coverage,
                'boundary_good_ratio': self.boundary_good_ratio,
                'boundary_min_spanning': self.boundary_min_spanning,
                'n_threshold': 0.05,
                'unmapped_gaps': self.unmapped_gaps,
                'mapped_with_n_gaps': self.mapped_with_n_gaps,
                'primary_alignments_only': True,
                'allele_specific_analysis': True,
                'effective_depth_used': True,
                'pipeline_mode': True,
                'downsampling_enabled': True,
                'max_pdf_length': self.max_pdf_length
            }
            with open(stats_file, 'w') as f:
                json.dump(stats, f, indent=2)
        
        if plot_coverage:
            plot_dir = os.path.join(output_dir, "coverage_plots")
            if os.path.exists(plot_dir):
                self.log(f"Plots saved to {plot_dir}", "info")
        
        return all_results


def main():
    parser = argparse.ArgumentParser(description='fill_quality_check.py - Gap filling quality assessment with allele-specific analysis')
    parser.add_argument('--before', required=True, help='Reference genome FASTA (with N gaps)')
    parser.add_argument('--after', required=True, help='Filled genome FASTA')
    parser.add_argument('-o', '--output', default='gap_mapping_results', help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Maximum threads for parallel processing')
    parser.add_argument('--threads-per-alignment', type=int, default=2, help='Threads per nucmer/minimap2')
    parser.add_argument('--min-gap-size', type=int, default=1, help='Minimum gap size')
    parser.add_argument('--keep-intermediate', action='store_true', help='Keep intermediate files')
    parser.add_argument('--quiet', action='store_true', help='Quiet mode')
    
    parser.add_argument('--hifi-reads', help='HiFi reads FASTQ file for quality assessment')
    parser.add_argument('--min-depth', type=int, default=5, help='Minimum depth for quality assessment')
    parser.add_argument('--min-coverage', type=float, default=0.5, help='Minimum coverage for quality assessment')
    parser.add_argument('--flank-size', type=int, default=5000,
                       help='Number of bases to expand on each side for quality assessment (default: 5000)')
    parser.add_argument('--min-zero-stretch', type=int, default=1,
                       help='Minimum consecutive zero coverage bases to flag as problem (default: 1)')
    parser.add_argument('--boundary-good-ratio', type=float, default=0.25,
                       help='Good boundary crossing ratio threshold (default: 0.25 = 25%%)')
    parser.add_argument('--boundary-min-spanning', type=int, default=5,
                       help='Minimum spanning reads for full score (default: 5)')
    parser.add_argument('--plot-coverage', action='store_true',
                       help='Generate PDF coverage plots (with automatic downsampling for large regions)')
    parser.add_argument('--min-alignment-coverage', type=float, default=0.7,
                       help='Minimum alignment coverage (currently disabled)')
    parser.add_argument('--low-coverage-threshold', type=int, default=3,
                       help='Coverage threshold for low coverage regions (default: 3, was 5)')
    parser.add_argument('--max-pdf-length', type=int, default=100000,
                       help='Maximum region length to generate PDF (default: 100000 bp)')
    
    args = parser.parse_args()
    
    mapper = SimpleGapMapperWithQuality(
        verbose=not args.quiet,
        threads_per_alignment=args.threads_per_alignment,
        keep_intermediate=args.keep_intermediate,
        hifi_reads=args.hifi_reads,
        min_depth=args.min_depth,
        min_coverage=args.min_coverage,
        flank_size=args.flank_size,
        min_zero_stretch_bp=args.min_zero_stretch,
        boundary_good_ratio=args.boundary_good_ratio,
        boundary_min_spanning=args.boundary_min_spanning,
        min_alignment_coverage=args.min_alignment_coverage,
        low_coverage_threshold=args.low_coverage_threshold,
        max_pdf_length=args.max_pdf_length
    )
    
    mapper.run(
        before_fasta=args.before,
        after_fasta=args.after,
        output_dir=args.output,
        max_workers=args.threads,
        min_gap_size=args.min_gap_size,
        plot_coverage=args.plot_coverage
    )


if __name__ == "__main__":
    main()
