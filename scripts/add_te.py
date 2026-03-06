#!/usr/bin/env python3
"""
Telomere Repair Sequence Merge Script (Telomere-First Version)
Purpose: Merge extracted telomere repair sequences to corresponding chromosome ends
Core Logic: PRIORITIZE SEQUENCES THAT RESTORE TELOMERES
Features: 
  1. Try all candidates, not just first success
  2. Validate telomeres BEFORE and AFTER repair
  3. ONLY accept sequences that restore/preserve telomeres
  4. Fall back to best available only if no telomere-restoring sequences found
"""

import sys
import os
import re
import argparse
import textwrap
import subprocess
import tempfile
import time
from typing import Dict, List, Tuple, Optional, Any, Set
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VERBOSE = False

def printv(message: str):
    if VERBOSE:
        print(message)

# ============================================================================
# Telomere detection module
# ============================================================================

class TelomereDetector:
    """Telomere detector: used to detect telomere repeat units in sequences"""
    
    # Standard telomere patterns
    TELOMERE_PATTERNS = {
        'vertebrate': 'TTAGGG',      # vertebrate telomere
        'vertebrate_rev': 'CCCTAA',  # reverse complement
        'alternative': ['TTAGGG', 'CCCTAA', 'TTAGGC', 'TTGGGG'],  # other variants
    }
    
    def __init__(self, min_repeats: int = 20, min_length: int = 500):
        """
        Initialize telomere detector
        
        Args:
            min_repeats: minimum number of repeats
            min_length: minimum telomere length
        """
        self.min_repeats = min_repeats
        self.min_length = min_length
        self.patterns = [
            self.TELOMERE_PATTERNS['vertebrate'],
            self.TELOMERE_PATTERNS['vertebrate_rev']
        ]
    
    def detect(self, sequence: str, region_name: str = "unknown") -> Dict[str, Any]:
        """
        Detect telomere regions in sequence
        
        Args:
            sequence: DNA sequence to detect
            region_name: region name (for logging)
            
        Returns:
            telomere detection result dictionary
        """
        if not sequence:
            return {
                'has_telomere': False,
                'telomere_count': 0,
                'telomeres': [],
                'total_telomere_length': 0,
                'max_repeats': 0,
                'max_repeat_region': None,
                'region_name': region_name,
                'sequence_length': 0
            }
        
        sequence = sequence.upper()
        seq_len = len(sequence)
        telomeres = []
        
        # detect forward telomeres (TTAGGG)
        forward_telomeres = self._detect_pattern(sequence, self.patterns[0], "forward")
        telomeres.extend(forward_telomeres)
        
        # detect reverse telomeres (CCCTAA)
        reverse_telomeres = self._detect_pattern(sequence, self.patterns[1], "reverse")
        telomeres.extend(reverse_telomeres)
        
        # sort by position
        telomeres.sort(key=lambda x: x['start'])
        
        # merge adjacent telomere regions
        merged_telomeres = self._merge_adjacent_telomeres(telomeres)
        
        # calculate statistics
        total_length = sum(t['length'] for t in merged_telomeres)
        max_repeats = max((t['repeats'] for t in merged_telomeres), default=0)
        max_repeat_region = next((t for t in merged_telomeres if t['repeats'] == max_repeats), None)
        
        result = {
            'has_telomere': len(merged_telomeres) > 0,
            'telomere_count': len(merged_telomeres),
            'telomeres': merged_telomeres,
            'total_telomere_length': total_length,
            'max_repeats': max_repeats,
            'max_repeat_region': max_repeat_region,
            'region_name': region_name,
            'sequence_length': seq_len,
            'telomere_density': total_length / seq_len if seq_len > 0 else 0
        }
        
        return result
    
    def _detect_pattern(self, sequence: str, pattern: str, 
                       pattern_type: str) -> List[Dict[str, Any]]:
        """
        Detect specific telomere repeat pattern
        
        Args:
            sequence: DNA sequence
            pattern: telomere pattern (e.g. 'TTAGGG')
            pattern_type: type (forward/reverse)
            
        Returns:
            list of telomere regions
        """
        telomeres = []
        pattern_len = len(pattern)
        seq_len = len(sequence)
        
        i = 0
        while i < seq_len - pattern_len:
            # check if current position matches pattern
            if sequence[i:i+pattern_len] == pattern:
                # calculate consecutive repeats
                repeat_count = 1
                j = i + pattern_len
                
                # extend forward to detect consecutive repeats
                while j <= seq_len - pattern_len:
                    if sequence[j:j+pattern_len] == pattern:
                        repeat_count += 1
                        j += pattern_len
                    else:
                        # check for possible mismatches (allow 1 mismatch)
                        mismatches = sum(1 for a, b in zip(sequence[j:j+pattern_len], pattern) if a != b)
                        if mismatches <= 1:  # allow 1 base mismatch
                            repeat_count += 1
                            j += pattern_len
                        else:
                            break
                
                # check if minimum repeat count is reached
                if repeat_count >= self.min_repeats:
                    end_pos = i + repeat_count * pattern_len
                    length = end_pos - i
                    
                    # check minimum length
                    if length >= self.min_length:
                        telomeres.append({
                            'start': i,
                            'end': end_pos,
                            'length': length,
                            'repeats': repeat_count,
                            'type': pattern_type,
                            'pattern': pattern,
                            'sequence': sequence[i:end_pos]
                        })
                
                i = j
            else:
                i += 1
        
        return telomeres
    
    def _merge_adjacent_telomeres(self, telomeres: List[Dict], 
                                  max_gap: int = 100) -> List[Dict]:
        """
        Merge adjacent telomere regions
        
        Args:
            telomeres: list of telomere regions
            max_gap: maximum allowed gap
            
        Returns:
            merged list of telomere regions
        """
        if not telomeres:
            return []
        
        merged = []
        current = telomeres[0].copy()
        
        for t in telomeres[1:]:
            # if current region and next region are adjacent or have small gap
            if t['start'] - current['end'] <= max_gap:
                # merge regions
                current['end'] = max(current['end'], t['end'])
                current['length'] = current['end'] - current['start']
                current['repeats'] = max(current['repeats'], t['repeats'])
                current['type'] = f"{current['type']}_merged"
                current['sequence'] = "merged_region"
            else:
                merged.append(current)
                current = t.copy()
        
        merged.append(current)
        return merged
    
    def compare_telomere_status(self, before: Dict[str, Any], 
                                after: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compare telomere status before and after repair
        
        Args:
            before: telomere detection result before repair
            after: telomere detection result after repair
            
        Returns:
            comparison result dictionary
        """
        comparison = {
            'status': 'UNKNOWN',
            'telomere_restored': False,
            'telomere_improved': False,
            'telomere_preserved': False,
            'telomere_damaged': False,
            'details': {}
        }
        
        # Case 1: no telomeres before, telomeres after (best)
        if not before['has_telomere'] and after['has_telomere']:
            comparison['status'] = 'RESTORED'
            comparison['telomere_restored'] = True
            comparison['details'] = {
                'new_telomeres': after['telomere_count'],
                'max_repeats': after['max_repeats'],
                'total_length': after['total_telomere_length'],
                'density': after['telomere_density']
            }
        
        # Case 2: telomeres before, telomeres after and better
        elif before['has_telomere'] and after['has_telomere']:
            length_diff = after['total_telomere_length'] - before['total_telomere_length']
            repeats_diff = after['max_repeats'] - before['max_repeats']
            density_diff = after['telomere_density'] - before['telomere_density']
            
            if length_diff > 0 or repeats_diff > 0 or density_diff > 0:
                comparison['status'] = 'IMPROVED'
                comparison['telomere_improved'] = True
                comparison['details'] = {
                    'length_increase': length_diff,
                    'repeats_increase': repeats_diff,
                    'density_increase': density_diff,
                    'before': {
                        'length': before['total_telomere_length'],
                        'repeats': before['max_repeats'],
                        'density': before['telomere_density']
                    },
                    'after': {
                        'length': after['total_telomere_length'],
                        'repeats': after['max_repeats'],
                        'density': after['telomere_density']
                    }
                }
            else:
                comparison['status'] = 'PRESERVED'
                comparison['telomere_preserved'] = True
                comparison['details'] = {
                    'length': before['total_telomere_length'],
                    'repeats': before['max_repeats']
                }
        
        # Case 3: telomeres before, no telomeres after (worst)
        elif before['has_telomere'] and not after['has_telomere']:
            comparison['status'] = 'DAMAGED'
            comparison['telomere_damaged'] = True
            comparison['details'] = {
                'lost_telomeres': before['telomere_count'],
                'lost_length': before['total_telomere_length'],
                'lost_repeats': before['max_repeats']
            }
        
        # Case 4: no telomeres in either
        elif not before['has_telomere'] and not after['has_telomere']:
            comparison['status'] = 'NO_TELOMERE'
            comparison['details'] = {
                'message': 'No telomeres detected before or after repair'
            }
        
        return comparison


# ============================================================================
# Minimap2 alignment module
# ============================================================================

class Minimap2Matcher:
    
    def __init__(self, mode: str = 'asm5', min_score: float = 0.7,
                 min_match_length: int = 100, min_mapq: int = 20):
        self.mode = mode
        self.min_score = min_score
        self.min_match_length = min_match_length
        self.min_mapq = min_mapq
        
        self._check_minimap2()
    
    def _check_minimap2(self):
        try:
            result = subprocess.run(['minimap2', '--version'], 
                                  capture_output=True, text=True, check=False)
            if result.returncode == 0:
                version = result.stdout.strip()
                printv(f"minimap2 installed: {version}")
            else:
                raise RuntimeError("minimap2 not correctly installed")
        except FileNotFoundError:
            raise RuntimeError("minimap2 not installed, please install: conda install -c bioconda minimap2")
    
    def find_best_match(self, query_seq: str, target_seq: str,
                       target_start_pos: int = 0) -> Optional[Dict]:
        
        if not query_seq or not target_seq:
            return None
        
        query_len = len(query_seq)
        target_len = len(target_seq)
        
        printv(f"    Aligning {query_len:,}bp vs {target_len:,}bp (mode: {self.mode})")
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_file:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as target_file:
                query_file.write(f">query\n{query_seq}\n")
                target_file.write(f">target\n{target_seq}\n")
                query_file.flush()
                target_file.flush()
                
                cmd = ['minimap2', '-x', self.mode, '-c', '--secondary=no',
                       target_file.name, query_file.name]
                
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, 
                                          check=True, timeout=300)
                    return self._parse_best_match(result.stdout, target_start_pos)
                except subprocess.CalledProcessError as e:
                    printv(f"    minimap2 execution failed: {e}")
                    return None
                except subprocess.TimeoutExpired:
                    printv(f"    minimap2 execution timeout")
                    return None
                finally:
                    try:
                        os.unlink(query_file.name)
                        os.unlink(target_file.name)
                    except:
                        pass
    
    def _parse_best_match(self, paf_output: str, target_start_pos: int) -> Optional[Dict]:
        best_match = None
        best_score = 0
        
        for line in paf_output.strip().split('\n'):
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 12:
                continue
            
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target_start = int(fields[7]) + target_start_pos
            target_end = int(fields[8]) + target_start_pos
            matches = int(fields[9])
            alignment_len = int(fields[10])
            mapq = int(fields[11])
            
            if alignment_len > 0:
                identity = matches / alignment_len
                coverage = (query_end - query_start) / query_len if query_len > 0 else 0
                match_length = query_end - query_start
                
                length_factor = min(1.0, match_length / 5000.0)
                mapq_factor = min(1.0, mapq / 60.0)
                score = identity * (0.5 + 0.5 * length_factor) * (0.7 + 0.3 * mapq_factor)
                
                match_info = {
                    'query_start': query_start,
                    'query_end': query_end,
                    'target_start': target_start,
                    'target_end': target_end,
                    'identity': identity,
                    'coverage': coverage,
                    'mapq': mapq,
                    'strand': strand,
                    'matches': matches,
                    'alignment_len': alignment_len,
                    'match_length': match_length,
                    'score': score,
                    'query_len': query_len
                }
                
                if score > best_score:
                    best_score = score
                    best_match = match_info
        
        return best_match


# ============================================================================
# Telomere repair main class
# ============================================================================

class TelomereRepair:
    
    def __init__(self, config: Dict = None):
        self.config = config or {
            'minimap2_mode': 'asm5',
            'min_score': 0.7,
            'min_match_length': 100,
            'min_mapq': 20,
            'search_range': 5000000,
            'min_telomere_repeats': 20,
            'min_telomere_length': 500
        }
        
        self.matcher = Minimap2Matcher(
            mode=self.config['minimap2_mode'],
            min_score=self.config['min_score'],
            min_match_length=self.config['min_match_length'],
            min_mapq=self.config['min_mapq']
        )
        
        # initialize telomere detector
        self.telomere_detector = TelomereDetector(
            min_repeats=self.config.get('min_telomere_repeats', 20),
            min_length=self.config.get('min_telomere_length', 500)
        )
    
    def find_telomere_overlap(self, chr_seq: str, repair_seq: str, 
                             repair_type: str) -> Optional[Dict]:
        
        chr_len = len(chr_seq)
        repair_len = len(repair_seq)
        search_range = self.config['search_range']
        
        printv(f"    Chromosome length: {chr_len:,} bp")
        printv(f"    Repair sequence length: {repair_len:,} bp")
        printv(f"    Search range: {search_range:,} bp")
        
        if repair_type == '5prime':
            printv(f"    5' end repair strategy: repair sequence vs chromosome start {search_range:,}bp")
            
            chr_region_end = min(search_range, chr_len)
            chr_region = chr_seq[:chr_region_end]
            chr_region_start = 0
            
            match = self.matcher.find_best_match(repair_seq, chr_region, chr_region_start)
            
            if match:
                printv(f"    Found 5' end overlap:")
                printv(f"      Repair sequence position: {match['query_start']:,}-{match['query_end']:,} "
                      f"(length: {match['match_length']:,}bp)")
                printv(f"      Chromosome position: {match['target_start']:,}-{match['target_end']:,}")
                printv(f"      Similarity: {match['identity']:.3f}")
                printv(f"      MAPQ: {match['mapq']}")
                printv(f"      Strand: {match['strand']}")
            
            return match
            
        elif repair_type == '3prime':
            printv(f"    3' end repair strategy: repair sequence vs chromosome end {search_range:,}bp")
            
            chr_region_start = max(0, chr_len - search_range)
            chr_region = chr_seq[chr_region_start:chr_len]
            
            match = self.matcher.find_best_match(repair_seq, chr_region, chr_region_start)
            
            if match:
                printv(f"    Found 3' end overlap:")
                printv(f"      Repair sequence position: {match['query_start']:,}-{match['query_end']:,} "
                      f"(length: {match['match_length']:,}bp)")
                printv(f"      Chromosome position: {match['target_start']:,}-{match['target_end']:,}")
                printv(f"      Similarity: {match['identity']:.3f}")
                printv(f"      MAPQ: {match['mapq']}")
                printv(f"      Strand: {match['strand']}")
            
            return match
        
        else:
            printv(f"    Unknown repair type: {repair_type}")
            return None
    
    def validate_match(self, match: Dict, repair_type: str, 
                      chr_len: int, repair_len: int) -> Tuple[bool, str]:
        if not match:
            return False, "No match"
        
        if match['identity'] < self.config['min_score']:
            return False, f"Insufficient similarity: {match['identity']:.3f} < {self.config['min_score']}"
        
        if match['match_length'] < self.config['min_match_length']:
            return False, f"Insufficient match length: {match['match_length']:,} < {self.config['min_match_length']}"
        
        if match['mapq'] < self.config['min_mapq']:
            return False, f"Insufficient MAPQ: {match['mapq']} < {self.config['min_mapq']}"
        
        return True, "Validation passed"
    
    def merge_sequences(self, chr_seq: str, repair_seq: str, 
                       match: Dict, repair_type: str) -> Tuple[str, Dict]:
        
        chr_len = len(chr_seq)
        repair_len = len(repair_seq)
        
        printv(f"    Starting sequence merge...")
        printv(f"    Match region: repair sequence[{match['query_start']:,}-{match['query_end']:,}] "
              f"vs chromosome[{match['target_start']:,}-{match['target_end']:,}]")
        
        if match['strand'] == '-':
            printv(f"    Processing reverse strand match: reverse complement repair sequence")
            original_query_start = match['query_start']
            original_query_end = match['query_end']
            match['query_start'] = repair_len - original_query_end
            match['query_end'] = repair_len - original_query_start
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                         'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
            repair_seq = ''.join(complement.get(base, base) for base in reversed(repair_seq))
            printv(f"    Reverse complement coordinates: repair sequence[{match['query_start']:,}-{match['query_end']:,}]")
        
        if repair_type == '5prime':
            repair_unique = repair_seq[:match['query_start']]
            overlap_from_repair = repair_seq[match['query_start']:match['query_end']]
            chr_remainder_start = match['target_end']
            chr_remainder = chr_seq[chr_remainder_start:] if chr_remainder_start < chr_len else ""
            merged_seq = repair_unique + overlap_from_repair + chr_remainder
            
            stats = {
                'method': '5prime_overlap',
                'repair_unique_len': len(repair_unique),
                'overlap_len': len(overlap_from_repair),
                'chr_remainder_len': len(chr_remainder),
                'cut_positions': {
                    'repair_cut': match['query_start'],
                    'chr_cut': match['target_end']
                }
            }
        
        elif repair_type == '3prime':
            chr_unique_end = match['target_start']
            chr_unique = chr_seq[:chr_unique_end]
            overlap_from_chr = chr_seq[match['target_start']:match['target_end']]
            repair_unique_start = match['query_end']
            repair_unique = repair_seq[repair_unique_start:] if repair_unique_start < repair_len else ""
            merged_seq = chr_unique + overlap_from_chr + repair_unique
            
            stats = {
                'method': '3prime_overlap',
                'chr_unique_len': len(chr_unique),
                'overlap_len': len(overlap_from_chr),
                'repair_unique_len': len(repair_unique),
                'cut_positions': {
                    'chr_cut': match['target_start'],
                    'repair_cut': match['query_end']
                }
            }
        
        else:
            raise ValueError(f"Unknown repair type: {repair_type}")
        
        stats['original_chr_len'] = chr_len
        stats['original_repair_len'] = repair_len
        stats['merged_len'] = len(merged_seq)
        stats['extension'] = len(merged_seq) - chr_len
        stats['identity'] = match['identity']
        stats['coverage'] = match['coverage']
        stats['mapq'] = match['mapq']
        stats['strand'] = match['strand']
        stats['score'] = match['score']
        
        return merged_seq, stats
    
    def repair_with_verification(self, chr_id: str, chr_seq: str, 
                                repair_info: Dict) -> Tuple[Optional[str], Dict]:
        """
        Repair process with telomere verification
        
        Args:
            chr_id: chromosome ID
            chr_seq: chromosome sequence
            repair_info: repair sequence information
            
        Returns:
            (repaired sequence, repair statistics)
        """
        repair_seq = repair_info['sequence']
        repair_type = repair_info.get('repair_type', 'unknown')
        repair_id = repair_info['id']
        
        printv(f"\n  Attempting repair with sequence: {repair_id[:50]}...")
        printv(f"    Length: {len(repair_seq):,} bp")
        printv(f"    Type: {repair_type}")
        
        # 1. detect telomeres before repair
        before_telomeres = self.telomere_detector.detect(
            chr_seq, region_name=f"{chr_id}_before"
        )
        
        # 2. find overlap region
        match = self.find_telomere_overlap(chr_seq, repair_seq, repair_type)
        
        if not match:
            printv(f"    No overlap found")
            return None, {'error': 'no_overlap'}
        
        # 3. validate alignment quality
        is_valid, valid_msg = self.validate_match(
            match, repair_type, len(chr_seq), len(repair_seq)
        )
        
        if not is_valid:
            printv(f"    Validation failed: {valid_msg}")
            return None, {'error': 'validation_failed', 'message': valid_msg}
        
        # 4. perform merge
        merged_seq, stats = self.merge_sequences(chr_seq, repair_seq, match, repair_type)
        
        # 5. detect telomeres after repair
        after_telomeres = self.telomere_detector.detect(
            merged_seq, region_name=f"{chr_id}_after"
        )
        
        # 6. compare telomere status
        telomere_comparison = self.telomere_detector.compare_telomere_status(
            before_telomeres, after_telomeres
        )
        
        # 7. add telomere information to statistics
        stats['telomere_before'] = before_telomeres
        stats['telomere_after'] = after_telomeres
        stats['telomere_comparison'] = telomere_comparison
        
        return merged_seq, stats


# ============================================================================
# Helper functions
# ============================================================================

def read_fasta(filename: str) -> Dict[str, str]:
    sequences = {}
    current_header = ""
    current_seq = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header:
                        sequences[current_header] = ''.join(current_seq).upper()
                    current_header = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_header:
                sequences[current_header] = ''.join(current_seq).upper()
    except FileNotFoundError:
        print(f"Error: File {filename} does not exist")
        raise
    except Exception as e:
        print(f"Error: Failed to read FASTA file {filename}: {e}")
        raise
    
    return sequences


def parse_repair_sequence_id(seq_id: str) -> Dict[str, Any]:
    info = {}
    
    parts = seq_id.split('_')
    if len(parts) >= 4:
        chr_info = parts[3]
        
        chr_name = chr_info
        repair_end = None
        
        if '_5prime' in chr_info.lower() or '5prime' in chr_info.lower():
            repair_end = '5prime'
            chr_name = re.sub(r'[_-]?5prime', '', chr_name, flags=re.IGNORECASE)
        elif '_3prime' in chr_info.lower() or '3prime' in chr_info.lower():
            repair_end = '3prime'
            chr_name = re.sub(r'[_-]?3prime', '', chr_name, flags=re.IGNORECASE)
        
        info['target_chromosome'] = chr_name.strip('_')
        info['repair_type'] = repair_end or 'unknown'
    
    return info


def parse_repair_sequence_description(description: str) -> Dict[str, Any]:
    info = {}
    
    parts = description.split()
    for part in parts:
        if ':' in part:
            key, value = part.split(':', 1)
            info[key] = value
    
    if 'region_type' in info:
        if 'first_alignment_with_preceding' in info['region_type']:
            info['repair_type'] = '5prime'
        elif 'last_alignment_with_subsequent' in info['region_type']:
            info['repair_type'] = '3prime'
    
    if 'aligned_to' in info:
        target = info['aligned_to']
        if '_5prime' in target:
            info['repair_type'] = '5prime'
            target = target.replace('_5prime', '')
        elif '_3prime' in target:
            info['repair_type'] = '3prime'
            target = target.replace('_3prime', '')
        info['target_chromosome'] = target
    
    return info


def group_repair_sequences_by_chromosome(repair_records: List[Any]) -> Dict[str, List[Dict]]:
    chromosome_groups = defaultdict(list)
    
    for record in repair_records:
        try:
            id_info = parse_repair_sequence_id(record.id)
            desc_info = parse_repair_sequence_description(record.description)
            
            repair_info = {**id_info, **desc_info}
            repair_info['sequence'] = str(record.seq).upper()
            repair_info['id'] = record.id
            repair_info['length'] = len(record.seq)
            repair_info['record'] = record
            
            target_chr = repair_info.get('target_chromosome', '')
            if target_chr:
                chromosome_groups[target_chr].append(repair_info)
                
        except Exception as e:
            printv(f"Warning: Error parsing repair sequence {record.id}: {e}")
            continue
    
    for chr_name in chromosome_groups:
        chromosome_groups[chr_name].sort(key=lambda x: x['length'], reverse=True)
        printv(f"  Chromosome {chr_name}: Found {len(chromosome_groups[chr_name])} repair sequences")
    
    return chromosome_groups


def find_matching_chromosome(chromosomes: Dict[str, str], target_chr: str) -> Optional[str]:
    if not target_chr:
        return None
    
    target_clean = target_chr.lower().strip().replace('chr', '').replace('chromosome', '')
    
    for chr_id in chromosomes.keys():
        chr_simple = chr_id.split()[0].lower().replace('chr', '').replace('chromosome', '')
        if target_clean == chr_simple:
            return chr_id
    
    for chr_id in chromosomes.keys():
        chr_simple = chr_id.split()[0].lower()
        if target_clean in chr_simple or chr_simple in target_clean:
            return chr_id
    
    match = re.search(r'(\d+)', target_clean)
    if match:
        num = match.group(1)
        for chr_id in chromosomes.keys():
            chr_simple = chr_id.split()[0]
            if num in chr_simple:
                return chr_id
    
    return None


def calculate_selection_score(result_info: Dict, strategy: str = 'minimal_extension',
                             has_telomere: bool = False) -> float:
    """
    Calculate selection score (prioritize results with telomeres)
    """
    if strategy == 'first_success':
        return 1.0
    
    # base score (based on length change)
    extension_abs = abs(result_info.get('extension', 0))
    extension_norm = max(0, 1.0 - min(1.0, extension_abs / 100000.0))
    quality = result_info.get('minimap2_info', {}).get('identity', 0.5)
    similarity = result_info.get('similarity', 0.7)
    
    base_score = (
        extension_norm * 0.3 +
        quality * 0.3 +
        similarity * 0.2
    )
    
    if has_telomere:
        # results with telomeres get high score
        telomere_after = result_info.get('telomere_after', {})
        telomere_comparison = result_info.get('telomere_comparison', {})
        
        # bonus based on telomere status
        status = telomere_comparison.get('status', 'UNKNOWN')
        if status == 'RESTORED':
            telomere_bonus = 0.3
        elif status == 'IMPROVED':
            telomere_bonus = 0.2
        elif status == 'PRESERVED':
            telomere_bonus = 0.1
        else:
            telomere_bonus = 0.0
        
        # bonus based on telomere quality
        tel_length = telomere_after.get('total_telomere_length', 0)
        max_repeats = telomere_after.get('max_repeats', 0)
        
        length_factor = min(1.0, tel_length / 5000.0)
        repeats_factor = min(1.0, max_repeats / 100.0)
        telomere_quality = (length_factor * 0.5 + repeats_factor * 0.5)
        
        final_score = base_score * 0.4 + telomere_bonus * 0.3 + telomere_quality * 0.3
    else:
        # results without telomeres get discounted score
        final_score = base_score * 0.5
    
    return final_score


def try_repair_chromosome(chr_id: str, chr_seq: str, 
                         repair_sequences: List[Dict],
                         config: Dict,
                         selection_strategy: str = 'minimal_extension') -> Tuple[Optional[str], Dict, Dict]:
    
    printv(f"\n{'='*60}")
    printv(f"Processing chromosome: {chr_id}")
    printv(f"  Length: {len(chr_seq):,} bp")
    printv(f"  Candidates: {len(repair_sequences)} sequences")
    printv(f"  Strategy: {selection_strategy}")
    printv(f"  Goal: Find sequence that RESTORES telomeres")
    printv(f"{'='*60}")
    
    repairer = TelomereRepair(config)
    
    # store all successful results WITH telomeres
    successful_with_telomere = []
    # store successful results without telomeres (as fallback)
    successful_without_telomere = []
    # store results that damage telomeres (never use)
    damaged_results = []
    
    for i, repair_info in enumerate(repair_sequences):
        printv(f"\n  >>> Attempt {i+1}/{len(repair_sequences)} <<<")
        
        # use repair process with verification
        merged_seq, stats = repairer.repair_with_verification(
            chr_id, chr_seq, repair_info
        )
        
        if merged_seq:
            # build result information
            result_info = {
                'chromosome': chr_id,
                'repair_id': repair_info['id'],
                'repair_type': repair_info.get('repair_type', 'unknown'),
                'overlap_length': stats.get('overlap_len', 0),
                'similarity': stats.get('identity', 0),
                'original_length': len(chr_seq),
                'repaired_length': len(merged_seq),
                'extension': stats.get('extension', 0),
                'attempt_number': i + 1,
                'total_attempts': len(repair_sequences),
                'sequence_length': len(repair_info['sequence']),
                'minimap2_info': {
                    'identity': stats.get('identity', 0),
                    'coverage': stats.get('coverage', 0),
                    'mapq': stats.get('mapq', 0),
                    'strand': stats.get('strand', '+'),
                    'score': stats.get('score', 0),
                },
                'merge_method': stats.get('method', 'unknown'),
                'cut_positions': stats.get('cut_positions', {}),
                'telomere_before': stats.get('telomere_before', {}),
                'telomere_after': stats.get('telomere_after', {}),
                'telomere_comparison': stats.get('telomere_comparison', {}),
                'merged_sequence': merged_seq,
                'stats': stats
            }
            
            tel_status = result_info['telomere_comparison'].get('status', 'UNKNOWN')
            has_tel_after = result_info['telomere_after'].get('has_telomere', False)
            
            # classify based on telomere status
            if tel_status in ['RESTORED', 'IMPROVED', 'PRESERVED']:
                printv(f"    TELOMERE SUCCESS: {tel_status}")
                if has_tel_after:
                    printv(f"       Length: {result_info['telomere_after'].get('total_telomere_length', 0):,}bp")
                    printv(f"       Repeats: {result_info['telomere_after'].get('max_repeats', 0)}")
                successful_with_telomere.append({
                    'score': 0.0,
                    'result_info': result_info,
                    'repair_info': repair_info,
                    'merged_seq': merged_seq
                })
            elif tel_status == 'NO_TELOMERE':
                printv(f"    NO TELOMERES after repair")
                successful_without_telomere.append({
                    'score': 0.0,
                    'result_info': result_info,
                    'repair_info': repair_info,
                    'merged_seq': merged_seq
                })
            elif tel_status == 'DAMAGED':
                printv(f"    TELOMERE DAMAGED - REJECTED")
                damaged_results.append({
                    'result_info': result_info,
                    'repair_info': repair_info
                })
        else:
            error = stats.get('error', 'unknown_error')
            printv(f"    Repair failed: {error}")
    
    # prioritize results with telomeres
    if successful_with_telomere:
        printv(f"\n  Found {len(successful_with_telomere)} sequences that restore/preserve telomeres")
        
        # calculate scores
        for item in successful_with_telomere:
            score = calculate_selection_score(item['result_info'], selection_strategy, has_telomere=True)
            item['score'] = score
            item['result_info']['selection_score'] = score
        
        # sort by score
        successful_with_telomere.sort(key=lambda x: x['score'], reverse=True)
        
        # show candidates
        if VERBOSE:
            printv(f"\n  Candidates with telomeres:")
            for rank, item in enumerate(successful_with_telomere[:3]):
                result = item['result_info']
                tel_status = result['telomere_comparison'].get('status', 'UNKNOWN')
                tel_length = result['telomere_after'].get('total_telomere_length', 0)
                printv(f"    {rank+1}. Score: {item['score']:.3f}, "
                      f"Telomere: {tel_status} ({tel_length:,}bp), "
                      f"Extension: {result['extension']:,}bp")
        
        best_item = successful_with_telomere[0]
        best_item['result_info']['telomere_success'] = True
        
    elif successful_without_telomere:
        # no results with telomeres found, consider fallback without telomeres
        printv(f"\n  WARNING: No sequences restore telomeres!")
        printv(f"  Found {len(successful_without_telomere)} sequences WITHOUT telomeres")
        printv(f"  Using best available, but telomeres will NOT be restored")
        
        for item in successful_without_telomere:
            score = calculate_selection_score(item['result_info'], selection_strategy, has_telomere=False)
            item['score'] = score
            item['result_info']['selection_score'] = score
        
        successful_without_telomere.sort(key=lambda x: x['score'], reverse=True)
        
        best_item = successful_without_telomere[0]
        best_item['result_info']['telomere_success'] = False
        best_item['result_info']['warning'] = 'No telomere restoration achieved'
    else:
        printv(f"\n  No successful repairs at all")
        return None, {
            'chromosome': chr_id,
            'repair_id': 'none',
            'reason': 'all_attempts_failed',
            'attempts': len(repair_sequences)
        }, None
    
    best_result_info = best_item['result_info']
    best_repair_info = best_item['repair_info']
    best_merged_seq = best_item['merged_seq']
    
    printv(f"\n  {'='*40}")
    printv(f"  SELECTED OPTIMAL RESULT:")
    printv(f"    Sequence: {best_result_info['repair_id'][:80]}...")
    printv(f"    Telomere status: {best_result_info['telomere_comparison'].get('status', 'UNKNOWN')}")
    if best_result_info['telomere_after'].get('has_telomere', False):
        printv(f"    Telomere length: {best_result_info['telomere_after'].get('total_telomere_length', 0):,} bp")
        printv(f"    Max repeats: {best_result_info['telomere_after'].get('max_repeats', 0)}")
    printv(f"    Extension: {best_result_info['extension']:,} bp")
    printv(f"    Score: {best_item['score']:.3f}")
    printv(f"  {'='*40}")
    
    best_result_info['selection_method'] = selection_strategy
    best_result_info['alternatives_count'] = len(successful_with_telomere) + len(successful_without_telomere)
    
    # record alternative candidates information
    all_alternatives = []
    for item in successful_with_telomere:
        if item != best_item:
            alt = item['result_info']
            all_alternatives.append({
                'repair_id': alt['repair_id'],
                'score': item['score'],
                'extension': alt['extension'],
                'similarity': alt['similarity'],
                'telomere_status': alt['telomere_comparison'].get('status', 'UNKNOWN'),
                'telomere_length': alt['telomere_after'].get('total_telomere_length', 0)
            })
    for item in successful_without_telomere:
        if item != best_item:
            alt = item['result_info']
            all_alternatives.append({
                'repair_id': alt['repair_id'],
                'score': item['score'],
                'extension': alt['extension'],
                'similarity': alt['similarity'],
                'telomere_status': 'NO_TELOMERE',
                'telomere_length': 0
            })
    
    best_result_info['all_alternatives_scores'] = sorted(all_alternatives, 
                                                         key=lambda x: x['score'], 
                                                         reverse=True)[:10]
    
    if 'merged_sequence' in best_result_info:
        del best_result_info['merged_sequence']
    
    return best_merged_seq, best_result_info, best_repair_info


# ============================================================================
# Main function
# ============================================================================

def merge_telomere_sequences(
    genome_file: str,
    repair_file: str,
    output_dir: str,
    search_range: int = 5000000,
    min_similarity: float = 0.7,
    verbose: bool = False,
    minimap2_mode: str = 'asm5',
    min_match_length: int = 100,
    min_mapq: int = 20,
    selection_strategy: str = 'minimal_extension',
    min_telomere_repeats: int = 20,
    min_telomere_length: int = 500
) -> Dict[str, Any]:
    
    global VERBOSE
    VERBOSE = verbose
    
    start_time = time.time()
    
    valid_strategies = ['first_success', 'minimal_extension', 'balanced']
    if selection_strategy not in valid_strategies:
        return {
            'success': False,
            'error': f"Invalid selection strategy: {selection_strategy}, options: {', '.join(valid_strategies)}"
        }
    
    try:
        if not os.path.exists(genome_file):
            return {
                'success': False,
                'error': f"Genome file does not exist: {genome_file}"
            }
        
        if not os.path.exists(repair_file):
            return {
                'success': False,
                'error': f"Repair sequence file does not exist: {repair_file}"
            }
        
        os.makedirs(output_dir, exist_ok=True)
        
        if verbose:
            print("\n" + "="*70)
            print("Telomere Repair Sequence Merge Script (Telomere-First Version)")
            print("="*70)
            print(f"Genome file: {genome_file}")
            print(f"Repair sequence file: {repair_file}")
            print(f"Output directory: {output_dir}")
            print(f"Search range: {search_range:,} bp")
            print(f"Minimum similarity: {min_similarity}")
            print(f"minimap2 mode: {minimap2_mode}")
            print(f"Minimum match length: {min_match_length}")
            print(f"Minimum MAPQ: {min_mapq}")
            print(f"Selection strategy: {selection_strategy}")
            print(f"Minimum telomere repeats: {min_telomere_repeats}")
            print(f"Minimum telomere length: {min_telomere_length:,} bp")
            print("="*70)
        
        if verbose:
            print(f"\nReading genome file: {genome_file}")
        chromosomes = read_fasta(genome_file)
        if verbose:
            print(f"Found {len(chromosomes)} chromosomes")
        
        if verbose:
            print(f"\nReading repair sequence file: {repair_file}")
        try:
            repair_records = list(SeqIO.parse(repair_file, "fasta"))
        except Exception as e:
            return {
                'success': False,
                'error': f"Cannot read repair sequence file: {e}"
            }
        
        if verbose:
            print(f"Found {len(repair_records)} repair sequences")
        
        if verbose:
            print("\nGrouping repair sequences by chromosome:")
        chromosome_groups = group_repair_sequences_by_chromosome(repair_records)
        
        minimap2_config = {
            'minimap2_mode': minimap2_mode,
            'min_score': min_similarity,
            'min_match_length': min_match_length,
            'min_mapq': min_mapq,
            'search_range': search_range,
            'min_telomere_repeats': min_telomere_repeats,
            'min_telomere_length': min_telomere_length
        }
        
        results = {'success': [], 'failed': [], 'skipped': []}
        repaired_chromosomes = {}
        used_repair_ids = set()
        
        for target_chr, repair_sequences in chromosome_groups.items():
            if not repair_sequences:
                continue
            
            chr_id = find_matching_chromosome(chromosomes, target_chr)
            if not chr_id:
                if verbose:
                    print(f"\nCannot find matching chromosome in genome: '{target_chr}'")
                for repair_info in repair_sequences:
                    results['skipped'].append({
                        'repair_id': repair_info['id'],
                        'reason': 'chromosome_not_found',
                        'target': target_chr
                    })
                continue
            
            if chr_id in repaired_chromosomes:
                if verbose:
                    print(f"\nChromosome {chr_id} already repaired, skipping this group")
                for repair_info in repair_sequences:
                    results['skipped'].append({
                        'repair_id': repair_info['id'],
                        'reason': 'chromosome_already_repaired',
                        'chromosome': chr_id
                    })
                continue
            
            chr_seq = chromosomes[chr_id]
            
            merged_seq, result_info, used_repair_info = try_repair_chromosome(
                chr_id, chr_seq, repair_sequences, minimap2_config, selection_strategy
            )
            
            if merged_seq:
                repaired_chromosomes[chr_id] = merged_seq
                results['success'].append(result_info)
                
                if used_repair_info:
                    used_repair_ids.add(used_repair_info['id'])
                
                # mark skipped alternative sequences
                for repair_info in repair_sequences:
                    if repair_info['id'] != used_repair_info['id']:
                        results['skipped'].append({
                            'repair_id': repair_info['id'],
                            'reason': 'not_used_alternative',
                            'chromosome': chr_id,
                            'used_repair': used_repair_info['id'],
                            'selection_strategy': selection_strategy,
                            'sequence_length': repair_info['length']
                        })
            else:
                results['failed'].append(result_info)
        
        # process ungrouped repair sequences
        all_repair_ids = {record.id for record in repair_records}
        unused_repair_ids = all_repair_ids - used_repair_ids
        
        for repair_id in unused_repair_ids:
            already_processed = False
            for result_list in results.values():
                for result in result_list:
                    if isinstance(result, dict) and result.get('repair_id') == repair_id:
                        already_processed = True
                        break
                if already_processed:
                    break
            
            if not already_processed:
                results['skipped'].append({
                    'repair_id': repair_id,
                    'reason': 'unrecognized_chromosome_or_not_grouped'
                })
        
        if verbose:
            print(f"\n{'='*70}")
            print("Saving repair results...")
            print('='*70)
        
        repaired_file = os.path.join(output_dir, "repaired_genome.fasta")
        report_file = os.path.join(output_dir, "repair_report.txt")
        
        repaired_records = []
        
        for chr_id, seq in chromosomes.items():
            if chr_id in repaired_chromosomes:
                repaired_seq = repaired_chromosomes[chr_id]
                record = SeqRecord(
                    Seq(repaired_seq),
                    id=chr_id,
                    description=f"repaired|strategy:{selection_strategy}|original_{len(seq)}bp|repaired_{len(repaired_seq)}bp"
                )
                repaired_records.append(record)
            else:
                record = SeqRecord(
                    Seq(seq),
                    id=chr_id,
                    description=f"original_{len(seq)}bp"
                )
                repaired_records.append(record)
        
        try:
            SeqIO.write(repaired_records, repaired_file, "fasta")
            if verbose:
                print(f"Repaired genome saved to: {repaired_file}")
        except Exception as e:
            return {
                'success': False,
                'error': f"Cannot save repaired genome file: {e}"
            }
        
        # generate enhanced report
        try:
            with open(report_file, 'w') as f:
                f.write("Telomere Repair Sequence Merge Report (Telomere-First Version)\n")
                f.write("="*70 + "\n")
                f.write(f"Genome file: {genome_file}\n")
                f.write(f"Repair sequence file: {repair_file}\n")
                f.write(f"Output directory: {output_dir}\n")
                f.write(f"Processing time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("Parameter settings:\n")
                f.write(f"  Search range: {search_range:,} bp\n")
                f.write(f"  Minimum similarity: {min_similarity}\n")
                f.write(f"  minimap2 mode: {minimap2_mode}\n")
                f.write(f"  Minimum match length: {min_match_length}\n")
                f.write(f"  Minimum MAPQ: {min_mapq}\n")
                f.write(f"  Selection strategy: {selection_strategy}\n")
                f.write(f"  Minimum telomere repeats: {min_telomere_repeats}\n")
                f.write(f"  Minimum telomere length: {min_telomere_length:,} bp\n\n")
                
                total_success = len([r for r in results['success'] if 'attempt_number' in r])
                total_failed = len([r for r in results['failed'] if isinstance(r, dict) and r.get('reason') == 'all_attempts_failed'])
                total_skipped = len(results['skipped'])
                
                f.write("Processing results:\n")
                f.write(f"  Total repair sequences: {len(repair_records)}\n")
                f.write(f"  Successfully repaired chromosomes: {total_success}\n")
                f.write(f"  Failed repair sequences: {total_failed}\n")
                f.write(f"  Skipped sequences: {total_skipped}\n\n")
                
                if results['success']:
                    f.write("Successfully repaired chromosomes:\n")
                    f.write("-"*70 + "\n")
                    for result in results['success']:
                        if 'attempt_number' in result:
                            f.write(f"Chromosome: {result['chromosome']}\n")
                            f.write(f"  Used repair sequence: {result['repair_id']}\n")
                            f.write(f"  Repair type: {result['repair_type']}\n")
                            f.write(f"  Selection strategy: {result.get('selection_method', selection_strategy)}\n")
                            f.write(f"  Selection score: {result.get('selection_score', 0):.3f}\n")
                            f.write(f"  Merge method: {result.get('merge_method', 'unknown')}\n")
                            f.write(f"  Original length: {result['original_length']:,} bp\n")
                            f.write(f"  Repaired length: {result['repaired_length']:,} bp\n")
                            f.write(f"  Extension length: {result['extension']:,} bp\n")
                            f.write(f"  Overlap length: {result['overlap_length']:,} bp\n")
                            f.write(f"  Similarity: {result['similarity']:.4f}\n")
                            f.write(f"  MAPQ: {result['minimap2_info']['mapq']}\n")
                            f.write(f"  Strand: {result['minimap2_info']['strand']}\n")
                            
                            # telomere validation information
                            tel_comparison = result.get('telomere_comparison', {})
                            tel_after = result.get('telomere_after', {})
                            
                            f.write(f"\n  Telomere Validation:\n")
                            f.write(f"    Status: {tel_comparison.get('status', 'UNKNOWN')}\n")
                            
                            if tel_after.get('has_telomere', False):
                                f.write(f"    Telomeres PRESENT after repair\n")
                                f.write(f"      Total length: {tel_after.get('total_telomere_length', 0):,} bp\n")
                                f.write(f"      Max repeats: {tel_after.get('max_repeats', 0)}\n")
                                f.write(f"      Telomere count: {tel_after.get('telomere_count', 0)}\n")
                                f.write(f"      Density: {tel_after.get('telomere_density', 0):.2%}\n")
                                
                                if tel_comparison.get('status') == 'RESTORED':
                                    f.write(f"      TELOMERES RESTORED (were absent, now present)\n")
                                elif tel_comparison.get('status') == 'IMPROVED':
                                    details = tel_comparison.get('details', {})
                                    f.write(f"      TELOMERES IMPROVED: +{details.get('length_increase', 0):,}bp\n")
                                elif tel_comparison.get('status') == 'PRESERVED':
                                    f.write(f"      Telomeres preserved (unchanged)\n")
                            else:
                                f.write(f"    WARNING: No telomeres detected after repair\n")
                                if 'warning' in result:
                                    f.write(f"      {result['warning']}\n")
                            
                            f.write(f"\n  Cut positions: repair sequence{result['cut_positions'].get('repair_cut', 'N/A')}, "
                                  f"chromosome{result['cut_positions'].get('chr_cut', 'N/A')}\n")
                            f.write(f"  Attempt number: {result['attempt_number']}/{result['total_attempts']}\n")
                            f.write(f"  Alternative count: {result.get('alternatives_count', 'N/A')}\n")
                            f.write(f"  Sequence length: {result.get('sequence_length', 'N/A'):,} bp\n")
                            
                            if 'all_alternatives_scores' in result and result['all_alternatives_scores']:
                                f.write(f"  Alternative scores:\n")
                                for alt in result['all_alternatives_scores'][:5]:
                                    f.write(f"    - {alt['repair_id'][:50]}...: score {alt['score']:.3f}, "
                                           f"telomere {alt.get('telomere_status', 'UNKNOWN')}, "
                                           f"extension {alt['extension']:,}bp\n")
                                if len(result['all_alternatives_scores']) > 5:
                                    f.write(f"    ... {len(result['all_alternatives_scores'])-5} more alternatives\n")
                            
                            f.write("\n")
                
                if results['failed']:
                    f.write("\nFailed repair sequences:\n")
                    f.write("-"*50 + "\n")
                    for fail in results['failed']:
                        if isinstance(fail, dict) and fail.get('repair_id') != 'none':
                            f.write(f"Chromosome: {fail.get('chromosome', 'unknown')}\n")
                            f.write(f"  Failure reason: {fail.get('reason', 'unknown')}\n")
                            f.write(f"  Attempts: {fail.get('attempts', 0)}\n")
                            f.write("\n")
                
                if results['skipped']:
                    f.write("\nSkipped sequences:\n")
                    f.write("-"*50 + "\n")
                    for skip in results['skipped']:
                        f.write(f"Repair sequence: {skip['repair_id']}\n")
                        f.write(f"  Skip reason: {skip['reason']}\n")
                        if 'target' in skip:
                            f.write(f"  Target chromosome: {skip['target']}\n")
                        if 'chromosome' in skip:
                            f.write(f"  Corresponding chromosome: {skip['chromosome']}\n")
                        if 'used_repair' in skip:
                            f.write(f"  Used sequence: {skip['used_repair']}\n")
                        if 'sequence_length' in skip:
                            f.write(f"  Sequence length: {skip['sequence_length']:,} bp\n")
                        f.write("\n")
            
            if verbose:
                print(f"Detailed report saved to: {report_file}")
        except Exception as e:
            print(f"Warning: Cannot save text report: {e}")
        
        total_time = time.time() - start_time
        
        # statistics
        success_chromosomes = len([r for r in results['success'] if 'attempt_number' in r])
        
        total_extension = 0
        restored_count = 0
        improved_count = 0
        preserved_count = 0
        no_telomere_count = 0
        
        if results['success']:
            for r in results['success']:
                if 'extension' in r:
                    total_extension += r.get('extension', 0)
                
                tel_status = r.get('telomere_comparison', {}).get('status', '')
                if tel_status == 'RESTORED':
                    restored_count += 1
                elif tel_status == 'IMPROVED':
                    improved_count += 1
                elif tel_status == 'PRESERVED':
                    preserved_count += 1
                elif tel_status == 'NO_TELOMERE':
                    no_telomere_count += 1
        
        avg_extension = total_extension / success_chromosomes if success_chromosomes > 0 else 0
        
        result_dict = {
            'success': True,
            'repaired_genome': repaired_file,
            'report_file': report_file,
            'results': results,
            'statistics': {
                'total_repair_sequences': len(repair_records),
                'success_chromosomes': success_chromosomes,
                'failed_sequences': total_failed,
                'skipped_sequences': total_skipped,
                'total_extension': total_extension,
                'avg_extension': avg_extension,
                'selection_strategy': selection_strategy,
                'telomere_stats': {
                    'restored': restored_count,
                    'improved': improved_count,
                    'preserved': preserved_count,
                    'no_telomere': no_telomere_count
                },
                'processing_time': total_time
            },
            'output_dir': output_dir
        }
        
        if verbose:
            print("\n" + "="*70)
            print("Processing completed!")
            print("="*70)
            print(f"Successfully repaired: {success_chromosomes} chromosomes")
            print(f"  Telomere status:")
            print(f"    - RESTORED: {restored_count}")
            print(f"    - IMPROVED: {improved_count}")
            print(f"    - PRESERVED: {preserved_count}")
            print(f"    - NO_TELOMERE: {no_telomere_count}")
            print(f"Repair failed: {total_failed} sequences")
            print(f"Skipped: {total_skipped} sequences")
            print(f"Total time: {total_time:.2f} seconds")
            
            if success_chromosomes > 0:
                print(f"\nTotal extension length: {total_extension:,} bp")
                print(f"Average extension length: {avg_extension:,.0f} bp")
                print(f"Repaired genome: {repaired_file}")
            
            print(f"Detailed report: {report_file}")
            print("="*70)
        
        return result_dict
        
    except Exception as e:
        import traceback
        error_trace = traceback.format_exc()
        if verbose:
            print(f"Error during merge process: {e}")
            print(error_trace)
        
        return {
            'success': False,
            'error': str(e),
            'traceback': error_trace
        }


def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description="Telomere Repair Sequence Merge Script (Telomere-First Version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Core Logic: PRIORITIZE SEQUENCES THAT RESTORE TELOMERES
        
        Features:
          1. Try ALL candidates, not just first success
          2. Validate telomeres BEFORE and AFTER repair
          3. ONLY accept sequences that restore/preserve telomeres
          4. Fall back to best available only if no telomere-restoring sequences found
        
        Selection strategies:
          - first_success: Stop at first successful result (NOT recommended)
          - minimal_extension: Minimal change strategy (default)
          - balanced: Balanced strategy (consider both extension and quality)
        
        Usage examples:
          python add_te.py -q genome.fasta -i repair_sequences.fa -o output_dir
        """)
    )
    
    parser.add_argument("-q", "--query", required=True, help="Genome FASTA file")
    parser.add_argument("-i", "--input", required=True, help="Repair sequences FASTA file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    
    parser.add_argument("--search-range", type=int, default=5000000,
                       help="Alignment search range(bp) (default: 5,000,000)")
    parser.add_argument("--min-similarity", type=float, default=0.6,
                       help="Minimum alignment similarity (default: 0.6)")
    
    parser.add_argument("--minimap2-mode", default='asm5',
                       choices=['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb'],
                       help="minimap2 alignment mode (default: asm5)")
    parser.add_argument("--min-match-length", type=int, default=100,
                       help="Minimum match length (default: 100)")
    parser.add_argument("--min-mapq", type=int, default=20,
                       help="Minimum mapping quality (default: 20)")
    
    parser.add_argument("--strategy", default='minimal_extension',
                       choices=['first_success', 'minimal_extension', 'balanced'],
                       help="Selection strategy (default: minimal_extension)")
    
    parser.add_argument("--min-telomere-repeats", type=int, default=20,
                       help="Minimum telomere repeats for detection (default: 20)")
    parser.add_argument("--min-telomere-length", type=int, default=500,
                       help="Minimum telomere length in bp (default: 500)")
    
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed output")
    
    return parser.parse_args()


def main():
    args = parse_command_line_args()
    
    result = merge_telomere_sequences(
        genome_file=args.query,
        repair_file=args.input,
        output_dir=args.output_dir,
        search_range=args.search_range,
        min_similarity=args.min_similarity,
        verbose=args.verbose,
        minimap2_mode=args.minimap2_mode,
        min_match_length=args.min_match_length,
        min_mapq=args.min_mapq,
        selection_strategy=args.strategy,
        min_telomere_repeats=args.min_telomere_repeats,
        min_telomere_length=args.min_telomere_length
    )
    
    if not result['success']:
        print(f"\nError: {result['error']}")
        if args.verbose and 'traceback' in result:
            print(f"Detailed error information:\n{result['traceback']}")
        sys.exit(1)
    else:
        stats = result['statistics']
        tel_stats = stats['telomere_stats']
        
        print(f"\n{'='*50}")
        print("PROCESSING COMPLETED")
        print(f"{'='*50}")
        print(f"Successfully repaired: {stats['success_chromosomes']} chromosomes")
        print(f"\nTelomere Status:")
        print(f"  RESTORED:  {tel_stats['restored']}  (was absent, now present)")
        print(f"  IMPROVED:  {tel_stats['improved']}  (longer or more repeats)")
        print(f"  PRESERVED: {tel_stats['preserved']} (unchanged)")
        print(f"  NO TELOMERE: {tel_stats['no_telomere']} (merged but no telomeres)")
        
        if tel_stats['restored'] > 0 or tel_stats['improved'] > 0:
            print(f"\nSUCCESS: Telomeres were restored/improved!")
        elif tel_stats['preserved'] > 0:
            print(f"\nExisting telomeres were preserved")
        else:
            print(f"\nWARNING: No telomeres were restored")
            print(f"   Check if repair sequences actually contain telomere repeats")
        
        print(f"\nTotal extension length: {stats.get('total_extension', 0):,} bp")
        print(f"Detailed report: {result['report_file']}")
        print(f"{'='*50}")


if __name__ == "__main__":
    main()