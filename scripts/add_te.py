#!/usr/bin/env python3
"""
Telomere Repair Sequence Merge Script (Telomere-First Version)
Modified for Plant Telomeres with Tabular Report
Purpose: Merge extracted telomere repair sequences to corresponding chromosome ends
Core Logic: PRIORITIZE SEQUENCES THAT RESTORE TELOMERES
Features: 
  1. Try all candidates, not just first success
  2. Validate telomeres BEFORE and AFTER repair
  3. ONLY accept sequences that restore/preserve telomeres
  4. Fall back to best available only if no telomere-restoring sequences found
  5. Support for both plant (TTTAGGG) and vertebrate (TTAGGG) telomeres
  6. Overlap length considered in scoring for reliable fusion
  7. Tabular report format for clear visualization
"""

import sys
import os
import re
import argparse
import textwrap
import subprocess
import tempfile
import time
from typing import Dict, List, Tuple, Optional, Any, Set, Union
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VERBOSE = False

def printv(message: str):
    if VERBOSE:
        print(message)

# ============================================================================
# Custom Exceptions for Error Handling
# ============================================================================

class RepairError(Exception):
    """Base class for repair errors"""
    pass

class NoOverlapError(RepairError):
    """No overlap found between sequences"""
    pass

class LowQualityError(RepairError):
    """Alignment quality too low"""
    pass

class TelomereNotRestoredError(RepairError):
    """Failed to restore telomeres"""
    pass

# ============================================================================
# Telomere detection module - Modified for Plant Telomeres
# ============================================================================

class TelomereDetector:
    """Telomere detector: used to detect telomere repeat units in sequences"""
    
    # Standard telomere patterns - supporting both plant and vertebrate
    TELOMERE_PATTERNS = {
        'vertebrate': 'TTAGGG',      # vertebrate telomere
        'vertebrate_rev': 'CCCTAA',  # reverse complement
        'plant': 'TTTAGGG',          # plant telomere (Arabidopsis, etc.)
        'plant_rev': 'CCCTAAA',      # plant telomere reverse complement
        'alternative': ['TTAGGG', 'CCCTAA', 'TTAGGC', 'TTGGGG', 'TTTAGGG', 'CCCTAAA'],  # all variants
    }
    
    def __init__(self, min_repeats: int = 20, organism: str = 'plant'):
        """
        Initialize telomere detector
        
        Args:
            min_repeats: minimum number of repeats
            organism: organism type ('plant', 'vertebrate', 'auto')
        """
        self.min_repeats = min_repeats
        self.organism = organism
        
        # Set detection patterns based on organism
        if organism == 'plant':
            self.patterns = [
                self.TELOMERE_PATTERNS['plant'],      # TTTAGGG
                self.TELOMERE_PATTERNS['plant_rev'],  # CCCTAAA
            ]
        elif organism == 'vertebrate':
            self.patterns = [
                self.TELOMERE_PATTERNS['vertebrate'],
                self.TELOMERE_PATTERNS['vertebrate_rev']
            ]
        else:  # auto - detect all patterns
            self.patterns = [
                self.TELOMERE_PATTERNS['plant'],
                self.TELOMERE_PATTERNS['plant_rev'],
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
                'sequence_length': 0,
                'telomere_patterns_found': []
            }
        
        sequence = sequence.upper()
        seq_len = len(sequence)
        telomeres = []
        patterns_found = set()
        
        # detect all telomere patterns
        for pattern in self.patterns:
            pattern_telomeres = self._detect_pattern(sequence, pattern, 
                                                     f"pattern_{pattern}")
            telomeres.extend(pattern_telomeres)
            if pattern_telomeres:
                patterns_found.add(pattern)
        
        # sort by position
        telomeres.sort(key=lambda x: x['start'])
        
        # merge adjacent telomere regions
        merged_telomeres = self._merge_adjacent_telomeres(telomeres)
        
        # calculate statistics
        total_length = sum(t['length'] for t in merged_telomeres)
        max_repeats = max((t['repeats'] for t in merged_telomeres), default=0)
        max_repeat_region = next((t for t in merged_telomeres if t['repeats'] == max_repeats), None)
        
        # determine primary telomere type
        telomere_type = 'unknown'
        if merged_telomeres:
            # count pattern frequencies
            type_count = defaultdict(int)
            for t in telomeres:
                type_count[t['pattern']] += 1
            # most common pattern as primary type
            if type_count:
                most_common = max(type_count.items(), key=lambda x: x[1])[0]
                if most_common in [self.TELOMERE_PATTERNS['plant'], self.TELOMERE_PATTERNS['plant_rev']]:
                    telomere_type = 'plant'
                elif most_common in [self.TELOMERE_PATTERNS['vertebrate'], self.TELOMERE_PATTERNS['vertebrate_rev']]:
                    telomere_type = 'vertebrate'
        
        result = {
            'has_telomere': len(merged_telomeres) > 0,
            'telomere_count': len(merged_telomeres),
            'telomeres': merged_telomeres,
            'total_telomere_length': total_length,
            'max_repeats': max_repeats,
            'max_repeat_region': max_repeat_region,
            'region_name': region_name,
            'sequence_length': seq_len,
            'telomere_density': total_length / seq_len if seq_len > 0 else 0,
            'telomere_type': telomere_type,
            'patterns_found': list(patterns_found)
        }
        
        return result
    
    def _detect_pattern(self, sequence: str, pattern: str, 
                       pattern_type: str) -> List[Dict[str, Any]]:
        """
        Detect specific telomere repeat pattern with strict matching
        
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
                
                # extend forward to detect consecutive repeats - strict matching
                while j <= seq_len - pattern_len:
                    if sequence[j:j+pattern_len] == pattern:
                        repeat_count += 1
                        j += pattern_len
                    else:
                        break
                
                # check if minimum repeat count is reached
                if repeat_count >= self.min_repeats:
                    end_pos = i + repeat_count * pattern_len
                    length = end_pos - i
                    
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
                                  max_gap: int = 50) -> List[Dict]:
        """
        Merge adjacent telomere regions with proper repeat counting
        
        Args:
            telomeres: list of telomere regions
            max_gap: maximum allowed gap (reduced to 50bp for stricter merging)
            
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
                # merge regions - repeats add up (they are non-overlapping)
                current['repeats'] += t['repeats']
                current['end'] = max(current['end'], t['end'])
                current['length'] = current['end'] - current['start']
                current['type'] = f"{current['type']}_merged"
                current['sequence'] = "merged_region"
            else:
                merged.append(current)
                current = t.copy()
        
        merged.append(current)
        return merged
    
    def validate_telomere_quality(self, telomere_region: Dict) -> Dict[str, Any]:
        """
        Validate quality of telomere region
        
        Args:
            telomere_region: telomere region dictionary
            
        Returns:
            quality assessment dictionary
        """
        sequence = telomere_region.get('sequence', '')
        pattern = telomere_region.get('pattern', '')
        pattern_len = len(pattern)
        
        if not sequence or not pattern or sequence == "merged_region":
            return {'valid': False, 'reason': 'missing_data'}
        
        # check repeat consistency
        repeats = len(sequence) // pattern_len
        mismatches = 0
        for i in range(repeats):
            start = i * pattern_len
            unit = sequence[start:start+pattern_len]
            if unit != pattern:
                # count mismatches
                diff = sum(1 for a, b in zip(unit, pattern) if a != b)
                mismatches += diff
        
        mismatch_rate = mismatches / (repeats * pattern_len) if repeats > 0 else 1.0
        
        # check GC content (telomeres are typically AT-rich)
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence) if len(sequence) > 0 else 0
        
        # plant telomere TTTAGGG has ~28.6% GC
        expected_gc_range = (0.20, 0.40)  # 20-40%
        
        return {
            'valid': mismatch_rate < 0.05,  # mismatch rate < 5%
            'mismatch_rate': mismatch_rate,
            'gc_content': gc_content,
            'gc_within_range': expected_gc_range[0] <= gc_content <= expected_gc_range[1],
            'repeats': repeats,
            'quality_score': 1.0 - min(1.0, mismatch_rate * 2)  # 0-1 score
        }
    
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
            'telomere_type_changed': False,
            'details': {}
        }
        
        # check if telomere type changed
        if before.get('telomere_type') != after.get('telomere_type'):
            comparison['telomere_type_changed'] = True
            comparison['telomere_type_before'] = before.get('telomere_type')
            comparison['telomere_type_after'] = after.get('telomere_type')
        
        # Since we know chromosomes lack telomeres before repair, we only care about:
        # Case 1: no telomeres before, telomeres after (best - RESTORED)
        if not before['has_telomere'] and after['has_telomere']:
            comparison['status'] = 'RESTORED'
            comparison['telomere_restored'] = True
            comparison['details'] = {
                'new_telomeres': after['telomere_count'],
                'max_repeats': after['max_repeats'],
                'total_length': after['total_telomere_length'],
                'density': after['telomere_density'],
                'telomere_type': after.get('telomere_type')
            }
        
        # Case 2: no telomeres before, no telomeres after (failure)
        elif not before['has_telomere'] and not after['has_telomere']:
            comparison['status'] = 'NO_TELOMERE'
            comparison['details'] = {
                'message': 'No telomeres detected after repair'
            }
        
        # Case 3: telomeres before (should not happen with our assumption)
        elif before['has_telomere']:
            comparison['status'] = 'PREEXISTING_TELOMERE'
            comparison['details'] = {
                'message': 'Chromosome already had telomeres before repair'
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
# Telomere repair main class - Modified with organism parameter
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
            'organism': 'plant'  # default to plant
        }
        
        self.matcher = Minimap2Matcher(
            mode=self.config['minimap2_mode'],
            min_score=self.config['min_score'],
            min_match_length=self.config['min_match_length'],
            min_mapq=self.config['min_mapq']
        )
        
        # initialize telomere detector with organism
        self.telomere_detector = TelomereDetector(
            min_repeats=self.config.get('min_telomere_repeats', 20),
            organism=self.config.get('organism', 'plant')
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
                                repair_info: Dict) -> Tuple[str, Dict]:
        """
        Repair process with telomere verification - uses exceptions for error cases
        
        Args:
            chr_id: chromosome ID
            chr_seq: chromosome sequence
            repair_info: repair sequence information
            
        Returns:
            (repaired sequence, repair statistics)
            
        Raises:
            NoOverlapError: if no overlap found
            LowQualityError: if alignment quality too low
            TelomereNotRestoredError: if telomeres not restored after merge
        """
        repair_seq = repair_info['sequence']
        repair_type = repair_info.get('repair_type', 'unknown')
        repair_id = repair_info['id']
        
        printv(f"\n  Attempting repair with sequence: {repair_id[:50]}...")
        printv(f"    Length: {len(repair_seq):,} bp")
        printv(f"    Type: {repair_type}")
        
        # 1. detect telomeres before repair (for logging only)
        before_telomeres = self.telomere_detector.detect(
            chr_seq, region_name=f"{chr_id}_before"
        )
        
        # 2. find overlap region
        match = self.find_telomere_overlap(chr_seq, repair_seq, repair_type)
        
        if not match:
            raise NoOverlapError("No overlap found between sequences")
        
        # 3. validate alignment quality
        is_valid, valid_msg = self.validate_match(
            match, repair_type, len(chr_seq), len(repair_seq)
        )
        
        if not is_valid:
            raise LowQualityError(f"Low quality alignment: {valid_msg}")
        
        # 4. perform merge
        merged_seq, stats = self.merge_sequences(chr_seq, repair_seq, match, repair_type)
        
        # 5. detect telomeres after repair
        after_telomeres = self.telomere_detector.detect(
            merged_seq, region_name=f"{chr_id}_after"
        )
        
        # 6. verify telomeres were restored
        if not after_telomeres.get('has_telomere', False):
            raise TelomereNotRestoredError("Failed to restore telomeres")
        
        # 7. verify telomere quality
        if after_telomeres.get('max_repeats', 0) < self.config.get('min_telomere_repeats', 20):
            raise TelomereNotRestoredError(
                f"Insufficient telomere repeats: {after_telomeres.get('max_repeats')} < {self.config.get('min_telomere_repeats')}"
            )
        
        # 8. add telomere information to statistics
        stats['telomere_before'] = before_telomeres
        stats['telomere_after'] = after_telomeres
        stats['telomere_comparison'] = self.telomere_detector.compare_telomere_status(
            before_telomeres, after_telomeres
        )
        
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
    """
    Robust parsing of repair sequence ID
    
    Supports formats:
    - prefix_chr1_5prime
    - prefix_chr01_3prime
    - prefix_chromosome_1_5p
    - prefix_scaffold_5_3prime
    - etc.
    """
    info = {
        'target_chromosome': None,
        'repair_type': 'unknown'
    }
    
    original_id = seq_id
    
    # 1. identify repair type
    type_patterns = [
        (r'[_-]?5prime', '5prime'),
        (r'[_-]?3prime', '3prime'),
        (r'[_-]?5p', '5prime'),
        (r'[_-]?3p', '3prime'),
        (r'[_-]?five_prime', '5prime'),
        (r'[_-]?three_prime', '3prime')
    ]
    
    for pattern, repair_type in type_patterns:
        if re.search(pattern, seq_id, re.IGNORECASE):
            info['repair_type'] = repair_type
            # remove type marker for chromosome extraction
            seq_id = re.sub(pattern, '', seq_id, flags=re.IGNORECASE)
            break
    
    # 2. identify chromosome name
    # match common chromosome naming patterns
    chr_patterns = [
        r'(?:chr|chromosome|scaffold|contig)[_\s]?([0-9]+|[A-Za-z])',
        r'([0-9]+)[_\s]?(?:chr|chromosome)',
        r'^([0-9]+)$',
        r'chr([0-9]+|[A-Za-z])',
        r'chromosome_?([0-9]+)',
        r'scaffold_?([0-9]+)',
        r'contig_?([0-9]+)'
    ]
    
    for pattern in chr_patterns:
        match = re.search(pattern, seq_id, re.IGNORECASE)
        if match:
            # reconstruct full chromosome name
            chr_num = match.group(1)
            # try to preserve original naming style
            if 'chr' in seq_id.lower():
                info['target_chromosome'] = f"chr{chr_num}"
            elif 'chromosome' in seq_id.lower():
                info['target_chromosome'] = f"chromosome_{chr_num}"
            elif 'scaffold' in seq_id.lower():
                info['target_chromosome'] = f"scaffold_{chr_num}"
            elif 'contig' in seq_id.lower():
                info['target_chromosome'] = f"contig_{chr_num}"
            else:
                info['target_chromosome'] = chr_num
            break
    
    # if still no match, try to extract any number as chromosome
    if not info['target_chromosome']:
        num_match = re.search(r'(\d+)', original_id)
        if num_match:
            info['target_chromosome'] = num_match.group(1)
    
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
    """
    Improved chromosome matching logic
    
    Args:
        chromosomes: dictionary of chromosome sequences keyed by ID
        target_chr: target chromosome name from repair sequence
        
    Returns:
        matching chromosome ID or None
    """
    if not target_chr:
        return None
    
    # normalize target chromosome name
    target_norm = re.sub(r'^(chr|chromosome|scaffold|contig)[_\s]?', '', 
                         target_chr, flags=re.IGNORECASE)
    target_norm = target_norm.lower().strip()
    
    # try exact match
    for chr_id in chromosomes.keys():
        chr_clean = re.sub(r'^(chr|chromosome|scaffold|contig)[_\s]?', '', 
                          chr_id.split()[0], flags=re.IGNORECASE)
        chr_clean = chr_clean.lower().strip()
        
        if target_norm == chr_clean:
            return chr_id
    
    # try contains match
    for chr_id in chromosomes.keys():
        chr_clean = re.sub(r'^(chr|chromosome|scaffold|contig)[_\s]?', '', 
                          chr_id.split()[0], flags=re.IGNORECASE)
        chr_clean = chr_clean.lower().strip()
        
        if target_norm in chr_clean or chr_clean in target_norm:
            return chr_id
    
    # try number match
    target_num = re.search(r'\d+', target_norm)
    if target_num:
        target_num = target_num.group()
        for chr_id in chromosomes.keys():
            chr_num = re.search(r'\d+', chr_id)
            if chr_num and chr_num.group() == target_num:
                return chr_id
    
    return None


def calculate_selection_score(result_info: Dict, has_telomere: bool) -> float:
    """
    Calculate selection score - ensures telomere-restoring sequences always win
    MODIFIED: Increased overlap length weight for better fusion
    
    Args:
        result_info: result information
        has_telomere: whether result has telomeres
    
    Returns:
        score (higher is better)
    """
    if has_telomere:
        # Results with telomeres: score based on telomere quality and overlap
        tel_after = result_info.get('telomere_after', {})
        tel_length = tel_after.get('total_telomere_length', 0)
        max_repeats = tel_after.get('max_repeats', 0)
        
        # Telomere quality score (0-1)
        length_score = min(1.0, tel_length / 5000.0)  # 5000bp = full score
        repeats_score = min(1.0, max_repeats / 100.0)  # 100 repeats = full score
        telomere_score = length_score * 0.6 + repeats_score * 0.4
        
        # Fusion quality score - INCREASED WEIGHT
        overlap = result_info.get('overlap_length', 0)
        repair_len = result_info.get('sequence_length', 1)
        overlap_score = min(1.0, overlap / min(repair_len, 500))  # 500bp overlap = full score
        
        # ALSO consider relative overlap (percentage of repair sequence)
        relative_overlap = overlap / repair_len if repair_len > 0 else 0
        relative_overlap_score = min(1.0, relative_overlap * 2)  # 50% overlap = full score
        
        identity = result_info.get('similarity', 0)
        
        # NEW WEIGHTS: Overlap 50%, Telomere 40%, Identity 10%
        final_score = (telomere_score * 0.4 + 
                      (overlap_score * 0.3 + relative_overlap_score * 0.2) +  # 50% total for overlap
                      identity * 0.1)
        
        if VERBOSE:
            print(f"    Score calculation (with telomere):")
            print(f"      Telomere length: {tel_length}bp -> {length_score:.3f}")
            print(f"      Telomere repeats: {max_repeats} -> {repeats_score:.3f}")
            print(f"      Telomere score: {telomere_score:.3f}")
            print(f"      Absolute overlap: {overlap}bp -> {overlap_score:.3f}")
            print(f"      Relative overlap: {relative_overlap:.1%} -> {relative_overlap_score:.3f}")
            print(f"      Identity: {identity:.3f}")
            print(f"      Final score: {final_score:.3f}")
        
        return final_score
    else:
        # Results without telomeres: low base score, always < any telomere result
        overlap = result_info.get('overlap_length', 0)
        repair_len = result_info.get('sequence_length', 1)
        overlap_score = min(1.0, overlap / min(repair_len, 500))
        relative_overlap = overlap / repair_len if repair_len > 0 else 0
        relative_overlap_score = min(1.0, relative_overlap * 2)
        identity = result_info.get('similarity', 0)
        
        # Ensure this is always less than the minimum telomere score (0.4)
        base_score = ((overlap_score * 0.6 + relative_overlap_score * 0.3 + identity * 0.1)) * 0.3
        
        if VERBOSE:
            print(f"    Score calculation (no telomere, discounted):")
            print(f"      Base score: {base_score:.3f}")
        
        return base_score


def try_repair_chromosome(chr_id: str, chr_seq: str, 
                         repair_sequences: List[Dict],
                         config: Dict,
                         selection_strategy: str = 'minimal_extension') -> Tuple[Optional[str], Dict, Dict, List[Dict]]:
    """
    Try to repair chromosome with all candidate sequences
    Strictly prioritizes sequences that restore telomeres
    
    Returns:
        (merged_seq, best_result_info, used_repair_info, all_candidates_info)
    """
    printv(f"\n{'='*60}")
    printv(f"Processing chromosome: {chr_id}")
    printv(f"  Length: {len(chr_seq):,} bp")
    printv(f"  Candidates: {len(repair_sequences)} sequences")
    printv(f"  Goal: Find sequence that RESTORES telomeres")
    printv(f"{'='*60}")
    
    repairer = TelomereRepair(config)
    
    # store results by category
    telomere_restored = []      # successfully restored telomeres
    failed_attempts = []        # failed attempts (no telomere or other errors)
    
    # store all candidates for tabular display
    all_candidates = []
    
    for i, repair_info in enumerate(repair_sequences):
        printv(f"\n  >>> Attempt {i+1}/{len(repair_sequences)} <<<")
        
        try:
            # use repair process with verification
            merged_seq, stats = repairer.repair_with_verification(
                chr_id, chr_seq, repair_info
            )
            
            # If we get here, telomeres were successfully restored
            tel_after = stats.get('telomere_after', {})
            
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
                'telomere_after': tel_after,
                'telomere_comparison': stats.get('telomere_comparison', {}),
                'merged_sequence': merged_seq,
                'stats': stats
            }
            
            # add to restored list
            telomere_restored.append({
                'merged_seq': merged_seq,
                'stats': stats,
                'result_info': result_info,
                'repair_info': repair_info,
                'rank': i + 1,
                'telomere_length': tel_after.get('total_telomere_length', 0),
                'max_repeats': tel_after.get('max_repeats', 0)
            })
            
            # add to candidates list for display
            all_candidates.append({
                'rank': i + 1,
                'repair_id_full': repair_info['id'],
                'repair_id': repair_info['id'][:40] + "..." if len(repair_info['id']) > 40 else repair_info['id'],
                'repair_type': repair_info.get('repair_type', 'unknown'),
                'overlap_length': stats.get('overlap_len', 0),
                'similarity': f"{stats.get('identity', 0):.3f}",
                'mapq': stats.get('mapq', 0),
                'extension': stats.get('extension', 0),
                'telomere_before': 'NO',  # by assumption
                'telomere_after': 'YES',
                'telomere_status': 'RESTORED',
                'telomere_length': tel_after.get('total_telomere_length', 0),
                'telomere_type': tel_after.get('telomere_type', 'unknown'),
                'selected': False
            })
            
            printv(f"    ✓ TELOMERES RESTORED: {tel_after.get('total_telomere_length', 0):,}bp, {tel_after.get('max_repeats', 0)} repeats")
            
        except TelomereNotRestoredError as e:
            # merged but no telomeres
            printv(f"    ✗ Failed: {e}")
            failed_attempts.append({
                'repair_info': repair_info,
                'reason': str(e),
                'rank': i + 1
            })
            
            all_candidates.append({
                'rank': i + 1,
                'repair_id': repair_info['id'][:40] + "..." if len(repair_info['id']) > 40 else repair_info['id'],
                'repair_type': repair_info.get('repair_type', 'unknown'),
                'overlap_length': 0,
                'similarity': '0.000',
                'mapq': 0,
                'extension': 0,
                'telomere_before': 'NO',
                'telomere_after': 'NO',
                'telomere_status': 'FAILED',
                'telomere_length': 0,
                'telomere_type': 'unknown',
                'selected': False,
                'error': str(e)
            })
            
        except (NoOverlapError, LowQualityError) as e:
            # couldn't even merge
            printv(f"    ✗ Failed: {e}")
            failed_attempts.append({
                'repair_info': repair_info,
                'reason': str(e),
                'rank': i + 1
            })
            
            all_candidates.append({
                'rank': i + 1,
                'repair_id': repair_info['id'][:40] + "..." if len(repair_info['id']) > 40 else repair_info['id'],
                'repair_type': repair_info.get('repair_type', 'unknown'),
                'overlap_length': 0,
                'similarity': '0.000',
                'mapq': 0,
                'extension': 0,
                'telomere_before': 'NO',
                'telomere_after': 'NO',
                'telomere_status': 'FAILED',
                'telomere_length': 0,
                'telomere_type': 'unknown',
                'selected': False,
                'error': str(e)
            })
    
    # sort candidates by rank
    all_candidates.sort(key=lambda x: x['rank'])
    
    # PRIORITY: choose from telomere-restored sequences
    if telomere_restored:
        printv(f"\n  Found {len(telomere_restored)} sequences that restore telomeres")
        
        # This ensures longer overlaps are preferred for more reliable fusion
        telomere_restored.sort(key=lambda x: (
            x['result_info'].get('overlap_length', 0),  # 优先重叠长度
            x['telomere_length'],                         # 其次端粒长度
            x['max_repeats']                               # 最后端粒重复数
        ), reverse=True)


        # show top candidates
        for rank, item in enumerate(telomere_restored[:3]):
            printv(f"    {rank+1}. Length: {item['telomere_length']:,}bp, Repeats: {item['max_repeats']}")
        
        # select best
        best_item = telomere_restored[0]
        best_item['result_info']['telomere_success'] = True
        

        best_repair_id = best_item['repair_info']['id']
        for candidate in all_candidates:
            candidate_id = candidate.get('repair_id_full', candidate['repair_id'])
            if best_repair_id in candidate_id or candidate_id in best_repair_id:
                candidate['selected'] = True
                printv(f"    Selected candidate: {best_repair_id[:60]}...")
                break
        
        return best_item['merged_seq'], best_item['result_info'], best_item['repair_info'], all_candidates
    
    else:
        # No telomeres restored
        printv(f"\n  WARNING: No sequences restore telomeres!")
        printv(f"  Chromosome {chr_id} will not be repaired")
        
        result_info = {
            'chromosome': chr_id,
            'repair_id': 'none',
            'reason': 'no_telomere_restored',
            'attempts': len(repair_sequences),
            'failed_attempts': failed_attempts
        }
        
        return None, result_info, None, all_candidates


def write_tabular_report(report_file: str, results: Dict, all_candidates_by_chr: Dict[str, List[Dict]], 
                         genome_file: str, repair_file: str, params: Dict):
    """
    Write report in tabular format for clear visualization
    FIXED: Properly handles results dictionary structure
    """
    # Validate and fix results structure
    if not isinstance(results, dict):
        print(f"Warning: results is not a dict, it's {type(results)}")
        results = {}
    
    # Ensure results contains the expected structure
    if 'results' not in results:
        results['results'] = {'success': [], 'failed': [], 'skipped': []}
    
    # Ensure success/failed/skipped are lists
    for key in ['success', 'failed', 'skipped']:
        if key not in results['results']:
            results['results'][key] = []
        elif not isinstance(results['results'][key], list):
            print(f"Warning: results['results']['{key}'] is not a list, it's {type(results['results'][key])}")
            results['results'][key] = []
    
    with open(report_file, 'w') as f:
        # Header
        f.write("=" * 100 + "\n")
        f.write("TELOMERE REPAIR SEQUENCE MERGE REPORT (TELOMERE-FIRST VERSION)\n")
        f.write("=" * 100 + "\n\n")
        
        # Parameters
        f.write("PARAMETERS:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Genome file:            {genome_file}\n")
        f.write(f"Repair sequence file:   {repair_file}\n")
        f.write(f"Output directory:       {params['output_dir']}\n")
        f.write(f"Search range:           {params['search_range']:,} bp\n")
        f.write(f"Minimum similarity:     {params['min_similarity']}\n")
        f.write(f"Minimap2 mode:          {params['minimap2_mode']}\n")
        f.write(f"Min match length:       {params['min_match_length']} bp\n")
        f.write(f"Min MAPQ:               {params['min_mapq']}\n")
        f.write(f"Min telomere repeats:   {params['min_telomere_repeats']}\n")
        f.write(f"Organism:               {params['organism']}\n")
        f.write(f"Processing time:        {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Get statistics safely
        stats = results.get('statistics', {})
        tel_stats = stats.get('telomere_stats', {})
        
        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total repair sequences:     {stats.get('total_repair_sequences', 0)}\n")
        f.write(f"Successfully repaired chr:  {stats.get('success_chromosomes', 0)}\n")
        f.write(f"Failed sequences:           {stats.get('failed_sequences', 0)}\n")
        f.write(f"Skipped sequences:          {stats.get('skipped_sequences', 0)}\n\n")
        
        f.write("TELOMERE STATUS:\n")
        f.write(f"  RESTORED:  {tel_stats.get('restored', 0):3d}  (was absent, now present)\n")
        f.write(f"  FAILED:    {stats.get('failed_sequences', 0):3d}  (no telomeres restored)\n\n")
        
        f.write(f"Total extension length: {stats.get('total_extension', 0):,} bp\n")
        f.write(f"Average extension:      {stats.get('avg_extension', 0):,.0f} bp\n")
        f.write(f"Total overlap length:   {stats.get('total_overlap', 0):,} bp\n")
        f.write(f"Average overlap:        {stats.get('avg_overlap', 0):,.0f} bp\n\n")
        
        # Detailed results for each chromosome
        for chr_name, candidates in all_candidates_by_chr.items():
            f.write("=" * 100 + "\n")
            f.write(f"CHROMOSOME: {chr_name}\n")
            f.write("=" * 100 + "\n")
            
            # Find which candidate was selected
            selected_id = None
            for candidate in candidates:
                if candidate.get('selected', False):
                    selected_id = candidate['repair_id']
                    break
            
            # Table header
            f.write("\n")
            f.write("{:<4} {:<45} {:<8} {:>10} {:>8} {:>6} {:>10} {:>12} {:>12} {:>12} {:>10} {:<8}\n".format(
                "Rank", "Repair Sequence ID", "Type", "Overlap", "Sim", "MAPQ", 
                "Extension", "Tel Before", "Tel After", "Tel Status", "Tel Length", "Selected"
            ))
            f.write("-" * 160 + "\n")
            
            # Table rows
            for candidate in sorted(candidates, key=lambda x: x['rank']):
                selected_mark = "→ SELECTED" if candidate.get('selected', False) else ""
                
                tel_before = candidate.get('telomere_before', 'NO')
                tel_after = candidate.get('telomere_after', 'NO')
                tel_status = candidate.get('telomere_status', 'UNKNOWN')
                tel_length = candidate.get('telomere_length', 0)
                tel_type = candidate.get('telomere_type', '')
                
                # Add telomere type to status if present
                if tel_after == 'YES' and tel_type != 'unknown':
                    tel_display = f"{tel_after} ({tel_type})"
                else:
                    tel_display = tel_after
                
                # Truncate long IDs
                repair_id = candidate['repair_id']
                if len(repair_id) > 42:
                    repair_id = repair_id[:39] + "..."
                
                f.write("{:<4} {:<45} {:<8} {:>10} {:>8} {:>6} {:>10} {:>12} {:>12} {:>12} {:>10} {:<8}\n".format(
                    candidate['rank'],
                    repair_id,
                    candidate.get('repair_type', 'unknown'),
                    f"{candidate.get('overlap_length', 0):,}",
                    candidate.get('similarity', '0.000'),
                    candidate.get('mapq', 0),
                    f"{candidate.get('extension', 0):+,}",
                    tel_before,
                    tel_display,
                    tel_status,
                    f"{tel_length:,}",
                    selected_mark
                ))
                
                # Show error if failed
                if 'error' in candidate:
                    f.write(f"      └─ ERROR: {candidate['error']}\n")
            
            f.write("\n")
            
            # Additional details for selected candidate
            if selected_id:
                # Find the selected result in results['success']
                success_list = results.get('results', {}).get('success', [])
                if isinstance(success_list, list):
                    for result in success_list:
                        if isinstance(result, dict) and result.get('repair_id', '') in selected_id:
                            f.write("SELECTED CANDIDATE DETAILS:\n")
                            f.write("-" * 50 + "\n")
                            f.write(f"  Merge method:          {result.get('merge_method', 'unknown')}\n")
                            f.write(f"  Original length:       {result.get('original_length', 0):,} bp\n")
                            f.write(f"  Repaired length:       {result.get('repaired_length', 0):,} bp\n")
                            
                            # Telomere details
                            tel_after = result.get('telomere_after', {})
                            if tel_after.get('has_telomere', False):
                                f.write(f"  Telomere details:\n")
                                f.write(f"    - Total length:      {tel_after.get('total_telomere_length', 0):,} bp\n")
                                f.write(f"    - Max repeats:       {tel_after.get('max_repeats', 0)}\n")
                                f.write(f"    - Telomere count:    {tel_after.get('telomere_count', 0)}\n")
                                f.write(f"    - Telomere type:     {tel_after.get('telomere_type', 'unknown')}\n")
                            
                            break
                
                f.write("\n")
        
        # Failed sequences summary
        failed_list = results.get('results', {}).get('failed', [])
        if failed_list and isinstance(failed_list, list):
            f.write("=" * 100 + "\n")
            f.write("FAILED SEQUENCES\n")
            f.write("=" * 100 + "\n\n")
            
            f.write("FAILED SEQUENCES:\n")
            f.write("-" * 50 + "\n")
            for fail in failed_list:
                if isinstance(fail, dict) and fail.get('repair_id') != 'none':
                    f.write(f"  {fail.get('repair_id', 'unknown')}: {fail.get('reason', 'unknown')}\n")
            f.write("\n")
        
        # Footer
        f.write("=" * 100 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 100 + "\n")


# ============================================================================
# Main function - Modified with organism parameter
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
    min_telomere_repeats: int = 20,
    organism: str = 'plant'  # default to plant
) -> Dict[str, Any]:
    
    global VERBOSE
    VERBOSE = verbose
    
    start_time = time.time()
    
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
            print(f"Minimum telomere repeats: {min_telomere_repeats}")
            print(f"Organism: {organism} (telomere pattern: {'plant (TTTAGGG)' if organism == 'plant' else 'vertebrate (TTAGGG)' if organism == 'vertebrate' else 'auto-detect'})")
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
            'organism': organism
        }
        
        results = {'success': [], 'failed': [], 'skipped': []}
        repaired_chromosomes = {}
        used_repair_ids = set()
        
        # Store all candidates for tabular report
        all_candidates_by_chr = {}
        
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
            
            merged_seq, result_info, used_repair_info, candidates = try_repair_chromosome(
                chr_id, chr_seq, repair_sequences, minimap2_config
            )
            
            # Store candidates for tabular report
            all_candidates_by_chr[chr_id] = candidates
            
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
                    description=f"repaired|organism:{organism}|original_{len(seq)}bp|repaired_{len(repaired_seq)}bp"
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
        
        # generate enhanced tabular report
        try:
            params = {
                'output_dir': output_dir,
                'search_range': search_range,
                'min_similarity': min_similarity,
                'minimap2_mode': minimap2_mode,
                'min_match_length': min_match_length,
                'min_mapq': min_mapq,
                'min_telomere_repeats': min_telomere_repeats,
                'organism': organism
            }
            
            # Calculate statistics
            success_list = results.get('success', [])
            if not isinstance(success_list, list):
                success_list = []
            
            failed_list = results.get('failed', [])
            if not isinstance(failed_list, list):
                failed_list = []
            
            skipped_list = results.get('skipped', [])
            if not isinstance(skipped_list, list):
                skipped_list = []
            
            success_chromosomes = len([r for r in success_list if isinstance(r, dict) and 'attempt_number' in r])
            
            total_extension = 0
            total_overlap = 0
            restored_count = 0
            
            if success_list and isinstance(success_list, list):
                for r in success_list:
                    if not isinstance(r, dict):
                        continue
                    if 'extension' in r:
                        total_extension += r.get('extension', 0)
                    if 'overlap_length' in r:
                        total_overlap += r.get('overlap_length', 0)
                    
                    tel_status = r.get('telomere_comparison', {}).get('status', '')
                    if tel_status == 'RESTORED':
                        restored_count += 1
            
            avg_extension = total_extension / success_chromosomes if success_chromosomes > 0 else 0
            avg_overlap = total_overlap / success_chromosomes if success_chromosomes > 0 else 0
            
            result_dict = {
                'success': True,
                'repaired_genome': repaired_file,
                'report_file': report_file,
                'results': results,  # contains 'success', 'failed', 'skipped'
                'statistics': {
                    'total_repair_sequences': len(repair_records),
                    'success_chromosomes': success_chromosomes,
                    'failed_sequences': len([f for f in failed_list if isinstance(f, dict) and f.get('reason') == 'no_telomere_restored']),
                    'skipped_sequences': len(skipped_list),
                    'total_extension': total_extension,
                    'avg_extension': avg_extension,
                    'total_overlap': total_overlap,
                    'avg_overlap': avg_overlap,
                    'organism': organism,
                    'telomere_stats': {
                        'restored': restored_count
                    },
                    'processing_time': time.time() - start_time
                },
                'output_dir': output_dir
            }
            
            # Write tabular report
            write_tabular_report(report_file, result_dict, all_candidates_by_chr, 
                                genome_file, repair_file, params)
            
            if verbose:
                print(f"Tabular report saved to: {report_file}")
        except Exception as e:
            print(f"Warning: Cannot save tabular report: {e}")
            import traceback
            traceback.print_exc()
        
        total_time = time.time() - start_time
        
        if verbose:
            print("\n" + "="*70)
            print("Processing completed!")
            print("="*70)
            print(f"Successfully repaired: {success_chromosomes} chromosomes")
            print(f"  Telomere status: {restored_count} RESTORED")
            print(f"Repair failed: {len(failed_list)} sequences")
            print(f"Skipped: {len(skipped_list)} sequences")
            print(f"Total time: {total_time:.2f} seconds")
            
            if success_chromosomes > 0:
                print(f"\nTotal extension length: {total_extension:,} bp")
                print(f"Average extension length: {avg_extension:,.0f} bp")
                print(f"Total overlap length: {total_overlap:,} bp")
                print(f"Average overlap length: {avg_overlap:,.0f} bp")
                print(f"Repaired genome: {repaired_file}")
            
            print(f"Tabular report: {report_file}")
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
        description="Telomere Repair Sequence Merge Script (Telomere-First Version) - Plant Telomere Support with Tabular Report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Core Logic: PRIORITIZE SEQUENCES THAT RESTORE TELOMERES
        
        Features:
          1. Try ALL candidates, not just first success
          2. Validate telomeres BEFORE and AFTER repair
          3. ONLY accept sequences that restore telomeres
          4. Support for both plant (TTTAGGG) and vertebrate (TTAGGG) telomeres
          5. Overlap length considered for reliable fusion
          6. Tabular report format for clear visualization
        
        Usage examples:
          python add_te.py -q genome.fasta -i repair_sequences.fa -o output_dir --organism plant
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
    
    parser.add_argument("--min-telomere-repeats", type=int, default=20,
                       help="Minimum telomere repeats for detection (default: 20)")
    
    # organism parameter - default to plant
    parser.add_argument("--organism", default='plant',
                       choices=['plant', 'vertebrate', 'auto'],
                       help="Organism type for telomere pattern detection: plant (TTTAGGG), vertebrate (TTAGGG), or auto-detect (default: plant)")
    
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
        min_telomere_repeats=args.min_telomere_repeats,
        organism=args.organism
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
        print(f"Organism: {stats.get('organism', 'plant')}")
        print(f"Successfully repaired: {stats['success_chromosomes']} chromosomes")
        print(f"\nTelomere Status:")
        print(f"  RESTORED:  {tel_stats.get('restored', 0):3d}  (was absent, now present)")
        
        if tel_stats.get('restored', 0) > 0:
            print(f"\n✓ SUCCESS: Telomeres were restored!")
        else:
            print(f"\n⚠ WARNING: No telomeres were restored")
            print(f"   Check if repair sequences actually contain telomere repeats")
        
        print(f"\nTotal extension length: {stats.get('total_extension', 0):,} bp")
        print(f"Average extension: {stats.get('avg_extension', 0):,.0f} bp")
        print(f"Total overlap length: {stats.get('total_overlap', 0):,} bp")
        print(f"Average overlap: {stats.get('avg_overlap', 0):,.0f} bp")
        print(f"\nDetailed tabular report: {result['report_file']}")
        print(f"{'='*50}")


if __name__ == "__main__":
    main()
