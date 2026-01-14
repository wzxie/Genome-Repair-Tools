#!/usr/bin/env python3
"""
Telomere Repair Sequence Merge Script (Optimized Validation Version)
Purpose: Merge extracted telomere repair sequences to corresponding chromosome ends
Features: 1. Intelligent selection (try from longest to shortest) 2. minimap2 unilateral alignment 3. Lenient validation logic 4. Support minimal change strategy
"""

import sys
import os
import re
import argparse
import textwrap
import subprocess
import tempfile
import json
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
                printv(f"✓ minimap2 installed: {version}")
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
            
            cigar = ""
            for i in range(12, len(fields)):
                if fields[i].startswith('cg:Z:'):
                    cigar = fields[i][5:]
                    break
            
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
                    'cigar': cigar,
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

class TelomereRepair:
    
    def __init__(self, config: Dict = None):
        self.config = config or {
            'minimap2_mode': 'asm5',
            'min_score': 0.7,
            'min_match_length': 100,
            'min_mapq': 20,
            'search_range': 5000000
        }
        
        self.matcher = Minimap2Matcher(
            mode=self.config['minimap2_mode'],
            min_score=self.config['min_score'],
            min_match_length=self.config['min_match_length'],
            min_mapq=self.config['min_mapq']
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
                printv(f"      Coverage: {match['coverage']:.3f}")
                printv(f"      MAPQ: {match['mapq']}")
                printv(f"      Strand: {match['strand']}")
                printv(f"      Combined score: {match['score']:.3f}")
                
                distance_to_start = match['target_start']
                printv(f"      Distance to chromosome start: {distance_to_start:,} bp")
            
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
                printv(f"      Coverage: {match['coverage']:.3f}")
                printv(f"      MAPQ: {match['mapq']}")
                printv(f"      Strand: {match['strand']}")
                printv(f"      Combined score: {match['score']:.3f}")
                
                distance_to_end = chr_len - match['target_start']
                printv(f"      Distance to chromosome end: {distance_to_end:,} bp")
            
            return match
        
        else:
            printv(f"    ⚠️  Unknown repair type: {repair_type}")
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
        
        search_range = self.config['search_range']
        
        if repair_type == '5prime':
            if match['target_start'] > search_range:
                printv(f"    ⚠️  Match position far from chromosome start: {match['target_start']:,} (search range: {search_range:,})")
            
            if match['query_end'] < repair_len * 0.3:
                printv(f"    ⚠️  Match may be at start of repair sequence: {match['query_end']:,}/{repair_len:,}")
            else:
                printv(f"    Match position reasonable in repair sequence: end {match['query_end']:,}/{repair_len:,}")
        
        elif repair_type == '3prime':
            expected_min_start = chr_len - search_range
            if match['target_start'] < expected_min_start:
                printv(f"    ⚠️  Match position may exceed search range: {match['target_start']:,} < {expected_min_start:,}")
            else:
                distance_to_end = chr_len - match['target_start']
                printv(f"    Match distance to chromosome end: {distance_to_end:,} bp (within {search_range:,}bp search range)")
            
            if match['query_start'] > repair_len * 0.5:
                printv(f"    ⚠️  Match may be at end of repair sequence: {match['query_start']:,}/{repair_len:,}")
            else:
                printv(f"    Match position reasonable in repair sequence: start {match['query_start']:,}/{repair_len:,}")
        
        if match['strand'] == '-':
            printv(f"    ⚠️  Match on reverse strand, will perform reverse complement")
        
        quality_score = self._calculate_quality_score(match)
        printv(f"    Combined quality score: {quality_score:.3f}")
        
        if quality_score < 0.5:
            return False, f"Combined quality score too low: {quality_score:.3f}"
        elif quality_score < 0.7:
            printv(f"    ⚠️  Medium quality score, but accepting this match")
        
        return True, "Validation passed"
    
    def _calculate_quality_score(self, match: Dict) -> float:
        identity_weight = 0.4
        length_weight = 0.3
        mapq_weight = 0.2
        coverage_weight = 0.1
        
        identity_norm = match['identity']
        length_norm = min(1.0, match['match_length'] / 10000.0)
        mapq_norm = min(1.0, match['mapq'] / 60.0)
        coverage_norm = match['coverage']
        
        score = (
            identity_norm * identity_weight +
            length_norm * length_weight +
            mapq_norm * mapq_weight +
            coverage_norm * coverage_weight
        )
        
        return score
    
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
            
            printv(f"    5' end merge:")
            printv(f"      Repair unique part: {len(repair_unique):,} bp")
            printv(f"      Overlap region: {len(overlap_from_repair):,} bp")
            printv(f"      Chromosome remainder: {len(chr_remainder):,} bp")
        
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
            
            printv(f"    3' end merge:")
            printv(f"      Chromosome unique part: {len(chr_unique):,} bp")
            printv(f"      Overlap region: {len(overlap_from_chr):,} bp")
            printv(f"      Repair unique part: {len(repair_unique):,} bp")
        
        else:
            raise ValueError(f"Unknown repair type: {repair_type}")
        
        merged_len = len(merged_seq)
        expected_len = stats.get('repair_unique_len', 0) + stats['overlap_len'] + stats.get('chr_remainder_len', stats.get('chr_unique_len', 0))
        
        if merged_len != expected_len:
            printv(f"    ⚠️  Merge length mismatch: actual {merged_len:,} vs expected {expected_len:,}")
        
        stats['original_chr_len'] = chr_len
        stats['original_repair_len'] = repair_len
        stats['merged_len'] = merged_len
        stats['extension'] = merged_len - chr_len
        stats['identity'] = match['identity']
        stats['coverage'] = match['coverage']
        stats['mapq'] = match['mapq']
        stats['strand'] = match['strand']
        stats['score'] = match['score']
        stats['quality_score'] = self._calculate_quality_score(match)
        
        printv(f"    Merge statistics:")
        printv(f"      Original chromosome: {chr_len:,} bp")
        printv(f"      After merge: {merged_len:,} bp")
        printv(f"      Extension: {stats['extension']:,} bp")
        printv(f"      Similarity: {match['identity']:.3f}")
        printv(f"      Quality score: {stats['quality_score']:.3f}")
        
        return merged_seq, stats

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
        for i, seq in enumerate(chromosome_groups[chr_name][:3]):
            printv(f"    Length {i+1}: {seq['length']:,} bp ({seq.get('repair_type', 'unknown')} end)")
        if len(chromosome_groups[chr_name]) > 3:
            printv(f"    ... {len(chromosome_groups[chr_name])-3} more shorter sequences")
    
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

def calculate_selection_score(result_info: Dict, strategy: str = 'minimal_extension') -> float:
    
    if strategy == 'first_success':
        return 1.0
    
    elif strategy == 'minimal_extension':
        extension_abs = abs(result_info.get('extension', 0))
        extension_norm = max(0, 1.0 - min(1.0, extension_abs / 100000.0))
        quality = result_info.get('minimap2_info', {}).get('quality_score', 0.5)
        similarity = result_info.get('similarity', 0.7)
        
        score = (
            extension_norm * 0.6 +
            quality * 0.3 +
            similarity * 0.1
        )
        
        return score
    
    elif strategy == 'balanced':
        extension_abs = abs(result_info.get('extension', 0))
        extension_norm = max(0, 1.0 - min(1.0, extension_abs / 100000.0))
        quality = result_info.get('minimap2_info', {}).get('quality_score', 0.5)
        similarity = result_info.get('similarity', 0.7)
        mapq = min(1.0, result_info.get('minimap2_info', {}).get('mapq', 20) / 60.0)
        
        score = (
            extension_norm * 0.4 +
            quality * 0.3 +
            similarity * 0.2 +
            mapq * 0.1
        )
        
        return score
    
    else:
        raise ValueError(f"Unknown selection strategy: {strategy}")

def try_repair_chromosome(chr_id: str, chr_seq: str, 
                         repair_sequences: List[Dict],
                         config: Dict,
                         selection_strategy: str = 'minimal_extension') -> Tuple[Optional[str], Dict, Dict]:
    
    printv(f"\nProcessing chromosome: {chr_id} ({len(chr_seq):,} bp)")
    printv(f"  {len(repair_sequences)} candidate repair sequences")
    printv(f"  Selection strategy: {selection_strategy}")
    
    repairer = TelomereRepair(config)
    
    all_successful_results = []
    
    for i, repair_info in enumerate(repair_sequences):
        repair_seq = repair_info['sequence']
        repair_type = repair_info.get('repair_type', 'unknown')
        repair_id = repair_info['id']
        
        printv(f"\n  Trying sequence {i+1}/{len(repair_sequences)}: {repair_id}")
        printv(f"    Length: {len(repair_seq):,} bp")
        printv(f"    Type: {repair_type}")
        
        if repair_type not in ['5prime', '3prime']:
            printv(f"    ⚠️  Unknown repair type, skipping")
            continue
        
        try:
            match = repairer.find_telomere_overlap(chr_seq, repair_seq, repair_type)
            
            if not match:
                printv(f"    ✗ No overlap found")
                continue
            
            is_valid, valid_msg = repairer.validate_match(
                match, repair_type, len(chr_seq), len(repair_seq)
            )
            
            if not is_valid:
                printv(f"    ✗ Validation failed: {valid_msg}")
                continue
            
            merged_seq, stats = repairer.merge_sequences(
                chr_seq, repair_seq, match, repair_type
            )
            
            extension = stats['extension']
            
            result_info = {
                'chromosome': chr_id,
                'repair_id': repair_id,
                'repair_type': repair_type,
                'overlap_length': stats['overlap_len'],
                'similarity': stats['identity'],
                'original_length': len(chr_seq),
                'repaired_length': len(merged_seq),
                'extension': extension,
                'extension_abs': abs(extension),
                'attempt_number': i + 1,
                'total_attempts': len(repair_sequences),
                'sequence_length': len(repair_seq),
                'minimap2_info': {
                    'identity': stats['identity'],
                    'coverage': stats['coverage'],
                    'mapq': stats['mapq'],
                    'strand': stats['strand'],
                    'score': stats['score'],
                    'quality_score': stats['quality_score']
                },
                'merge_method': stats['method'],
                'cut_positions': stats['cut_positions'],
                'merged_sequence': merged_seq,
                'stats': stats
            }
            
            printv(f"    ✓ Sequence {i+1} repair successful!")
            printv(f"      Overlap length: {stats['overlap_len']:,} bp")
            printv(f"      Extension length: {extension:,} bp")
            printv(f"      Similarity: {stats['identity']:.3f}")
            printv(f"      Quality score: {stats['quality_score']:.3f}")
            
            if selection_strategy == 'first_success':
                printv(f"    → Using first successful result (first_success strategy)")
                del result_info['merged_sequence']
                return merged_seq, result_info, repair_info
            
            all_successful_results.append({
                'score': 0.0,
                'result_info': result_info,
                'repair_info': repair_info,
                'merged_seq': merged_seq
            })
            
        except Exception as e:
            printv(f"    ✗ Processing error: {e}")
            import traceback
            traceback.print_exc()
    
    if not all_successful_results:
        printv(f"  ✗ All {len(repair_sequences)} sequences failed")
        return None, {
            'chromosome': chr_id,
            'repair_id': 'none',
            'reason': 'all_attempts_failed',
            'attempts': len(repair_sequences)
        }, None
    
    printv(f"\n  ✓ Found {len(all_successful_results)} successful merge results")
    printv(f"    Selecting optimal result based on '{selection_strategy}' strategy...")
    
    for item in all_successful_results:
        score = calculate_selection_score(item['result_info'], selection_strategy)
        item['score'] = score
    
    all_successful_results.sort(key=lambda x: x['score'], reverse=True)
    
    if VERBOSE:
        printv(f"\n  All candidate results:")
        for rank, item in enumerate(all_successful_results[:5]):
            result_info = item['result_info']
            printv(f"    {rank+1}. Score: {item['score']:.3f}, "
                  f"Repair sequence: {result_info['repair_id'][:50]}..., "
                  f"Extension: {result_info['extension']:,} bp, "
                  f"Similarity: {result_info['similarity']:.3f}")
    
    best_item = all_successful_results[0]
    best_score = best_item['score']
    best_result_info = best_item['result_info']
    best_repair_info = best_item['repair_info']
    best_merged_seq = best_item['merged_seq']
    
    if len(all_successful_results) > 1:
        second_score = all_successful_results[1]['score']
        if abs(best_score - second_score) < 0.1:
            printv(f"    ⚠️  Warning: Best result and second best result have close scores ({best_score:.3f} vs {second_score:.3f})")
            printv(f"    Manual inspection may be needed")
    
    printv(f"\n  ✓ Selected optimal result:")
    printv(f"    Using sequence: {best_result_info['repair_id']}")
    printv(f"    Selection strategy: {selection_strategy}")
    printv(f"    Combined score: {best_score:.3f}")
    printv(f"    Extension length: {best_result_info['extension']:,} bp")
    printv(f"    Similarity: {best_result_info['similarity']:.3f}")
    printv(f"    Quality score: {best_result_info['minimap2_info']['quality_score']:.3f}")
    
    best_result_info['selection_method'] = selection_strategy
    best_result_info['selection_score'] = best_score
    best_result_info['alternatives_count'] = len(all_successful_results)
    
    best_result_info['all_alternatives_scores'] = []
    for item in all_successful_results:
        alt_info = item['result_info']
        best_result_info['all_alternatives_scores'].append({
            'repair_id': alt_info['repair_id'],
            'score': item['score'],
            'extension': alt_info['extension'],
            'similarity': alt_info['similarity'],
            'quality_score': alt_info['minimap2_info']['quality_score']
        })
    
    if 'merged_sequence' in best_result_info:
        del best_result_info['merged_sequence']
    
    return best_merged_seq, best_result_info, best_repair_info

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
    selection_strategy: str = 'minimal_extension'
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
            print("Telomere Repair Sequence Merge Script (Optimized Validation Version)")
            print("="*70)
            print("Strategy: Intelligent selection + minimap2 unilateral alignment + lenient validation")
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
            print("="*70)
        
        if verbose:
            print(f"Reading genome file: {genome_file}")
        chromosomes = read_fasta(genome_file)
        if verbose:
            print(f"Found {len(chromosomes)} chromosomes")
        
        if verbose:
            print("\nChromosome list:")
            for i, (chr_id, seq) in enumerate(list(chromosomes.items())[:5]):
                print(f"  {i+1:2d}. {chr_id:30s} ({len(seq):,} bp)")
            if len(chromosomes) > 5:
                print(f"  ... {len(chromosomes)-5} more")
        
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
            'search_range': search_range
        }
        
        results = {'success': [], 'failed': [], 'skipped': []}
        repaired_chromosomes = {}
        used_repair_ids = set()
        
        for target_chr, repair_sequences in chromosome_groups.items():
            if not repair_sequences:
                continue
            
            if verbose:
                print(f"\n{'='*70}")
                print(f"Processing chromosome group: {target_chr} ({len(repair_sequences)} candidate sequences)")
                print('='*70)
            
            chr_id = find_matching_chromosome(chromosomes, target_chr)
            if not chr_id:
                if verbose:
                    print(f"✗ Cannot find matching chromosome in genome: '{target_chr}'")
                for repair_info in repair_sequences:
                    results['skipped'].append({
                        'repair_id': repair_info['id'],
                        'reason': 'chromosome_not_found',
                        'target': target_chr
                    })
                continue
            
            if chr_id in repaired_chromosomes:
                if verbose:
                    print(f"⚠️  Chromosome {chr_id} already repaired, skipping this group")
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
                
                for repair_info in repair_sequences:
                    results['failed'].append({
                        'repair_id': repair_info['id'],
                        'chromosome': chr_id,
                        'reason': 'repair_attempt_failed',
                        'target': target_chr,
                        'sequence_length': repair_info['length']
                    })
        
        all_repair_ids = {record.id for record in repair_records}
        unused_repair_ids = all_repair_ids - used_repair_ids
        
        for repair_id in unused_repair_ids:
            already_processed = False
            for result_list in results.values():
                for result in result_list:
                    if result.get('repair_id') == repair_id:
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
        
        try:
            with open(report_file, 'w') as f:
                f.write("Telomere Repair Sequence Merge Report (Optimized Validation Version)\n")
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
                f.write(f"  Selection strategy: {selection_strategy}\n\n")
                
                total_success = len([r for r in results['success'] if 'attempt_number' in r])
                total_failed = len([r for r in results['failed'] if r.get('reason') == 'repair_attempt_failed'])
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
                            f.write(f"  Selection score: {result.get('selection_score', 'N/A')}\n")
                            f.write(f"  Merge method: {result.get('merge_method', 'unknown')}\n")
                            f.write(f"  Original length: {result['original_length']:,} bp\n")
                            f.write(f"  Repaired length: {result['repaired_length']:,} bp\n")
                            f.write(f"  Extension length: {result['extension']:,} bp\n")
                            f.write(f"  Overlap length: {result['overlap_length']:,} bp\n")
                            f.write(f"  Similarity: {result['similarity']:.4f}\n")
                            f.write(f"  Quality score: {result['minimap2_info']['quality_score']:.3f}\n")
                            f.write(f"  MAPQ: {result['minimap2_info']['mapq']}\n")
                            f.write(f"  Strand: {result['minimap2_info']['strand']}\n")
                            f.write(f"  Cut positions: repair sequence{result['cut_positions'].get('repair_cut', 'N/A')}, "
                                  f"chromosome{result['cut_positions'].get('chr_cut', 'N/A')}\n")
                            f.write(f"  Attempt number: {result['attempt_number']}/{result['total_attempts']}\n")
                            f.write(f"  Alternative count: {result.get('alternatives_count', 'N/A')}\n")
                            f.write(f"  Sequence length: {result.get('sequence_length', 'N/A'):,} bp\n")
                            
                            if selection_strategy != 'first_success' and 'all_alternatives_scores' in result:
                                f.write(f"  Alternative scores:\n")
                                for alt in result['all_alternatives_scores'][:3]:
                                    f.write(f"    - {alt['repair_id'][:50]}...: score {alt['score']:.3f}, "
                                           f"extension {alt['extension']:,}bp, similarity {alt['similarity']:.3f}\n")
                                if len(result['all_alternatives_scores']) > 3:
                                    f.write(f"    ... {len(result['all_alternatives_scores'])-3} more alternatives\n")
                            
                            f.write("\n")
                
                if results['failed']:
                    f.write("\nFailed repair sequences:\n")
                    f.write("-"*50 + "\n")
                    for fail in results['failed']:
                        if fail.get('repair_id') != 'none':
                            f.write(f"Repair sequence: {fail.get('repair_id', 'unknown')}\n")
                            f.write(f"  Target chromosome: {fail.get('chromosome', 'unknown')}\n")
                            f.write(f"  Failure reason: {fail.get('reason', 'unknown')}\n")
                            if 'sequence_length' in fail:
                                f.write(f"  Sequence length: {fail['sequence_length']:,} bp\n")
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
                        if 'selection_strategy' in skip:
                            f.write(f"  Selection strategy: {skip['selection_strategy']}\n")
                        if 'sequence_length' in skip:
                            f.write(f"  Sequence length: {skip['sequence_length']:,} bp\n")
                        f.write("\n")
            
            if verbose:
                print(f"Detailed report saved to: {report_file}")
        except Exception as e:
            print(f"Warning: Cannot save text report: {e}")
        
        total_time = time.time() - start_time
        
        success_chromosomes = len([r for r in results['success'] if 'attempt_number' in r])
        
        total_extension = 0
        if results['success']:
            total_extension = sum(r.get('extension', 0) for r in results['success'] if 'extension' in r)
        
        avg_extension = total_extension / success_chromosomes if success_chromosomes > 0 else 0
        
        if selection_strategy in ['minimal_extension', 'balanced'] and results['success']:
            min_extensions = []
            for result in results['success']:
                if 'all_alternatives_scores' in result:
                    alternatives = result['all_alternatives_scores']
                    min_ext = min(alt['extension'] for alt in alternatives)
                    actual_ext = result['extension']
                    min_extensions.append((actual_ext, min_ext))
            
            if min_extensions:
                avg_actual_ext = sum(ext[0] for ext in min_extensions) / len(min_extensions)
                avg_min_ext = sum(ext[1] for ext in min_extensions) / len(min_extensions)
                selection_efficiency = (avg_min_ext / avg_actual_ext) if avg_actual_ext != 0 else 1.0
            else:
                selection_efficiency = 1.0
        else:
            selection_efficiency = 1.0
        
        result_dict = {
            'success': True,
            'repaired_genome': repaired_file,
            'report_file': report_file,
            'results': results,
            'statistics': {
                'total_repair_sequences': len(repair_records),
                'success_chromosomes': success_chromosomes,
                'failed_sequences': len([r for r in results['failed'] if r.get('reason') == 'repair_attempt_failed']),
                'skipped_sequences': len(results['skipped']),
                'total_extension': total_extension,
                'avg_extension': avg_extension,
                'selection_strategy': selection_strategy,
                'selection_efficiency': selection_efficiency if selection_strategy in ['minimal_extension', 'balanced'] else 'N/A',
                'processing_time': total_time
            },
            'output_dir': output_dir
        }
        
        if verbose:
            print("\n" + "="*70)
            print("Processing completed!")
            print("="*70)
            print(f"Successfully repaired: {success_chromosomes} chromosomes")
            print(f"Repair failed: {len([r for r in results['failed'] if r.get('reason') == 'repair_attempt_failed'])} sequences")
            print(f"Skipped: {len(results['skipped'])} sequences")
            print(f"Total time: {total_time:.2f} seconds")
            print(f"Selection strategy: {selection_strategy}")
            
            if success_chromosomes > 0:
                print(f"\nTotal extension length: {total_extension:,} bp")
                print(f"Average extension length: {avg_extension:,.0f} bp")
                
                if selection_strategy in ['minimal_extension', 'balanced'] and selection_efficiency != 1.0:
                    print(f"Selection efficiency: {selection_efficiency:.2%} (closer to 100% means smaller changes selected)")
                
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
        description="Telomere Repair Sequence Merge Script (Optimized Validation Version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Optimization strategies:
          1. Intelligent selection (try from longest to shortest)
          2. minimap2 unilateral alignment (repair sequence vs chromosome end region)
          3. Lenient validation logic, focus on match quality rather than strict position limits
          4. Support reverse strand matching and reverse complement processing
          5. Support multiple selection strategies, including minimal change strategy
        
        Selection strategy description:
          - first_success: Stop at first successful result (original strategy)
          - minimal_extension: Minimal change strategy (evaluate all possibilities, select smallest extension)
          - balanced: Balanced strategy (consider both extension and match quality)
        
        Usage examples:
          python telomere_repair_optimized.py -q genome.fasta -i repair_sequences.fasta -o output_dir
          
          Optional parameters:
          --search-range: Search range(bp), default 5,000,000
          --min-similarity: Minimum alignment similarity, default 0.7
          --minimap2-mode: minimap2 alignment mode, default asm5
          --min-match-length: Minimum match length, default 100
          --min-mapq: Minimum mapping quality, default 20
          --strategy: Selection strategy, default minimal_extension
        """)
    )
    
    parser.add_argument("-q", "--query", required=True, help="Genome FASTA file")
    parser.add_argument("-i", "--input", required=True, help="Repair sequences FASTA file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    
    parser.add_argument("--search-range", type=int, default=5000000,
                       help="Alignment search range(bp) (default: 5,000,000)")
    parser.add_argument("--min-similarity", type=float, default=0.7,
                       help="Minimum alignment similarity (default: 0.7)")
    
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
        selection_strategy=args.strategy
    )
    
    if not result['success']:
        print(f"\nError: {result['error']}")
        if args.verbose and 'traceback' in result:
            print(f"Detailed error information:\n{result['traceback']}")
        sys.exit(1)
    else:
        stats = result['statistics']
        print(f"\nProcessing completed!")
        print(f"Successfully repaired {stats['success_chromosomes']} chromosomes")
        print(f"Total extension length: {stats.get('total_extension', 0):,} bp")
        print(f"Average extension length: {stats.get('avg_extension', 0):,.0f} bp")
        print(f"Selection strategy: {stats['selection_strategy']}")
        
        if stats['selection_strategy'] in ['minimal_extension', 'balanced']:
            print(f"Selection efficiency: {stats.get('selection_efficiency', 'N/A'):.2%}")
        
        print(f"Detailed report: {result['report_file']}")

if __name__ == "__main__":
    main()