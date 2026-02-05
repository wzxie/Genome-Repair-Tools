#!/usr/bin/env python3
"""
Telomere Sequence Extraction Script (Enhanced Version) - With length limit
Purpose: Analyze NUCMER coords files, extract correct collinear regions and their surrounding sequences
Enhanced features:
1. Analyze complete alignment direction distribution
2. Intelligent extraction strategy determination
3. Verify telomere position correctness
4. Add length limit to prevent extracting overly long sequences
"""

import sys
import os
import re
import argparse
import textwrap
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Dict, List, Tuple, Optional, Any, Set
import json

VERBOSE = False

def printv(message: str):
    if VERBOSE:
        print(message)

def find_telomere_sequences(sequence, min_repeats=20, min_telomere_length=500):
    telomere_patterns = [
        (re.compile(r'(TTTAGGG){%d,}' % min_repeats), "Plant telomere forward"),
        (re.compile(r'(CCCTAAA){%d,}' % min_repeats), "Plant telomere reverse"),
        (re.compile(r'(TTAGGG){%d,}' % min_repeats), "Universal telomere forward"),
        (re.compile(r'(CCCTAA){%d,}' % min_repeats), "Universal telomere reverse"),
    ]
    
    matches = []
    sequence_upper = sequence.upper()
    
    for pattern, pattern_name in telomere_patterns:
        for match in pattern.finditer(sequence_upper):
            if len(match.group()) >= min_telomere_length:
                matches.append({
                    'pattern': pattern_name,
                    'sequence': match.group(),
                    'start': match.start(),
                    'end': match.end(),
                    'length': len(match.group())
                })
    
    return matches

def parse_data_line_enhanced(line: str) -> Optional[Tuple]:
    pattern = r'(\d+\.?\d*|\w+[\w\-\.]*)'
    matches = re.findall(pattern, line)
    
    if len(matches) < 11:
        return None
    
    try:
        s1 = int(matches[0])
        e1 = int(matches[1])
        s2 = int(matches[2])
        e2 = int(matches[3])
        len1 = int(matches[4])
        len2 = int(matches[5])
        identity = float(matches[6])
        len_r = int(matches[7])
        len_q = int(matches[8])
        
        tags = matches[9:]
        
        if len(tags) >= 2:
            ref_name = tags[-2]
            query_name = tags[-1]
        else:
            ref_name = tags[0] if tags else "unknown"
            query_name = ""
        
        is_reverse = (s2 > e2)
        
        if s1 > e1:
            s1, e1 = e1, s1
        
        return s1, e1, s2, e2, len1, len2, identity, len_r, len_q, ref_name, query_name, is_reverse
        
    except (ValueError, IndexError) as e:
        printv(f"Failed to parse data line: {e} - {line}")
        return None

def parse_coords_file_enhanced(coords_file: str) -> Dict:
    alignments_by_pair = {}
    
    try:
        with open(coords_file, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error: Cannot read coords file {coords_file}: {e}")
        return {}
    
    data_started = False
    header_found = False
    
    for line_num, line in enumerate(lines):
        line = line.strip()
        
        if not line:
            continue
        
        if line.startswith('=' * 80) or line.startswith('=' * 60):
            header_found = True
            continue
        
        if '[S1]' in line and '[E1]' in line and '[S2]' in line:
            header_found = True
            continue
        
        if header_found and not data_started:
            data_started = True
        
        if data_started and not line.startswith('NUCMER'):
            parsed = parse_data_line_enhanced(line)
            if parsed:
                s1, e1, s2, e2, len1, len2, identity, len_r, len_q, ref_name, query_name, is_reverse = parsed
                
                key = (ref_name, query_name)
                if key not in alignments_by_pair:
                    alignments_by_pair[key] = {
                        'alignments': [],
                        'ref_len': len_r,
                        'query_len': len_q
                    }
                
                alignments_by_pair[key]['alignments'].append({
                    's1': s1,
                    'e1': e1,
                    's2': s2,
                    'e2': e2,
                    'len1': len1,
                    'len2': len2,
                    'identity': identity,
                    'ref_name': ref_name,
                    'query_name': query_name,
                    'ref_len': len_r,
                    'query_len': len_q,
                    'ref_start': min(s1, e1),
                    'ref_end': max(s1, e1),
                    'is_reverse': is_reverse
                })
    
    return alignments_by_pair

def analyze_contig_structure(contig_file: str, min_repeats: int = 20, min_telomere_length: int = 500) -> Dict[str, Any]:
    telomere_positions = {}
    
    for record in SeqIO.parse(contig_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        
        telomere_matches = find_telomere_sequences(sequence, min_repeats, min_telomere_length)
        
        telomere_regions = []
        for match in telomere_matches:
            telomere_regions.append({
                'start': match['start'] + 1,
                'end': match['end'],
                'length': match['length'],
                'pattern': match['pattern']
            })
        
        telomere_regions.sort(key=lambda x: x['start'])
        
        telomere_positions[seq_id] = {
            'sequence_length': len(record.seq),
            'telomeres': telomere_regions,
            'has_telomeres': len(telomere_regions) > 0,
            'sequence': sequence
        }
    
    return telomere_positions

def check_extracted_sequence_for_telomere(
    extracted_seq: str, 
    min_repeats: int = 20, 
    min_telomere_length: int = 500,
    check_end: str = 'both'
) -> Dict[str, Any]:
    matches = find_telomere_sequences(extracted_seq, min_repeats, min_telomere_length)
    
    has_5prime_telomere = False
    has_3prime_telomere = False
    
    for match in matches:
        if match['start'] < len(extracted_seq) * 0.1:
            has_5prime_telomere = True
        if match['end'] > len(extracted_seq) * 0.9:
            has_3prime_telomere = True
    
    if check_end == '5prime':
        telomeres_valid = has_5prime_telomere
    elif check_end == '3prime':
        telomeres_valid = has_3prime_telomere
    else:
        telomeres_valid = has_5prime_telomere or has_3prime_telomere
    
    return {
        'has_telomere': len(matches) > 0,
        'telomeres': matches,
        'has_5prime_telomere': has_5prime_telomere,
        'has_3prime_telomere': has_3prime_telomere,
        'telomeres_valid': telomeres_valid
    }

def analyze_alignment_direction(alignments: List[Dict]) -> Dict[str, Any]:
    forward_count = sum(1 for align in alignments if not align['is_reverse'])
    reverse_count = sum(1 for align in alignments if align['is_reverse'])
    total = forward_count + reverse_count
    
    if total == 0:
        return {
            'total_alignments': 0,
            'forward_count': 0,
            'reverse_count': 0,
            'forward_ratio': 0.0,
            'is_mostly_forward': False
        }
    
    forward_ratio = forward_count / total
    
    return {
        'total_alignments': total,
        'forward_count': forward_count,
        'reverse_count': reverse_count,
        'forward_ratio': forward_ratio,
        'is_mostly_forward': forward_ratio >= 0.5
    }

def smart_truncate_extraction_region(
    region_type: str,
    extract_start: int,
    extract_end: int,
    alignment_start: int,
    alignment_end: int,
    telomeres_in_region: List,
    max_total_length: int = 5000000,
    max_alignment_length: int = 10000
) -> Tuple[int, int]:
    total_length = extract_end - extract_start + 1
    
    if total_length <= max_total_length:
        return extract_start, extract_end
    
    printv(f"  Extracted region too long: {total_length:,}bp > {max_total_length:,}bp")
    
    if region_type == 'first_alignment_with_preceding':
        telomere_adjusted_end = extract_start
        if telomeres_in_region:
            last_telomere = max(telomeres_in_region, key=lambda x: x['end'])
            telomere_adjusted_end = last_telomere['end']
            printv(f"    Preserving telomere region to: {telomere_adjusted_end:,}bp")
        
        alignment_length = alignment_end - alignment_start + 1
        alignment_part_to_keep = min(max_alignment_length, alignment_length)
        desired_end = max(telomere_adjusted_end, alignment_start + alignment_part_to_keep - 1)
        
        final_end = min(extract_start + max_total_length - 1, desired_end)
        
        printv(f"    After truncation: {extract_start:,}-{final_end:,} ({final_end - extract_start + 1:,}bp)")
        return extract_start, final_end
        
    else:
        telomere_adjusted_start = extract_end
        if telomeres_in_region:
            first_telomere = min(telomeres_in_region, key=lambda x: x['start'])
            telomere_adjusted_start = first_telomere['start']
            printv(f"    Preserving telomere region from: {telomere_adjusted_start:,}bp")
        
        alignment_length = alignment_end - alignment_start + 1
        alignment_part_to_keep = min(max_alignment_length, alignment_length)
        desired_start = min(telomere_adjusted_start, alignment_end - alignment_part_to_keep + 1)
        
        final_start = max(extract_end - max_total_length + 1, desired_start)
        
        printv(f"    After truncation: {final_start:,}-{extract_end:,} ({extract_end - final_start + 1:,}bp)")
        return final_start, extract_end

def find_first_alignment_with_preceding_sequence(
    alignments_by_pair: Dict, 
    contig_structure: Dict,
    min_repeats: int = 20,
    min_telomere_length: int = 500,
    threshold: float = 95.0, 
    min_length: int = 100,
    max_extraction_length: int = 50000
) -> List[Dict]:
    regions_to_extract = []
    
    for (ref_name, query_name), data in alignments_by_pair.items():
        alignments = data['alignments']
        
        if not alignments:
            continue
        
        direction_info = analyze_alignment_direction(alignments)
        
        alignments.sort(key=lambda x: x['s1'])
        
        first_good_alignment = None
        for align in alignments:
            if align['identity'] >= threshold and align['len1'] >= min_length:
                first_good_alignment = align
                break
        
        if first_good_alignment:
            contig_id = first_good_alignment['ref_name']
            if contig_id not in contig_structure:
                continue
                
            contig_len = contig_structure[contig_id]['sequence_length']
            contig_telomeres = contig_structure[contig_id]['telomeres']
            
            extract_start = 1
            extract_end = first_good_alignment['e1']
            
            telomeres_in_region = []
            for telomere in contig_telomeres:
                if telomere['start'] >= extract_start and telomere['end'] <= extract_end:
                    telomeres_in_region.append(telomere)
            
            has_telomere_in_region = len(telomeres_in_region) > 0
            
            region = {
                'ref_name': ref_name,
                'query_name': query_name,
                'contig_id': contig_id,
                'alignment_start': first_good_alignment['s1'],
                'alignment_end': first_good_alignment['e1'],
                'extract_start': extract_start,
                'extract_end': extract_end,
                'type': 'first_alignment_with_preceding',
                'alignment_info': first_good_alignment,
                'alignment_length': first_good_alignment['len1'],
                'identity': first_good_alignment['identity'],
                'contig_length': contig_len,
                'telomeres_in_region': telomeres_in_region,
                'has_telomere_in_region': has_telomere_in_region,
                'extracted_region_length': extract_end - extract_start + 1,
                'is_reverse': first_good_alignment['is_reverse'],
                'direction_info': direction_info,
                'max_extraction_length': max_extraction_length
            }
            regions_to_extract.append(region)
    
    return regions_to_extract

def find_last_alignment_with_subsequent_sequence(
    alignments_by_pair: Dict, 
    contig_structure: Dict,
    min_repeats: int = 20,
    min_telomere_length: int = 500,
    threshold: float = 95.0, 
    min_length: int = 100,
    max_extraction_length: int = 50000
) -> List[Dict]:
    regions_to_extract = []
    
    for (ref_name, query_name), data in alignments_by_pair.items():
        alignments = data['alignments']
        
        if not alignments:
            continue
        
        direction_info = analyze_alignment_direction(alignments)
        
        alignments.sort(key=lambda x: x['e1'], reverse=True)
        
        last_good_alignment = None
        for align in alignments:
            if align['identity'] >= threshold and align['len1'] >= min_length:
                last_good_alignment = align
                break
        
        if last_good_alignment:
            contig_id = last_good_alignment['ref_name']
            if contig_id not in contig_structure:
                continue
                
            contig_len = contig_structure[contig_id]['sequence_length']
            contig_telomeres = contig_structure[contig_id]['telomeres']
            
            extract_start = last_good_alignment['s1']
            extract_end = contig_len
            
            telomeres_in_region = []
            for telomere in contig_telomeres:
                if telomere['start'] >= extract_start and telomere['end'] <= extract_end:
                    telomeres_in_region.append(telomere)
            
            has_telomere_in_region = len(telomeres_in_region) > 0
            
            region = {
                'ref_name': ref_name,
                'query_name': query_name,
                'contig_id': contig_id,
                'alignment_start': last_good_alignment['s1'],
                'alignment_end': last_good_alignment['e1'],
                'extract_start': extract_start,
                'extract_end': extract_end,
                'type': 'last_alignment_with_subsequent',
                'alignment_info': last_good_alignment,
                'alignment_length': last_good_alignment['len1'],
                'identity': last_good_alignment['identity'],
                'contig_length': contig_len,
                'telomeres_in_region': telomeres_in_region,
                'has_telomere_in_region': has_telomere_in_region,
                'extracted_region_length': extract_end - extract_start + 1,
                'is_reverse': last_good_alignment['is_reverse'],
                'direction_info': direction_info,
                'max_extraction_length': max_extraction_length
            }
            regions_to_extract.append(region)
    
    return regions_to_extract

def find_last_forward_alignment_for_5prime(
    alignments_by_pair: Dict,
    contig_structure: Dict,
    ref_name: str,
    query_name: str,
    threshold: float = 95.0,
    min_length: int = 100
) -> Optional[Dict]:
    if (ref_name, query_name) not in alignments_by_pair:
        return None
    
    alignments = alignments_by_pair[(ref_name, query_name)]['alignments']
    
    forward_alignments = [
        align for align in alignments 
        if not align['is_reverse'] and align['identity'] >= threshold and align['len1'] >= min_length
    ]
    
    if not forward_alignments:
        return None
    
    forward_alignments.sort(key=lambda x: x['s1'])
    
    return forward_alignments[0]

def find_last_forward_alignment_for_3prime(
    alignments_by_pair: Dict,
    contig_structure: Dict,
    ref_name: str,
    query_name: str,
    threshold: float = 95.0,
    min_length: int = 100
) -> Optional[Dict]:
    if (ref_name, query_name) not in alignments_by_pair:
        return None
    
    alignments = alignments_by_pair[(ref_name, query_name)]['alignments']
    
    forward_alignments = [
        align for align in alignments 
        if not align['is_reverse'] and align['identity'] >= threshold and align['len1'] >= min_length
    ]
    
    if not forward_alignments:
        return None
    
    forward_alignments.sort(key=lambda x: x['e1'], reverse=True)
    
    return forward_alignments[0]

def extract_sequence_region(
    alignments_by_pair: Dict,
    contig_structure: Dict,
    region: Dict, 
    output_dir: str,
    min_repeats: int = 20,
    min_telomere_length: int = 500,
    extend_before: int = 0,
    extend_after: int = 0,
    max_total_length: int = 5000000,
    max_alignment_keep: int = 10000
) -> Optional[str]:
    try:
        contig_id = region['contig_id']
        
        if contig_id not in contig_structure:
            printv(f"  Error: Cannot find contig {contig_id}")
            return None
        
        contig_seq = contig_structure[contig_id]['sequence']
        seq_len = contig_structure[contig_id]['sequence_length']
        
        direction_info = region.get('direction_info', {})
        is_mostly_forward = direction_info.get('is_mostly_forward', False)
        current_is_reverse = region.get('is_reverse', False)
        
        need_special_handling = False
        if is_mostly_forward and current_is_reverse:
            need_special_handling = True
            printv(f"  Note: Main alignment direction is forward({direction_info['forward_ratio']:.1%}), but current alignment is reverse, special handling needed")
        
        extract_start = region['extract_start']
        extract_end = region['extract_end']
        
        if need_special_handling and region['type'] == 'first_alignment_with_preceding':
            last_forward_align = find_last_forward_alignment_for_5prime(
                alignments_by_pair,
                contig_structure,
                region['ref_name'],
                region['query_name']
            )
            
            if last_forward_align:
                printv(f"  Found forward alignment region for 5' end extraction: {last_forward_align['s1']}-{last_forward_align['e1']}")
                extract_start = 1
                extract_end = last_forward_align['e1']
                current_is_reverse = False
            else:
                printv(f"  Warning: Cannot find suitable forward alignment region, abandoning this extraction")
                return None
        
        elif need_special_handling and region['type'] == 'last_alignment_with_subsequent':
            last_forward_align = find_last_forward_alignment_for_3prime(
                alignments_by_pair,
                contig_structure,
                region['ref_name'],
                region['query_name']
            )
            
            if last_forward_align:
                printv(f"  Found forward alignment region for 3' end extraction: {last_forward_align['s1']}-{last_forward_align['e1']}")
                extract_start = last_forward_align['s1']
                extract_end = seq_len
                current_is_reverse = False
            else:
                printv(f"  Warning: Cannot find suitable forward alignment region, abandoning this extraction")
                return None
        
        original_length = extract_end - extract_start + 1
        if original_length > max_total_length:
            printv(f"  Original extraction region too long: {original_length:,}bp > {max_total_length:,}bp")
            
            extract_start, extract_end = smart_truncate_extraction_region(
                region_type=region['type'],
                extract_start=extract_start,
                extract_end=extract_end,
                alignment_start=region['alignment_start'],
                alignment_end=region['alignment_end'],
                telomeres_in_region=region.get('telomeres_in_region', []),
                max_total_length=max_total_length,
                max_alignment_length=max_alignment_keep
            )
        
        start = extract_start - extend_before - 1
        end = extract_end + extend_after - 1
        
        if start < 0:
            start = 0
        if end >= seq_len:
            end = seq_len - 1
        if start > end:
            printv(f"  Warning: Invalid extraction region {start+1}-{end+1} (sequence length {seq_len})")
            return None
        
        extracted_seq = contig_seq[start:end+1]
        
        if len(extracted_seq) == 0:
            printv(f"  Warning: Extracted empty sequence {start+1}-{end+1}")
            return None
        
        if current_is_reverse and not is_mostly_forward:
            printv(f"  Note: Alignment is reverse, performing reverse complement")
            extracted_seq = str(Seq(extracted_seq).reverse_complement())
            is_reversed_for_extraction = True
        else:
            is_reversed_for_extraction = False
        
        check_end = '5prime' if region['type'] == 'first_alignment_with_preceding' else '3prime'
        telomere_check = check_extracted_sequence_for_telomere(
            extracted_seq, min_repeats, min_telomere_length, check_end
        )
        
        if not telomere_check['telomeres_valid']:
            printv(f"  Warning: Extracted {check_end} end sequence does not contain correct telomere positions, abandoning this region")
            if telomere_check['has_telomere']:
                printv(f"    Sequence contains telomeres but not in correct positions: 5' end telomere={telomere_check['has_5prime_telomere']}, 3' end telomere={telomere_check['has_3prime_telomere']}")
            return None
        
        query_name = region['query_name']
        ref_name = region['ref_name']
        
        if '_' in query_name:
            chromosome_part = query_name.split('_')[0]
        else:
            chromosome_part = query_name
        
        clean_contig_name = contig_id.split()[0] if ' ' in contig_id else contig_id
        clean_contig_name = clean_contig_name.replace('|', '_').replace(':', '_')
        
        alignment_start = region['alignment_start']
        alignment_end = region['alignment_end']
        alignment_length = region['alignment_length']
        identity = region['identity']
        telomere_count = len(telomere_check['telomeres'])
        region_type = region['type']
        extract_start_abs = start + 1
        extract_end_abs = end + 1
        extracted_length = len(extracted_seq)
        
        forward_ratio = direction_info.get('forward_ratio', 0.0)
        is_mostly_forward = direction_info.get('is_mostly_forward', False)
        
        extracted_id = f"{clean_contig_name}_{chromosome_part}_{extract_start_abs}_{extract_end_abs}"
        description_parts = [
            f"region_type:{region_type}",
            f"alignment:{alignment_start}-{alignment_end}",
            f"aligned_to:{query_name}",
            f"ref_seq:{ref_name}",
            f"identity:{identity:.1f}%",
            f"alignment_len:{alignment_length}",
            f"extracted_len:{extracted_length}",
            f"telomeres_in_seq:{telomere_count}",
            f"has_telomere:YES",
            f"current_is_reverse:{'YES' if region.get('is_reverse', False) else 'NO'}",
            f"mostly_forward:{'YES' if is_mostly_forward else 'NO'}",
            f"forward_ratio:{forward_ratio:.1%}",
            f"was_reversed:{'YES' if is_reversed_for_extraction else 'NO'}",
            f"truncated:{'YES' if original_length > extracted_length else 'NO'}"
        ]
        
        description = " ".join(description_parts)
        
        extracted_record = SeqRecord(
            Seq(extracted_seq),
            id=extracted_id,
            description=description
        )
        
        output_file = os.path.join(output_dir, f"{extracted_id}.fa")
        SeqIO.write([extracted_record], output_file, "fasta")
        
        printv(f"  Extracted: {extracted_id} ({len(extracted_seq):,} bp)")
        if original_length > extracted_length:
            printv(f"    Truncated: {original_length:,}bp → {extracted_length:,}bp")
        printv(f"    Direction: {'Reverse complemented' if is_reversed_for_extraction else 'Forward'}")
        printv(f"    Telomeres: 5' end={telomere_check['has_5prime_telomere']}, 3' end={telomere_check['has_3prime_telomere']}")
        
        return output_file
        
    except Exception as e:
        printv(f"  Extraction error: {e}")
        import traceback
        traceback.print_exc()
        return None

def extract_sequences_from_coords(
    coords_file: str,
    contig_file: str,
    output_dir: str,
    fragment_type: str = "both",
    min_repeats: int = 20,
    min_telomere_length: int = 500,
    similarity_threshold: float = 95.0,
    min_alignment_length: int = 100,
    extend_before: int = 0,
    extend_after: int = 0,
    max_total_length: int = 5000000,
    max_alignment_keep: int = 10000
) -> Dict[str, List[str]]:
    print(f"Analyzing coords file: {coords_file}")
    print(f"Using contig file: {contig_file}")
    print(f"Extraction type: {fragment_type}")
    print(f"Telomere parameters: {min_repeats} repeats, {min_telomere_length}bp")
    print(f"Alignment filtering: similarity≥{similarity_threshold}%, length≥{min_alignment_length}bp")
    print(f"Length limit: maximum extraction length={max_total_length:,}bp, maximum alignment kept={max_alignment_keep:,}bp")
    
    os.makedirs(output_dir, exist_ok=True)
    
    alignments_by_pair = parse_coords_file_enhanced(coords_file)
    
    if not alignments_by_pair:
        print("Warning: No valid alignment information in coords file")
        return {'5prime': [], '3prime': []}
    
    print(f"Found {len(alignments_by_pair)} sequence pair alignments")
    
    print("Analyzing contig telomere structure...")
    contig_structure = analyze_contig_structure(
        contig_file, min_repeats, min_telomere_length
    )
    
    extracted_files = {'5prime': [], '3prime': []}
    
    if fragment_type in ['5prime', 'both']:
        print("\nExtracting 5' end regions (first collinear region and preceding sequences)...")
        regions_5prime = find_first_alignment_with_preceding_sequence(
            alignments_by_pair,
            contig_structure,
            min_repeats,
            min_telomere_length,
            similarity_threshold,
            min_alignment_length,
            max_total_length
        )
        
        print(f"  Found {len(regions_5prime)} 5' end candidate regions")
        
        for i, region in enumerate(regions_5prime):
            printv(f"  Processing 5' end candidate region {i+1}/{len(regions_5prime)}: {region['contig_id']}")
            extracted_file = extract_sequence_region(
                alignments_by_pair,
                contig_structure,
                region,
                output_dir,
                min_repeats,
                min_telomere_length,
                extend_before,
                extend_after,
                max_total_length,
                max_alignment_keep
            )
            if extracted_file:
                extracted_files['5prime'].append(extracted_file)
    
    if fragment_type in ['3prime', 'both']:
        print("\nExtracting 3' end regions (last collinear region and subsequent sequences)...")
        regions_3prime = find_last_alignment_with_subsequent_sequence(
            alignments_by_pair,
            contig_structure,
            min_repeats,
            min_telomere_length,
            similarity_threshold,
            min_alignment_length,
            max_total_length
        )
        
        print(f"  Found {len(regions_3prime)} 3' end candidate regions")
        
        for i, region in enumerate(regions_3prime):
            printv(f"  Processing 3' end candidate region {i+1}/{len(regions_3prime)}: {region['contig_id']}")
            extracted_file = extract_sequence_region(
                alignments_by_pair,
                contig_structure,
                region,
                output_dir,
                min_repeats,
                min_telomere_length,
                extend_before,
                extend_after,
                max_total_length,
                max_alignment_keep
            )
            if extracted_file:
                extracted_files['3prime'].append(extracted_file)
    
    all_extracted_records = []
    for file_type in ['5prime', '3prime']:
        for file_path in extracted_files[file_type]:
            if os.path.exists(file_path):
                try:
                    for record in SeqIO.parse(file_path, "fasta"):
                        all_extracted_records.append(record)
                except Exception as e:
                    print(f"Warning: Cannot read file {file_path}: {e}")
    
    if all_extracted_records:
        all_extracted_file = os.path.join(output_dir, f"all_extracted_regions.fa")
        SeqIO.write(all_extracted_records, all_extracted_file, "fasta")
        print(f"\nMerged {len(all_extracted_records)} extracted regions to: {all_extracted_file}")
    
    return extracted_files

def generate_extraction_report(
    extracted_files: Dict[str, List[str]],
    output_dir: str,
    coords_file: str,
    contig_file: str,
    params: Dict[str, Any]
) -> None:
    report_file = os.path.join(output_dir, "extraction_report.txt")
    
    with open(report_file, 'w') as f:
        f.write("Telomere Sequence Extraction Report (Enhanced Version) - With Length Limit\n")
        f.write("="*60 + "\n")
        f.write(f"Input coords file: {coords_file}\n")
        f.write(f"Input contig file: {contig_file}\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Extraction time: {params.get('timestamp', 'unknown')}\n\n")
        
        f.write("Extraction parameters:\n")
        for key, value in params.items():
            if key != 'timestamp':
                f.write(f"  {key}: {value}\n")
        f.write("\n")
        
        total_5prime = len(extracted_files.get('5prime', []))
        total_3prime = len(extracted_files.get('3prime', []))
        total_extracted = total_5prime + total_3prime
        
        f.write("Extraction statistics:\n")
        f.write(f"  5' end extraction regions: {total_5prime}\n")
        f.write(f"  3' end extraction regions: {total_3prime}\n")
        f.write(f"  Total extraction regions: {total_extracted}\n\n")
        
        if total_5prime > 0:
            f.write("5' end extraction files:\n")
            for file_path in extracted_files['5prime'][:20]:
                file_name = os.path.basename(file_path)
                try:
                    record = next(SeqIO.parse(file_path, "fasta"))
                    seq_len = len(record.seq)
                    f.write(f"  {file_name}: {seq_len:,} bp\n")
                    f.write(f"      Description: {record.description}\n")
                except:
                    f.write(f"  {file_name}\n")
            
            if total_5prime > 20:
                f.write(f"  ... {total_5prime - 20} more files not shown\n")
            f.write("\n")
        
        if total_3prime > 0:
            f.write("3' end extraction files:\n")
            for file_path in extracted_files['3prime'][:20]:
                file_name = os.path.basename(file_path)
                try:
                    record = next(SeqIO.parse(file_path, "fasta"))
                    seq_len = len(record.seq)
                    f.write(f"  {file_name}: {seq_len:,} bp\n")
                    f.write(f"      Description: {record.description}\n")
                except:
                    f.write(f"  {file_name}\n")
            
            if total_3prime > 20:
                f.write(f"  ... {total_3prime - 20} more files not shown\n")
    
    print(f"Extraction report saved: {report_file}")

def extract_telomere_regions(
    coords_file: str,
    contig_file: str,
    output_dir: str,
    fragment_type: str = "both",
    min_repeats: int = 20,
    min_telomere_length: int = 500,
    similarity_threshold: float = 95.0,
    min_alignment_length: int = 100,
    extend_before: int = 0,
    extend_after: int = 0,
    max_total_length: int = 5000000,
    max_alignment_keep: int = 10000,
    verbose: bool = False
) -> Dict[str, Any]:
    
    global VERBOSE
    VERBOSE = verbose
    
    try:
        if not os.path.exists(coords_file):
            return {
                'success': False,
                'error': f"Coords file does not exist: {coords_file}"
            }
        
        if not os.path.exists(contig_file):
            return {
                'success': False,
                'error': f"Contig file does not exist: {contig_file}"
            }
        
        os.makedirs(output_dir, exist_ok=True)
        
        extracted_files = extract_sequences_from_coords(
            coords_file=coords_file,
            contig_file=contig_file,
            output_dir=output_dir,
            fragment_type=fragment_type,
            min_repeats=min_repeats,
            min_telomere_length=min_telomere_length,
            similarity_threshold=similarity_threshold,
            min_alignment_length=min_alignment_length,
            extend_before=extend_before,
            extend_after=extend_after,
            max_total_length=max_total_length,
            max_alignment_keep=max_alignment_keep
        )
        
        params = {
            'fragment_type': fragment_type,
            'min_repeats': min_repeats,
            'min_telomere_length': min_telomere_length,
            'similarity_threshold': similarity_threshold,
            'min_alignment_length': min_alignment_length,
            'extend_before': extend_before,
            'extend_after': extend_after,
            'max_total_length': max_total_length,
            'max_alignment_keep': max_alignment_keep,
            'timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        generate_extraction_report(
            extracted_files,
            output_dir,
            coords_file,
            contig_file,
            params
        )
        
        total_5prime = len(extracted_files.get('5prime', []))
        total_3prime = len(extracted_files.get('3prime', []))
        total_extracted = total_5prime + total_3prime
        
        all_extracted_file = os.path.join(output_dir, "all_extracted_regions.fa")
        report_file = os.path.join(output_dir, "extraction_report.txt")
        
        result = {
            'success': True,
            'extracted_files': extracted_files,
            'all_extracted_file': all_extracted_file if os.path.exists(all_extracted_file) else None,
            'report_file': report_file,
            'statistics': {
                '5prime_regions': total_5prime,
                '3prime_regions': total_3prime,
                'total_regions': total_extracted
            },
            'parameters': params
        }
        
        return result
        
    except Exception as e:
        import traceback
        error_trace = traceback.format_exc()
        if verbose:
            print(f"Error during extraction process: {e}")
            print(error_trace)
        
        return {
            'success': False,
            'error': str(e),
            'traceback': error_trace
        }

def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description="Telomere Sequence Extraction Script (Enhanced Version) - Intelligently extract telomere sequence regions from NUCMER alignment results (with length limit)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Usage examples:
          # Extract all regions (5' end and 3' end), using default length limits
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir
          
          # Extract only 5' end regions
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir --type 5prime
          
          # Extract only 3' end regions
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir --type 3prime
          
          # Custom parameters
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir \\
            --min-repeats 20 --min-length 500 --similarity 95.0 --min-alignment 100
          
          # Extend extraction regions
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir \\
            --extend-before 1000 --extend-after 1000
            
          # Custom length limits
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir \\
            --max-total-length 50000 --max-alignment-keep 10000
            
          # Stricter or looser length limits
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir \\
            --max-total-length 20000 --max-alignment-keep 5000  # Stricter
          python telomere_extract_regions_enhanced.py -coords matches.coords -c telomere_contigs.fa -o output_dir \\
            --max-total-length 100000 --max-alignment-keep 20000  # Looser
        """)
    )
    
    parser.add_argument("-coords", "--coords-file", required=True, help="NUCMER coords file")
    parser.add_argument("-c", "--contig-file", required=True, help="Telomere-containing contigs file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    
    parser.add_argument("--type", choices=['5prime', '3prime', 'both'], default='both',
                       help="Extraction type: 5' end, 3' end, or both (default: both)")
    parser.add_argument("--min-repeats", type=int, default=20,
                       help="Minimum telomere repeats (default: 20)")
    parser.add_argument("--min-length", type=int, default=500,
                       help="Minimum telomere length(bp) (default: 500)")
    parser.add_argument("--similarity", type=float, default=95.0,
                       help="Alignment similarity threshold(%%%) (default: 95.0)")
    parser.add_argument("--min-alignment", type=int, default=100,
                       help="Minimum alignment length(bp) (default: 100)")
    parser.add_argument("--extend-before", type=int, default=0,
                       help="Additional extension bp before alignment start (default: 0)")
    parser.add_argument("--extend-after", type=int, default=0,
                       help="Additional extension bp after alignment end (default: 0)")
    parser.add_argument("--max-total-length", type=int, default=5000000,
                       help="Maximum total extraction length(bp) (default: 5,000,000)")
    parser.add_argument("--max-alignment-keep", type=int, default=10000,
                       help="Maximum alignment length to keep(bp) (default: 10000)")
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed output")
    
    return parser.parse_args()

def main():
    args = parse_command_line_args()
    
    print("Telomere Sequence Extraction Script (Enhanced Version) - With Length Limit")
    print("="*60)
    print("Enhanced features:")
    print("  1. Analyze complete alignment direction distribution")
    print("  2. Intelligent extraction strategy determination")
    print("  3. Verify telomere position correctness")
    print("  4. Intelligent length limits to prevent extracting overly long sequences")
    print("="*60)
    print(f"Input coords file: {args.coords_file}")
    print(f"Input contig file: {args.contig_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Extraction type: {args.type}")
    print(f"Telomere parameters: {args.min_repeats} repeats, {args.min_length}bp")
    print(f"Alignment filtering: similarity≥{args.similarity}%, length≥{args.min_alignment}bp")
    print(f"Sequence extension: {args.extend_before}bp before alignment, {args.extend_after}bp after alignment")
    print(f"Length limits: total length≤{args.max_total_length:,}bp, alignment part≤{args.max_alignment_keep:,}bp")
    print("="*60)
    
    result = extract_telomere_regions(
        coords_file=args.coords_file,
        contig_file=args.contig_file,
        output_dir=args.output_dir,
        fragment_type=args.type,
        min_repeats=args.min_repeats,
        min_telomere_length=args.min_length,
        similarity_threshold=args.similarity,
        min_alignment_length=args.min_alignment,
        extend_before=args.extend_before,
        extend_after=args.extend_after,
        max_total_length=args.max_total_length,
        max_alignment_keep=args.max_alignment_keep,
        verbose=args.verbose
    )
    
    if result['success']:
        stats = result['statistics']
        print("\n" + "="*60)
        print("Extraction completed!")
        print(f"5' end extraction: {stats['5prime_regions']} regions")
        print(f"3' end extraction: {stats['3prime_regions']} regions")
        print(f"Total extraction: {stats['total_regions']} regions")
        
        if result['all_extracted_file']:
            total_length = 0
            record_count = 0
            try:
                for record in SeqIO.parse(result['all_extracted_file'], "fasta"):
                    total_length += len(record.seq)
                    record_count += 1
                if record_count > 0:
                    avg_length = total_length / record_count
                    print(f"Average sequence length: {avg_length:,.1f} bp")
                    print(f"Longest sequence: {max(len(record.seq) for record in SeqIO.parse(result['all_extracted_file'], 'fasta')):,} bp")
            except:
                pass
            
            print(f"All extracted sequences: {result['all_extracted_file']}")
        
        print(f"Detailed report: {result['report_file']}")
        print("="*60)
    else:
        print(f"\nError: {result['error']}")
        if args.verbose and 'traceback' in result:
            print(f"Detailed error information:\n{result['traceback']}")
        sys.exit(1)

if __name__ == "__main__":
    main()