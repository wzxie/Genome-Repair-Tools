#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import argparse
from pathlib import Path
import shutil
import gzip
from typing import List, Dict, Any, Optional, Tuple, Union
from dataclasses import dataclass, field
import json
import tempfile
import re
import time
import hashlib
import random


@dataclass
class AssemblyConfig:
    hifi_files: Optional[List[str]] = None
    ont_ul_files: Optional[List[str]] = None
    clr_files: Optional[List[str]] = None
    
    threads: int = 32
    output_prefix: str = "assembly"
    
    hifiasm_primary: bool = False
    hifiasm_merge_haplotypes: bool = True
    hifiasm_l: int = 3
    verkko_memory_gb: int = 64
    verkko_cleanup: bool = True
    nextdenovo_genome_size: str = "1g"
    flye_genome_size: str = "1g"
    flye_iterations: int = 1
    flye_nano_type: str = "hq"
    
    auto_estimate_genome: bool = True
    target_depth: int = 30
    
    hifiasm_threads: Optional[int] = None
    verkko_threads: Optional[int] = None
    nextdenovo_threads: Optional[int] = None
    flye_threads: Optional[int] = None
    
    tools_to_run: List[str] = field(default_factory=lambda: ['all'])
    
    estimated_genome_size: Optional[str] = None
    estimated_genome_size_bp: Optional[float] = None
    
    def __post_init__(self):
        self.validate_files()
    
    def validate_files(self):
        for file_list in [self.hifi_files, self.ont_ul_files, self.clr_files]:
            if file_list:
                for file_path in file_list:
                    if not Path(file_path).exists():
                        raise FileNotFoundError(f"File not found: {file_path}")
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            k: v for k, v in self.__dict__.items() 
            if not k.startswith('_')
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'AssemblyConfig':
        return cls(**config_dict)
    
    @classmethod
    def from_json(cls, json_file: str) -> 'AssemblyConfig':
        with open(json_file, 'r') as f:
            config_dict = json.load(f)
        return cls.from_dict(config_dict)


class FastGenomeSizeEstimator:
    """
    FAST genome size estimation using ONLY:
    1. File size estimation (metadata-based) - Method 1
    2. Read statistics (sampling-based) - Method 2
    NO k-mer analysis for maximum speed
    """
    
    @staticmethod
    def _format_size(size_bp: float) -> str:
        if size_bp >= 1e9:
            return f"{size_bp/1e9:.2f}g"
        elif size_bp >= 1e6:
            return f"{int(size_bp/1e6)}m"
        else:
            return f"{int(size_bp/1000)}k"
    
    @staticmethod
    def _parse_genome_size(size_str: str) -> float:
        size_str = size_str.lower().strip()
        try:
            if size_str.endswith('g'):
                return float(size_str[:-1]) * 1e9
            elif size_str.endswith('m'):
                return float(size_str[:-1]) * 1e6
            elif size_str.endswith('k'):
                return float(size_str[:-1]) * 1e3
            else:
                return float(size_str)
        except:
            return 130e6
    
    @staticmethod
    def get_read_stats_sampled(file_path: Path, max_reads: int = 100000) -> Dict[str, Any]:
        stats = {
            'total_bases': 0,
            'total_reads': 0,
            'long_reads_count': 0,
            'long_reads_bases': 0,
            'n50': 0,
            'mean_length': 0,
            'max_length': 0,
            'min_length': float('inf')
        }
        
        try:
            lengths = []
            read_count = 0
            
            if str(file_path).endswith('.gz'):
                cmd = f"zcat {file_path} | awk 'NR%4==2 {{print length($0)}}' | head -{max_reads}"
            else:
                cmd = f"awk 'NR%4==2 {{print length($0)}}' {file_path} | head -{max_reads}"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        length = int(line.strip())
                        lengths.append(length)
                        stats['total_bases'] += length
                        stats['total_reads'] += 1
                        read_count += 1
                        
                        if length >= 10000:
                            stats['long_reads_count'] += 1
                            stats['long_reads_bases'] += length
                        
                        if length > stats['max_length']:
                            stats['max_length'] = length
                        if length < stats['min_length']:
                            stats['min_length'] = length
                
                if lengths:
                    lengths.sort(reverse=True)
                    stats['mean_length'] = sum(lengths) / len(lengths)
                    
                    half_total = sum(lengths) / 2
                    cumulative = 0
                    for length in lengths:
                        cumulative += length
                        if cumulative >= half_total:
                            stats['n50'] = length
                            break
        
        except subprocess.TimeoutExpired:
            print(f"  Warning: Read stats timeout for {file_path}")
        except Exception as e:
            print(f"  Warning: Failed to get read stats for {file_path}: {e}")
        
        return stats
    
    @staticmethod
    def estimate_from_read_stats(reads_files: List[str], 
                                target_depth: int = 30,
                                min_length: int = 5000,
                                max_reads_per_file: int = 100000) -> Optional[float]:
        total_bases = 0
        total_reads = 0
        long_reads_bases = 0
        long_reads_count = 0
        total_files_processed = 0
        
        print("  Method 2: Estimating from read statistics (sampling)...")
        print(f"    Sampling {max_reads_per_file:,} reads per file")
        
        for file_path in reads_files:
            print(f"    Processing: {Path(file_path).name}")
            stats = FastGenomeSizeEstimator.get_read_stats_sampled(
                Path(file_path), max_reads_per_file
            )
            
            if stats['total_reads'] > 0:
                total_bases += stats['total_bases']
                total_reads += stats['total_reads']
                long_reads_bases += stats['long_reads_bases']
                long_reads_count += stats['long_reads_count']
                total_files_processed += 1
                print(f"      Sampled {stats['total_reads']:,} reads, avg length: {stats['mean_length']:.0f} bp")
        
        if total_files_processed == 0 or total_reads == 0:
            print("    No reads found for statistical estimation")
            return None
        
        avg_length = total_bases / total_reads
        long_read_ratio = long_reads_bases / total_bases if total_bases > 0 else 0
        long_read_percentage = long_read_ratio * 100
        
        print(f"    Total sampled: {total_reads:,} reads, avg length: {avg_length:.0f} bp")
        print(f"    Long reads (>=10kb): {long_reads_count:,} reads, {long_read_percentage:.1f}% of bases")
        
        if long_read_ratio > 0.5:
            effective_bases = long_reads_bases
            print(f"    Using long reads for estimation ({long_read_ratio:.1%} of data)")
        else:
            effective_bases = total_bases
            print(f"    Using all reads for estimation")
        
        estimated_bp = effective_bases / target_depth
        
        if avg_length > 10000:
            estimated_bp = estimated_bp * 0.85
            print(f"    Applied ONT correction factor (0.85)")
        elif avg_length > 5000:
            estimated_bp = estimated_bp * 0.95
            print(f"    Applied HiFi correction factor (0.95)")
        
        if estimated_bp > 0:
            print(f"    Statistical estimate: {FastGenomeSizeEstimator._format_size(estimated_bp)}")
            return estimated_bp
        
        return None
    
    @staticmethod
    def estimate_from_file_size(reads_files: List[str], 
                               target_depth: int = 30,
                               read_type: str = "hifi") -> Optional[float]:
        total_bases = 0
        total_raw_size = 0
        
        print("  Method 1: Estimating from file size (metadata only)...")
        
        if read_type == "hifi":
            compression_ratio = 3.5
            sequence_ratio = 0.8
            error_correction = 0.95
        elif read_type == "ont":
            compression_ratio = 2.5
            sequence_ratio = 0.7
            error_correction = 0.85
        else:
            compression_ratio = 3.0
            sequence_ratio = 0.75
            error_correction = 0.90
        
        for file_path in reads_files:
            file_path = Path(file_path)
            file_size = file_path.stat().st_size
            total_raw_size += file_size
            
            if str(file_path).endswith('.gz'):
                estimated_uncompressed = file_size * compression_ratio
                estimated_bases = estimated_uncompressed * sequence_ratio
                total_bases += estimated_bases
                print(f"    {file_path.name}: {file_size/1e9:.2f} GB (compressed)")
            else:
                try:
                    with open(file_path, 'r') as f:
                        lines = 0
                        bases = 0
                        for i, line in enumerate(f):
                            if i >= 1000:
                                break
                            if i % 4 == 1:
                                bases += len(line.strip())
                                lines += 1
                        
                        if lines > 0:
                            avg_length = bases / lines
                            avg_line_size = (avg_length * 2 + 80)
                            total_reads = file_size / avg_line_size
                            total_bases += total_reads * avg_length
                    print(f"    {file_path.name}: {file_size/1e9:.2f} GB (uncompressed)")
                except:
                    total_bases += file_size * sequence_ratio
        
        if total_bases == 0:
            print("    Could not estimate from file size")
            return None
        
        total_bases_corrected = total_bases * error_correction
        estimated_bp = total_bases_corrected / target_depth
        
        size_gb = estimated_bp / 1e9
        print(f"    Total raw file size: {total_raw_size / 1e9:.2f} GB")
        print(f"    Estimated total bases: {total_bases_corrected / 1e9:.2f} Gb")
        print(f"    File size estimate: {FastGenomeSizeEstimator._format_size(estimated_bp)}")
        
        return estimated_bp
    
    @staticmethod
    def estimate_genome_size(reads_files: List[str],
                            read_type: str = "hifi",
                            target_depth: int = 30,
                            user_provided: Optional[str] = None) -> Tuple[str, str]:
        print("\n" + "="*60)
        print("FAST GENOME SIZE ESTIMATION")
        print("Using: File size estimation (metadata) + Read statistics (sampling)")
        print("NO k-mer analysis - Maximum speed")
        print("="*60)
        
        if user_provided:
            print(f"\nUsing user-provided genome size: {user_provided}")
            return user_provided, "user-provided"
        
        estimations = []
        
        file_estimate_bp = FastGenomeSizeEstimator.estimate_from_file_size(
            reads_files, target_depth, read_type
        )
        
        if file_estimate_bp and file_estimate_bp > 0:
            estimations.append(('file_size', file_estimate_bp))
        
        stat_estimate_bp = FastGenomeSizeEstimator.estimate_from_read_stats(
            reads_files, target_depth, max_reads_per_file=100000
        )
        
        if stat_estimate_bp and stat_estimate_bp > 0:
            estimations.append(('read_stats', stat_estimate_bp))
        
        if not estimations:
            print("\nWARNING: Could not estimate genome size, using default 1g")
            return "1g", "default"
        
        if len(estimations) == 1:
            method, value = estimations[0]
            final_bp = value
            method_desc = method.replace('_', '-')
        else:
            file_bp = next(v for m, v in estimations if m == 'file_size')
            stat_bp = next(v for m, v in estimations if m == 'read_stats')
            
            ratio = max(stat_bp, file_bp) / min(stat_bp, file_bp)
            
            print(f"\nEstimation results:")
            print(f"  Method 1 (File size):     {FastGenomeSizeEstimator._format_size(file_bp)}")
            print(f"  Method 2 (Read stats):     {FastGenomeSizeEstimator._format_size(stat_bp)}")
            print(f"  Ratio: {ratio:.2f}x")
            
            if ratio <= 1.5:
                final_bp = (stat_bp + file_bp) / 2
                method_desc = "average (file-size + read-stats)"
                print(f"  Estimates are consistent, using average")
            else:
                if read_type == "hifi":
                    weight_stat = 0.7
                    weight_file = 0.3
                    final_bp = stat_bp * weight_stat + file_bp * weight_file
                    method_desc = "weighted average (70% read-stats, 30% file-size)"
                    print(f"  HiFi data: trusting read statistics more (70/30 weight)")
                elif read_type == "ont":
                    weight_stat = 0.5
                    weight_file = 0.5
                    final_bp = stat_bp * weight_stat + file_bp * weight_file
                    method_desc = "average (50/50 weight)"
                    print(f"  ONT data: using equal weights")
                else:
                    weight_stat = 0.65
                    weight_file = 0.35
                    final_bp = stat_bp * weight_stat + file_bp * weight_file
                    method_desc = "weighted average (65% read-stats, 35% file-size)"
                    print(f"  CLR data: trusting read statistics more")
        
        min_size = 10e6
        max_size = 200e9
        
        if final_bp < min_size:
            print(f"  Estimate too small ({final_bp/1e6:.1f} Mb), adjusting to minimum {min_size/1e6:.0f} Mb")
            final_bp = min_size
        elif final_bp > max_size:
            print(f"  Estimate too large ({final_bp/1e9:.1f} Gb), adjusting to maximum {max_size/1e9:.0f} Gb")
            final_bp = max_size
        
        formatted_size = FastGenomeSizeEstimator._format_size(final_bp)
        
        print(f"\nFinal genome size: {formatted_size}")
        print(f"  ({final_bp/1e6:.1f} Mb / {final_bp/1e9:.2f} Gb)")
        print(f"  Method: {method_desc}")
        
        return formatted_size, method_desc


def print_title():
    print("\n" + "="*70)
    print("      Genome Assembly Tool Control Center (HiFi / ONT / CLR)")
    print("      hifiasm | verkko | nextDenovo | flye")
    print("="*70)
    print("Supports: HiFi, ONT Ultra-Long, PacBio CLR")
    print("FAST genome size estimation (file size + sampling) - NO k-mer")
    print("hifiasm: Outputs both haplotypes (hap1 + hap2) and merges them")
    print("hifiasm ONT support: --ont for simplex, --ul for ultra-long+HiFi")
    print("="*70)


def ask(prompt: str, options: Optional[List[str]] = None, 
        default: Optional[str] = None, required: bool = True) -> str:
    while True:
        print(prompt)
        if options:
            print(f"Options: {', '.join(options)}")
        if default is not None:
            prompt_line = f"[Default: {default}] > "
        else:
            prompt_line = "> "

        choice = input(prompt_line).strip()

        if not choice and default is not None:
            return default
        if not choice and required:
            print("This field is required, please enter again.\n")
            continue
        if options and choice.lower() not in [opt.lower() for opt in options]:
            print(f"Invalid input, please choose from {options}.\n")
            continue
        return choice


def validate_files(file_list: str) -> Optional[List[str]]:
    if not file_list:
        return None
    files = file_list.split()
    for f in files:
        if not Path(f).exists():
            print(f"File not found: {f}")
            return None
    return files


def gfa_to_fasta_awk(gfa_file: Path, fasta_file: Path) -> bool:
    print(f"Converting {gfa_file} -> {fasta_file}")
    
    try:
        cmd = f"awk '/^S/{{print \">\"$2; print $3}}' {gfa_file} > {fasta_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            if os.path.getsize(fasta_file) > 0:
                with open(fasta_file, 'r') as f:
                    content = f.read(1000)
                    if '*' in content and len(content.strip()) < 100:
                        with open(fasta_file, 'r') as f_full:
                            seqs = []
                            current_seq = []
                            for line in f_full:
                                if line.startswith('>'):
                                    if current_seq:
                                        seqs.append(''.join(current_seq))
                                    current_seq = []
                                elif not line.startswith('>'):
                                    current_seq.append(line.strip())
                            if current_seq:
                                seqs.append(''.join(current_seq))
                            
                            all_asterisks = all(seq == '*' or all(c in '*Nn' for c in seq) for seq in seqs)
                            if all_asterisks:
                                print(f"Converted file contains only asterisks/Ns, removing")
                                os.remove(fasta_file)
                                return False
                
                print(f"Conversion complete: {fasta_file}")
                return True
            else:
                print(f"Converted file is empty: {fasta_file}")
                if fasta_file.exists():
                    os.remove(fasta_file)
                return False
        else:
            print(f"awk conversion failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"GFA to FASTA conversion failed: {e}")
        return False


def filter_empty_sequences(fasta_file: Path) -> bool:
    temp_file = fasta_file.with_suffix('.temp.fa')
    
    try:
        with open(fasta_file, 'r') as infile, open(temp_file, 'w') as outfile:
            current_header = None
            current_seq = []
            
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    if current_header and current_seq:
                        seq_str = ''.join(current_seq)
                        if seq_str and seq_str != '*' and not all(c in '*Nn' for c in seq_str):
                            outfile.write(f"{current_header}\n")
                            for i in range(0, len(seq_str), 80):
                                outfile.write(seq_str[i:i+80] + "\n")
                    
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_header and current_seq:
                seq_str = ''.join(current_seq)
                if seq_str and seq_str != '*' and not all(c in '*Nn' for c in seq_str):
                    outfile.write(f"{current_header}\n")
                    for i in range(0, len(seq_str), 80):
                        outfile.write(seq_str[i:i+80] + "\n")
        
        if temp_file.stat().st_size > 0:
            shutil.move(temp_file, fasta_file)
            return True
        else:
            temp_file.unlink()
            return False
        
    except Exception as e:
        print(f"Filtering empty sequences failed {fasta_file}: {e}")
        if temp_file.exists():
            temp_file.unlink()
        return False


def merge_fasta_files(fasta_files: List[Path], output_name: str = "all_assemblies.fasta") -> Optional[Path]:
    if not fasta_files:
        print("No FASTA files found to merge")
        return None
    
    unique_files = []
    seen_names = set()
    for fasta_file in fasta_files:
        if fasta_file.name not in seen_names:
            seen_names.add(fasta_file.name)
            unique_files.append(fasta_file)
    
    print(f"Found {len(fasta_files)} files, {len(unique_files)} after deduplication")
    
    print("\nFiltering empty sequences...")
    filtered_files = []
    for fasta_file in unique_files:
        if filter_empty_sequences(fasta_file):
            if fasta_file.stat().st_size > 0:
                filtered_files.append(fasta_file)
            else:
                print(f"File {fasta_file.name} is empty after filtering, skipped")
                fasta_file.unlink()
        else:
            filtered_files.append(fasta_file)
    
    if not filtered_files:
        print("All files are empty after filtering")
        return None
    
    output_path = Path.cwd() / output_name
    
    print(f"\nMerging {len(filtered_files)} FASTA files to {output_path}")
    
    with open(output_path, 'w') as outfile:
        for fasta_file in filtered_files:
            try:
                with open(fasta_file, 'r') as infile:
                    content = infile.read().strip()
                    if content:
                        outfile.write(content)
                        outfile.write("\n")
                    print(f"  Added: {fasta_file.name}")
            except Exception as e:
                print(f"  Processing {fasta_file} failed: {e}")
    
    if output_path.exists() and output_path.stat().st_size > 0:
        with open(output_path, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        print(f"Merge complete! Total {seq_count} valid sequences")
        print(f"Merged file: {output_path}")
        return output_path
    else:
        print("Merged file is empty or not found")
        if output_path.exists():
            output_path.unlink()
        return None


class AssemblyTool:
    
    def __init__(self, name: str, supported_types: List[str]):
        self.name = name
        self.supported_types = supported_types
        self.description = ""
        self.config = None
    
    def validate_input(self, data: Dict[str, Any]) -> bool:
        for data_type in self.supported_types:
            if data.get(data_type):
                return True
        return False
    
    def get_threads(self, data: Dict[str, Any]) -> int:
        tool_thread_key = f"{self.name}_threads"
        if tool_thread_key in data and data[tool_thread_key]:
            return data[tool_thread_key]
        return data.get('threads', 32)
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        raise NotImplementedError
    
    def process_output(self, output_path: str) -> List[Path]:
        raise NotImplementedError
    
    def run(self, data: Dict[str, Any], config: AssemblyConfig) -> Optional[List[Path]]:
        if not self.validate_input(data):
            print(f"{self.name} missing required input data type")
            return None
        
        try:
            cmd, output_path = self.build_command(data, config)
            if not cmd:
                return None
            
            print(f"Running {self.name}: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"{self.name} completed successfully")
                return self.process_output(output_path)
            else:
                print(f"{self.name} failed: {result.stderr}")
                return None
                
        except FileNotFoundError:
            print(f"{self.name} command not found, please ensure it's installed and in PATH")
            return None
        except Exception as e:
            print(f"{self.name} execution error: {e}")
            return None


class HifiasmTool(AssemblyTool):
    
    def __init__(self):
        super().__init__("hifiasm", ['hifi', 'clr', 'ont_ul'])
        self.description = "PacBio HiFi/CLR/ONT assembly tool - outputs both haplotypes"
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        prefix = f"{data['output_prefix']}.hifiasm"
        threads = self.get_threads(data)
        
        cmd = ["hifiasm", "-o", prefix, "-t", str(threads)]
        
        has_hifi = bool(data.get('hifi'))
        has_clr = bool(data.get('clr'))
        has_ont = bool(data.get('ont_ul'))
        
        if has_hifi:
            cmd += data['hifi']
            if has_ont:
                cmd += ["--ul"] + data['ont_ul']
                cmd += ["-l", "0"]
                print("hifiasm: HiFi + Ultra-Long hybrid mode (--ul -l 0)")
        elif has_clr:
            cmd += ["--pacbio"] + data['clr']
            if has_ont:
                cmd += ["--ul"] + data['ont_ul']
                cmd += ["-l", "0"]
                print("hifiasm: CLR + Ultra-Long hybrid mode (--ul -l 0)")
        elif has_ont:
            cmd += ["--ont"] + data['ont_ul']
            print("hifiasm: ONT simplex mode (--ont)")
        else:
            print("hifiasm requires HiFi, CLR, or ONT input.")
            return [], ""
        
        if hasattr(config, 'hifiasm_primary') and config.hifiasm_primary:
            cmd += ["--primary"]
            print("hifiasm: Primary assembly mode (single assembly)")
        else:
            print("hifiasm: Haplotype-resolved mode (outputs hap1 and hap2)")
        
        return cmd, prefix
    
    def process_output(self, output_path: str) -> List[Path]:
        fasta_files = []
        base_dir = Path.cwd()
        
        output_name = Path(output_path).name
        
        found_gfa_files = []
        
        patterns = [
            f"{output_name}.bp.p_ctg.gfa",
            f"{output_name}*.bp.p_ctg.gfa",
            "*.bp.p_ctg.gfa"
        ]
        
        for pattern in patterns:
            try:
                matching_files = list(base_dir.glob(pattern))
                for gfa_file in matching_files:
                    if gfa_file.exists() and gfa_file.stat().st_size > 0:
                        if gfa_file not in found_gfa_files:
                            found_gfa_files.append(gfa_file)
                            print(f"Found GFA file: {gfa_file}")
            except Exception as e:
                continue
        
        if not found_gfa_files:
            for pattern in ["**/*.bp.p_ctg.gfa"]:
                try:
                    matching_files = list(base_dir.glob(pattern))
                    for gfa_file in matching_files:
                        if gfa_file.exists() and gfa_file.stat().st_size > 0:
                            if "hifiasm" in gfa_file.name.lower() and gfa_file not in found_gfa_files:
                                found_gfa_files.append(gfa_file)
                                print(f"Found GFA file recursively: {gfa_file}")
                except Exception as e:
                    continue
        
        if not found_gfa_files:
            print(f"No .bp.p_ctg.gfa files found for hifiasm")
            print(f"Current directory: {base_dir}")
            print(f"Output name: {output_name}")
            print("Note: .bp.p_ctg.gfa is produced when using ONT ultra-long data with HiFi")
            return fasta_files
        
        merge_haplotypes = getattr(self.config, 'hifiasm_merge_haplotypes', True) if hasattr(self, 'config') else True
        
        if merge_haplotypes and len(found_gfa_files) > 1:
            print(f"\nFound {len(found_gfa_files)} .bp.p_ctg.gfa files, merging...")
            all_contigs = []
            
            for gfa_file in found_gfa_files:
                temp_fasta = base_dir / f"temp_{gfa_file.stem}.fa"
                
                if gfa_to_fasta_awk(gfa_file, temp_fasta):
                    with open(temp_fasta, 'r') as f:
                        lines = f.readlines()
                    
                    prefix = gfa_file.stem.replace(f"{output_name}.", "").replace(".bp.p_ctg", "")
                    if prefix:
                        modified_lines = []
                        for line in lines:
                            if line.startswith('>'):
                                header = line[1:].strip()
                                modified_lines.append(f">{prefix}_{header}\n")
                            else:
                                modified_lines.append(line)
                        all_contigs.extend(modified_lines)
                    else:
                        all_contigs.extend(lines)
                    
                    if temp_fasta.exists():
                        temp_fasta.unlink()
            
            if all_contigs:
                merged_fasta_path = base_dir / f"{self.name}_bp_p_ctg.fa"
                with open(merged_fasta_path, 'w') as out_f:
                    out_f.writelines(all_contigs)
                
                total_count = sum(1 for line in all_contigs if line.startswith('>'))
                print(f"\nMerged: {merged_fasta_path}")
                print(f"  Total sequences: {total_count}")
                
                fasta_files.append(merged_fasta_path)
                return fasta_files
        
        elif len(found_gfa_files) == 1:
            gfa_file = found_gfa_files[0]
            fasta_path = base_dir / f"{self.name}_bp_p_ctg.fa"
            
            if gfa_to_fasta_awk(gfa_file, fasta_path):
                fasta_files.append(fasta_path)
                return fasta_files
        else:
            for gfa_file in found_gfa_files:
                output_name = f"{self.name}_{gfa_file.stem}.fa"
                fasta_path = base_dir / output_name
                
                if gfa_to_fasta_awk(gfa_file, fasta_path):
                    fasta_files.append(fasta_path)
            
            if fasta_files:
                return fasta_files
        
        return fasta_files


class VerkkoTool(AssemblyTool):
    
    def __init__(self):
        super().__init__("verkko", ['hifi'])
        self.description = "HiFi + ONT hybrid assembly tool"
    
    def _calculate_recommended_memory(self, genome_size_bp: float, 
                                       has_ont: bool = True) -> int:
        genome_size_gb = genome_size_bp / 1e9
        base_memory = 32
        
        if has_ont:
            if genome_size_gb <= 0.1:
                multiplier = 1.5
            elif genome_size_gb <= 0.5:
                multiplier = 2.0
            elif genome_size_gb <= 1.0:
                multiplier = 2.5
            elif genome_size_gb <= 3.0:
                multiplier = 3.0
            else:
                multiplier = 3.5
        else:
            if genome_size_gb <= 0.1:
                multiplier = 1.0
            elif genome_size_gb <= 0.5:
                multiplier = 1.2
            elif genome_size_gb <= 1.0:
                multiplier = 1.5
            elif genome_size_gb <= 3.0:
                multiplier = 2.0
            else:
                multiplier = 2.5
        
        recommended = int(base_memory * multiplier * (genome_size_gb ** 0.7))
        min_memory = 32
        max_memory = 1024
        
        recommended = max(min_memory, min(recommended, max_memory))
        return recommended
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        out_dir = f"{data['output_prefix']}.verkko"
        Path(out_dir).parent.mkdir(parents=True, exist_ok=True)
        
        threads = self.get_threads(data)
        cmd = ["verkko", "-d", out_dir]
        
        has_hifi = bool(data.get('hifi'))
        has_ont = bool(data.get('ont_ul'))
        
        if has_hifi:
            cmd += ["--hifi"] + data['hifi']
        
        if has_ont:
            cmd += ["--nano"] + data['ont_ul']
        else:
            cmd += ["--no-nano"]
        
        if config.verkko_cleanup:
            cmd += ["--clean"]
        
        cmd += ["--local-cpus", str(threads)]
        
        if config.estimated_genome_size_bp and config.estimated_genome_size_bp > 0:
            genome_size_gb = config.estimated_genome_size_bp / 1e9
            mem_gb = self._calculate_recommended_memory(config.estimated_genome_size_bp, has_ont)
            
            print(f"verkko: Estimated genome size: {genome_size_gb:.2f} Gb")
            print(f"verkko: Data type: {'HiFi+ONT' if has_ont else 'HiFi-only'}")
            print(f"verkko: Recommended memory: {mem_gb} GB")
            
            if config.verkko_memory_gb > mem_gb:
                print(f"verkko: Using user-specified memory: {config.verkko_memory_gb} GB")
                mem_gb = config.verkko_memory_gb
        else:
            mem_gb = config.verkko_memory_gb
            print(f"verkko: No genome size estimate, using default memory: {mem_gb} GB")
        
        mem_mb = mem_gb * 1024
        cmd += ["--local-memory", str(mem_mb)]
        
        return cmd, out_dir
    
    def process_output(self, output_path: str) -> List[Path]:
        fasta_files = []
        output_dir = Path(output_path)
        
        if not output_dir.exists():
            print(f"verkko output directory not found: {output_path}")
            return fasta_files
        
        priority_files = ["assembly.fasta", "assembly.homopolymer-compressed.fasta"]
        
        for filename in priority_files:
            fasta_file = output_dir / filename
            if fasta_file.exists() and fasta_file.stat().st_size > 0:
                with open(fasta_file, 'r') as f:
                    content = f.read(1000)
                    if '*' not in content or len(content.strip()) > 10:
                        new_name = f"{self.name}.fa"
                        new_path = Path.cwd() / new_name
                        
                        print(f"Copying {fasta_file} -> {new_path}")
                        shutil.copy2(fasta_file, new_path)
                        fasta_files.append(new_path)
                        return fasta_files
        
        backup_patterns = ["**/*.fasta", "**/*.fa"]
        for pattern in backup_patterns:
            for fasta_file in output_dir.glob(pattern):
                if fasta_file.stat().st_size > 0:
                    with open(fasta_file, 'r') as f:
                        content = f.read(1000)
                        if '*' not in content or len(content.strip()) > 10:
                            new_name = f"{self.name}.fa"
                            new_path = Path.cwd() / new_name
                            
                            print(f"Copying {fasta_file} -> {new_path}")
                            shutil.copy2(fasta_file, new_path)
                            fasta_files.append(new_path)
                            return fasta_files
        
        print(f"No valid FASTA files found in {output_path}")
        return fasta_files


class NextDenovoTool(AssemblyTool):
    
    def __init__(self):
        super().__init__("nextdenovo", ['hifi', 'ont_ul', 'clr'])
        self.description = "ONT/HiFi/CLR assembly tool"
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        workdir = f"{data['output_prefix']}.nextDenovo"
        Path(workdir).mkdir(parents=True, exist_ok=True)
        cfg_file = f"{workdir}/run.cfg"
        
        read_type = "ont"
        if data.get('ont_ul'):
            input_files = data['ont_ul']
            read_type = "ont"
            print("nextDenovo: Using ONT data")
        elif data.get('hifi'):
            input_files = data['hifi']
            read_type = "hifi"
            print("nextDenovo: Using HiFi data")
        elif data.get('clr'):
            input_files = data['clr']
            read_type = "clr"
            print("nextDenovo: Using CLR data")
        else:
            print("nextDenovo requires ONT, HiFi, or CLR input.")
            return [], ""
        
        with open(f"{workdir}/input.fofn", 'w') as f:
            for read_file in input_files:
                abs_path = Path(read_file).resolve()
                f.write(f"{abs_path}\n")
        
        threads = self.get_threads(data)
        
        if threads >= 64:
            parallel_jobs = 8
            per_task_threads = 8
        elif threads >= 32:
            parallel_jobs = 4
            per_task_threads = 8
        elif threads >= 16:
            parallel_jobs = 2
            per_task_threads = 8
        else:
            parallel_jobs = 1
            per_task_threads = threads
        
        pa_correction = min(5, max(1, parallel_jobs // 2))
        
        original_genome_size = config.nextdenovo_genome_size
        user_specified_size = (config.nextdenovo_genome_size != "1g")
        
        if config.auto_estimate_genome and not user_specified_size and config.estimated_genome_size:
            config.nextdenovo_genome_size = config.estimated_genome_size
            print(f"nextDenovo using estimated genome size: {config.nextdenovo_genome_size}")
        elif user_specified_size:
            print(f"nextDenovo using user-specified genome size: {config.nextdenovo_genome_size}")
        else:
            print(f"nextDenovo using default genome size: {config.nextdenovo_genome_size}")
        
        read_cutoff = self.adjust_read_cutoff(config.nextdenovo_genome_size)
        
        if read_type == "ont":
            minimap2_raw = f"-t {per_task_threads} -x ava-ont"
            minimap2_cns = f"-t {per_task_threads} -x ava-ont -k 17 -w 17 --minlen 200 --maxhan1 1000"
        elif read_type == "hifi":
            minimap2_raw = f"-t {per_task_threads} -x ava-pb"
            minimap2_cns = f"-t {per_task_threads} -x ava-pb -k 17 -w 17 --minlen 200 --maxhan1 1000"
        else:
            minimap2_raw = f"-t {per_task_threads} -x ava-pb"
            minimap2_cns = f"-t {per_task_threads} -x ava-pb -k 17 -w 17 --minlen 200 --maxhan1 1000"
        
        config_content = f"""
[General]
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
parallel_jobs = {parallel_jobs}
input_type = raw
read_type = {read_type}
input_fofn = input.fofn
workdir = {workdir}

[correct_option]
read_cutoff = {read_cutoff}
genome_size = {config.nextdenovo_genome_size}
sort_options = -m 20g -t {per_task_threads}
minimap2_options_raw = {minimap2_raw}
pa_correction = {pa_correction}
correction_options = -p {per_task_threads}

[assemble_option]
minimap2_options_cns = {minimap2_cns}
nextgraph_options = -a 1
"""
        
        with open(cfg_file, 'w') as f:
            f.write(config_content.strip())
        
        print(f"nextDenovo configuration file generated: {cfg_file}")
        print(f"Thread allocation: {threads} total threads -> {parallel_jobs} parallel tasks x {per_task_threads} threads/task")
        print(f"Read type: {read_type}")
        
        config.nextdenovo_genome_size = original_genome_size
        
        return ["nextDenovo", cfg_file], workdir
    
    def adjust_read_cutoff(self, genome_size: str) -> str:
        genome_size_lower = genome_size.lower()
        if genome_size_lower.endswith('g'):
            size_value = float(genome_size_lower[:-1]) * 1e9
        elif genome_size_lower.endswith('m'):
            size_value = float(genome_size_lower[:-1]) * 1e6
        elif genome_size_lower.endswith('k'):
            size_value = float(genome_size_lower[:-1]) * 1e3
        else:
            size_value = float(genome_size_lower)
        
        if size_value >= 1e9:
            return "2k"
        elif size_value >= 100e6:
            return "1k"
        else:
            return "500"
    
    def process_output(self, output_path: str) -> List[Path]:
        fasta_files = []
        output_dir = Path(output_path)
        
        if not output_dir.exists():
            print(f"nextdenovo output directory not found: {output_path}")
            return fasta_files
        
        target_file = output_dir / "03.ctg_graph" / "nd.asm.fasta"
        
        if target_file.exists() and target_file.stat().st_size > 0:
            new_name = f"{self.name}.fa"
            new_path = Path.cwd() / new_name
            
            file_size_mb = target_file.stat().st_size / (1024 * 1024)
            print(f"Found nextDenovo assembly: {target_file} ({file_size_mb:.2f} MB)")
            print(f"Copying to: {new_path}")
            shutil.copy2(target_file, new_path)
            fasta_files.append(new_path)
        else:
            print(f"nextDenovo assembly file not found at expected location: {target_file}")
            alt_files = list(output_dir.glob("**/*.fasta")) + list(output_dir.glob("**/*.fa"))
            if alt_files:
                print(f"Found alternative FASTA files:")
                for f in alt_files[:3]:
                    print(f"  {f}")
        
        return fasta_files


class FlyeTool(AssemblyTool):
    
    def __init__(self):
        super().__init__("flye", ['hifi', 'ont_ul', 'clr'])
        self.description = "Long-read assembly tool for multiple data types"
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        out_dir = f"{data['output_prefix']}.flye"
        Path(out_dir).parent.mkdir(parents=True, exist_ok=True)
        
        threads = self.get_threads(data)
        
        original_genome_size = config.flye_genome_size
        user_specified_size = (config.flye_genome_size != "1g")
        
        if config.auto_estimate_genome and not user_specified_size and config.estimated_genome_size:
            config.flye_genome_size = config.estimated_genome_size
            print(f"Flye using estimated genome size: {config.flye_genome_size}")
        elif user_specified_size:
            print(f"Flye using user-specified genome size: {config.flye_genome_size}")
        else:
            print(f"Flye using default genome size: {config.flye_genome_size}")
        
        cmd = ["flye", "--out-dir", out_dir, "--threads", str(threads)]
        
        if data.get('hifi'):
            cmd += ["--pacbio-hifi"] + data['hifi']
        elif data.get('ont_ul'):
            cmd += [f"--nano-{config.flye_nano_type}"] + data['ont_ul']
        elif data.get('clr'):
            cmd += ["--pacbio-raw"] + data['clr']
        else:
            print("Flye requires HiFi / ONT / CLR input.")
            return [], ""
        
        cmd += ["--genome-size", config.flye_genome_size]
        cmd += ["--iterations", str(config.flye_iterations)]
        
        config.flye_genome_size = original_genome_size
        
        return cmd, out_dir
    
    def process_output(self, output_path: str) -> List[Path]:
        fasta_files = []
        output_dir = Path(output_path)
        
        if not output_dir.exists():
            print(f"flye output directory not found: {output_path}")
            return fasta_files
        
        priority_files = ["assembly.fasta", "assembly.fa", "contigs.fasta", "contigs.fa"]
        
        for filename in priority_files:
            fasta_file = output_dir / filename
            if fasta_file.exists() and fasta_file.stat().st_size > 0:
                new_name = f"{self.name}.fa"
                new_path = Path.cwd() / new_name
                
                print(f"Copying {fasta_file} -> {new_path}")
                shutil.copy2(fasta_file, new_path)
                fasta_files.append(new_path)
                return fasta_files
        
        backup_patterns = ["*.fasta", "*.fa"]
        for pattern in backup_patterns:
            for fasta_file in output_dir.glob(pattern):
                if fasta_file.stat().st_size > 0:
                    new_name = f"{self.name}.fa"
                    new_path = Path.cwd() / new_name
                    
                    print(f"Copying {fasta_file} -> {new_path}")
                    shutil.copy2(fasta_file, new_path)
                    fasta_files.append(new_path)
                    return fasta_files
        
        print(f"No FASTA files found in {output_path}")
        return fasta_files


class AssemblyPipeline:
    
    def __init__(self, config: Optional[AssemblyConfig] = None):
        self.config = config or AssemblyConfig()
        self.tools = {
            'hifiasm': HifiasmTool(),
            'verkko': VerkkoTool(),
            'nextdenovo': NextDenovoTool(),
            'flye': FlyeTool()
        }
        for tool in self.tools.values():
            tool.config = self.config
        self.results = {}
        self.completed_tools = []
    
    def prepare_data(self) -> Dict[str, Any]:
        data = self.config.to_dict()
        
        for key in ['hifi_files', 'ont_ul_files', 'clr_files']:
            if data[key]:
                data[key.replace('_files', '')] = data[key]
                del data[key]
        
        return data
    
    def get_available_tools(self, data: Dict[str, Any]) -> List[str]:
        available = []
        for tool_name, tool in self.tools.items():
            if tool.validate_input(data):
                available.append(tool_name)
        return available
    
    def _parse_genome_size(self, size_str: str) -> float:
        size_str = size_str.lower().strip()
        try:
            if size_str.endswith('g'):
                return float(size_str[:-1]) * 1e9
            elif size_str.endswith('m'):
                return float(size_str[:-1]) * 1e6
            elif size_str.endswith('k'):
                return float(size_str[:-1]) * 1e3
            else:
                return float(size_str)
        except:
            return 130e6
    
    def run(self, selected_tools: Optional[List[str]] = None) -> Dict[str, Any]:
        print_title()
        
        data = self.prepare_data()
        
        if self.config.auto_estimate_genome:
            print("\n" + "="*60)
            print("STEP 1: Fast genome size estimation")
            print("="*60)
            
            all_reads = []
            has_hifi = False
            read_type = "hifi"
            
            if data.get('hifi'):
                all_reads.extend(data['hifi'])
                has_hifi = True
                read_type = "hifi"
            
            if data.get('ont_ul'):
                all_reads.extend(data['ont_ul'])
                if not has_hifi:
                    read_type = "ont"
            
            if data.get('clr'):
                all_reads.extend(data['clr'])
                if not has_hifi and not data.get('ont_ul'):
                    read_type = "clr"
            
            if not all_reads:
                print("Error: No input files found for estimation")
                return {}
            
            print(f"Input data type: {read_type.upper()}")
            print(f"Target depth: {self.config.target_depth}X")
            
            estimated_size, method = FastGenomeSizeEstimator.estimate_genome_size(
                reads_files=all_reads,
                read_type=read_type,
                target_depth=self.config.target_depth
            )
            
            if estimated_size:
                self.config.estimated_genome_size = estimated_size
                self.config.estimated_genome_size_bp = self._parse_genome_size(estimated_size)
            else:
                print("Warning: Failed to estimate genome size, using default 1g")
                self.config.estimated_genome_size = "1g"
                self.config.estimated_genome_size_bp = 1e9
        
        available_tools = self.get_available_tools(data)
        
        if not available_tools:
            print("No assembly tools suitable for current data types")
            return {}
        
        if selected_tools is None:
            if 'all' in self.config.tools_to_run:
                selected_tools = available_tools
            else:
                selected_tools = [t for t in self.config.tools_to_run if t in available_tools]
        else:
            selected_tools = [t for t in selected_tools if t in available_tools]
        
        if not selected_tools:
            print("No available tools selected")
            return {}
        
        print(f"\n{'='*60}")
        print("STEP 2: Running selected assembly tools...")
        print("="*60)
        print(f"Will run the following tools: {', '.join(selected_tools)}")
        print(f"Output prefix: {self.config.output_prefix}")
        print(f"Thread count: {self.config.threads}")
        if self.config.estimated_genome_size:
            print(f"Estimated genome size: {self.config.estimated_genome_size}")
        
        self.results = {}
        self.completed_tools = []
        
        for tool_name in selected_tools:
            tool = self.tools[tool_name]
            print(f"\n{'='*60}")
            print(f"Running {tool_name.upper()} ...")
            print("="*60)
            
            result = tool.run(data, self.config)
            if result:
                self.completed_tools.append((tool_name, tool_name))
                self.results[tool_name] = {
                    'files': result,
                    'status': 'success'
                }
                print(f"✓ {tool_name} completed successfully")
            else:
                self.results[tool_name] = {
                    'files': [],
                    'status': 'failed'
                }
                print(f"✗ {tool_name} failed")
        
        if self.completed_tools:
            print(f"\n{'='*60}")
            print("STEP 3: Processing assembly result files...")
            print("="*60)
            
            all_fasta_files = []
            for tool_name, _ in self.completed_tools:
                if tool_name in self.results and self.results[tool_name]['files']:
                    all_fasta_files.extend(self.results[tool_name]['files'])
            
            if all_fasta_files:
                merge_output_name = f"{self.config.output_prefix}_all_assemblies.fasta"
                merged_result = merge_fasta_files(all_fasta_files, merge_output_name)
                
                if merged_result:
                    self.results['merged'] = {
                        'file': merged_result,
                        'status': 'success'
                    }
        
        return {
            'completed_tools': self.completed_tools,
            'results': self.results,
            'config': self.config.to_dict()
        }
    
    def run_tool(self, tool_name: str) -> Optional[List[Path]]:
        if tool_name not in self.tools:
            print(f"Unknown tool: {tool_name}")
            return None
        
        data = self.prepare_data()
        tool = self.tools[tool_name]
        
        if not tool.validate_input(data):
            print(f"{tool_name} cannot handle current data type")
            return None
        
        print(f"\n{'='*60}")
        print(f"Running {tool_name.upper()} ...")
        print("="*60)
        
        result = tool.run(data, self.config)
        if result:
            self.results[tool_name] = {
                'files': result,
                'status': 'success'
            }
            print(f"✓ {tool_name} completed successfully")
        else:
            self.results[tool_name] = {
                'files': [],
                'status': 'failed'
            }
            print(f"✗ {tool_name} failed")
        
        return result
    
    def get_summary(self) -> str:
        summary = []
        summary.append("="*60)
        summary.append("Assembly Pipeline Run Summary")
        summary.append("="*60)
        
        if not self.results:
            summary.append("No tools have been run yet")
            return "\n".join(summary)
        
        summary.append(f"Output prefix: {self.config.output_prefix}")
        summary.append(f"Thread count: {self.config.threads}")
        summary.append(f"Automatic genome estimation: {'ON' if self.config.auto_estimate_genome else 'OFF'}")
        if self.config.estimated_genome_size:
            summary.append(f"Estimated genome size: {self.config.estimated_genome_size} "
                          f"({self.config.estimated_genome_size_bp/1e9:.2f} Gb)")
        summary.append("")
        
        summary.append("Input data:")
        for data_type, files in [('HiFi', self.config.hifi_files),
                               ('ONT UL', self.config.ont_ul_files),
                               ('CLR', self.config.clr_files)]:
            if files:
                summary.append(f"  - {data_type}: {len(files)} files")
        
        summary.append("\nRun results:")
        for tool_name, result in self.results.items():
            if tool_name != 'merged':
                status_text = "✓ SUCCESS" if result['status'] == 'success' else "✗ FAILED"
                file_count = len(result['files'])
                summary.append(f"  - {tool_name}: {status_text} ({file_count} files)")
        
        if 'merged' in self.results:
            merged_file = self.results['merged']['file']
            if merged_file and merged_file.exists():
                size_mb = merged_file.stat().st_size / (1024 * 1024)
                size_gb = size_mb / 1024
                summary.append(f"\nMerged assembly file: {merged_file.name}")
                summary.append(f"  Size: {size_gb:.2f} GB ({size_mb:.1f} MB)")
                
                with open(merged_file, 'r') as f:
                    seq_count = sum(1 for line in f if line.startswith('>'))
                summary.append(f"  Sequences: {seq_count}")
        
        return "\n".join(summary)


def collect_input_interactive() -> Dict[str, Any]:
    print("Step 1: Please provide your sequencing data paths\n")
    
    hifi_input = ask("PacBio HiFi reads (FASTQ/GZ):", required=False)
    hifi_files = validate_files(hifi_input) if hifi_input else None
    
    ont_input = ask("Oxford Nanopore Ultra-Long reads (FASTQ/GZ):", required=False)
    ont_files = validate_files(ont_input) if ont_input else None
    
    clr_input = ask("PacBio CLR reads (FASTQ/GZ):", required=False)
    clr_files = validate_files(clr_input) if clr_input else None
    
    if not (hifi_files or ont_files or clr_files):
        print("At least one type of sequencing data is required (HiFi / ONT / CLR).")
        sys.exit(1)
    
    threads = int(ask("Thread count:", default="32"))
    output_prefix = ask("Output prefix/project name:", default="assembly")
    
    return {
        'hifi_files': hifi_files,
        'ont_ul_files': ont_files,
        'clr_files': clr_files,
        'threads': threads,
        'output_prefix': output_prefix
    }


def collect_tool_config_interactive() -> Dict[str, Any]:
    print("\nStep 2: Tool-specific configuration (press Enter for defaults)")
    
    print("\nGenome size automatic estimation configuration:")
    auto_estimate = ask("Automatically estimate genome size?", options=['y','n'], default='y')
    auto_estimate_genome = (auto_estimate == 'y')
    
    target_depth = 30
    if auto_estimate_genome:
        target_depth = int(ask("Target sequencing depth (X):", default="30"))
    
    print("\nHifiasm configuration:")
    print("  Note: --ont for ONT simplex, --ul for HiFi+ONT hybrid")
    hifiasm_primary = ask("Run hifiasm in primary mode (single assembly, not haplotype-resolved)?", 
                          options=['y','n'], default='n') == 'y'
    hifiasm_merge_haplotypes = ask("Merge hap1 and hap2 into single FASTA file?", 
                                    options=['y','n'], default='y') == 'y'
    
    verkko_memory_gb = int(ask("verkko: Maximum memory (GB):", default="64"))
    nextdenovo_genome_size = ask("nextDenovo: Genome size (e.g., 3g, 1.5m):", default="1g")
    flye_genome_size = ask("Flye: Genome size (e.g., 3g):", default="1g")
    flye_iterations = int(ask("Flye: Polishing iterations:", default="1"))
    flye_nano_type = ask("Flye: ONT data type?", options=['raw', 'hq'], default='hq')
    
    return {
        'auto_estimate_genome': auto_estimate_genome,
        'target_depth': target_depth,
        'hifiasm_primary': hifiasm_primary,
        'hifiasm_merge_haplotypes': hifiasm_merge_haplotypes,
        'verkko_memory_gb': verkko_memory_gb,
        'nextdenovo_genome_size': nextdenovo_genome_size,
        'flye_genome_size': flye_genome_size,
        'flye_iterations': flye_iterations,
        'flye_nano_type': flye_nano_type
    }


def interactive_mode():
    print_title()
    
    input_data = collect_input_interactive()
    tool_config = collect_tool_config_interactive()
    
    print("\nStep 3: Select assembly tools to run")
    
    temp_config = AssemblyConfig(**input_data, **tool_config)
    pipeline = AssemblyPipeline(temp_config)
    data = pipeline.prepare_data()
    available_tools = pipeline.get_available_tools(data)
    
    if not available_tools:
        print("No assembly tools suitable for current data types")
        sys.exit(1)
    
    print(f"Available tools: {', '.join(available_tools)}")
    
    choice = ask(
        "\nSelect assembly tools to run:",
        options=['all', 'custom'] + available_tools,
        default='all'
    )
    
    tools_to_run = []
    if choice == 'all':
        tools_to_run = available_tools
    elif choice == 'custom':
        for tool in available_tools:
            if ask(f"Run {tool}?", options=['y','n'], default='y') == 'y':
                tools_to_run.append(tool)
    else:
        tools_to_run = [choice]
    
    if not tools_to_run:
        print("No tools selected.")
        sys.exit(1)
    
    config = AssemblyConfig(
        **input_data,
        **tool_config,
        tools_to_run=tools_to_run
    )
    
    pipeline = AssemblyPipeline(config)
    
    confirm = ask("\nRun now?", options=['y','n'], default='n')
    if confirm != 'y':
        print("Cancelled.")
        sys.exit(0)
    
    results = pipeline.run()
    
    print(pipeline.get_summary())


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Genome Assembly Tool Control Center (HiFi/ONT/CLR) - FAST estimation (NO k-mer)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
 %(prog)s -hifi reads.fastq.gz -t 64 -o my_assembly --run-hifiasm --run-flye
 %(prog)s -hifi hifi.fastq.gz -t 64 --run-all
 %(prog)s -hifi hifi.fastq.gz -t 64 --run-nextdenovo --no-auto-estimate --nextdenovo-genome-size 100m
 %(prog)s -hifi hifi.fastq.gz -t 64 --run-hifiasm --hifiasm-primary --no-hifiasm-merge
 %(prog)s -ont ont.fq.gz -t 64 --run-hifiasm --run-flye

NOTE: Genome size estimation uses ONLY file size + sampling methods (NO k-mer)
      This makes it extremely fast (10-30 seconds for 80GB ONT data)
      Hifiasm ONT support: --ont for simplex reads, --ul for ultra-long+HiFi hybrid
        """
    )
    
    parser.add_argument('-hifi', '--hifi-files', nargs='+', help='HiFi reads (FASTQ/GZ format)')
    parser.add_argument('-ont', '--ont-ul-files', nargs='+', help='ONT Ultra-Long reads (FASTQ/GZ format)')
    parser.add_argument('-clr', '--clr-files', nargs='+', help='PacBio CLR reads (FASTQ/GZ format)')
    
    parser.add_argument('-t', '--threads', type=int, default=32, help='Default thread count')
    parser.add_argument('-o', '--output', default='assembly', help='Output prefix/project name')
    
    parser.add_argument('--auto-estimate', action='store_true', default=True, 
                       help='Automatically estimate genome size (default ON)')
    parser.add_argument('--no-auto-estimate', action='store_false', dest='auto_estimate',
                       help='Disable automatic genome size estimation')
    parser.add_argument('--target-depth', type=int, default=30, 
                       help='Target sequencing depth (for estimation)')
    
    parser.add_argument('--run-all', action='store_true', help='Run all supported assembly software')
    parser.add_argument('--run-hifiasm', action='store_true', help='Run hifiasm')
    parser.add_argument('--run-verkko', action='store_true', help='Run verkko')
    parser.add_argument('--run-nextdenovo', action='store_true', help='Run nextDenovo')
    parser.add_argument('--run-flye', action='store_true', help='Run flye')
    
    parser.add_argument('--hifiasm-threads', type=int, help='hifiasm-specific thread count')
    parser.add_argument('--verkko-threads', type=int, help='verkko-specific thread count')
    parser.add_argument('--nextdenovo-threads', type=int, help='nextDenovo-specific thread count')
    parser.add_argument('--flye-threads', type=int, help='flye-specific thread count')
    
    parser.add_argument('--hifiasm-primary', action='store_true', 
                       help='Run hifiasm in primary mode (single assembly, not haplotype-resolved)')
    parser.add_argument('--no-hifiasm-merge', action='store_false', dest='hifiasm_merge_haplotypes',
                       help='Do not merge hap1 and hap2 into single FASTA file')
    
    parser.add_argument('--verkko-memory', type=int, default=64, help='verkko maximum memory (GB)')
    parser.add_argument('--verkko-no-cleanup', action='store_false', dest='verkko_cleanup',
                       help='Keep intermediate files for verkko')
    parser.add_argument('--flye-genome-size', default='1g', help='Flye genome size estimate')
    parser.add_argument('--nextdenovo-genome-size', default='1g', help='nextDenovo genome size estimate')
    parser.add_argument('--flye-iterations', type=int, default=1, help='Flye polishing iterations')
    parser.add_argument('--flye-nano-type', choices=['raw', 'hq'], default='hq', help='Flye ONT data type')
    
    parser.add_argument('-i', '--interactive', action='store_true', help='Interactive mode')
    parser.add_argument('--config', help='JSON configuration file path')
    
    parser.set_defaults(hifiasm_merge_haplotypes=True)
    
    return parser.parse_args()


def cli_mode(args):
    tools_to_run = []
    if args.run_all:
        tools_to_run = ['all']
    else:
        tool_mapping = {
            'run_hifiasm': 'hifiasm',
            'run_verkko': 'verkko',
            'run_nextdenovo': 'nextdenovo',
            'run_flye': 'flye'
        }
        
        for arg_name, tool_name in tool_mapping.items():
            if getattr(args, arg_name, False):
                tools_to_run.append(tool_name)
    
    if not tools_to_run:
        print("No assembly tools selected.")
        print("   Use --run-all or specify specific software switches (e.g., --run-hifiasm)")
        print("   Use --help to see all available options")
        sys.exit(1)
    
    config = AssemblyConfig(
        hifi_files=args.hifi_files,
        ont_ul_files=args.ont_ul_files,
        clr_files=args.clr_files,
        threads=args.threads,
        output_prefix=args.output,
        auto_estimate_genome=args.auto_estimate,
        target_depth=args.target_depth,
        hifiasm_primary=args.hifiasm_primary,
        hifiasm_merge_haplotypes=args.hifiasm_merge_haplotypes,
        verkko_memory_gb=args.verkko_memory,
        verkko_cleanup=args.verkko_cleanup,
        nextdenovo_genome_size=args.nextdenovo_genome_size,
        flye_genome_size=args.flye_genome_size,
        flye_iterations=args.flye_iterations,
        flye_nano_type=args.flye_nano_type,
        hifiasm_threads=args.hifiasm_threads,
        verkko_threads=args.verkko_threads,
        nextdenovo_threads=args.nextdenovo_threads,
        flye_threads=args.flye_threads,
        tools_to_run=tools_to_run
    )
    
    for file_list in [config.hifi_files, config.ont_ul_files, config.clr_files]:
        if file_list:
            for f in file_list:
                if not Path(f).exists():
                    print(f"File not found: {f}")
                    sys.exit(1)
    
    pipeline = AssemblyPipeline(config)
    results = pipeline.run()
    
    print(pipeline.get_summary())


def main():
    args = parse_arguments()
    
    if args.config:
        pipeline = AssemblyPipeline(AssemblyConfig.from_json(args.config))
        results = pipeline.run()
        print(pipeline.get_summary())
    elif len(sys.argv) > 1 and not args.interactive:
        cli_mode(args)
    else:
        interactive_mode()


def run_pipeline_from_dict(config_dict: Dict[str, Any]) -> Dict[str, Any]:
    try:
        config = AssemblyConfig.from_dict(config_dict)
        pipeline = AssemblyPipeline(config)
        results = pipeline.run()
        return results
    except Exception as e:
        print(f"Pipeline execution failed: {e}")
        return {}


def run_pipeline_from_json(json_file: str) -> Dict[str, Any]:
    try:
        config = AssemblyConfig.from_json(json_file)
        pipeline = AssemblyPipeline(config)
        results = pipeline.run()
        return results
    except Exception as e:
        print(f"Pipeline execution failed: {e}")
        return {}


def run_specific_tools(config_dict: Dict[str, Any], tools: List[str]) -> Dict[str, Any]:
    try:
        config = AssemblyConfig.from_dict(config_dict)
        config.tools_to_run = tools
        pipeline = AssemblyPipeline(config)
        results = pipeline.run()
        return results
    except Exception as e:
        print(f"Running specific tools failed: {e}")
        return {}


if __name__ == "__main__":
    main()
