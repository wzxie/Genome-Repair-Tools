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


@dataclass
class AssemblyConfig:
    hifi_files: Optional[List[str]] = None
    ont_ul_files: Optional[List[str]] = None
    clr_files: Optional[List[str]] = None
    
    threads: int = 32
    output_prefix: str = "assembly"
    
    hifiasm_primary: bool = True
    hifiasm_l: int = 3
    verkko_memory_gb: int = 64
    verkko_cleanup: bool = True
    nextdenovo_genome_size: str = "1g"
    nextdenovo_read_cutoff: str = "1k"
    flye_genome_size: str = "1g"
    flye_iterations: int = 1
    flye_nano_type: str = "hq"
    
    auto_estimate_genome: bool = True
    estimate_method: str = "kmer"
    target_depth: int = 30
    kmer_size: int = 21
    
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


class ONTGenomeSizeEstimator:
    
    @staticmethod
    def get_read_stats(file_path: Path) -> Dict[str, Any]:
        stats = {
            'total_bases': 0,
            'total_reads': 0,
            'long_reads_count': 0,
            'long_reads_bases': 0,
            'n50': 0,
            'mean_length': 0
        }
        
        try:
            lengths = []
            
            if str(file_path).endswith('.gz'):
                cmd = f"zcat {file_path} | awk 'NR%4==2 {{print length($0)}}' | head -100000"
            else:
                cmd = f"awk 'NR%4==2 {{print length($0)}}' {file_path} | head -100000"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        length = int(line.strip())
                        lengths.append(length)
                        stats['total_bases'] += length
                        stats['total_reads'] += 1
                        
                        if length >= 10000:
                            stats['long_reads_count'] += 1
                            stats['long_reads_bases'] += length
            
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
        
        except Exception as e:
            print(f"Warning: Failed to get read stats for {file_path}: {e}")
        
        return stats
    
    @staticmethod
    def estimate_from_long_reads_only(reads_files: List[str], 
                                      target_depth: int = 30,
                                      min_length: int = 10000) -> Optional[str]:
        total_long_bases = 0
        total_reads = 0
        long_reads_ratio = 0
        
        for file_path in reads_files:
            stats = ONTGenomeSizeEstimator.get_read_stats(Path(file_path))
            total_long_bases += stats['long_reads_bases']
            total_reads += stats['total_reads']
            
            if stats['total_bases'] > 0:
                file_long_ratio = stats['long_reads_bases'] / stats['total_bases']
                long_reads_ratio = max(long_reads_ratio, file_long_ratio)
        
        if total_long_bases == 0:
            print("   Warning: No long reads found (>=10kb)")
            return None
        
        if long_reads_ratio == 0:
            long_reads_ratio = 0.3
        
        estimated_total_bases = total_long_bases / long_reads_ratio
        estimated_size_bp = estimated_total_bases / target_depth
        
        if estimated_size_bp >= 1e9:
            return f"{estimated_size_bp/1e9:.1f}g"
        elif estimated_size_bp >= 1e6:
            return f"{int(estimated_size_bp/1e6)}m"
        else:
            return f"{int(estimated_size_bp/1000)}k"
    
    @staticmethod
    def estimate_ont_with_meryl(reads_files: List[str], 
                                k: int = 31, 
                                threads: int = 32) -> Optional[str]:
        try:
            temp_dir = Path(tempfile.mkdtemp(prefix="ont_kmer_est_"))
            
            filtered_reads = temp_dir / "filtered_reads.fastq"
            
            print(f"  Extracting long reads for k-mer analysis...")
            for file_path in reads_files:
                if str(file_path).endswith('.gz'):
                    cmd = f"zcat {file_path} | awk 'BEGIN{{RS=\"@\";FS=\"\\n\"}} length($2)>=10000 {{print \"@\"$0}}' >> {filtered_reads}"
                else:
                    cmd = f"awk 'BEGIN{{RS=\"@\";FS=\"\\n\"}} length($2)>=10000 {{print \"@\"$0}}' {file_path} >> {filtered_reads}"
                
                subprocess.run(cmd, shell=True, check=True, capture_output=True)
            
            if not filtered_reads.exists() or filtered_reads.stat().st_size == 0:
                print("  No long reads available for k-mer analysis")
                return None
            
            meryl_db = temp_dir / "meryl_db"
            
            cmd = [
                "meryl", f"k={k}", "count",
                f"memory={threads*2}",
                f"threads={threads}",
                "output", str(meryl_db),
                str(filtered_reads)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  meryl count failed: {result.stderr[:200]}")
                return None
            
            kmer_hist = temp_dir / "kmer.hist"
            cmd = f"meryl histogram {meryl_db} > {kmer_hist}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0 or not kmer_hist.exists():
                print(f"  meryl histogram failed")
                return None
            
            with open(kmer_hist, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 10:
                print("  Insufficient k-mer spectrum data")
                return None
            
            total_kmers = 0
            kmer_counts = []
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        coverage = int(parts[0])
                        count = int(parts[1])
                        if coverage > 5:
                            total_kmers += count
                            kmer_counts.append((coverage, count))
                    except:
                        continue
            
            if not kmer_counts:
                return None
            
            kmer_counts.sort(key=lambda x: x[1], reverse=True)
            main_coverage = kmer_counts[0][0]
            
            genome_size_bp = total_kmers / main_coverage if main_coverage > 0 else 0
            
            if genome_size_bp <= 0:
                return None
            
            if genome_size_bp >= 1e9:
                return f"{genome_size_bp/1e9:.1f}g"
            elif genome_size_bp >= 1e6:
                return f"{int(genome_size_bp/1e6)}m"
            else:
                return f"{int(genome_size_bp/1000)}k"
            
        except Exception as e:
            print(f"  ONT k-mer analysis error: {e}")
            return None
        finally:
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
    
    @staticmethod
    def cross_validate(results: List[str]) -> str:
        if not results:
            return "1g"
        
        numeric_sizes = []
        for size_str in results:
            try:
                size_str = size_str.lower()
                if size_str.endswith('g'):
                    numeric_sizes.append(float(size_str[:-1]) * 1e9)
                elif size_str.endswith('m'):
                    numeric_sizes.append(float(size_str[:-1]) * 1e6)
                elif size_str.endswith('k'):
                    numeric_sizes.append(float(size_str[:-1]) * 1e3)
                else:
                    numeric_sizes.append(float(size_str))
            except:
                continue
        
        if not numeric_sizes:
            return results[0]
        
        numeric_sizes.sort()
        median_size = numeric_sizes[len(numeric_sizes)//2]
        
        if median_size >= 1e9:
            return f"{median_size/1e9:.1f}g"
        elif median_size >= 1e6:
            return f"{int(median_size/1e6)}m"
        else:
            return f"{int(median_size/1000)}k"


class GenomeSizeEstimator:
    
    @staticmethod
    def check_tool_availability(tool_name: str) -> bool:
        try:
            subprocess.run([tool_name, "--version" if tool_name != "meryl" else "k=1"], 
                          capture_output=True, check=True)
            return True
        except:
            try:
                subprocess.run([tool_name, "--help"], capture_output=True, check=True)
                return True
            except:
                return False
    
    @staticmethod
    def estimate_from_file_size(reads_files: List[str], target_depth: int = 30) -> str:
        total_bases = 0
        for file_path in reads_files:
            file_path = Path(file_path)
            if file_path.suffix == '.gz':
                with gzip.open(file_path, 'rb') as f:
                    compressed_size = file_path.stat().st_size
                    estimated_bases = compressed_size * 3 * 0.75
                    total_bases += int(estimated_bases)
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
                            file_size = file_path.stat().st_size
                            est_lines = (file_size / (avg_length * 4 + 100)) * 4
                            total_bases += int(est_lines / 4 * avg_length)
                except:
                    file_size = file_path.stat().st_size
                    total_bases += int(file_size * 0.75)
        
        if total_bases == 0:
            return "1g"
        
        estimated_size_bp = total_bases / target_depth
        
        if estimated_size_bp >= 1e9:
            return f"{estimated_size_bp/1e9:.1f}g"
        elif estimated_size_bp >= 1e6:
            return f"{int(estimated_size_bp/1e6)}m"
        else:
            return f"{int(estimated_size_bp/1000)}k"
    
    @staticmethod
    def estimate_using_meryl(reads_files: List[str], k: int = 21, 
                           threads: int = 32, temp_dir: Optional[Path] = None) -> Optional[str]:
        try:
            if temp_dir is None:
                temp_dir = Path(tempfile.mkdtemp(prefix="kmer_est_"))
            else:
                temp_dir.mkdir(parents=True, exist_ok=True)
            
            meryl_db = temp_dir / "meryl_db"
            kmer_hist = temp_dir / "kmer.hist"
            
            print(f"Running meryl for k-mer analysis (k={k})...")
            
            cmd = [
                "meryl", f"k={k}", "count",
                f"memory={threads*2}",
                f"threads={threads}",
                "output", str(meryl_db)
            ] + reads_files
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"meryl count failed: {result.stderr[:200]}")
                return None
            
            cmd = [
                "meryl", "histogram",
                str(meryl_db),
                ">", str(kmer_hist)
            ]
            
            result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"meryl histogram failed: {result.stderr[:200]}")
                return None
            
            if not kmer_hist.exists() or kmer_hist.stat().st_size == 0:
                print("k-mer spectrum file is empty")
                return None
            
            with open(kmer_hist, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 10:
                print("Insufficient k-mer spectrum data")
                return None
            
            total_kmers = 0
            kmer_counts = []
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        coverage = int(parts[0])
                        count = int(parts[1])
                        total_kmers += count
                        kmer_counts.append((coverage, count))
                    except:
                        continue
            
            if total_kmers == 0:
                print("No valid k-mer data found")
                return None
            
            kmer_counts.sort(key=lambda x: x[1], reverse=True)
            if len(kmer_counts) > 0:
                main_coverage = kmer_counts[0][0]
                
                genome_size_bp = total_kmers / main_coverage if main_coverage > 0 else 0
                genome_size_bp = genome_size_bp + (k - 1)
                
                if genome_size_bp <= 0:
                    return None
                
                if genome_size_bp >= 1e9:
                    return f"{genome_size_bp/1e9:.1f}g"
                elif genome_size_bp >= 1e6:
                    return f"{int(genome_size_bp/1e6)}m"
                else:
                    return f"{int(genome_size_bp/1000)}k"
            
            return None
            
        except Exception as e:
            print(f"meryl analysis error: {e}")
            return None
        finally:
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
    
    @staticmethod
    def estimate_using_bbtools(reads_files: List[str], k: int = 31,
                             threads: int = 32, temp_dir: Optional[Path] = None) -> Optional[str]:
        try:
            try:
                subprocess.run(["kmercountexact.sh", "--version"], 
                             capture_output=True, check=True)
            except:
                print("BBTools not installed, skipping k-mer analysis")
                return None
            
            if temp_dir is None:
                temp_dir = Path(tempfile.mkdtemp(prefix="bbtools_kmer_"))
            else:
                temp_dir.mkdir(parents=True, exist_ok=True)
            
            output_prefix = temp_dir / "kmer_counts"
            
            print(f"Running BBTools for k-mer analysis (k={k})...")
            
            cmd = [
                "kmercountexact.sh", f"k={k}",
                f"threads={threads}",
                f"-Xmx{threads*2}g",
                f"in={' '.join(reads_files)}",
                f"out={output_prefix}"
            ]
            
            result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                output = result.stdout + result.stderr
                
                patterns = [
                    r"Genome size:\s*(\d+)\s*bp",
                    r"Unique Kmers:\s*(\d+)",
                    r"Estimated genome size:\s*(\d+)\s*bp"
                ]
                
                for pattern in patterns:
                    match = re.search(pattern, output, re.IGNORECASE)
                    if match:
                        genome_size_bp = int(match.group(1))
                        
                        if genome_size_bp >= 1e9:
                            return f"{genome_size_bp/1e9:.1f}g"
                        elif genome_size_bp >= 1e6:
                            return f"{int(genome_size_bp/1e6)}m"
                        else:
                            return f"{int(genome_size_bp/1000)}k"
            
            print("Could not extract genome size from BBTools output")
            return None
            
        except Exception as e:
            print(f"BBTools analysis error: {e}")
            return None
        finally:
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
    
    @staticmethod
    def hybrid_estimation(reads_files: List[str], method: str = "hybrid",
                         k: int = 21, threads: int = 32, 
                         target_depth: int = 30,
                         read_type: str = "hifi") -> str:
        estimated_size = None
        
        if read_type == "ont":
            print("ONT data: Using specialized estimation methods...")
            
            print("  Strategy 1: File size based estimation with long reads...")
            estimated_size = ONTGenomeSizeEstimator.estimate_from_long_reads_only(
                reads_files, target_depth=target_depth, min_length=10000
            )
            if estimated_size:
                print(f"  Strategy 1 result: {estimated_size}")
            
            if not estimated_size and method in ["kmer", "hybrid"]:
                print("  Strategy 2: Large k-mer analysis (k=31)...")
                estimated_size = ONTGenomeSizeEstimator.estimate_ont_with_meryl(
                    reads_files, k=31, threads=threads
                )
                if estimated_size:
                    print(f"  Strategy 2 result: {estimated_size}")
            
            if not estimated_size:
                print("  Strategy 3: Quick file size estimation...")
                estimated_size = GenomeSizeEstimator.estimate_from_file_size(
                    reads_files, target_depth
                )
                print(f"  Strategy 3 result: {estimated_size}")
        
        else:
            if method == "kmer" or method == "hybrid":
                print("Trying meryl for k-mer analysis...")
                if GenomeSizeEstimator.check_tool_availability("meryl"):
                    estimated_size = GenomeSizeEstimator.estimate_using_meryl(
                        reads_files, k, threads
                    )
                
                if not estimated_size:
                    print("meryl unavailable or failed, trying BBTools...")
                    estimated_size = GenomeSizeEstimator.estimate_using_bbtools(
                        reads_files, k, threads
                    )
            
            if not estimated_size or method == "quick":
                print("Using quick file size estimation...")
                estimated_size = GenomeSizeEstimator.estimate_from_file_size(
                    reads_files, target_depth
                )
        
        return estimated_size or "1g"


def print_title():
    print("\n" + "="*70)
    print("      Genome Assembly Tool Control Center (HiFi / ONT / CLR)")
    print("      hifiasm | verkko | nextDenovo | flye")
    print("="*70)
    print("Supports: HiFi, ONT Ultra-Long, PacBio CLR")
    print("Integrated k-mer analysis for automatic genome size estimation")
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


def decompress_gz_files(input_files: List[str], temp_dir: Path, 
                       threads: int) -> List[str]:
    new_files = []
    
    try:
        subprocess.run(["pigz", "--version"], capture_output=True, check=True)
        use_pigz = True
        print("pigz detected, using multi-thread decompression")
    except:
        use_pigz = False
        print("pigz not detected, using single-thread gunzip")
    
    for f in input_files:
        f_path = Path(f)
        if f.endswith('.gz'):
            base_name = f_path.stem
            target = temp_dir / base_name
            print(f"Decompressing {f} -> {target}")
            
            if use_pigz:
                pigz_threads = min(threads, 8)
                try:
                    subprocess.run(
                        ["pigz", "-d", "-c", "-p", str(pigz_threads), f],
                        stdout=open(target, 'wb'), check=True
                    )
                except Exception as e:
                    print(f"pigz decompression failed, trying gunzip: {e}")
                    with open(target, 'wb') as out_f:
                        subprocess.run(["gunzip", "-c", f], stdout=out_f, check=True)
            else:
                with open(target, 'wb') as out_f:
                    subprocess.run(["gunzip", "-c", f], stdout=out_f, check=True)
                    
            new_files.append(str(target))
        else:
            new_files.append(f)
    return new_files


def gfa_to_fasta_awk(gfa_file: Path, fasta_file: Path) -> bool:
    print(f"Converting {gfa_file} -> {fasta_file}")
    
    try:
        cmd = f"awk '/^S/{{print \">\"$2; print $3}}' {gfa_file} > {fasta_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            if os.path.getsize(fasta_file) > 0:
                print(f"Conversion complete: {fasta_file}")
                return True
            else:
                print(f"Converted file is empty: {fasta_file}")
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
            write_sequence = True
            
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    if current_header and write_sequence and current_seq:
                        seq_str = ''.join(current_seq)
                        if seq_str and seq_str != '*' and not all(c in '*Nn' for c in seq_str):
                            outfile.write(f"{current_header}\n")
                            for i in range(0, len(seq_str), 80):
                                outfile.write(seq_str[i:i+80] + "\n")
                    
                    current_header = line
                    current_seq = []
                    write_sequence = True
                elif write_sequence:
                    current_seq.append(line)
            
            if current_header and write_sequence and current_seq:
                seq_str = ''.join(current_seq)
                if seq_str and seq_str != '*' and not all(c in '*Nn' for c in seq_str):
                    outfile.write(f"{current_header}\n")
                    for i in range(0, len(seq_str), 80):
                        outfile.write(seq_str[i:i+80] + "\n")
        
        shutil.move(temp_file, fasta_file)
        return True
        
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
        super().__init__("hifiasm", ['hifi', 'clr'])
        self.description = "PacBio HiFi/CLR assembly tool"
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        prefix = f"{data['output_prefix']}.hifiasm"
        threads = self.get_threads(data)
        
        cmd = ["hifiasm", "-o", prefix, "-t", str(threads)]
        
        if data.get('hifi'):
            cmd += data['hifi']
        elif data.get('clr'):
            if not data.get('ont_ul'):
                cmd += ["--ul"] + data['clr']
            else:
                cmd += ["--pacbio"] + data['clr']
        else:
            print("hifiasm requires HiFi or CLR input.")
            return [], ""
        
        if data.get('ont_ul'):
            cmd += ["--ul"] + data['ont_ul']
        
        cmd += ["--primary"]
        
        return cmd, prefix
    
    def process_output(self, output_path: str) -> List[Path]:
        fasta_files = []
        base_dir = Path.cwd()
        
        output_name = Path(output_path).name
        
        search_patterns = [
            f"{output_name}.p_ctg.gfa",
            f"{output_name}*.p_ctg.gfa",
            f"{output_name}.gfa",
            f"{output_name}*.gfa",
            "*.p_ctg.gfa",
            "*.gfa"
        ]
        
        found_gfa = None
        
        for pattern in search_patterns:
            try:
                matching_files = list(base_dir.glob(pattern))
                for gfa_file in matching_files:
                    if gfa_file.exists() and gfa_file.stat().st_size > 0:
                        found_gfa = gfa_file
                        print(f"Found GFA file in current directory: {gfa_file}")
                        break
                if found_gfa:
                    break
            except Exception as e:
                print(f"Pattern search error for {pattern}: {e}")
                continue
        
        if not found_gfa:
            output_path_obj = Path(output_path)
            if output_path_obj.exists() and output_path_obj.is_dir():
                for pattern in ["*.p_ctg.gfa", "*.gfa"]:
                    try:
                        matching_files = list(output_path_obj.glob(pattern))
                        for gfa_file in matching_files:
                            if gfa_file.exists() and gfa_file.stat().st_size > 0:
                                found_gfa = gfa_file
                                print(f"Found GFA file in output directory: {gfa_file}")
                                break
                        if found_gfa:
                            break
                    except Exception as e:
                        continue
            else:
                parent_dir = output_path_obj.parent
                if parent_dir.exists():
                    for pattern in [f"*{output_name}*.p_ctg.gfa", f"*{output_name}*.gfa"]:
                        try:
                            matching_files = list(parent_dir.glob(pattern))
                            for gfa_file in matching_files:
                                if gfa_file.exists() and gfa_file.stat().st_size > 0:
                                    found_gfa = gfa_file
                                    print(f"Found GFA file in parent directory: {gfa_file}")
                                    break
                            if found_gfa:
                                break
                        except Exception as e:
                            continue
        
        if not found_gfa:
            for pattern in ["**/*.p_ctg.gfa", "**/*.gfa"]:
                try:
                    matching_files = list(base_dir.glob(pattern))
                    for gfa_file in matching_files:
                        if gfa_file.exists() and gfa_file.stat().st_size > 0:
                            if "hifiasm" in gfa_file.name.lower() or "p_ctg" in gfa_file.name:
                                found_gfa = gfa_file
                                print(f"Found GFA file recursively: {gfa_file}")
                                break
                    if found_gfa:
                        break
                except Exception as e:
                    continue
        
        if found_gfa:
            fasta_name = f"{self.name}.fa"
            fasta_path = base_dir / fasta_name
            
            if gfa_to_fasta_awk(found_gfa, fasta_path):
                fasta_files.append(fasta_path)
                return fasta_files
        else:
            print(f"No hifiasm GFA files found")
            print(f"Current directory: {base_dir}")
            print(f"Output name: {output_name}")
            
            all_gfa_files = list(base_dir.glob("**/*.gfa"))
            if all_gfa_files:
                print(f"All GFA files found:")
                for gfa_file in all_gfa_files[:10]:
                    print(f"  {gfa_file}")
        
        return fasta_files


class VerkkoTool(AssemblyTool):
    
    def __init__(self):
        super().__init__("verkko", ['hifi'])
        self.description = "HiFi + ONT hybrid assembly tool"
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        out_dir = f"{data['output_prefix']}.verkko"
        Path(out_dir).parent.mkdir(parents=True, exist_ok=True)
        
        threads = self.get_threads(data)
        cmd = ["verkko", "-d", out_dir]
        
        if data.get('hifi'):
            cmd += ["--hifi"] + data['hifi']
        else:
            print("verkko requires HiFi data.")
            return [], ""
        
        if data.get('ont_ul'):
            cmd += ["--nano"] + data['ont_ul']
        else:
            cmd += ["--no-nano"]
        
        if config.verkko_cleanup:
            cmd += ["--cleanup"]
        
        cmd += ["--local-cpus", str(threads)]
        
        if config.estimated_genome_size_bp:
            genome_size_mb = config.estimated_genome_size_bp / 1e6
            
            base_memory_gb = 128
            base_genome_mb = 130
            
            recommended_memory = int(base_memory_gb * (genome_size_mb / base_genome_mb))
            
            min_memory = 64
            max_memory = 512
            recommended_memory = max(min_memory, min(recommended_memory, max_memory))
            
            mem_gb = recommended_memory
            print(f"verkko: Based on estimated genome size {genome_size_mb:.0f}Mb, using {mem_gb}GB memory")
        else:
            mem_gb = config.verkko_memory_gb
            print(f"verkko: Using configured memory: {mem_gb}GB")
        
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
        self.description = "ONT/HiFi/CLR assembly tool - automatic genome size estimation"
    
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
        need_estimation = config.auto_estimate_genome and not user_specified_size
        
        if need_estimation:
            print(f"nextDenovo: Using pre-estimated genome size...")
            if config.estimated_genome_size:
                config.nextdenovo_genome_size = config.estimated_genome_size
                print(f"nextDenovo using pre-estimated genome size: {config.nextdenovo_genome_size}")
            else:
                print(f"nextDenovo: Estimating genome size...")
                estimated_genome_size = GenomeSizeEstimator.hybrid_estimation(
                    reads_files=input_files,
                    method=config.estimate_method,
                    k=config.kmer_size if read_type != "ont" else 31,
                    threads=config.threads,
                    target_depth=config.target_depth,
                    read_type=read_type
                )
                config.nextdenovo_genome_size = estimated_genome_size
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
            print("Expected path: [output_prefix].nextDenovo/03.ctg_graph/nd.asm.fasta")
            
            alt_files = list(output_dir.glob("**/*.fasta")) + list(output_dir.glob("**/*.fa"))
            if alt_files:
                print(f"Found alternative FASTA files:")
                for f in alt_files[:3]:
                    print(f"  {f}")
        
        return fasta_files


class FlyeTool(AssemblyTool):
    
    def __init__(self):
        super().__init__("flye", ['hifi', 'ont_ul', 'clr'])
        self.description = "Long-read assembly tool for multiple data types - automatic genome size estimation"
    
    def build_command(self, data: Dict[str, Any], config: AssemblyConfig) -> Tuple[List[str], str]:
        out_dir = f"{data['output_prefix']}.flye"
        Path(out_dir).parent.mkdir(parents=True, exist_ok=True)
        
        threads = self.get_threads(data)
        
        original_genome_size = config.flye_genome_size
        
        user_specified_size = (config.flye_genome_size != "1g")
        need_estimation = config.auto_estimate_genome and not user_specified_size
        
        if need_estimation:
            print(f"Flye: Using pre-estimated genome size...")
            if config.estimated_genome_size:
                config.flye_genome_size = config.estimated_genome_size
                print(f"Flye using pre-estimated genome size: {config.flye_genome_size}")
            else:
                print(f"Flye: Estimating genome size...")
                
                input_files = []
                read_type = "hifi"
                if data.get('hifi'):
                    input_files.extend(data['hifi'])
                    read_type = "hifi"
                if data.get('ont_ul'):
                    input_files.extend(data['ont_ul'])
                    read_type = "ont"
                if data.get('clr'):
                    input_files.extend(data['clr'])
                    read_type = "clr"
                
                estimated_genome_size = GenomeSizeEstimator.hybrid_estimation(
                    reads_files=input_files,
                    method=config.estimate_method,
                    k=config.kmer_size if read_type != "ont" else 31,
                    threads=config.threads,
                    target_depth=config.target_depth,
                    read_type=read_type
                )
                
                print(f"Flye using estimated genome size: {estimated_genome_size}")
                config.flye_genome_size = estimated_genome_size
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
            print("STEP 1: Estimating genome size for all tools...")
            print("="*60)
            
            all_reads = []
            if data.get('hifi'):
                all_reads.extend(data['hifi'])
            if data.get('ont_ul'):
                all_reads.extend(data['ont_ul'])
            if data.get('clr'):
                all_reads.extend(data['clr'])
            
            if not all_reads:
                print("Error: No input files found for estimation")
                return {}
            
            read_type = "hifi"
            if data.get('ont_ul'):
                read_type = "ont"
            elif data.get('clr'):
                read_type = "clr"
            
            print(f"Input data type: {read_type.upper()}")
            print(f"Estimation method: {self.config.estimate_method}")
            print(f"Target depth: {self.config.target_depth}x")
            
            estimated_size = GenomeSizeEstimator.hybrid_estimation(
                reads_files=all_reads,
                method=self.config.estimate_method,
                k=self.config.kmer_size if read_type != "ont" else 31,
                threads=self.config.threads,
                target_depth=self.config.target_depth,
                read_type=read_type
            )
            
            self.config.estimated_genome_size = estimated_size
            self.config.estimated_genome_size_bp = self._parse_genome_size(estimated_size)
            
            print(f"\nEstimated genome size: {estimated_size}")
            print(f"   ({self.config.estimated_genome_size_bp/1e6:.1f} Mb)")
            
            if 'verkko' in (selected_tools or self.config.tools_to_run):
                genome_size_mb = self.config.estimated_genome_size_bp / 1e6
                
                base_memory_gb = 128
                base_genome_mb = 130
                
                recommended_memory = int(base_memory_gb * (genome_size_mb / base_genome_mb))
                
                min_memory = 64
                max_memory = 512
                recommended_memory = max(min_memory, min(recommended_memory, max_memory))
                
                print(f"\nBased on estimated genome size {genome_size_mb:.0f}Mb, recommended verkko memory: {recommended_memory}GB")
                
                if self.config.verkko_memory_gb == 64:
                    self.config.verkko_memory_gb = recommended_memory
                    print(f"Automatically setting verkko_memory_gb = {recommended_memory}")
                else:
                    print(f"Using user-specified verkko_memory_gb = {self.config.verkko_memory_gb}")
        
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
            
            result = tool.run(data, self.config)
            if result:
                self.completed_tools.append((tool_name, tool_name))
                self.results[tool_name] = {
                    'files': result,
                    'status': 'success'
                }
                print(f"{tool_name} completed")
            else:
                self.results[tool_name] = {
                    'files': [],
                    'status': 'failed'
                }
                print(f"{tool_name} failed")
        
        if self.completed_tools:
            print(f"\n{'='*60}")
            print("Processing assembly result files...")
            
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
        
        result = tool.run(data, self.config)
        if result:
            self.results[tool_name] = {
                'files': result,
                'status': 'success'
            }
            print(f"{tool_name} completed")
        else:
            self.results[tool_name] = {
                'files': [],
                'status': 'failed'
            }
            print(f"{tool_name} failed")
        
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
            summary.append(f"Estimated genome size: {self.config.estimated_genome_size}")
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
                status_text = "success" if result['status'] == 'success' else "failed"
                file_count = len(result['files'])
                summary.append(f"  - {tool_name}: {status_text} ({file_count} files)")
        
        if 'merged' in self.results:
            merged_file = self.results['merged']['file']
            if merged_file and merged_file.exists():
                size_mb = merged_file.stat().st_size / (1024 * 1024)
                summary.append(f"\nMerged file: {merged_file.name} ({size_mb:.1f} MB)")
        
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
    
    estimate_method = "quick"
    if auto_estimate_genome:
        estimate_method = ask("Estimation method?", options=['quick','kmer','hybrid'], default='hybrid')
    
    target_depth = 30
    if auto_estimate_genome:
        target_depth = int(ask("Target sequencing depth (X):", default="30"))
    
    verkko_memory_gb = int(ask("verkko: Maximum memory (GB) [Auto estimation will override this]:", default="64"))
    nextdenovo_genome_size = ask("nextDenovo: Estimated genome size (e.g., 3g, 1.5m) [Auto estimation will override]:", default="1g")
    flye_genome_size = ask("Flye: Estimated genome size (e.g., 3g) [Auto estimation will override]:", default="1g")
    flye_iterations = int(ask("Flye: Polishing iterations:", default="1"))
    flye_nano_type = ask("Flye: ONT data type?", options=['raw', 'hq'], default='hq')
    
    return {
        'auto_estimate_genome': auto_estimate_genome,
        'estimate_method': estimate_method,
        'target_depth': target_depth,
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
        description="Genome Assembly Tool Control Center (HiFi/ONT/CLR) - Integrated automatic genome size estimation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -hifi reads.fastq.gz -t 64 -o my_assembly --run-hifiasm --run-flye
  %(prog)s -hifi hifi.fastq.gz -t 64 --run-all --estimate-method kmer
  %(prog)s -hifi hifi.fastq.gz -t 64 --run-nextdenovo --no-auto-estimate --nextdenovo-genome-size 100m
  %(prog)s -hifi hifi.fastq.gz -t 64 --run-all --estimate-method hybrid
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
    parser.add_argument('--estimate-method', choices=['quick', 'kmer', 'hybrid'], 
                       default='hybrid', help='Genome estimation method')
    parser.add_argument('--target-depth', type=int, default=30, 
                       help='Target sequencing depth (for estimation)')
    parser.add_argument('--kmer-size', type=int, default=21, 
                       help='k-mer analysis size (for k-mer method)')
    
    parser.add_argument('--run-all', action='store_true', help='Run all supported assembly software')
    parser.add_argument('--run-hifiasm', action='store_true', help='Run hifiasm')
    parser.add_argument('--run-verkko', action='store_true', help='Run verkko')
    parser.add_argument('--run-nextdenovo', action='store_true', help='Run nextDenovo')
    parser.add_argument('--run-flye', action='store_true', help='Run flye')
    
    parser.add_argument('--hifiasm-threads', type=int, help='hifiasm-specific thread count')
    parser.add_argument('--verkko-threads', type=int, help='verkko-specific thread count')
    parser.add_argument('--nextdenovo-threads', type=int, help='nextDenovo-specific thread count')
    parser.add_argument('--flye-threads', type=int, help='flye-specific thread count')
    
    parser.add_argument('--verkko-memory', type=int, default=64, help='verkko maximum memory (GB) [Auto estimation will override]')
    parser.add_argument('--flye-genome-size', default='1g', help='Flye genome size estimate (overrides automatic estimation)')
    parser.add_argument('--nextdenovo-genome-size', default='1g', help='nextDenovo genome size estimate (overrides automatic estimation)')
    parser.add_argument('--flye-iterations', type=int, default=1, help='Flye polishing iterations')
    parser.add_argument('--flye-nano_type', choices=['raw', 'hq'], default='hq', help='Flye ONT data type')
    
    parser.add_argument('-i', '--interactive', action='store_true', help='Interactive mode')
    parser.add_argument('--config', help='JSON configuration file path')
    
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
        estimate_method=args.estimate_method,
        target_depth=args.target_depth,
        kmer_size=args.kmer_size,
        verkko_memory_gb=args.verkko_memory,
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
    
    if len(sys.argv) > 1 and not args.interactive:
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
