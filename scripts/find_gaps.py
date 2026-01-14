#!/usr/bin/env python3
"""
Genome Gap Finder Tool - Unified API Version
Find gap positions in FASTA files for each chromosome, support parallel processing
"""

import re
import sys
import os
import json
import logging
import multiprocessing
import time
import csv
from concurrent.futures import ProcessPoolExecutor
from typing import List, Dict, Tuple, Optional, Any, Union
from collections import defaultdict
from pathlib import Path

class GapFinderAPI:
    """Unified API for gap finding"""
    
    def __init__(self, verbose: bool = True, log_file: Optional[str] = None):
        """
        Initialize GapFinder API
        
        Args:
            verbose: Show detailed output
            log_file: Log file path (optional)
        """
        self.verbose = verbose
        self.log_file = log_file
        self._setup_logging()
        
        multiprocessing.set_start_method('spawn', force=True)
    
    def _setup_logging(self):
        """Configure logging system"""
        self.logger = logging.getLogger('GapFinder')
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
    
    def find_gaps(
        self,
        fasta_file: str,
        output_dir: str = ".",
        output_prefix: str = "gap_analysis",
        gap_char: str = "N",
        min_gap_size: int = 1,
        parallel_threshold_mb: int = 100,
        max_workers: Optional[int] = None,
        save_json: bool = True,
        save_bed: bool = False,
        save_summary: bool = True
    ) -> Dict[str, Any]:
        """
        Find gap positions in FASTA file
        
        Args:
            fasta_file: Input FASTA file path
            output_dir: Output directory
            output_prefix: Output file prefix
            gap_char: Character representing gaps
            min_gap_size: Minimum gap size
            parallel_threshold_mb: Threshold for parallel processing (MB)
            max_workers: Maximum worker processes (None=auto-detect)
            save_json: Save results in JSON format
            save_bed: Save results in BED format
            save_summary: Save text summary
            
        Returns:
            Dictionary containing gap information and statistics
        """
        self.log(f"Starting gap position search", "info")
        self.log(f"Input file: {fasta_file}", "info")
        self.log(f"Output directory: {output_dir}", "info")
        
        start_time = time.time()
        
        os.makedirs(output_dir, exist_ok=True)
        
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
        
        file_size_mb = os.path.getsize(fasta_file) / (1024 * 1024)
        self.log(f"File size: {file_size_mb:.2f} MB", "info")
        
        result = self._analyze_fasta_file(
            fasta_file=fasta_file,
            gap_char=gap_char,
            min_gap_size=min_gap_size,
            parallel_threshold_mb=parallel_threshold_mb,
            max_workers=max_workers
        )
        
        result["api_info"] = {
            "execution_time": time.time() - start_time,
            "input_file": fasta_file,
            "parameters": {
                "gap_char": gap_char,
                "min_gap_size": min_gap_size,
                "parallel_threshold_mb": parallel_threshold_mb,
                "max_workers": max_workers
            },
            "output_files": {}
        }
        
        output_files = self._save_results(
            result=result,
            output_dir=output_dir,
            output_prefix=output_prefix,
            save_json=save_json,
            save_bed=save_bed,
            save_summary=save_summary
        )
        
        result["api_info"]["output_files"] = output_files
        
        self._print_summary(result)
        
        return result
    
    def _analyze_fasta_file(
        self,
        fasta_file: str,
        gap_char: str = "N",
        min_gap_size: int = 1,
        parallel_threshold_mb: int = 100,
        max_workers: Optional[int] = None
    ) -> Dict[str, Any]:
        """Analyze gap positions in FASTA file"""
        try:
            from Bio import SeqIO
        except ImportError:
            raise ImportError("BioPython not installed. Install with: conda install -c bioconda biopython")
        
        self.log("Reading FASTA file...", "info")
        
        records_data = []
        total_genome_length = 0
        chromosomes_info = []
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id
            sequence = str(record.seq)
            seq_length = len(sequence)
            
            total_genome_length += seq_length
            records_data.append((seq_id, sequence))
            
            chromosomes_info.append({
                "id": seq_id,
                "length": seq_length
            })
        
        if not records_data:
            raise ValueError("No valid sequences found in FASTA file")
        
        self.log(f"Found {len(records_data)} sequences", "info")
        self.log(f"Total genome length: {total_genome_length:,} bp", "info")
        
        file_size_mb = os.path.getsize(fasta_file) / (1024 * 1024)
        use_parallel = file_size_mb > parallel_threshold_mb
        
        if use_parallel:
            self.log("Using parallel processing...", "info")
            if max_workers is None:
                max_workers = min(multiprocessing.cpu_count(), len(records_data))
            gaps = self._find_gaps_parallel(
                records_data=records_data,
                gap_char=gap_char,
                min_gap_size=min_gap_size,
                max_workers=max_workers
            )
        else:
            self.log("Using sequential processing...", "info")
            gaps = self._find_gaps_sequential(
                records_data=records_data,
                gap_char=gap_char,
                min_gap_size=min_gap_size
            )
        
        self.log(f"Found {len(gaps)} gap regions", "info")
        
        stats_by_chrom = {}
        for gap in gaps:
            chrom = gap['chrom']
            if chrom not in stats_by_chrom:
                stats_by_chrom[chrom] = {
                    'count': 0,
                    'total_length': 0,
                    'max_gap': 0,
                    'min_gap': float('inf'),
                    'avg_gap': 0
                }
            
            stats_by_chrom[chrom]['count'] += 1
            stats_by_chrom[chrom]['total_length'] += gap['length']
            stats_by_chrom[chrom]['max_gap'] = max(stats_by_chrom[chrom]['max_gap'], gap['length'])
            stats_by_chrom[chrom]['min_gap'] = min(stats_by_chrom[chrom]['min_gap'], gap['length'])
        
        for chrom in stats_by_chrom:
            if stats_by_chrom[chrom]['count'] > 0:
                stats_by_chrom[chrom]['avg_gap'] = stats_by_chrom[chrom]['total_length'] / stats_by_chrom[chrom]['count']
        
        total_gaps = len(gaps)
        total_gap_length = sum(gap['length'] for gap in gaps)
        
        if total_gaps > 0:
            max_gap = max(gaps, key=lambda x: x['length'])
            min_gap = min(gaps, key=lambda x: x['length'])
            avg_gap_length = total_gap_length / total_gaps
            gap_coverage = total_gap_length / total_genome_length * 100
        else:
            max_gap = min_gap = None
            avg_gap_length = gap_coverage = 0
        
        result = {
            'metadata': {
                'input_file': fasta_file,
                'file_size_mb': file_size_mb,
                'total_sequences': len(records_data),
                'chromosomes': chromosomes_info,
                'total_genome_length': total_genome_length,
                'gap_char': gap_char,
                'min_gap_size': min_gap_size,
                'processing_mode': 'parallel' if use_parallel else 'sequential',
                'parallel_threshold_mb': parallel_threshold_mb,
                'max_workers': max_workers if use_parallel else 1
            },
            'statistics': {
                'total_gaps': total_gaps,
                'total_gap_length': total_gap_length,
                'gap_coverage_percent': round(gap_coverage, 4),
                'average_gap_length': round(avg_gap_length, 2) if total_gaps > 0 else 0,
                'max_gap': max_gap,
                'min_gap': min_gap,
                'by_chromosome': stats_by_chrom
            },
            'gaps': gaps
        }
        
        return result
    
    def _find_gaps_parallel(
        self,
        records_data: List[Tuple[str, str]],
        gap_char: str = "N",
        min_gap_size: int = 1,
        max_workers: Optional[int] = None
    ) -> List[Dict[str, Any]]:
        """Find gaps in parallel for all chromosomes"""
        if max_workers is None:
            max_workers = min(multiprocessing.cpu_count(), len(records_data))
        
        self.log(f"Using {max_workers} worker processes for parallel processing", "info")
        
        args_list = [(seq_id, sequence, gap_char, min_gap_size) for seq_id, sequence in records_data]
        
        all_gaps = []
        
        try:
            with ProcessPoolExecutor(
                max_workers=max_workers,
                mp_context=multiprocessing.get_context('spawn')
            ) as executor:
                futures = [executor.submit(self._process_chromosome, args) for args in args_list]
                
                for i, future in enumerate(futures):
                    try:
                        chromosome_gaps = future.result(timeout=3600)
                        all_gaps.extend(chromosome_gaps)
                        
                        if self.verbose and (i + 1) % max(1, len(futures) // 10) == 0:
                            progress = (i + 1) / len(futures) * 100
                            self.log(f"  Progress: {progress:.1f}% ({i + 1}/{len(futures)})", "info")
                            
                    except Exception as e:
                        seq_id = args_list[i][0]
                        self.log(f"Error processing sequence {seq_id}: {e}", "error")
        
        except Exception as e:
            self.log(f"Parallel processing failed: {e}", "error")
            self.log("Falling back to sequential processing...", "warning")
            all_gaps = self._find_gaps_sequential(records_data, gap_char, min_gap_size)
        
        return all_gaps
    
    def _find_gaps_sequential(
        self,
        records_data: List[Tuple[str, str]],
        gap_char: str = "N",
        min_gap_size: int = 1
    ) -> List[Dict[str, Any]]:
        """Find gaps sequentially for all chromosomes"""
        all_gaps = []
        
        for seq_id, sequence in records_data:
            gap_pattern = re.compile(f'{re.escape(gap_char)}+')
            
            for match in gap_pattern.finditer(sequence):
                start = match.start()
                end = match.end()
                gap_length = end - start
                
                if gap_length >= min_gap_size:
                    all_gaps.append({
                        'chrom': seq_id,
                        'start': start,
                        'end': end,
                        'length': gap_length,
                        'position_human': f"{start + 1}-{end}"
                    })
        
        return all_gaps
    
    def _process_chromosome(self, args: tuple) -> List[Dict[str, Any]]:
        """Process gap finding for single chromosome (multiprocessing worker function)"""
        seq_id, sequence, gap_char, min_gap_size = args
        gaps = []
        
        gap_pattern = re.compile(f'{re.escape(gap_char)}+')
        
        for match in gap_pattern.finditer(sequence):
            start = match.start()
            end = match.end()
            gap_length = end - start
            
            if gap_length >= min_gap_size:
                gaps.append({
                    'chrom': seq_id,
                    'start': start,
                    'end': end,
                    'length': gap_length,
                    'position_human': f"{start + 1}-{end}"
                })
        
        return gaps
    
    def _save_results(
        self,
        result: Dict[str, Any],
        output_dir: str,
        output_prefix: str,
        save_json: bool = True,
        save_bed: bool = False,
        save_summary: bool = True
    ) -> Dict[str, str]:
        """Save analysis results to files"""
        output_files = {}
        
        if save_json:
            json_file = os.path.join(output_dir, f"{output_prefix}.json")
            with open(json_file, 'w') as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
            output_files['json'] = json_file
            self.log(f"JSON results saved: {json_file}", "info")
        
        if save_bed:
            bed_file = os.path.join(output_dir, f"{output_prefix}.bed")
            self._save_bed_file(result['gaps'], bed_file)
            output_files['bed'] = bed_file
            self.log(f"BED format results saved: {bed_file}", "info")
        
        if save_summary:
            summary_file = os.path.join(output_dir, f"{output_prefix}_summary.txt")
            self._save_text_summary(result, summary_file)
            output_files['summary'] = summary_file
            self.log(f"Text summary saved: {summary_file}", "info")
        
        csv_file = os.path.join(output_dir, f"{output_prefix}.csv")
        self._save_csv_file(result, csv_file)
        output_files['csv'] = csv_file
        self.log(f"CSV format results saved: {csv_file}", "info")
        
        return output_files
    
    def _save_bed_file(self, gaps: List[Dict[str, Any]], bed_file: str):
        """Save BED format file"""
        with open(bed_file, 'w') as f:
            for gap in gaps:
                chrom = gap['chrom']
                start = gap['start']
                end = gap['end']
                name = f"gap_{chrom}:{start+1}-{end}"
                score = gap['length']
                
                f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\n")
    
    def _save_text_summary(self, result: Dict[str, Any], summary_file: str):
        """Save text summary"""
        metadata = result['metadata']
        stats = result['statistics']
        
        with open(summary_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("Genome Gap Analysis Summary\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Input file: {metadata['input_file']}\n")
            f.write(f"File size: {metadata['file_size_mb']:.2f} MB\n")
            f.write(f"Number of sequences: {metadata['total_sequences']}\n")
            f.write(f"Total genome length: {metadata['total_genome_length']:,} bp\n")
            f.write(f"Gap character: '{metadata['gap_char']}'\n")
            f.write(f"Minimum gap size: {metadata['min_gap_size']} bp\n")
            f.write(f"Processing mode: {metadata['processing_mode']}\n\n")
            
            f.write("-" * 60 + "\n")
            f.write("Overall Statistics\n")
            f.write("-" * 60 + "\n")
            f.write(f"Total gaps: {stats['total_gaps']}\n")
            f.write(f"Total gap length: {stats['total_gap_length']:,} bp\n")
            f.write(f"Genome coverage: {stats['gap_coverage_percent']:.4f}%\n")
            
            if stats['total_gaps'] > 0:
                f.write(f"Average gap length: {stats['average_gap_length']:,.1f} bp\n")
                f.write(f"Largest gap: {stats['max_gap']['chrom']}:{stats['max_gap']['start']+1:,}-{stats['max_gap']['end']:,} "
                       f"({stats['max_gap']['length']:,} bp)\n")
                f.write(f"Smallest gap: {stats['min_gap']['chrom']}:{stats['min_gap']['start']+1:,}-{stats['min_gap']['end']:,} "
                       f"({stats['min_gap']['length']:,} bp)\n")
            
            f.write("\n" + "-" * 60 + "\n")
            f.write("Statistics by Chromosome\n")
            f.write("-" * 60 + "\n")
            
            for chrom, chr_stats in stats['by_chromosome'].items():
                f.write(f"\nChromosome: {chrom}\n")
                f.write(f"  Gap count: {chr_stats['count']}\n")
                f.write(f"  Total gap length: {chr_stats['total_length']:,} bp\n")
                if chr_stats['count'] > 0:
                    f.write(f"  Largest gap: {chr_stats['max_gap']:,} bp\n")
                    f.write(f"  Smallest gap: {chr_stats['min_gap']:,} bp\n")
                    f.write(f"  Average gap: {chr_stats['avg_gap']:,.1f} bp\n")
    
    def _save_csv_file(self, result: Dict[str, Any], csv_file: str):
        """Save CSV format results"""
        gaps = result['gaps']
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            writer.writerow(['Chromosome', 'Start_0based', 'End_1based', 
                           'Start_1based', 'End_1based_human', 'Length_bp'])
            
            for gap in gaps:
                chrom = gap['chrom']
                start_0based = gap['start']
                end_1based = gap['end']
                start_1based = start_0based + 1
                length = gap['length']
                
                writer.writerow([chrom, start_0based, end_1based, 
                               start_1based, f"{start_1based}-{end_1based}", length])
    
    def _print_summary(self, result: Dict[str, Any]):
        """Print analysis summary"""
        stats = result['statistics']
        
        self.log("\n" + "="*60, "info")
        self.log("Gap finding completed!", "info")
        self.log("="*60, "info")
        
        self.log(f"Total gaps found: {stats['total_gaps']}", "info")
        self.log(f"Total gap length: {stats['total_gap_length']:,} bp", "info")
        self.log(f"Genome coverage: {stats['gap_coverage_percent']:.4f}%", "info")
        
        if stats['total_gaps'] > 0:
            max_gap = stats['max_gap']
            min_gap = stats['min_gap']
            
            self.log(f"Largest gap: {max_gap['chrom']}:{max_gap['start']+1:,}-{max_gap['end']:,} "
                    f"({max_gap['length']:,} bp)", "info")
            self.log(f"Smallest gap: {min_gap['chrom']}:{min_gap['start']+1:,}-{min_gap['end']:,} "
                    f"({min_gap['length']:,} bp)", "info")
            self.log(f"Average gap length: {stats['average_gap_length']:,.1f} bp", "info")
        
        self.log("\nMajor chromosome statistics:", "info")
        chr_stats = stats['by_chromosome']
        
        sorted_chroms = sorted(chr_stats.items(), key=lambda x: x[1]['count'], reverse=True)
        
        for chrom, chr_stat in sorted_chroms[:10]:
            self.log(f"  {chrom}: {chr_stat['count']} gaps, "
                    f"Total length: {chr_stat['total_length']:,} bp", "info")
        
        if len(sorted_chroms) > 10:
            self.log(f"  ... {len(sorted_chroms) - 10} more chromosomes", "info")
    
    def batch_find_gaps(
        self,
        fasta_files: List[str],
        output_dir: str = ".",
        batch_prefix: str = "batch",
        **kwargs
    ) -> Dict[str, Any]:
        """
        Find gaps in multiple FASTA files in batch
        
        Args:
            fasta_files: List of FASTA file paths
            output_dir: Output directory
            batch_prefix: Batch processing prefix
            **kwargs: Other parameters passed to find_gaps
            
        Returns:
            Batch processing summary
        """
        self.log(f"Starting batch gap finding for {len(fasta_files)} FASTA files", "info")
        
        results = {}
        for i, fasta_file in enumerate(fasta_files):
            self.log(f"Processing file {i+1}/{len(fasta_files)}: {os.path.basename(fasta_file)}", "info")
            
            file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
            output_prefix = f"{batch_prefix}_{i+1:03d}_{file_prefix}"
            
            try:
                result = self.find_gaps(
                    fasta_file=fasta_file,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    **kwargs
                )
                results[os.path.basename(fasta_file)] = {
                    "success": True,
                    "result": result,
                    "output_prefix": output_prefix
                }
            except Exception as e:
                results[os.path.basename(fasta_file)] = {
                    "success": False,
                    "error": str(e),
                    "output_prefix": output_prefix
                }
                self.log(f"Processing failed: {e}", "warning")
        
        summary = {
            "total": len(fasta_files),
            "successful": sum(1 for r in results.values() if r["success"]),
            "failed": sum(1 for r in results.values() if not r["success"]),
            "results": results,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        summary_file = os.path.join(output_dir, f"{batch_prefix}_summary.json")
        with open(summary_file, "w") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        self.log(f"Batch gap finding completed, successful: {summary['successful']}/{summary['total']}", "info")
        self.log(f"Summary report saved: {summary_file}", "info")
        
        return summary
    
    def find_gaps_in_directory(
        self,
        directory: str,
        pattern: str = "*.fasta",
        recursive: bool = False,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Find gaps in all FASTA files in directory
        
        Args:
            directory: Directory path
            pattern: File pattern match
            recursive: Recursively search subdirectories
            **kwargs: Other parameters passed to find_gaps
            
        Returns:
            Batch processing summary
        """
        import glob
        
        if not os.path.exists(directory):
            raise FileNotFoundError(f"Directory not found: {directory}")
        
        if recursive:
            pattern_str = os.path.join(directory, "**", pattern)
        else:
            pattern_str = os.path.join(directory, pattern)
        
        fasta_files = glob.glob(pattern_str, recursive=recursive)
        
        if not fasta_files:
            raise FileNotFoundError(f"No files matching {pattern} found in directory {directory}")
        
        self.log(f"Found {len(fasta_files)} FASTA files in directory", "info")
        
        return self.batch_find_gaps(fasta_files, **kwargs)
    
    def get_gap_positions(
        self,
        gap_results: Dict[str, Any],
        chromosome: Optional[str] = None,
        min_gap_length: int = 100,
        max_gap_length: Optional[int] = None
    ) -> List[int]:
        """
        Extract gap positions from gap results
        
        Args:
            gap_results: Results dictionary returned by find_gaps
            chromosome: Specify chromosome (None=include all)
            min_gap_length: Minimum gap length
            max_gap_length: Maximum gap length (optional)
            
        Returns:
            List of gap positions (start positions, 0-based)
        """
        gaps = gap_results.get("gaps", [])
        
        filtered_gaps = []
        for gap in gaps:
            if chromosome and gap["chrom"] != chromosome:
                continue
            
            if gap["length"] < min_gap_length:
                continue
            
            if max_gap_length and gap["length"] > max_gap_length:
                continue
            
            gap_mid = gap["start"] + gap["length"] // 2
            filtered_gaps.append(gap_mid)
        
        return sorted(filtered_gaps)
    
    def get_gap_statistics(
        self,
        gap_results: Dict[str, Any],
        chromosome: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Get detailed gap statistics
        
        Args:
            gap_results: Results dictionary returned by find_gaps
            chromosome: Specify chromosome (optional)
            
        Returns:
            Detailed statistics dictionary
        """
        gaps = gap_results.get("gaps", [])
        
        if chromosome:
            gaps = [gap for gap in gaps if gap["chrom"] == chromosome]
        
        if not gaps:
            return {
                "total_gaps": 0,
                "total_length": 0,
                "length_distribution": {},
                "chromosome": chromosome or "all"
            }
        
        total_gaps = len(gaps)
        total_length = sum(gap["length"] for gap in gaps)
        lengths = [gap["length"] for gap in gaps]
        
        length_bins = {
            "0-100bp": 0,
            "101-500bp": 0,
            "501-1000bp": 0,
            "1001-5000bp": 0,
            "5001-10000bp": 0,
            "10001-50000bp": 0,
            "50001-100000bp": 0,
            "100001+": 0
        }
        
        for length in lengths:
            if length <= 100:
                length_bins["0-100bp"] += 1
            elif length <= 500:
                length_bins["101-500bp"] += 1
            elif length <= 1000:
                length_bins["501-1000bp"] += 1
            elif length <= 5000:
                length_bins["1001-5000bp"] += 1
            elif length <= 10000:
                length_bins["5001-10000bp"] += 1
            elif length <= 50000:
                length_bins["10001-50000bp"] += 1
            elif length <= 100000:
                length_bins["50001-100000bp"] += 1
            else:
                length_bins["100001+"] += 1
        
        stats = {
            "chromosome": chromosome or "all",
            "total_gaps": total_gaps,
            "total_length": total_length,
            "average_length": total_length / total_gaps if total_gaps > 0 else 0,
            "max_length": max(lengths) if lengths else 0,
            "min_length": min(lengths) if lengths else 0,
            "length_distribution": length_bins,
            "median_length": sorted(lengths)[len(lengths) // 2] if lengths else 0,
            "gap_density": total_gaps / gap_results["metadata"]["total_genome_length"] * 1e6 if gap_results["metadata"]["total_genome_length"] > 0 else 0
        }
        
        return stats

def find_gaps_in_fasta(
    fasta_file: str,
    output_json: Optional[str] = None,
    gap_char: str = "N",
    min_gap_size: int = 1,
    parallel_threshold_mb: int = 100,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Simplified interface for finding gaps in FASTA file
    
    Example:
        >>> from find_gaps import find_gaps_in_fasta
        >>> result = find_gaps_in_fasta(
        ...     fasta_file="genome.fasta",
        ...     output_json="gaps.json",
        ...     min_gap_size=100,
        ...     verbose=False
        ... )
    """
    finder = GapFinderAPI(verbose=verbose)
    
    if output_json:
        output_dir = os.path.dirname(output_json) or "."
        output_prefix = os.path.splitext(os.path.basename(output_json))[0]
    else:
        output_dir = "."
        output_prefix = "gap_analysis"
    
    return finder.find_gaps(
        fasta_file=fasta_file,
        output_dir=output_dir,
        output_prefix=output_prefix,
        gap_char=gap_char,
        min_gap_size=min_gap_size,
        parallel_threshold_mb=parallel_threshold_mb,
        save_json=output_json is not None
    )

def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Find gap positions in FASTA files for each chromosome (Unified API Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Basic usage
  python find_gaps.py genome.fasta
  
  # Specify output file, use API mode
  python find_gaps.py genome.fasta -o gaps.json --api-mode --verbose
  
  # Batch process all FASTA files in directory
  python find_gaps.py --directory genomes/ --pattern "*.fa" --api-mode --output-dir results
  
  # Python module call
  from find_gaps import GapFinderAPI
  finder = GapFinderAPI(verbose=True)
  result = finder.find_gaps(
      fasta_file="genome.fasta",
      output_dir="results",
      output_prefix="analysis",
      min_gap_size=100
  )
        """
    )
    
    import argparse
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('input_file', nargs='?', help='Input FASTA file path')
    input_group.add_argument('--directory', help='Directory path containing FASTA files')
    input_group.add_argument('--batch', nargs='+', help='Multiple FASTA file paths')
    
    parser.add_argument('-o', '--output', help='Output JSON file path')
    parser.add_argument('--output-dir', default='.', help='Output directory, default current directory')
    parser.add_argument('--prefix', default='gap_analysis', help='Output file prefix')
    parser.add_argument('--gap-char', default='N', help='Character representing gaps, default N')
    parser.add_argument('--min-gap-size', type=int, default=1, help='Minimum gap size, default 1bp')
    parser.add_argument('--parallel-threshold', type=int, default=100, 
                       help='Parallel processing threshold (MB), default 100MB')
    parser.add_argument('--max-workers', type=int, help='Maximum worker processes (auto-detect)')
    parser.add_argument('--save-bed', action='store_true', help='Save BED format results')
    parser.add_argument('--save-csv', action='store_true', help='Save CSV format results')
    parser.add_argument('--api-mode', action='store_true', help='Use API mode')
    parser.add_argument('--verbose', action='store_true', help='Show detailed output')
    parser.add_argument('--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--log-file', help='Log file path')
    parser.add_argument('--recursive', action='store_true', help='Recursively search subdirectories')
    parser.add_argument('--pattern', default='*.fasta', help='File pattern match, default *.fasta')
    
    args = parser.parse_args()
    
    verbose = args.verbose and not args.quiet
    
    try:
        if args.api_mode:
            finder = GapFinderAPI(verbose=verbose, log_file=args.log_file)
            
            if args.directory:
                result = finder.find_gaps_in_directory(
                    directory=args.directory,
                    pattern=args.pattern,
                    recursive=args.recursive,
                    output_dir=args.output_dir,
                    output_prefix=args.prefix,
                    gap_char=args.gap_char,
                    min_gap_size=args.min_gap_size,
                    parallel_threshold_mb=args.parallel_threshold,
                    max_workers=args.max_workers,
                    save_json=True,
                    save_bed=args.save_bed,
                    save_summary=True
                )
                
            elif args.batch:
                result = finder.batch_find_gaps(
                    fasta_files=args.batch,
                    output_dir=args.output_dir,
                    batch_prefix=args.prefix,
                    gap_char=args.gap_char,
                    min_gap_size=args.min_gap_size,
                    parallel_threshold_mb=args.parallel_threshold,
                    max_workers=args.max_workers
                )
                
            else:
                result = finder.find_gaps(
                    fasta_file=args.input_file,
                    output_dir=args.output_dir,
                    output_prefix=args.prefix,
                    gap_char=args.gap_char,
                    min_gap_size=args.min_gap_size,
                    parallel_threshold_mb=args.parallel_threshold,
                    max_workers=args.max_workers,
                    save_json=True,
                    save_bed=args.save_bed,
                    save_summary=True
                )
            
            if args.output:
                output_json = args.output
                with open(output_json, 'w') as f:
                    json.dump(result, f, indent=2, ensure_ascii=False)
                if verbose:
                    print(f"Results saved to: {output_json}")
            
        else:
            from Bio import SeqIO
            
            if args.input_file:
                fasta_file = args.input_file
            else:
                print("Error: Original mode requires input file specification", file=sys.stderr)
                sys.exit(1)
            
            if not os.path.exists(fasta_file):
                print(f"Error: File '{fasta_file}' not found", file=sys.stderr)
                sys.exit(1)
            
            output_json = args.output
            if output_json is None:
                base_name = os.path.splitext(os.path.basename(fasta_file))[0]
                output_json = f"{base_name}_gaps.json"
            
            result = find_gaps_in_fasta(
                fasta_file=fasta_file,
                output_json=output_json,
                gap_char=args.gap_char,
                min_gap_size=args.min_gap_size,
                parallel_threshold_mb=args.parallel_threshold,
                verbose=verbose
            )
            
            stats = result['statistics']
            print(f"\nAnalysis completed:", file=sys.stderr)
            print(f"Total gaps found: {stats['total_gaps']}", file=sys.stderr)
            print(f"Total gap length: {stats['total_gap_length']:,} bp", file=sys.stderr)
            print(f"Genome coverage: {stats['gap_coverage_percent']:.4f}%", file=sys.stderr)
            
            if stats['total_gaps'] > 0:
                print(f"Largest gap: {stats['max_gap']['chrom']}:{stats['max_gap']['start']+1:,}-{stats['max_gap']['end']:,} "
                      f"({stats['max_gap']['length']:,} bp)", file=sys.stderr)
                print(f"Smallest gap: {stats['min_gap']['chrom']}:{stats['min_gap']['start']+1:,}-{stats['min_gap']['end']:,} "
                      f"({stats['min_gap']['length']:,} bp)", file=sys.stderr)
                print(f"Average gap length: {stats['average_gap_length']:,.1f} bp", file=sys.stderr)
            
            print(f"Results saved to: {output_json}", file=sys.stderr)
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ImportError as e:
        print(f"Error: Missing dependency - {e}", file=sys.stderr)
        print("Install required dependencies: conda install -c bioconda biopython", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Exception while processing file - {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)
    main()