#!/usr/bin/env python3

import sys
import argparse
import subprocess
import os
import re
import time
import math
import shutil
import json
from typing import List, Dict, Optional, Tuple, Any, Union
from pathlib import Path

class GapFillerAPI:
    
    def __init__(self, verbose: bool = True, log_file: Optional[str] = None):
        self.verbose = verbose
        self.log_file = log_file
        self.temp_dirs = []
        self._setup_logging()
    
    def _setup_logging(self):
        import logging
        
        self.logger = logging.getLogger('GapFiller')
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
        if level == "info":
            self.logger.info(message)
        elif level == "warning":
            self.logger.warning(message)
        elif level == "error":
            self.logger.error(message)
        elif level == "debug":
            self.logger.debug(message)
    
    def preprocess_contigs_files(
        self,
        input_files: List[str],
        output_dir: str,
        min_gap_length: int = 100,
        min_contig_length: int = 1000
    ) -> List[str]:
        self.log(f"Starting contig file preprocessing, splitting contigs with {min_gap_length}N or longer gaps", "info")
        processed_files = []
        os.makedirs(output_dir, exist_ok=True)
        
        for i, input_file in enumerate(input_files):
            self.log(f"Processing file: {input_file}", "info")
            
            contig_dict = _read_fasta_dict(input_file)
            processed_dict = {}
            split_count = 0
            total_contigs = len(contig_dict)
            total_length_before = sum(len(seq) for seq in contig_dict.values())
            
            for contig_id, sequence in contig_dict.items():
                gap_pattern = 'N' * min_gap_length
                if gap_pattern in sequence:
                    gap_matches = list(re.finditer(f'N{{{min_gap_length},}}', sequence))
                    
                    if gap_matches:
                        last_end = 0
                        for gap_idx, match in enumerate(gap_matches, 1):
                            gap_start = match.start()
                            gap_end = match.end()
                            
                            if gap_start > last_end:
                                part_seq = sequence[last_end:gap_start]
                                if len(part_seq) >= min_contig_length:
                                    new_id = f"{contig_id}_p{gap_idx}"
                                    processed_dict[new_id] = part_seq
                                    split_count += 1
                            
                            last_end = gap_end
                        
                        if last_end < len(sequence):
                            part_seq = sequence[last_end:]
                            if len(part_seq) >= min_contig_length:
                                new_id = f"{contig_id}_p{len(gap_matches)+1}"
                                processed_dict[new_id] = part_seq
                                split_count += 1
                else:
                    if len(sequence) >= min_contig_length:
                        processed_dict[contig_id] = sequence
            
            output_file = os.path.join(output_dir, f"reads_{i}.fa")
            with open(output_file, 'w') as f:
                for seq_id, seq in processed_dict.items():
                    f.write(f'>{seq_id}\n{seq}\n')
            
            processed_files.append(output_file)
            
            kept_count = len(processed_dict)
            total_length_after = sum(len(seq) for seq in processed_dict.values())
            
            self.log(f"  Original: {total_contigs} contigs, {total_length_before:,} bp", "info")
            self.log(f"  Processed: {kept_count} contigs, {total_length_after:,} bp", "info")
            self.log(f"  Split operations: {split_count}", "info")
            
            if total_contigs > 0:
                self.log(f"  Retention ratio: {kept_count/total_contigs:.1%} (contigs), {total_length_after/total_length_before:.1%} (bases)", "info")
        
        self.log(f"Preprocessing completed, generated {len(processed_files)} files", "info")
        return processed_files
    
    def run_tgsgapcloser(
        self,
        draft_fasta: str,
        reads_fasta: Union[str, List[str]],
        output_dir: str,
        prefix: str = "tgs_stage",
        minmap_arg: str = "-x asm5",
        threads: int = 32,
        tgstype: str = "ont",
        keep_tgs_temp: bool = False,
        preprocess_contigs: bool = True,
        min_gap_length: int = 100,
        min_contig_length: int = 1000
    ) -> Dict[str, Any]:
        self.log(f"Starting TGS-GapCloser pre-filling", "info")
        self.log(f"Input genome: {draft_fasta}", "info")
        self.log(f"Output directory: {output_dir}/{prefix}", "info")
        
        start_time = time.time()
        
        tgs_output_dir = os.path.join(output_dir, prefix)
        os.makedirs(tgs_output_dir, exist_ok=True)
        
        original_cwd = os.getcwd()
        
        draft_fasta = os.path.abspath(draft_fasta)
        if isinstance(reads_fasta, str):
            reads_files = [os.path.abspath(reads_fasta)]
        else:
            reads_files = [os.path.abspath(f) for f in reads_fasta]
        
        if not os.path.exists(draft_fasta):
            raise FileNotFoundError(f"Genome file not found: {draft_fasta}")
        
        for read_file in reads_files:
            if not os.path.exists(read_file):
                raise FileNotFoundError(f"Read file not found: {read_file}")
        
        try:
            os.chdir(tgs_output_dir)
            self.log(f"Changed to working directory: {tgs_output_dir}", "info")
            
            if preprocess_contigs:
                self.log(f"Preprocessing contigs: splitting contigs with {min_gap_length}N or longer gaps", "info")
                
                processed_reads_files = []
                for i, read_file in enumerate(reads_files):
                    self.log(f"Processing file {i+1}/{len(reads_files)}: {read_file}", "info")
                    
                    try:
                        contig_dict = _read_fasta_dict(read_file)
                    except Exception as e:
                        self.log(f"Cannot read file {read_file}: {e}", "error")
                        raise
                    
                    processed_dict = {}
                    split_count = 0
                    
                    for contig_id, sequence in contig_dict.items():
                        gap_pattern = 'N' * min_gap_length
                        if gap_pattern in sequence:
                            gap_matches = list(re.finditer(f'N{{{min_gap_length},}}', sequence))
                            
                            if gap_matches:
                                last_end = 0
                                for gap_idx, match in enumerate(gap_matches, 1):
                                    gap_start = match.start()
                                    gap_end = match.end()
                                    
                                    if gap_start > last_end:
                                        part_seq = sequence[last_end:gap_start]
                                        if len(part_seq) >= min_contig_length:
                                            new_id = f"{contig_id}_p{gap_idx}"
                                            processed_dict[new_id] = part_seq
                                            split_count += 1
                                    
                                    last_end = gap_end
                                
                                if last_end < len(sequence):
                                    part_seq = sequence[last_end:]
                                    if len(part_seq) >= min_contig_length:
                                        new_id = f"{contig_id}_p{len(gap_matches)+1}"
                                        processed_dict[new_id] = part_seq
                                        split_count += 1
                        else:
                            if len(sequence) >= min_contig_length:
                                processed_dict[contig_id] = sequence
                    
                    output_file = f"reads_{i}.fa"
                    with open(output_file, 'w') as f:
                        for seq_id, seq in processed_dict.items():
                            f.write(f'>{seq_id}\n{seq}\n')
                    
                    processed_reads_files.append(output_file)
                    self.log(f"  Generated file: {output_file}, splits: {split_count}", "info")
                
                reads_files = processed_reads_files
                self.log(f"Preprocessing completed, using {len(reads_files)} files", "info")
            else:
                temp_reads_files = []
                for i, read_file in enumerate(reads_files):
                    temp_file = f"input_{i}.fa"
                    shutil.copy2(read_file, temp_file)
                    temp_reads_files.append(temp_file)
                reads_files = temp_reads_files
            
            draft_basename = os.path.basename(draft_fasta)
            if not os.path.exists(draft_basename):
                self.log(f"Copying genome file to working directory: {draft_fasta} -> {draft_basename}", "info")
                shutil.copy2(draft_fasta, draft_basename)
            
            draft_fasta_local = draft_basename
            
            cmd = [
                "tgsgapcloser",
                "--scaff", draft_fasta_local,
                "--reads", ",".join(reads_files),
                "--output", "tgs_filled",
                "--minmap_arg", minmap_arg,
                "--thread", str(threads),
                "--tgstype", tgstype,
                "--ne"
            ]
            
            cmd_str = " ".join(cmd)
            self.log(f"TGS-GapCloser command: {cmd_str}", "info")
            
            stdout_file = "tgs_stdout.log"
            stderr_file = "tgs_stderr.log"
            
            try:
                with open(stdout_file, "w") as stdout_f, open(stderr_file, "w") as stderr_f:
                    self.log(f"Starting TGS-GapCloser...", "info")
                    process = subprocess.Popen(
                        cmd,
                        stdout=stdout_f,
                        stderr=stderr_f,
                        universal_newlines=True
                    )
                    return_code = process.wait()
                
                if return_code != 0:
                    self.log(f"TGS-GapCloser failed, return code: {return_code}", "error")
                    if os.path.exists(stderr_file):
                        with open(stderr_file, "r") as f:
                            error_msg = f.read()
                            self.log(f"Error message: {error_msg[:500]}", "error")
                    raise RuntimeError(f"TGS-GapCloser execution failed")
                
                tgs_output_file = "tgs_filled.scaff_seqs"
                if not os.path.exists(tgs_output_file):
                    possible_files = [
                        "tgs_filled.fa",
                        "tgs_filled.fasta",
                        "tgs_filled.scaffold.fa",
                        "tgs_filled.scaffold.fasta",
                        "tgs_filled_scaffolds.fasta",
                    ]
                    
                    found = False
                    for file_name in possible_files:
                        if os.path.exists(file_name):
                            tgs_output_file = file_name
                            found = True
                            self.log(f"Found output file: {tgs_output_file}", "info")
                            break
                    
                    if not found:
                        current_files = os.listdir('.')
                        self.log(f"Current directory files: {current_files}", "warning")
                        raise FileNotFoundError(f"TGS-GapCloser output file not found")
                
                tgs_stats = self._analyze_gap_filling(draft_fasta_local, tgs_output_file)
                
                result = {
                    "success": True,
                    "tgs_output_file": os.path.join(tgs_output_dir, tgs_output_file),
                    "tgs_output_dir": tgs_output_dir,
                    "stdout_file": os.path.join(tgs_output_dir, stdout_file),
                    "stderr_file": os.path.join(tgs_output_dir, stderr_file),
                    "command": cmd_str,
                    "execution_time": time.time() - start_time,
                    "gaps_before": tgs_stats["gaps_before"],
                    "gaps_after": tgs_stats["gaps_after"],
                    "gaps_closed_by_tgs": tgs_stats["gaps_before"] - tgs_stats["gaps_after"],
                    "tgs_stats": tgs_stats,
                    "preprocess_contigs": preprocess_contigs
                }
                
                self.log(f"TGS-GapCloser completed, filled {result['gaps_closed_by_tgs']} gaps", "info")
                self.log(f"Output file: {tgs_output_file}", "info")
                
                return result
                
            except Exception as e:
                self.log(f"TGS-GapCloser execution error: {str(e)}", "error")
                raise
                
        finally:
            os.chdir(original_cwd)
            self.log(f"Restored working directory: {original_cwd}", "info")
    
    def _analyze_gap_filling(self, original_fasta: str, filled_fasta: str) -> Dict[str, Any]:
        def count_gaps(fasta_file: str) -> int:
            gap_count = 0
            with open(fasta_file, 'r') as f:
                seq_id = ""
                current_seq = []
                for line in f:
                    if line.startswith('>'):
                        if seq_id and current_seq:
                            seq = ''.join(current_seq)
                            gap_count += len(re.findall(r'N{100,}', seq))
                        seq_id = line[1:].strip().split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                
                if seq_id and current_seq:
                    seq = ''.join(current_seq)
                    gap_count += len(re.findall(r'N{100,}', seq))
            return gap_count
        
        gaps_before = count_gaps(original_fasta)
        gaps_after = count_gaps(filled_fasta)
        
        return {
            "gaps_before": gaps_before,
            "gaps_after": gaps_after,
            "closure_rate": (gaps_before - gaps_after) / gaps_before if gaps_before > 0 else 0
        }
    
    def fill_gaps(
        self,
        draft_fasta: str,
        gap_patches: Union[str, List[str]],
        output_dir: str = ".",
        output_prefix: str = "gap_filled",
        flank_size: int = 5000,
        min_alignment_length: int = 1000,
        min_identity: float = 40.0,
        max_fill_length: int = 1000000,
        aligner: str = "minimap2",
        threads: int = 1,
        overwrite: bool = False,
        keep_temp: bool = False,
        return_dict: bool = True,
        use_tgs: bool = True,
        tgs_minmap: str = "-x asm5",
        tgs_threads: int = 32,
        tgs_type: str = "ont",
        keep_tgs_temp: bool = False,
        preprocess_contigs: bool = True,
        min_gap_length: int = 100,
        min_contig_length: int = 1000
    ) -> Dict[str, Any]:
        self.log(f"Starting two-stage gap filling analysis", "info")
        start_time = time.time()
        
        if isinstance(gap_patches, str):
            gap_patches = [gap_patches]
        
        os.makedirs(output_dir, exist_ok=True)
        
        tgs_result = None
        stage1_output = draft_fasta
        
        if use_tgs:
            try:
                self.log(f"=== Stage 1: TGS-GapCloser pre-filling ===", "info")
                tgs_result = self.run_tgsgapcloser(
                    draft_fasta=draft_fasta,
                    reads_fasta=gap_patches,
                    output_dir=output_dir,
                    prefix=f"{output_prefix}_tgs",
                    minmap_arg=tgs_minmap,
                    threads=tgs_threads,
                    tgstype=tgs_type,
                    keep_tgs_temp=keep_tgs_temp,
                    preprocess_contigs=preprocess_contigs,
                    min_gap_length=min_gap_length,
                    min_contig_length=min_contig_length
                )
                stage1_output = tgs_result["tgs_output_file"]
                
            except Exception as e:
                if not os.path.exists(draft_fasta):
                    raise RuntimeError(f"TGS stage failed and original file unavailable: {str(e)}")
                self.log(f"TGS-GapCloser stage failed, will continue with original file: {str(e)}", "warning")
                stage1_output = draft_fasta
        
        self.log(f"=== Stage 2: Flanking alignment fine filling ===", "info")
        
        try:
            result = gap_filler(
                draft_genome=stage1_output,
                gapcloser_contig=gap_patches,
                flanking_len=flank_size,
                min_alignment_length=min_alignment_length,
                min_alignment_identity=min_identity,
                max_filling_len=max_fill_length,
                aligner=aligner,
                prefix=output_prefix,
                output_dir=output_dir,
                threads=threads,
                overwrite=overwrite,
                minimapoption="-x asm5",
                noplot=True,
                keep_temp=keep_temp
            )
            
            if return_dict:
                result["api_info"] = {
                    "execution_time": time.time() - start_time,
                    "input_files": {
                        "draft_fasta": draft_fasta,
                        "gap_patches": gap_patches
                    },
                    "parameters": {
                        "flank_size": flank_size,
                        "min_alignment_length": min_alignment_length,
                        "min_identity": min_identity,
                        "max_fill_length": max_fill_length,
                        "aligner": aligner,
                        "threads": threads,
                        "output_dir": output_dir,
                        "output_prefix": output_prefix,
                        "use_tgs": use_tgs,
                        "tgs_minmap": tgs_minmap,
                        "tgs_threads": tgs_threads,
                        "tgs_type": tgs_type,
                        "preprocess_contigs": preprocess_contigs,
                        "min_gap_length": min_gap_length,
                        "min_contig_length": min_contig_length
                    }
                }
                
                if tgs_result:
                    result["tgs_stage"] = tgs_result
                    result["total_gaps_closed"] = (
                        tgs_result.get("gaps_closed_by_tgs", 0) + 
                        result.get("gaps_closed", 0)
                    )
                else:
                    result["tgs_stage"] = None
                    result["total_gaps_closed"] = result.get("gaps_closed", 0)
                
                combined_report = self._generate_combined_report(result, tgs_result)
                report_file = os.path.join(output_dir, f"{output_prefix}.combined_report.txt")
                with open(report_file, "w") as f:
                    f.write(combined_report)
                result["combined_report_file"] = report_file
                
                return result
            else:
                return {
                    "filled_fasta": result.get("filled_fasta"),
                    "detail_file": result.get("detail_file"),
                    "stat_file": result.get("stat_file"),
                    "agp_file": result.get("modified_agp"),
                    "tgs_output": tgs_result["tgs_output_file"] if tgs_result else None
                }
                
        except Exception as e:
            self.log(f"Flanking alignment filling stage failed: {str(e)}", "error")
            raise
    
    def _generate_combined_report(self, flank_result: Dict, tgs_result: Optional[Dict]) -> str:
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("             Two-Stage Gap Filling Comprehensive Report")
        report_lines.append("=" * 80)
        report_lines.append("")
        
        report_lines.append("[Overall Statistics]")
        report_lines.append(f"Total execution time: {flank_result['api_info']['execution_time']:.2f} seconds")
        report_lines.append(f"Output directory: {flank_result['api_info']['parameters']['output_dir']}")
        report_lines.append("")
        
        if tgs_result:
            report_lines.append("[Stage 1: TGS-GapCloser]")
            report_lines.append(f"Status: Completed successfully")
            report_lines.append(f"Execution time: {tgs_result['execution_time']:.2f} seconds")
            
            if tgs_result.get("preprocess_contigs", False):
                report_lines.append(f"Contig preprocessing: Yes (split contigs with ≥{flank_result['api_info']['parameters'].get('min_gap_length', 100)}N gaps)")
            else:
                report_lines.append(f"Contig preprocessing: No")
            
            report_lines.append(f"Original gap count: {tgs_result['gaps_before']}")
            report_lines.append(f"Gaps after TGS: {tgs_result['gaps_after']}")
            report_lines.append(f"Gaps filled by TGS: {tgs_result['gaps_closed_by_tgs']}")
            report_lines.append(f"Filling rate: {tgs_result['tgs_stats']['closure_rate']:.2%}")
            report_lines.append(f"Output file: {tgs_result['tgs_output_file']}")
            report_lines.append("")
        else:
            report_lines.append("[Stage 1: TGS-GapCloser]")
            report_lines.append("Status: Not executed or execution failed")
            report_lines.append("")
        
        report_lines.append("[Stage 2: Flanking alignment filling]")
        report_lines.append(f"Input gap count: {flank_result.get('gaps_remaining', 'N/A') + flank_result.get('gaps_closed', 0)}")
        report_lines.append(f"Filled gaps: {flank_result.get('gaps_closed', 'N/A')}")
        report_lines.append(f"Remaining gaps: {flank_result.get('gaps_remaining', 'N/A')}")
        report_lines.append(f"Total filled length: {flank_result.get('total_filled_length', 'N/A'):,} bp")
        report_lines.append(f"Final genome length: {flank_result.get('total_genome_length', 'N/A'):,} bp")
        report_lines.append(f"GC content: {flank_result.get('gc_content', 'N/A'):.2%}")
        report_lines.append("")
        
        report_lines.append("[Output files]")
        report_lines.append(f"1. Final filled genome: {flank_result.get('filled_fasta', 'N/A')}")
        report_lines.append(f"2. Filling details: {flank_result.get('detail_file', 'N/A')}")
        report_lines.append(f"3. Genome statistics: {flank_result.get('stat_file', 'N/A')}")
        report_lines.append(f"4. AGP file: {flank_result.get('modified_agp', 'N/A')}")
        if tgs_result:
            report_lines.append(f"5. TGS output genome: {tgs_result.get('tgs_output_file', 'N/A')}")
            report_lines.append(f"6. TGS log: {tgs_result.get('stdout_file', 'N/A')}")
        report_lines.append(f"7. This report: {flank_result.get('combined_report_file', 'N/A')}")
        report_lines.append("")
        
        report_lines.append("[Parameter settings]")
        params = flank_result['api_info']['parameters']
        for key, value in params.items():
            report_lines.append(f"  {key}: {value}")
        
        return "\n".join(report_lines)
    
    def cleanup(self):
        for temp_dir in self.temp_dirs:
            if os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir)
                    self.log(f"Cleaned temporary directory: {temp_dir}", "info")
                except Exception as e:
                    self.log(f"Cannot clean temporary directory {temp_dir}: {e}", "warning")
        self.temp_dirs.clear()
    
    def __del__(self):
        self.cleanup()

def _decompress(file):
    if 'gzip compressed data' in subprocess.run(f'file {file}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout.decode():
        subprocess.run(f'gzip -d {file}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        file = '.'.join(file.split('.')[:-1])
    return file

def _read_fasta_dict(fastafile):
    if not os.path.exists(fastafile):
        current_dir = os.getcwd()
        possible_path = os.path.join(current_dir, fastafile)
        if os.path.exists(possible_path):
            fastafile = possible_path
        else:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            possible_path = os.path.join(script_dir, fastafile)
            if os.path.exists(possible_path):
                fastafile = possible_path
    
    fastaDict = {}
    with open(fastafile, 'r') as fil:
        allline = fil.read()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict

def _reverse_complement(seq: str):
    trans = str.maketrans('ATCG', 'TAGC')
    return seq[::-1].translate(trans)

def _minimap(reffasta, qryfasta, prefix, suffix, minimapoption, tmp_dir, aligner, output_dir):
    print(f'[Info] Running {aligner}...')
    
    coordfile = os.path.join(tmp_dir, f'{prefix}.{suffix}.paf')
    
    if not os.path.exists(coordfile):
        cmdr = subprocess.run(
            f'{aligner} {minimapoption} -c -o {coordfile} {reffasta} {qryfasta}',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
        )
        if '[morecore]' in cmdr.stderr.decode("utf-8") or cmdr.returncode < 0:
            print('[Error] Memory insufficient.')
            sys.exit(1)
        elif cmdr.returncode != 0:
            print(f'[Error] Unexpected error in {aligner}:')
            print(f'cmd: {cmdr.args}')
            print(f'returncode: {cmdr.returncode}')
            print('stdout:', cmdr.stdout.decode("utf-8"))
            print('stderr:', cmdr.stderr.decode("utf-8"))
            sys.exit(1)
    
    if os.path.getsize(coordfile) == 0:
        print(f'[Error] No alignment found.')
        sys.exit(1)
    
    return coordfile

def _prepare_temp_directory(prefix, overwrite, output_dir=None):
    base_prefix = prefix
    
    if output_dir:
        tmp_dir = os.path.join(output_dir, f'tmp_{base_prefix}')
    else:
        tmp_dir = f'tmp_{base_prefix}'
    
    if os.path.exists(tmp_dir):
        if overwrite:
            print(f'[Info] Overwriting existing temporary directory: {tmp_dir}')
            shutil.rmtree(tmp_dir)
        else:
            print(f'[Info] Using existing temporary directory: {tmp_dir}')
    
    os.makedirs(tmp_dir, exist_ok=True)
    
    print(f'[Info] Temporary directory created: {tmp_dir}')
    return tmp_dir

def _check_prerequisite(prerequisitelist: list):
    prerequisitenotfound = []
    for prerequisite in prerequisitelist:
        cmd = subprocess.run(f'which {prerequisite}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if cmd.stdout == b'':
            prerequisitenotfound.append(prerequisite)
    if prerequisitenotfound:
        for prerequisite in prerequisitenotfound:
            print(f'[Error] prerequisite not found: {prerequisite}')
        print('[Error] Please make sure these software have been installed, exported to $PATH, and authorized executable.')
        sys.exit(1)

def gap_filler(
    draft_genome: str,
    gapcloser_contig: list,
    flanking_len: int = 5000,
    min_alignment_length: int = 1000,
    min_alignment_identity: float = 40.0,
    max_filling_len: int = 1000000,
    aligner: str = 'minimap2',
    prefix: str = 'quarTeT',
    output_dir: str = ".",
    threads: int = 1,
    overwrite: bool = False,
    minimapoption: str = '-x asm5',
    noplot: bool = True,
    keep_temp: bool = False
) -> dict:
    
    if min_alignment_length > flanking_len:
        raise ValueError('min_alignment_length must be <= flanking_len')
    if not (0 <= min_alignment_identity <= 100):
        raise ValueError('min_alignment_identity must be between 0 and 100')
    
    os.makedirs(output_dir, exist_ok=True)
    
    tmp_dir = _prepare_temp_directory(prefix, overwrite, output_dir)
    
    draftgenomefile = _decompress(draft_genome)
    gapclosercontigfilelist = [_decompress(f) for f in gapcloser_contig]
    flanking = flanking_len
    minalignmentlength2 = min_alignment_length
    minalignmentidentity2 = min_alignment_identity / 100.0
    maxfillinglen = max_filling_len
    minimapoption = f"{minimapoption} -t {threads}"
    
    print(f'[Info] Starting GapFiller with parameters:')
    print(f'  Draft genome: {draftgenomefile}')
    print(f'  Contig files: {gapclosercontigfilelist}')
    print(f'  Flanking length: {flanking}')
    print(f'  Min alignment length: {minalignmentlength2}')
    print(f'  Min alignment identity: {minalignmentidentity2:.2%}')
    print(f'  Max filling length: {maxfillinglen}')
    print(f'  Aligner: {aligner}')
    print(f'  Prefix: {prefix}')
    print(f'  Output directory: {output_dir}')
    print(f'  Threads: {threads}')
    print(f'  Temporary directory: {tmp_dir}')
    print(f'[Note] Temporary files will NOT be automatically cleaned up.')
    
    _check_prerequisite([aligner])
    
    print('[Info] Getting gaps flanking sequence...')
    draftgenomedict = _read_fasta_dict(draftgenomefile)
    flankingdict = {}
    gapdict = {}
    for sid, seq in draftgenomedict.items():
        if 'N' * 100 in seq:
            i = 1
            gapsitelist = [r.span() for r in re.finditer(r'N{100,}', seq)]
            for start_n, end_n in gapsitelist:
                start = max(start_n - flanking, 0)
                end = min(end_n + flanking, len(seq))
                leftseq = seq[start:start_n]
                rightseq = seq[end_n:end]
                if 'N' * 100 in leftseq or 'N' * 100 in rightseq:
                    print(f'[Warning] Flanking of gap {sid}.{i} contains another gap.')
                else:
                    flankingdict[f'{sid}.{i}.L'] = leftseq
                    flankingdict[f'{sid}.{i}.R'] = rightseq
                gapdict[f'{sid}.{i}'] = seq[start_n:end_n]
                i += 1
    if not flankingdict:
        raise ValueError('Input genome does not have valid gap (≥100 Ns).')
    
    flankingfastafile = os.path.join(tmp_dir, f'{prefix}.gap.flanking.fasta')
    with open(flankingfastafile, 'w') as f:
        for sid, seq in flankingdict.items():
            f.write(f'>{sid}\n{seq}\n')
    del draftgenomedict
    
    print(f'[Info] Flanking sequences saved to: {flankingfastafile}')
    
    gapcloserdict = {}
    for gapfillfile in gapclosercontigfilelist:
        gapfiller = os.path.basename(gapfillfile)
        print(f'[Info] Gap-filling with {gapfiller}...')
        
        gapfillfasta = _read_fasta_dict(gapfillfile)
        gapfilldict = {}
        needsplit = False
        for sid, seq in gapfillfasta.items():
            if 'N' * 100 in seq:
                needsplit = True
                parts = re.split(r'N{100,}', seq)
                for idx, part in enumerate(parts, 1):
                    if part:
                        gapfilldict[f'{sid}_tig{idx}'] = part
            else:
                gapfilldict[sid] = seq
        
        if needsplit:
            print(f'[Info] Gaps found in {gapfiller}. Splitting into contigs.')
            split_file = os.path.join(tmp_dir, f'{gapfiller}.splitcontig.fasta')
            with open(split_file, 'w') as c:
                for tigid, seq in gapfilldict.items():
                    c.write(f'>{tigid}\n{seq}\n')
            gapfillfile = split_file
            gapfillfasta = gapfilldict
        
        tmp_gapfill_fasta = os.path.join(tmp_dir, f'{prefix}.gapfillfasta.fasta')
        with open(tmp_gapfill_fasta, 'w') as tmpf:
            for sid, seq in gapfillfasta.items():
                tmpf.write(f'>{sid}\n{seq}\n')
        del gapfillfasta, gapfilldict
        
        paf_file = _minimap(gapfillfile, flankingfastafile, prefix, 
                         f'flank_map_{os.path.basename(gapfillfile)}', 
                         minimapoption, tmp_dir, aligner, output_dir)
        
        print(f'[Info] Alignment results saved to: {paf_file}')
        
        allalignment = []
        with open(paf_file) as paf:
            for line in paf:
                fields = line.strip().split()
                if len(fields) < 11:
                    continue
                qryid, qrylen, qrystart, qryend, strand, refid, reflen, refstart, refend, match, alignlen = fields[:11]
                alignlen = int(alignlen)
                match = int(match)
                if alignlen >= minalignmentlength2 and match / alignlen >= minalignmentidentity2:
                    allalignment.append({
                        'qryid': qryid,
                        'qrylen': int(qrylen),
                        'qrystart': int(qrystart) + 1,
                        'qryend': int(qryend),
                        'strand': strand,
                        'refid': refid,
                        'reflen': int(reflen),
                        'refstart': int(refstart) + 1,
                        'refend': int(refend),
                        'match': match,
                        'alignlen': alignlen,
                        'identity': match / alignlen,
                        'gapid': '.'.join(qryid.split('.')[:-1]),
                        'LR': qryid.split('.')[-1]
                    })
        
        if not allalignment:
            print(f'[Info] No valid alignment for {gapfiller}.')
            continue
        
        gapfillfasta = _read_fasta_dict(tmp_gapfill_fasta)
        
        print(f'[Info] Analyzing alignments for {gapfiller}...')
        for gapid in gapdict:
            Leftanchor = [a for a in allalignment if a['gapid'] == gapid and a['LR'] == 'L']
            Rightanchor = [a for a in allalignment if a['gapid'] == gapid and a['LR'] == 'R']
            if not Leftanchor or not Rightanchor:
                continue
            
            best_score = gapcloserdict.get(gapid, {}).get('score', 0)
            for Laln in Leftanchor:
                for Raln in Rightanchor:
                    if Laln['refid'] != Raln['refid'] or Laln['strand'] != Raln['strand']:
                        continue
                    
                    current_strand = Laln['strand']
                    L = Laln
                    R = Raln
                    
                    if current_strand == '-':
                        L, R = R, L
                    
                    if L['refend'] < R['refstart']:
                        score = (L['identity'] + R['identity']) / 2
                        if score <= best_score:
                            continue
                        
                        fill_start = L['refend'] + 1
                        fill_end = R['refstart'] - 1
                        if fill_start > fill_end:
                            fill_seq = ''
                        else:
                            fill_seq = gapfillfasta[L['refid']][fill_start - 1:fill_end]
                        
                        if not fill_seq or len(fill_seq) > maxfillinglen:
                            continue
                        
                        if current_strand == '-':
                            fill_seq = _reverse_complement(fill_seq)
                        
                        gapcloserdict[gapid] = {
                            'sid': f'{gapfillfile}@{L["refid"]}',
                            'range': f'{fill_start}-{fill_end}',
                            'seq': fill_seq,
                            'strand': current_strand,
                            'score': score
                        }
                        best_score = score
    
    if not gapcloserdict:
        print('[Info] No gap can be closed. Will output original genome with empty filling records.')
        filledfastafile = os.path.join(output_dir, f'{prefix}.genome.filled.fasta')
        shutil.copy2(draftgenomefile, filledfastafile)
    
        return {
            'detail_file': os.path.join(output_dir, f'{prefix}.genome.filled.detail'),
            'filled_fasta': filledfastafile,
            'modified_agp': os.path.join(output_dir, f'{prefix}.genome.filled.modified.agp'),
            'stat_file': os.path.join(output_dir, f'{prefix}.genome.filled.stat'),
            'gaps_closed': 0,
            'gaps_remaining': len(gapdict),
            'total_filled_length': 0,
            'total_genome_length': sum(len(seq) for seq in _read_fasta_dict(draftgenomefile).values()),
            'gc_content': 0,
            'temp_directory': tmp_dir,
            'output_directory': output_dir
         }
    
    print('[Info] Generating filled FASTA and reports...')
    
    filldetailfile = os.path.join(output_dir, f'{prefix}.genome.filled.detail')
    totalfilledlen = sum(len(v['seq']) for v in gapcloserdict.values())
    with open(filldetailfile, 'w') as de:
        de.write(f'# Gap Closed: {len(gapcloserdict)}\n')
        de.write(f'# Total Filled length: {totalfilledlen}\n')
        de.write(f'# Gap Remains: {len(gapdict) - len(gapcloserdict)}\n')
        de.write('# Seqid\tGap_identifier\tStatus\tCloserTigID\tCloserRange\tCloserLength\tCloserStrand\tCloserIdentity\n')
        for gapid in gapdict:
            chrom = '.'.join(gapid.split('.')[:-1])
            gap_num = gapid.split('.')[-1]
            if gapid in gapcloserdict:
                info = gapcloserdict[gapid]
                de.write(f'{chrom}\t{gap_num}\tClosed\t{info["sid"]}\t{info["range"]}\t{len(info["seq"])}\t{info["strand"]}\t{info["score"]:.4f}\n')
            else:
                de.write(f'{chrom}\t{gap_num}\tNot_closed\n')
    print(f'[Output] Filling detail written to: {filldetailfile}')
    
    filledfastafile = os.path.join(output_dir, f'{prefix}.genome.filled.fasta')
    filledragpfile = os.path.join(output_dir, f'{prefix}.genome.filled.modified.agp')
    draftgenomedict = _read_fasta_dict(draftgenomefile)
    with open(filledfastafile, 'w') as w, open(filledragpfile, 'w') as ragp:
        for sid, seq in draftgenomedict.items():
            if 'N' * 100 not in seq:
                w.write(f'>{sid}\n{seq}\n')
                ragp.write(f'{sid}\t1\t{len(seq)}\t1\tW\t{sid}\t1\t{len(seq)}\t+\n')
                continue
            
            matches = list(re.finditer(r'N{100,}', seq))
            parts = []
            last_end = 0
            agp_idx = 1
            
            for match in matches:
                if match.start() > last_end:
                    subseq = seq[last_end:match.start()]
                    parts.append(subseq)
                    ragp.write(f'{sid}\t{last_end+1}\t{match.start()}\t{agp_idx}\tW\t{sid}\t{last_end+1}\t{match.start()}\t+\n')
                    agp_idx += 1
                gapid = f'{sid}.{len(parts)//2 + 1}'
                if gapid in gapcloserdict:
                    fill_info = gapcloserdict[gapid]
                    fill_seq = fill_info['seq']
                    parts.append(fill_seq)
                    ragp.write(f'{sid}\t{match.start()+1}\t{match.start()+len(fill_seq)}\t{agp_idx}\tW\t{fill_info["sid"]}\t{fill_info["range"].split("-")[0]}\t{fill_info["range"].split("-")[1]}\t{fill_info["strand"]}\n')
                    agp_idx += 1
                else:
                    gap_seq = seq[match.start():match.end()]
                    parts.append(gap_seq)
                    ragp.write(f'{sid}\t{match.start()+1}\t{match.end()}\t{agp_idx}\tN\t{len(gap_seq)}\tscaffold\tyes\tunspecified\n')
                    agp_idx += 1
                last_end = match.end()
            
            if last_end < len(seq):
                subseq = seq[last_end:]
                parts.append(subseq)
                ragp.write(f'{sid}\t{last_end+1}\t{len(seq)}\t{agp_idx}\tW\t{sid}\t{last_end+1}\t{len(seq)}\t+\n')
            
            newseq = ''.join(parts)
            w.write(f'>{sid}\n{newseq}\n')
    
    print(f'[Output] Filled genome FASTA written to: {filledfastafile}')
    print(f'[Output] Modified AGP written to: {filledragpfile}')
    
    filledstatfile = os.path.join(output_dir, f'{prefix}.genome.filled.stat')
    chrfastadict = _read_fasta_dict(filledfastafile)
    totallen = sum(len(seq) for seq in chrfastadict.values())
    gc = sum(seq.count('G') + seq.count('C') for seq in chrfastadict.values())
    gccontent = gc / totallen if totallen else 0
    
    with open(filledstatfile, 'w') as info:
        info.write(f'# Total Size: {totallen}\n')
        info.write(f'# GC content: {gccontent:.6f}\n')
        info.write('# ChrID\tLength\tGapcount\tGaplocus\n')
        for sid, seq in chrfastadict.items():
            gaps = [r.span() for r in re.finditer(r'N{100,}', seq)]
            if gaps:
                loci = '\t'.join(f'{s+1}-{e}' for s, e in gaps)
                info.write(f'{sid}\t{len(seq)}\t{len(gaps)}\t{loci}\n')
            else:
                info.write(f'{sid}\t{len(seq)}\t0\n')
    print(f'[Output] Genome statistics written to: {filledstatfile}')
    
    if not noplot:
        print('[Info] Plotting functionality has been removed from this version.')
    
    result = {
        'detail_file': filldetailfile,
        'filled_fasta': filledfastafile,
        'modified_agp': filledragpfile,
        'stat_file': filledstatfile,
        'gaps_closed': len(gapcloserdict),
        'gaps_remaining': len(gapdict) - len(gapcloserdict),
        'total_filled_length': totalfilledlen,
        'total_genome_length': totallen,
        'gc_content': gccontent,
        'temp_directory': tmp_dir,
        'output_directory': output_dir
    }
    
    if not keep_temp and os.path.exists(tmp_dir):
        print(f'[Info] Cleaning up temporary directory: {tmp_dir}')
        shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        print(f'[Info] Temporary files are kept in: {tmp_dir}')
    
    print(f'[Info] Script completed successfully!')
    
    return result

def main():
    parser = argparse.ArgumentParser(
        description='Two-stage gap filling tool: contig preprocessing + TGS-GapCloser pre-filling + flanking alignment fine filling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # Complete two-stage process (includes contig preprocessing by default)
  python gap_filler.py -d draft.fasta -g contigs.fasta -o results
  
  # Skip contig preprocessing
  python gap_filler.py -d draft.fasta -g contigs.fasta --no-preprocess
  
  # Only use flanking alignment (skip TGS stage)
  python gap_filler.py -d draft.fasta -g contigs.fasta --no-tgs
        """
    )
    
    parser.add_argument('-d', '--draft-fasta', dest='draft_genome', required=True, 
                       help='Draft genome file (FASTA).')
    parser.add_argument('-g', '--gap-patches', dest='gapcloser_contig', nargs='+', required=True, 
                       help='Contig files for gap closing (FASTA).')
    
    parser.add_argument('--no-tgs', dest='use_tgs', action='store_false', default=True,
                       help='Skip TGS-GapCloser stage, only use flanking alignment.')
    parser.add_argument('--tgs-minmap', dest='tgs_minmap', type=str, default='-x asm5',
                       help='Minimap2 arguments for TGS-GapCloser, default: -x asm5')
    parser.add_argument('--tgs-threads', dest='tgs_threads', type=int, default=32,
                       help='Threads for TGS-GapCloser, default: 32')
    parser.add_argument('--tgs-type', dest='tgs_type', choices=['ont', 'pacbio'], default='ont',
                       help='Read type for TGS-GapCloser, default: ont')
    parser.add_argument('--keep-tgs-temp', dest='keep_tgs_temp', action='store_true',
                       help='Keep TGS-GapCloser temporary files.')
    
    parser.add_argument('--no-preprocess', dest='preprocess_contigs', action='store_false', default=True,
                       help='Skip contig preprocessing (do not split contigs with gaps).')
    parser.add_argument('--min-gap-length', dest='min_gap_length', type=int, default=100,
                       help='Minimum gap length to trigger contig splitting, default: 100')
    parser.add_argument('--min-contig-length', dest='min_contig_length', type=int, default=1000,
                       help='Minimum contig length to keep after preprocessing, default: 1000')
    
    parser.add_argument('-f', '--flank-size', dest='flanking_len', type=int, default=5000, 
                       help='Flanking length (bp), default: 5000')
    parser.add_argument('-l', '--min-length', dest='min_alignment_length', type=int, default=1000, 
                       help='Min alignment length (bp), default: 1000')
    parser.add_argument('-i', '--min-identity', dest='min_alignment_identity', type=float, default=40, 
                       help='Min alignment identity (%%), default: 40')
    parser.add_argument('-m', '--max-fill', dest='max_filling_len', type=int, default=1000000, 
                       help='Max fill length, default: 1000000')
    parser.add_argument('-a', '--aligner', choices=['minimap2', 'unimap'], default='minimap2', 
                       help='Aligner, default: minimap2')
    parser.add_argument('-p', '--prefix', default='quarTeT', 
                       help='Output prefix, default: quarTeT')
    parser.add_argument('-o', '--output-dir', default='.', 
                       help='Output directory, default: current directory')
    parser.add_argument('-t', '--threads', type=int, default=1, 
                       help='Threads for flanking alignment, default: 1')
    
    parser.add_argument('--overwrite', action='store_true', 
                       help='Overwrite existing temporary directory.')
    parser.add_argument('--minimap-options', default='-x asm5', 
                       help='Extra minimap2 options for flanking alignment, default: -x asm5')
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files for debugging.')
    parser.add_argument('--api-mode', action='store_true',
                       help='Use API mode for better error handling and logging.')
    parser.add_argument('--version', action='version', version='GapFiller V2.3.0')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args()
    
    try:
        if args.api_mode:
            filler = GapFillerAPI(verbose=True)
            result = filler.fill_gaps(
                draft_fasta=args.draft_genome,
                gap_patches=args.gapcloser_contig,
                output_dir=args.output_dir,
                output_prefix=args.prefix,
                flank_size=args.flanking_len,
                min_alignment_length=args.min_alignment_length,
                min_identity=args.min_alignment_identity,
                max_fill_length=args.max_filling_len,
                aligner=args.aligner,
                threads=args.threads,
                overwrite=args.overwrite,
                keep_temp=args.keep_temp,
                use_tgs=args.use_tgs,
                tgs_minmap=args.tgs_minmap,
                tgs_threads=args.tgs_threads,
                tgs_type=args.tgs_type,
                keep_tgs_temp=args.keep_tgs_temp,
                preprocess_contigs=args.preprocess_contigs,
                min_gap_length=args.min_gap_length,
                min_contig_length=args.min_contig_length
            )
        else:
            if args.use_tgs:
                print(f"\n{'='*60}")
                print("Stage 1: TGS-GapCloser pre-filling (includes contig preprocessing)")
                print(f"{'='*60}")
                
                filler = GapFillerAPI(verbose=False)
                tgs_result = filler.run_tgsgapcloser(
                    draft_fasta=args.draft_genome,
                    reads_fasta=args.gapcloser_contig,
                    output_dir=args.output_dir,
                    prefix=f"{args.prefix}_tgs",
                    minmap_arg=args.tgs_minmap,
                    threads=args.tgs_threads,
                    tgstype=args.tgs_type,
                    keep_tgs_temp=args.keep_tgs_temp,
                    preprocess_contigs=args.preprocess_contigs,
                    min_gap_length=args.min_gap_length,
                    min_contig_length=args.min_contig_length
                )
                
                print(f"\n{'='*60}")
                print("Stage 2: Flanking alignment fine filling")
                print(f"{'='*60}")
                
                draft_for_flanking = tgs_result["tgs_output_file"]
            else:
                draft_for_flanking = args.draft_genome
            
            result = gap_filler(
                draft_genome=draft_for_flanking,
                gapcloser_contig=args.gapcloser_contig,
                flanking_len=args.flanking_len,
                min_alignment_length=args.min_alignment_length,
                min_alignment_identity=args.min_alignment_identity,
                max_filling_len=args.max_filling_len,
                aligner=args.aligner,
                prefix=args.prefix,
                output_dir=args.output_dir,
                threads=args.threads,
                overwrite=args.overwrite,
                minimapoption=args.minimap_options,
                noplot=True,
                keep_temp=args.keep_temp
            )
            
            if args.use_tgs:
                result["tgs_stage"] = {
                    "tgs_output_file": tgs_result["tgs_output_file"],
                    "gaps_closed_by_tgs": tgs_result["gaps_closed_by_tgs"],
                    "tgs_stats": tgs_result["tgs_stats"],
                    "preprocess_contigs": tgs_result.get("preprocess_contigs", False)
                }
                result["total_gaps_closed"] = tgs_result["gaps_closed_by_tgs"] + result["gaps_closed"]
            else:
                result["tgs_stage"] = None
                result["total_gaps_closed"] = result["gaps_closed"]
        
        print("\n" + "="*60)
        print("Two-stage gap filling completed!")
        print("="*60)
        
        if args.use_tgs and "tgs_stage" in result and result["tgs_stage"]:
            print(f"Stage 1 (TGS-GapCloser):")
            print(f"  Contig preprocessing: {'Yes' if result['tgs_stage'].get('preprocess_contigs', False) else 'No'}")
            print(f"  Filled gaps: {result['tgs_stage']['gaps_closed_by_tgs']}")
            print(f"  Output file: {os.path.basename(result['tgs_stage']['tgs_output_file'])}")
        
        print(f"\nStage 2 (Flanking alignment):")
        print(f"  Filled gaps: {result['gaps_closed']}")
        print(f"  Remaining gaps: {result['gaps_remaining']}")
        print(f"  Total filled length: {result['total_filled_length']:,} bp")
        print(f"  Final genome length: {result['total_genome_length']:,} bp")
        print(f"  GC content: {result['gc_content']:.2%}")
        
        if args.use_tgs and "tgs_stage" in result and result["tgs_stage"]:
            print(f"\nTotal filled gaps: {result['total_gaps_closed']}")
        
        print(f"\nOutput files (in {result['output_directory']}):")
        print(f"  Final filled genome: {os.path.basename(result['filled_fasta'])}")
        print(f"  Filling detail report: {os.path.basename(result['detail_file'])}")
        print(f"  Genome statistics: {os.path.basename(result['stat_file'])}")
        print(f"  AGP file: {os.path.basename(result['modified_agp'])}")
        
        if args.use_tgs and "tgs_stage" in result and result["tgs_stage"]:
            print(f"  TGS output genome: {os.path.basename(result['tgs_stage']['tgs_output_file'])}")
        
        print(f"\n{'='*60}")
        print("Completed!")
        print(f"{'='*60}")
        
    except Exception as e:
        print(f'[Error] {e}', file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()