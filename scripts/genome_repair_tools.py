#!/usr/bin/env python3

import os
import sys
import json
import time
import argparse
import logging
import traceback
import shutil
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from enum import Enum

def import_all_modules():
    modules_status = {
        'assembly_gap_fill': False,
        'gap_complete_enhanced': False,
        'gap_fixer': False,
        'telomere_recovery': False,
        'biopython': False
    }
    
    imported_modules = {}
    
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        if current_dir not in sys.path:
            sys.path.insert(0, current_dir)
        
        try:
            from assembly_gap_fill import assembly_gap_fill_from_dict as run_assembly_gap_fill
            modules_status['assembly_gap_fill'] = True
            imported_modules['run_assembly_gap_fill'] = run_assembly_gap_fill
        except ImportError as e:
            print(f"✗ Cannot import assembly_gap_fill module: {e}")
            imported_modules['run_assembly_gap_fill'] = None
        
        try:
            from gap_complete_controller_module_enhanced import (
                run_gap_complete,
                analyze_strategies,
                fill_specific_chromosome,
                create_controller
            )
            modules_status['gap_complete_enhanced'] = True
            imported_modules['run_gap_complete'] = run_gap_complete
            imported_modules['analyze_strategies'] = analyze_strategies
            imported_modules['fill_specific_chromosome'] = fill_specific_chromosome
            imported_modules['create_controller'] = create_controller
        except ImportError as e:
            print(f"✗ Cannot import gap_complete_controller_module_enhanced module: {e}")
            print("❌ Required module missing: gap_complete_controller_module_enhanced.py")
            imported_modules['run_gap_complete'] = None
            imported_modules['analyze_strategies'] = None
            imported_modules['fill_specific_chromosome'] = None
            imported_modules['create_controller'] = None
        
        try:
            from gap_fixer import GenomeGapFixerAdapter
            modules_status['gap_fixer'] = True
            imported_modules['GenomeGapFixerAdapter'] = GenomeGapFixerAdapter
        except ImportError as e:
            print(f"✗ Cannot import gap_fixer module: {e}")
            imported_modules['GenomeGapFixerAdapter'] = None
        
        try:
            from telomere_recovery_controller import TelomereRecoveryAdapter
            modules_status['telomere_recovery'] = True
            imported_modules['TelomereRecoveryAdapter'] = TelomereRecoveryAdapter
        except ImportError as e:
            print(f"✗ Cannot import telomere_recovery_controller module: {e}")
            imported_modules['TelomereRecoveryAdapter'] = None
        
        try:
            from Bio import SeqIO
            modules_status['biopython'] = True
            imported_modules['SeqIO'] = SeqIO
        except ImportError as e:
            print(f"✗ Cannot import BioPython module: {e}")
            imported_modules['SeqIO'] = None
        
        return imported_modules, modules_status
        
    except Exception as e:
        print(f"Error during module import: {e}")
        return {}, modules_status

MODULES, MODULES_STATUS = import_all_modules()

class Stage(Enum):
    ASSEMBLE_FILL = "assemble_fill"
    PATCH_REPAIR = "patch_repair"
    CORRECT_REFILL = "correct_refill"
    TELOMERE_RECOVER = "telomere_recover"
    
    @classmethod
    def from_number(cls, num: int):
        stages = {
            1: cls.ASSEMBLE_FILL,
            2: cls.PATCH_REPAIR,
            3: cls.CORRECT_REFILL,
            4: cls.TELOMERE_RECOVER
        }
        return stages.get(num, None)

@dataclass
class StageResult:
    stage: Stage
    success: bool
    output_files: Dict[str, str] = field(default_factory=dict)
    statistics: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    start_time: float = 0
    end_time: float = 0
    
    @property
    def duration(self) -> float:
        return self.end_time - self.start_time
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'stage': self.stage.value,
            'success': self.success,
            'output_files': self.output_files,
            'statistics': self.statistics,
            'error': self.error,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'duration': self.duration
        }

def extract_original_chrom_name(seq_id: str) -> str:
    if seq_id.endswith('.fa'):
        seq_id = seq_id[:-3]
    elif seq_id.endswith('.fasta'):
        seq_id = seq_id[:-6]
    
    first_column = seq_id.split()[0] if seq_id.split() else seq_id
    return first_column.strip()

def rename_sequences_in_fasta(input_file: str, output_file: str) -> bool:
    if MODULES.get('SeqIO') is None:
        return False
    
    try:
        records = []
        for record in MODULES['SeqIO'].parse(input_file, "fasta"):
            original_name = extract_original_chrom_name(record.id)
            record.id = original_name
            record.description = ""
            records.append(record)
        
        if records:
            MODULES['SeqIO'].write(records, output_file, "fasta")
            return True
        return False
    except Exception as e:
        print(f"Failed to rename sequences: {e}")
        return False

def setup_logger(name: str, verbose: bool = True, log_file: str = None) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    if logger.hasHandlers():
        logger.handlers.clear()
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def check_file_exists(filepath: str, description: str = "") -> Tuple[bool, str]:
    if not filepath:
        return False, f"{description}: File path is empty"
    
    if not os.path.exists(filepath):
        return False, f"{description}: File does not exist - {filepath}"
    
    if os.path.getsize(filepath) == 0:
        return False, f"{description}: File is empty - {filepath}"
    
    return True, f"{description}: File exists ({os.path.getsize(filepath)} bytes)"

def merge_fasta_files(input_files: List[str], output_file: str, rename_to_original: bool = True) -> bool:
    if MODULES.get('SeqIO') is None:
        return False
    
    try:
        all_records = []
        seen_ids = set()
        
        for input_file in input_files:
            if os.path.exists(input_file):
                for record in MODULES['SeqIO'].parse(input_file, "fasta"):
                    if rename_to_original:
                        original_name = extract_original_chrom_name(record.id)
                        record.id = original_name
                        record.description = ""
                    
                    base_id = record.id
                    counter = 1
                    while record.id in seen_ids:
                        record.id = f"{base_id}_{counter}"
                        counter += 1
                    
                    seen_ids.add(record.id)
                    all_records.append(record)
        
        if all_records:
            MODULES['SeqIO'].write(all_records, output_file, "fasta")
            return True
        
        return False
    except Exception as e:
        print(f"Failed to merge FASTA files: {e}")
        return False

def analyze_genome_stats(fasta_file: str) -> Dict[str, Any]:
    if MODULES.get('SeqIO') is None or not os.path.exists(fasta_file):
        return {}
    
    try:
        stats = {
            'total_chromosomes': 0,
            'total_length': 0,
            'chromosomes': {},
            'has_gaps': [],
            'no_gaps': [],
            'max_length': 0,
            'min_length': float('inf'),
            'avg_length': 0
        }
        
        for record in MODULES['SeqIO'].parse(fasta_file, "fasta"):
            chrom_id = record.id
            seq = str(record.seq)
            length = len(seq)
            
            stats['total_chromosomes'] += 1
            stats['total_length'] += length
            stats['max_length'] = max(stats['max_length'], length)
            stats['min_length'] = min(stats['min_length'], length)
            
            chrom_info = {
                'id': chrom_id,
                'length': length,
                'has_gap': 'N' * 10 in seq.upper(),
                'gap_count': seq.upper().count('N' * 10)
            }
            
            stats['chromosomes'][chrom_id] = chrom_info
            
            if chrom_info['has_gap']:
                stats['has_gaps'].append(chrom_id)
            else:
                stats['no_gaps'].append(chrom_id)
        
        if stats['total_chromosomes'] > 0:
            stats['avg_length'] = stats['total_length'] / stats['total_chromosomes']
        
        return stats
    except Exception as e:
        print(f"Failed to analyze genome statistics: {e}")
        return {}

def extract_single_chromosome(genome_file: str, chrom_id: str, output_file: str) -> bool:
    if MODULES.get('SeqIO') is None:
        return False
    
    try:
        for record in MODULES['SeqIO'].parse(genome_file, "fasta"):
            if record.id == chrom_id:
                MODULES['SeqIO'].write([record], output_file, "fasta")
                return True
        
        for record in MODULES['SeqIO'].parse(genome_file, "fasta"):
            if chrom_id.lower() in record.id.lower():
                MODULES['SeqIO'].write([record], output_file, "fasta")
                return True
        
        return False
        
    except Exception as e:
        print(f"Failed to extract chromosome {chrom_id}: {e}")
        return False

def find_contigs_files_in_directory(directory: str) -> List[str]:
    contigs_files = []
    directory = os.path.abspath(directory)
    
    if not os.path.exists(directory):
        return contigs_files
    
    contigs_patterns = [
        "*contig*.fasta", "*contig*.fa", "*scaffold*.fasta", "*scaffold*.fa",
        "*assembly*.fasta", "*assembly*.fa","assemblies.fasta","shasta.fa","flye.fa","nextdenovo.fa","verkko.fa"
    ]
    
    for pattern in contigs_patterns:
        for file_path in Path(directory).rglob(pattern):
            file_str = str(file_path)
            if os.path.getsize(file_str) > 1000:
                contigs_files.append(file_str)
    
    return list(set(contigs_files))

def prepare_final_contigs_files(generated_contigs: List[str], external_contigs: List[str], 
                               output_dir: str, logger: logging.Logger) -> Optional[str]:
    all_contigs = []
    
    if external_contigs:
        for contig in external_contigs:
            if os.path.exists(contig):
                all_contigs.append(contig)
                logger.info(f"  Added external contigs: {os.path.basename(contig)}")
    
    if generated_contigs:
        for contig in generated_contigs:
            if os.path.exists(contig):
                all_contigs.append(contig)
                logger.info(f"  Added generated contigs: {os.path.basename(contig)}")
    
    merged_dir = os.path.join(output_dir, "merged_contigs")
    os.makedirs(merged_dir, exist_ok=True)
    
    final_contigs_file = os.path.join(merged_dir, "final_contigs.fasta")
    
    # 处理 0、1、多个文件的情况
    if len(all_contigs) == 0:
        logger.warning("⚠️  No contigs files available, creating empty final_contigs.fasta")
        # 创建空的 contigs 文件
        with open(final_contigs_file, 'w') as f:
            f.write("# No contigs available - created as placeholder\n")
        logger.info(f"  Created empty contigs file: {final_contigs_file}")
        
    elif len(all_contigs) == 1:
        shutil.copy2(all_contigs[0], final_contigs_file)
        source_type = "external" if all_contigs[0] in external_contigs else "generated"
        logger.info(f"  Single {source_type} contigs file copied: {os.path.basename(all_contigs[0])}")
    else:
        logger.info(f"  Merging {len(all_contigs)} contigs files...")
        logger.info(f"    - External: {len([c for c in all_contigs if c in external_contigs])} files")
        logger.info(f"    - Generated: {len([c for c in all_contigs if c in generated_contigs])} files")
        
        if merge_fasta_files(all_contigs, final_contigs_file):
            logger.info(f"  ✓ Successfully merged {len(all_contigs)} contigs files")
            logger.info(f"  ✓ Merged contigs will be used for subsequent stages (highest quality)")
        else:
            logger.warning("  Failed to merge contigs files, using first available file")
            shutil.copy2(all_contigs[0], final_contigs_file)

    # 总是返回路径，即使文件可能是空的
    if os.path.exists(final_contigs_file):
        file_size = os.path.getsize(final_contigs_file)
        file_size_mb = file_size / (1024 * 1024)
        status = "empty" if file_size < 100 else "valid"
        logger.info(f"  Final contigs file: {final_contigs_file} ({file_size_mb:.2f} MB, {status})")
        return final_contigs_file
    else:
        logger.error("  Failed to create final contigs file")
        return None

def run_stage1_independently(args):
    logger = setup_logger('Stage1_Independent', verbose=not args.quiet, log_file=getattr(args, 'log_file', None))
    valid, msg = check_file_exists(args.query_fasta, "Input genome file")
    if not valid:
        logger.error(msg)
        return {"status": "error", "error": msg}
    
    has_data = False
    if args.hifi:
        for f in args.hifi:
            if os.path.exists(f):
                has_data = True
                break
    if args.ont:
        for f in args.ont:
            if os.path.exists(f):
                has_data = True
                break
    if args.clr:
        for f in args.clr:
            if os.path.exists(f):
                has_data = True
                break
    if args.contigs_file and os.path.exists(args.contigs_file):
        has_data = True
    
    if not has_data:
        logger.error("No data files provided")
        return {"status": "error", "error": "Please provide at least one data file (--hifi/--ont/--clr/-c)"}
    
    if MODULES.get('run_assembly_gap_fill') is None:
        error = "assembly_gap_fill module not available"
        logger.error(error)
        return {"status": "error", "error": error}
    
    output_dir = os.path.join(args.output_dir, "assemble_fill")
    os.makedirs(output_dir, exist_ok=True)
    
    assembly_tools = []
    if args.run_hifiasm:
        assembly_tools.append('hifiasm')
    if args.run_verkko:
        assembly_tools.append('verkko')
    if args.run_nextdenovo:
        assembly_tools.append('nextdenovo')
    if args.run_flye:
        assembly_tools.append('flye')
    if args.run_shasta:
        assembly_tools.append('shasta')
    
    if not assembly_tools:
        assembly_tools = None
        logger.info("No assembly tools specified, using all available tools")
    else:
        logger.info(f"Using specified assembly tools: {', '.join(assembly_tools)}")
    
    config_dict = {
        'query_fasta': args.query_fasta,
        'output_prefix': os.path.join(output_dir, "assembled_filled"),
        'threads': args.threads,
        'verbose': not args.quiet,
        'use_tgs': args.tgs,
        'assembly_tools': assembly_tools
    }
    
    if args.hifi:
        config_dict['hifi_files'] = [f for f in args.hifi if os.path.exists(f)]
    
    if args.ont:
        config_dict['ont_files'] = [f for f in args.ont if os.path.exists(f)]
    
    if args.clr:
        config_dict['clr_files'] = [f for f in args.clr if os.path.exists(f)]
    
    external_contigs = []
    if args.contigs_file and os.path.exists(args.contigs_file):
        config_dict['contigs_files'] = [args.contigs_file]
        external_contigs.append(args.contigs_file)
    
    logger.info(f"Data usage: HiFi={len(config_dict.get('hifi_files', []))}, "
                f"ONT={len(config_dict.get('ont_files', []))}, "
                f"CLR={len(config_dict.get('clr_files', []))}, "
                f"Contigs={len(external_contigs)}")
    
    try:
        assembly_result = MODULES['run_assembly_gap_fill'](config_dict)
        
        if assembly_result.get('status') != 'success':
            error = assembly_result.get('message', 'Unknown error')
            logger.error(f"Stage 1 failed: {error}")
            return {"status": "error", "error": error}
        
        filled_genome = assembly_result.get('output_file', '')
        if not filled_genome or not os.path.exists(filled_genome):
            error = "Filled genome file not found"
            logger.error(error)
            return {"status": "error", "error": error}
        
        generated_contigs = []
        if 'contigs_file' in assembly_result and assembly_result['contigs_file']:
            if os.path.exists(assembly_result['contigs_file']):
                generated_contigs.append(assembly_result['contigs_file'])
        
        additional_contigs = find_contigs_files_in_directory(output_dir)
        if additional_contigs:
            generated_contigs.extend(additional_contigs)
        
        generated_contigs = list(set(generated_contigs))
        
        logger.info("Preparing final contigs file...")
        final_contigs_file = prepare_final_contigs_files(
            generated_contigs, 
            external_contigs, 
            args.output_dir, 
            logger
        )
        
        contigs_info = {
            "generated_contigs": generated_contigs,
            "external_contigs": external_contigs,
            "final_contigs_file": final_contigs_file
        }
        
        result = {
            "status": "success",
            "output_file": filled_genome,
            "output_dir": output_dir,
            "summary": assembly_result.get('summary', {}),
            "contigs_info": contigs_info,
            "assembly_tools_used": assembly_tools,
            "final_contigs_file": final_contigs_file
        }
        
        logger.info(f"✓ Stage 1 (AssembleFill) completed successfully")
        logger.info(f"   Output file: {filled_genome}")
        logger.info(f"   Assembly tools: {assembly_tools if assembly_tools else 'All available'}")
        
        if generated_contigs:
            logger.info(f"   Generated contigs: {len(generated_contigs)} files")
        if final_contigs_file:
            logger.info(f"   Final contigs file: {final_contigs_file}")
        
        return result
        
    except Exception as e:
        logger.error(f"Stage 1 execution exception: {e}")
        logger.error(traceback.format_exc())
        return {"status": "error", "error": str(e)}

def run_stage2_independently(args):
    logger = setup_logger('Stage2_Independent', verbose=not args.quiet, log_file=getattr(args, 'log_file', None))
    valid, msg = check_file_exists(args.query_fasta, "Input genome file")
    if not valid:
        logger.error(msg)
        return {"status": "error", "error": msg}
    
    final_contigs_file = None
    
    # 1. 最高优先级：merged_contigs/final_contigs.fasta（合并后的最佳版本）
    merged_dir = os.path.join(args.output_dir, "merged_contigs")
    final_contigs_candidate = os.path.join(merged_dir, "final_contigs.fasta")
    
    if os.path.exists(final_contigs_candidate):
        file_size = os.path.getsize(final_contigs_candidate)
        if file_size > 100:  # 有效文件
            final_contigs_file = final_contigs_candidate
            logger.info(f"✓ Using merged contigs file (highest quality): {final_contigs_file}")
    
    # 2. 其次：用户直接提供的 -c 参数文件
    if not final_contigs_file and args.contigs_file and os.path.exists(args.contigs_file):
        final_contigs_file = args.contigs_file
        logger.info(f"✓ Using user-provided contigs file: {final_contigs_file}")
    
    # 3. 两种方式都没有找到
    if not final_contigs_file:
        logger.error("❌ No contigs file found for Stage 2")
        logger.error("Please provide contigs file with -c/--contigs or run Stage 1 first")
        return {"status": "error", "error": "No contigs available for Stage 2"}
    
    # 检查文件是否为空
    file_size = os.path.getsize(final_contigs_file)
    if file_size < 100:  # 小于 100 bytes 认为是空文件
        logger.warning(f"⚠️  Contigs file is very small ({file_size} bytes), may be empty")
        logger.warning("Stage 2 requires valid contigs for repair")
        return {"status": "error", "error": "Contigs file is empty or too small"}
    
    logger.info(f"  File size: {file_size/1024/1024:.2f} MB")
    logger.info(f"  Source: {'merged contigs' if final_contigs_file == final_contigs_candidate else 'user-provided'}")
    
    if MODULES.get('run_gap_complete') is None:
        error = "gap_complete_controller_module_enhanced module not available"
        logger.error(error)
        return {"status": "error", "error": error}
    
    output_dir = os.path.join(args.output_dir, "patch_repair")
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info("Running enhanced GAP repair module...")
    
    try:
        gap_complete_result = MODULES['run_gap_complete'](
            query_fasta=args.query_fasta,
            reference_fasta=final_contigs_file,
            threads=args.threads,
            output_dir=output_dir,
            repair_mode=args.repair_mode,
            min_gap_size=100,
            gap_length=100,
            max_search_distance=500000,
            search_step=100000,
            quiet=args.quiet
        )
        
        if gap_complete_result.get('status') != 'success':
            error = gap_complete_result.get('error', 'Unknown error')
            logger.error(f"Stage 2 failed: {error}")
            return {"status": "error", "error": error}
        
        final_genome = None
        output_dir_path = Path(output_dir)
        possible_files = [
            output_dir_path / "final_genome_module_enhanced.fasta",
            output_dir_path / "final_genome.fasta",
            output_dir_path / "final_genome_module.fasta"
        ]
        
        for file_path in possible_files:
            if file_path.exists():
                final_genome = str(file_path)
                break
        
        if not final_genome:
            for file_path in output_dir_path.glob("**/*.fasta"):
                if "final" in file_path.name.lower() or "complete" in file_path.name.lower():
                    final_genome = str(file_path)
                    break
        
        if not final_genome:
            error = "Repaired genome file not found"
            logger.error(error)
            return {"status": "error", "error": error}
        
        result = {
            "status": "success",
            "output_file": final_genome,
            "output_dir": output_dir,
            "summary": gap_complete_result,
            "method": "enhanced_gap_complete",
            "strategies": gap_complete_result.get('strategies', {}),
            "statistics": gap_complete_result.get('statistics', {}),
            "contigs_file": final_contigs_file
        }
        
        logger.info(f"✓ Stage 2 (PatchRepair) completed successfully")
        logger.info(f"   Output file: {final_genome}")
        logger.info(f"   Contigs used: {final_contigs_file}")
        
        if 'strategies' in gap_complete_result:
            strategies = gap_complete_result['strategies']
            logger.info(f"   Chromosome strategy statistics:")
            for chrom, strategy in strategies.items():
                logger.info(f"     {chrom}: {strategy}")
        
        return result
        
    except Exception as e:
        logger.error(f"Stage 2 execution exception: {e}")
        logger.error(traceback.format_exc())
        return {"status": "error", "error": str(e)}

def run_stage3_independently(args):
    logger = setup_logger('Stage3_Independent', verbose=not args.quiet, log_file=getattr(args, 'log_file', None))
    valid, msg = check_file_exists(args.query_fasta, "Input genome file")
    if not valid:
        logger.error(msg)
        return {"status": "error", "error": msg}
    
    final_contigs_file = None
    
    # 1. 最高优先级：merged_contigs/final_contigs.fasta（合并后的最佳版本）
    merged_dir = os.path.join(args.output_dir, "merged_contigs")
    final_contigs_candidate = os.path.join(merged_dir, "final_contigs.fasta")
    
    if os.path.exists(final_contigs_candidate):
        file_size = os.path.getsize(final_contigs_candidate)
        if file_size > 100:  # 有效文件
            final_contigs_file = final_contigs_candidate
            logger.info(f"✓ Using merged contigs file (highest quality): {final_contigs_file}")
    
    # 2. 其次：用户直接提供的 -c 参数文件
    if not final_contigs_file and args.contigs_file and os.path.exists(args.contigs_file):
        final_contigs_file = args.contigs_file
        logger.info(f"✓ Using user-provided contigs file: {final_contigs_file}")
    
    # 3. 两种方式都没有找到
    if not final_contigs_file:
        logger.error("❌ No contigs file found for Stage 3")
        logger.error("Please provide contigs file with -c/--contigs or run Stage 1 first")
        return {"status": "error", "error": "No contigs available for Stage 3"}
    
    # 检查文件是否为空
    file_size = os.path.getsize(final_contigs_file)
    if file_size < 100:  # 小于 100 bytes 认为是空文件
        logger.warning(f"⚠️  Contigs file is very small ({file_size} bytes), may be empty")
        logger.warning("Stage 3 requires valid contigs for error correction")
    
    logger.info(f"  File size: {file_size/1024/1024:.2f} MB")
    logger.info(f"  Source: {'merged contigs' if final_contigs_file == final_contigs_candidate else 'user-provided'}")

    
    if file_size < 100 and MODULES.get('GenomeGapFixerAdapter') is not None:
        logger.warning("No contigs available for error correction, proceeding to refill only")
        return run_final_fill_independently(args, logger, args.query_fasta, final_contigs_file)
    
    genome_stats = analyze_genome_stats(args.query_fasta)
    chroms_with_gaps = genome_stats.get('has_gaps', [])
    
    if not chroms_with_gaps:
        logger.info("No chromosomes with GAPs, skipping error correction")
        return run_final_fill_independently(args, logger, args.query_fasta, final_contigs_file)
    
    if MODULES.get('GenomeGapFixerAdapter') is None:
        logger.warning("gap_fixer module not available, skipping error correction")
        return run_final_fill_independently(args, logger, args.query_fasta, final_contigs_file)
    
    output_dir = os.path.join(args.output_dir, "correct_refill")
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"Starting error correction for {len(chroms_with_gaps)} chromosomes with GAPs")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Using contigs: {final_contigs_file}")
    
    fixed_chrom_files = []
    
    for chrom_id in chroms_with_gaps:
        logger.info(f"Processing chromosome: {chrom_id}")
        
        chrom_file = os.path.join(output_dir, f"{chrom_id}_input.fa")
        if not extract_single_chromosome(args.query_fasta, chrom_id, chrom_file):
            logger.warning(f"Cannot extract chromosome {chrom_id}, skipping")
            continue
        
        config_dict = {
            'query_fasta': chrom_file,
            'reference_fasta': final_contigs_file,
            'threads': args.threads,
            'output_fasta': os.path.join(output_dir, f"{chrom_id}_fixed.fa"),
            'work_dir': os.path.join(output_dir, f"work_{chrom_id}"),
            'repair_mode': args.repair_mode,
            'verbose': not args.quiet
        }
        
        logger.debug(f"Chromosome {chrom_id} config: threads={config_dict['threads']}")
        
        try:
            fixer = MODULES['GenomeGapFixerAdapter'](config_dict)
            fixer_result = fixer.run()
            
            if fixer_result.get('status') == 'success':
                fixed_file = fixer_result.get('output_files', {}).get('final_genome', '')
                if fixed_file and os.path.exists(fixed_file):
                    fixed_chrom_files.append(fixed_file)
                    logger.info(f"  ✓ Chromosome {chrom_id} correction successful")
                else:
                    if os.path.exists(chrom_file):
                        fixed_chrom_files.append(chrom_file)
                        logger.warning(f"  ⚠ Fixed file not found for {chrom_id}, using original")
            else:
                if os.path.exists(chrom_file):
                    fixed_chrom_files.append(chrom_file)
                logger.warning(f"  ⚠ Chromosome {chrom_id} correction failed: {fixer_result.get('error')}")
                    
        except Exception as e:
            if os.path.exists(chrom_file):
                fixed_chrom_files.append(chrom_file)
            logger.warning(f"  ⚠ Chromosome {chrom_id} correction exception: {e}")
    
    if not fixed_chrom_files:
        logger.warning("No chromosomes successfully corrected, using original file")
        corrected_genome = args.query_fasta
    else:
        corrected_genome = os.path.join(output_dir, "corrected_genome.fasta")
        all_chrom_files = []
        
        all_chrom_files.extend(fixed_chrom_files)
        
        for chrom_id in genome_stats.get('no_gaps', []):
            chrom_file = os.path.join(output_dir, f"{chrom_id}_no_gap.fa")
            if extract_single_chromosome(args.query_fasta, chrom_id, chrom_file):
                all_chrom_files.append(chrom_file)
        
        if not merge_fasta_files(all_chrom_files, corrected_genome):
            logger.warning("Failed to merge corrected chromosomes, using original file")
            corrected_genome = args.query_fasta
    
    logger.info("Starting final refill...")
    return run_final_fill_independently(args, logger, corrected_genome, final_contigs_file)

def run_final_fill_independently(args, logger, input_genome: str, contigs_file: Optional[str] = None):
    if MODULES.get('run_assembly_gap_fill') is None:
        logger.warning("assembly_gap_fill module not available, skipping fill")
        return {
            "status": "success",
            "output_file": input_genome,
            "message": "Fill module not available, using original file",
            "skipped": True
        }
    
    output_dir = os.path.join(args.output_dir, "correct_refill", "final_fill")
    os.makedirs(output_dir, exist_ok=True)
    
    config_dict = {
        'query_fasta': input_genome,
        'output_prefix': os.path.join(output_dir, "final_filled"),
        'threads': args.threads,
        'verbose': not args.quiet,
        'use_tgs': args.tgs,
        'assembly_tools': None,
        'skip_assembly': True
    }
    
    if contigs_file and os.path.exists(contigs_file):
        config_dict['contigs_files'] = [contigs_file]
        logger.info(f"Final refill using contigs file: {contigs_file}")
    else:
        logger.warning("No contigs available for final refill, using input file")
        return {
            "status": "success",
            "output_file": input_genome,
            "message": "No contigs available, using original file",
            "skipped": True
        }
    
    logger.info("Final refill: Using existing contigs file (no re-assembly)")
    
    try:
        fill_result = MODULES['run_assembly_gap_fill'](config_dict)
        
        if fill_result.get('status') != 'success':
            logger.warning(f"Final refill failed, using corrected file: {fill_result.get('message')}")
            return {
                "status": "success",
                "output_file": input_genome,
                "message": "Refill failed, using corrected file",
                "skipped": True
            }
        
        filled_genome = fill_result.get('output_file', '')
        if not filled_genome or not os.path.exists(filled_genome):
            logger.warning("Final refill file not found, using corrected file")
            return {
                "status": "success",
                "output_file": input_genome,
                "message": "Refill file not found, using corrected file",
                "skipped": True
            }
        
        renamed_genome = os.path.join(output_dir, "final_genome_renamed.fasta")
        if rename_sequences_in_fasta(filled_genome, renamed_genome):
            filled_genome = renamed_genome
            logger.info("✓ Sequences renamed to original chromosome names")
        
        result = {
            "status": "success",
            "output_file": filled_genome,
            "output_dir": output_dir,
            "summary": {
                "correction_statistics": {
                    "chromosomes_with_gaps": len(analyze_genome_stats(args.query_fasta).get('has_gaps', [])),
                    "chromosomes_processed": len(analyze_genome_stats(input_genome).get('has_gaps', []))
                },
                "fill_statistics": fill_result.get('summary', {}),
                "assembly_tools_used": None,
                "fill_mode": "contigs_only"
            },
            "contigs_file": contigs_file,
            "assembly_tools_used": None,
            "fill_mode": "contigs_only"
        }
        
        logger.info(f"✓ Stage 3 (CorrectRefill) completed successfully")
        logger.info(f"   Output file: {filled_genome}")
        logger.info(f"   Fill mode: Contigs only (no re-assembly)")
        return result
        
    except Exception as e:
        logger.error(f"Final refill execution exception: {e}")
        logger.error(traceback.format_exc())
        return {
            "status": "success",
            "output_file": input_genome,
            "message": f"Refill exception, using corrected file: {str(e)}",
            "skipped": True
        }

def run_stage4_independently(args):
    logger = setup_logger('Stage4_Independent', verbose=not args.quiet, log_file=getattr(args, 'log_file', None))
    valid, msg = check_file_exists(args.query_fasta, "Input genome file")
    if not valid:
        logger.error(msg)
        return {"status": "error", "error": msg}
    
    final_contigs_file = None
    
    # 1. 最高优先级：merged_contigs/final_contigs.fasta（合并后的最佳版本）
    merged_dir = os.path.join(args.output_dir, "merged_contigs")
    final_contigs_candidate = os.path.join(merged_dir, "final_contigs.fasta")
    
    if os.path.exists(final_contigs_candidate):
        file_size = os.path.getsize(final_contigs_candidate)
        if file_size > 100:  # 有效文件
            final_contigs_file = final_contigs_candidate
            logger.info(f"✓ Using merged contigs file (highest quality): {final_contigs_file}")
    
    # 2. 其次：用户直接提供的 -c 参数文件
    if not final_contigs_file and args.contigs_file and os.path.exists(args.contigs_file):
        final_contigs_file = args.contigs_file
        logger.info(f"✓ Using user-provided contigs file: {final_contigs_file}")
    
    # 3. 两种方式都没有找到
    if not final_contigs_file:
        logger.error("❌ No contigs file found for Stage 4")
        logger.error("Please provide contigs file with -c/--contigs or run Stage 1 first")
        return {"status": "error", "error": "No contigs available for Stage 4"}
    
    # 检查文件是否为空
    file_size = os.path.getsize(final_contigs_file)
    if file_size < 100:  # 小于 100 bytes 认为是空文件
        logger.warning(f"⚠️  Contigs file is very small ({file_size} bytes), may be empty")
        logger.warning("Stage 4 requires valid contigs for telomere recovery")
        return {"status": "error", "error": "Contigs file is empty or too small"}
    
    logger.info(f"  File size: {file_size/1024/1024:.2f} MB")
    logger.info(f"  Source: {'merged contigs' if final_contigs_file == final_contigs_candidate else 'user-provided'}")    

    if MODULES.get('TelomereRecoveryAdapter') is None:
        error = "telomere_recovery_controller module not available"
        logger.error(error)
        return {"status": "error", "error": error}
    
    output_dir = os.path.join(args.output_dir, "telomere_recover")
    input_dir = os.path.join(output_dir, "input")
    intermediate_dir = os.path.join(output_dir, "intermediate")
    output_file_dir = os.path.join(output_dir, "output")
    
    for dir_path in [output_dir, input_dir, intermediate_dir, output_file_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    input_copy = os.path.join(input_dir, "complete_genome.fasta")
    shutil.copy2(args.query_fasta, input_copy)
    
    config_dict = {
        'query_fasta': input_copy,
        'contigs_file': final_contigs_file,
        'output_file': os.path.join(output_file_dir, "final_genome.fasta"),
        'threads': args.threads,
        'verbose': not args.quiet
    }
    
    logger.info("Running telomere recovery module...")
    
    try:
        recovery = MODULES['TelomereRecoveryAdapter'](config_dict)
        recovery_result = recovery.run()
        
        if recovery_result.get('status') != 'success':
            error = recovery_result.get('error', recovery_result.get('message', 'Unknown error'))
            logger.error(f"Stage 4 failed: {error}")
            return {"status": "error", "error": error}
        
        final_genome = recovery_result.get('output_files', {}).get('final_genome', '')
        if not final_genome or not os.path.exists(final_genome):
            error = "Telomere recovered genome file not found"
            logger.error(error)
            return {"status": "error", "error": error}
        
        renamed_genome = os.path.join(output_file_dir, "final_genome_renamed.fasta")
        if rename_sequences_in_fasta(final_genome, renamed_genome):
            final_genome = renamed_genome
            logger.info("✓ Sequences renamed to original chromosome names")
        
        result = {
            "status": "success",
            "output_file": final_genome,
            "output_dir": output_file_dir,
            "summary": recovery_result,
            "contigs_file": final_contigs_file
        }
        
        logger.info(f"✓ Stage 4 (TelomereRecover) completed successfully")
        logger.info(f"   Output file: {final_genome}")
        logger.info(f"   Contigs used: {final_contigs_file}")
        return result
        
    except Exception as e:
        logger.error(f"Stage 4 execution exception: {e}")
        logger.error(traceback.format_exc())
        return {"status": "error", "error": str(e)}

class GenomeRepairController:
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.logger = setup_logger(
            'GenomeRepairTools',
            verbose=config.get('verbose', True),
            log_file=config.get('log_file')
        )
        
        self.work_dir = os.path.abspath(config.get('work_dir', 'genomerepair_results'))
        self.stage_results: Dict[Stage, StageResult] = {}
        self.current_stage: Optional[Stage] = None
        
        self.external_contigs: List[str] = []
        self.final_contigs_file: Optional[str] = None
        
    # 保存用户直接提供的 contigs 文件路径
        if 'contigs_file' in self.config and self.config['contigs_file']:
            if os.path.exists(self.config['contigs_file']):
                self.user_contigs_file = self.config['contigs_file']  # 新增
                self.external_contigs.append(self.config['contigs_file'])
                self.logger.info(f"User-provided contigs file: {self.config['contigs_file']}")
        
        os.makedirs(self.work_dir, exist_ok=True)
        
        self.logger.info("=" * 80)
        self.logger.info("Genome Repair Tools v2.5 with Simplified Contigs Management")
        self.logger.info("=" * 80)
        
        if not MODULES_STATUS.get('gap_complete_enhanced'):
            self.logger.error("❌ Required module gap_complete_controller_module_enhanced not available")
            self.logger.error("Please ensure gap_complete_controller_module_enhanced.py is in the same directory")
            raise ImportError("Required module not available")
    
    def prepare_final_contigs_file(self, generated_contigs: List[str]) -> Optional[str]:
        all_contigs = []
        
        if self.external_contigs:
            all_contigs.extend(self.external_contigs)
            self.logger.info(f"Adding {len(self.external_contigs)} external contigs files")
        
        if generated_contigs:
            all_contigs.extend(generated_contigs)
            self.logger.info(f"Adding {len(generated_contigs)} generated contigs files")
        
        merged_dir = os.path.join(self.work_dir, "merged_contigs")
        os.makedirs(merged_dir, exist_ok=True)
        
        final_contigs_file = os.path.join(merged_dir, "final_contigs.fasta")
        
        # 处理 0、1、多个文件的情况
        if len(all_contigs) == 0:
            self.logger.warning("⚠️  No contigs files available, creating empty final_contigs.fasta")
            # 创建空的 contigs 文件
            with open(final_contigs_file, 'w') as f:
                f.write("# No contigs available - created as placeholder\n")
            self.logger.info(f"  Created empty contigs file: {final_contigs_file}")
            
        elif len(all_contigs) == 1:
            shutil.copy2(all_contigs[0], final_contigs_file)
            self.logger.info(f"  Single contigs file copied: {os.path.basename(all_contigs[0])}")
        else:
            self.logger.info(f"  Merging {len(all_contigs)} contigs files...")
            if merge_fasta_files(all_contigs, final_contigs_file):
                self.logger.info(f"  ✓ Successfully merged {len(all_contigs)} contigs files")
            else:
                self.logger.warning("  Failed to merge contigs files, using first available file")
                shutil.copy2(all_contigs[0], final_contigs_file)
        
        if os.path.exists(final_contigs_file) and os.path.getsize(final_contigs_file) > 0:
            file_size_mb = os.path.getsize(final_contigs_file) / (1024 * 1024)
            status = "empty" if os.path.getsize(final_contigs_file) < 100 else "valid"
            self.logger.info(f"  Final contigs file: {final_contigs_file} ({file_size_mb:.2f} MB, {status})")
            self.final_contigs_file = final_contigs_file
            return final_contigs_file
        else:
            self.logger.error("  Failed to create final contigs file")
            return None
    
    def get_contigs_for_stage(self, stage_num: int) -> Optional[str]:
        if stage_num == 1:
            return None
    
    # 1. 最高优先级：Stage 1 合并生成的 final_contigs.fasta
        if self.final_contigs_file and os.path.exists(self.final_contigs_file):
            file_size = os.path.getsize(self.final_contigs_file)
            if file_size > 100:  # 有效文件
                self.logger.info(f"Stage {stage_num}: Using merged contigs file (highest quality): {self.final_contigs_file}")
                return self.final_contigs_file
            else:
                self.logger.warning(f"Stage {stage_num}: Merged contigs file is empty/small ({file_size} bytes)")
    
    # 2. 其次：merged_contigs 目录下的 final_contigs.fasta
        merged_dir = os.path.join(self.work_dir, "merged_contigs")
        final_contigs_candidate = os.path.join(merged_dir, "final_contigs.fasta")
    
        if os.path.exists(final_contigs_candidate):
            file_size = os.path.getsize(final_contigs_candidate)
            if file_size > 100:  # 有效文件
                self.logger.info(f"Stage {stage_num}: Using contigs from merged_contigs directory")
                return final_contigs_candidate
    
    # 3. 最后备选：用户直接提供的原始 contigs 文件
        if self.user_contigs_file and os.path.exists(self.user_contigs_file):
            file_size = os.path.getsize(self.user_contigs_file)
            if file_size > 100:  # 有效文件
                self.logger.info(f"Stage {stage_num}: Using user-provided contigs file (fallback): {self.user_contigs_file}")
                return self.user_contigs_file
            else:
                self.logger.warning(f"Stage {stage_num}: User-provided contigs file is empty/small ({file_size} bytes)")
    
    # 所有方式都失败了
        self.logger.warning(f"Stage {stage_num}: No valid contigs file available")
        return None
    
    def validate_config(self) -> Tuple[bool, str]:
        self.logger.info("Validating configuration parameters...")
        
        if 'query_fasta' not in self.config:
            return False, "Input genome file must be specified"
        
        valid, msg = check_file_exists(self.config['query_fasta'], "Input genome file")
        if not valid:
            return False, msg
        self.logger.info(f"✓ {msg}")
        
        has_data = False
        data_types = []
        
        if 'hifi_files' in self.config and self.config['hifi_files']:
            has_data = True
            data_types.append(f"HiFi({len(self.config['hifi_files'])} files)")
        
        if 'ont_files' in self.config and self.config['ont_files']:
            has_data = True
            data_types.append(f"ONT({len(self.config['ont_files'])} files)")
        
        if 'clr_files' in self.config and self.config['clr_files']:
            has_data = True
            data_types.append(f"CLR({len(self.config['clr_files'])} files)")
        
        if self.external_contigs:
            has_data = True
            data_types.append(f"Contigs({len(self.external_contigs)} files)")
        
        if not has_data:
            return False, "Please provide at least one data file (sequencing data or contigs file)"
        
        self.logger.info(f"✓ Data sources: {', '.join(data_types)}")
        
        threads = self.config.get('threads', 16)
        if threads <= 0:
            return False, "Thread count must be greater than 0"
        
        self.logger.info(f"✓ Configuration validated - using {threads} threads")
        return True, ""
    
    def run_stage(self, stage_num: int, input_file: str = None) -> StageResult:
        stage = Stage.from_number(stage_num)
        if not stage:
            return StageResult(
                stage=Stage(f"stage{stage_num}"),
                success=False,
                error=f"Invalid stage number: {stage_num}"
            )
        
        self.current_stage = stage
        stage_func = getattr(self, f"_run_{stage.value}", None)
        
        if not stage_func:
            return StageResult(
                stage=stage,
                success=False,
                error=f"Stage not implemented: {stage.value}"
            )
        
        self.logger.info(f"\n{'='*80}")
        self.logger.info(f"Running Stage {stage_num}: {stage.value}")
        self.logger.info(f"{'='*80}")
        
        if input_file is None:
            prev_stage = Stage.from_number(stage_num - 1)
            if prev_stage and prev_stage in self.stage_results:
                prev_result = self.stage_results[prev_stage]
                if 'output_genome' in prev_result.output_files:
                    input_file = prev_result.output_files['output_genome']
            else:
                input_file = self.config.get('query_fasta')
        
        start_time = time.time()
        result = stage_func(input_file)
        result.start_time = start_time
        result.end_time = time.time()
        
        self.stage_results[stage] = result
        
        if result.success:
            self.logger.info(f"✓ Stage {stage_num} completed successfully (duration: {result.duration:.1f}s)")
        else:
            self.logger.error(f"✗ Stage {stage_num} failed: {result.error}")
        
        return result
    
    def _run_assemble_fill(self, input_file: str) -> StageResult:
        result = StageResult(stage=Stage.ASSEMBLE_FILL, success=False)
        
        if MODULES.get('run_assembly_gap_fill') is None:
            result.error = "assembly_gap_fill module not available, skipping stage 1"
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {'skipped': True, 'reason': 'Module not available'}
            return result
        
        output_dir = os.path.join(self.work_dir, "assemble_fill")
        os.makedirs(output_dir, exist_ok=True)
        
        config_dict = {
            'query_fasta': input_file,
            'output_prefix': os.path.join(output_dir, "assembled_filled"),
            'threads': self.config.get('threads', 16),
            'verbose': self.config.get('verbose', True),
            'use_tgs': self.config.get('use_tgs', False),
            'assembly_tools': self.config.get('assembly_tools', None)
        }
        
        if 'hifi_files' in self.config and self.config['hifi_files']:
            config_dict['hifi_files'] = self.config['hifi_files']
        
        if 'ont_files' in self.config and self.config['ont_files']:
            config_dict['ont_files'] = self.config['ont_files']
        
        if 'clr_files' in self.config and self.config['clr_files']:
            config_dict['clr_files'] = self.config['clr_files']
        
        if self.external_contigs:
            config_dict['contigs_files'] = self.external_contigs
        
        data_sources = []
        if 'hifi_files' in config_dict:
            data_sources.append(f"HiFi({len(config_dict['hifi_files'])} files)")
        if 'ont_files' in config_dict:
            data_sources.append(f"ONT({len(config_dict['ont_files'])} files)")
        if 'clr_files' in config_dict:
            data_sources.append(f"CLR({len(config_dict['clr_files'])} files)")
        if 'contigs_files' in config_dict:
            data_sources.append(f"Contigs({len(config_dict['contigs_files'])} files)")
        
        assembly_tools = config_dict.get('assembly_tools', None)
        if assembly_tools:
            self.logger.info(f"Stage 1 assembly tools: {', '.join(assembly_tools)}")
        else:
            self.logger.info("Stage 1 assembly tools: All available tools (default)")
        
        self.logger.info(f"Stage 1 data sources: {', '.join(data_sources)}")
        
        try:
            assembly_result = MODULES['run_assembly_gap_fill'](config_dict)
            
            if assembly_result.get('status') != 'success':
                self.logger.warning(f"Stage 1 failed: {assembly_result.get('message')}, using original file")
                result.success = True
                result.output_files = {'output_genome': input_file}
                result.statistics = {'skipped': True, 'reason': 'Assembly failed'}
                return result
            
            filled_genome = assembly_result.get('output_file', '')
            if not filled_genome or not os.path.exists(filled_genome):
                self.logger.warning("Filled genome file not found, using original file")
                result.success = True
                result.output_files = {'output_genome': input_file}
                result.statistics = {'skipped': True, 'reason': 'Output file not found'}
                return result
            
            generated_contigs = []
            
            if 'contigs_file' in assembly_result and assembly_result['contigs_file']:
                if os.path.exists(assembly_result['contigs_file']):
                    generated_contigs.append(assembly_result['contigs_file'])
            
            output_contigs = find_contigs_files_in_directory(output_dir)
            if output_contigs:
                generated_contigs.extend(output_contigs)
            
            generated_contigs = list(set(generated_contigs))
            
            self.logger.info("Creating final contigs file...")
            final_contigs_file = self.prepare_final_contigs_file(generated_contigs)
            
            renamed_genome = os.path.join(output_dir, "assembled_filled_renamed.fasta")
            if rename_sequences_in_fasta(filled_genome, renamed_genome):
                filled_genome = renamed_genome
                self.logger.info("✓ Sequences renamed to original chromosome names")
            
            result.success = True
            result.output_files = {
                'output_genome': filled_genome,
                'log_file': os.path.join(output_dir, "assemble_fill.log")
            }
            
            summary = assembly_result.get('summary', {})
            summary['generated_contigs_count'] = len(generated_contigs)
            summary['assembly_tools_used'] = assembly_tools
            summary['final_contigs_file'] = final_contigs_file
            result.statistics.update(summary)
            
            if generated_contigs:
                self.logger.info(f"✓ Generated {len(generated_contigs)} contigs files")
            if final_contigs_file:
                self.logger.info(f"✓ Final contigs file created: {final_contigs_file}")
            
        except Exception as e:
            self.logger.error(f"Stage 1 execution exception: {str(e)}, using original file")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {'skipped': True, 'reason': f'Exception: {str(e)}'}
        
        return result
    
    def _run_patch_repair(self, input_file: str) -> StageResult:
        result = StageResult(stage=Stage.PATCH_REPAIR, success=False)
        
        contigs_file = self.get_contigs_for_stage(2)
        
        if not contigs_file:
            self.logger.warning("Stage 2 requires contigs file, skipping")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {'skipped': True, 'reason': 'No contigs available'}
            return result
        
        if MODULES.get('run_gap_complete') is None:
            result.error = "gap_complete_controller_module_enhanced module not available"
            return result
        
        output_dir = os.path.join(self.work_dir, "patch_repair")
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info("Running enhanced GAP repair...")
        self.logger.info("Using intelligent chromosome strategies")
        self.logger.info(f"Contigs file: {contigs_file}")
        
        try:
            gap_complete_result = MODULES['run_gap_complete'](
                query_fasta=input_file,
                reference_fasta=contigs_file,
                threads=self.config.get('threads', 16),
                output_dir=output_dir,
                repair_mode=self.config.get('repair_mode', 'aggressive'),
                min_gap_size=self.config.get('min_gap_size', 100),
                quiet=not self.config.get('verbose', True)
            )
            
            if gap_complete_result.get('status') != 'success':
                result.error = f"Patch repair failed: {gap_complete_result.get('error', 'Unknown error')}"
                return result
            
            final_genome = None
            output_dir_path = Path(output_dir)
            
            possible_files = [
                output_dir_path / "final_genome_module_enhanced.fasta",
                output_dir_path / "final_genome.fasta",
                output_dir_path / "final_genome_module.fasta"
            ]
            
            for file_path in possible_files:
                if file_path.exists():
                    final_genome = str(file_path)
                    break
            
            if not final_genome:
                for file_path in output_dir_path.glob("**/*.fasta"):
                    if "final" in file_path.name.lower() or "complete" in file_path.name.lower():
                        final_genome = str(file_path)
                        break
            
            if not final_genome:
                result.error = "Repaired genome file not found"
                return result
            
            renamed_genome = os.path.join(output_dir, "final_genome_renamed.fasta")
            if rename_sequences_in_fasta(final_genome, renamed_genome):
                final_genome = renamed_genome
                self.logger.info("✓ Sequences renamed to original chromosome names")
            
            result.success = True
            result.output_files = {
                'output_genome': final_genome,
                'log_file': os.path.join(output_dir, "patch_repair.log")
            }
            
            stats = {
                'method': 'enhanced_gap_complete',
                'strategies': gap_complete_result.get('strategies', {}),
                'categories': gap_complete_result.get('categories', {}),
                'output_dir': output_dir,
                'contigs_file': contigs_file,
                'contigs_source': 'final_contigs' if contigs_file == self.final_contigs_file else 'external_contigs'
            }
            
            if 'strategies' in gap_complete_result:
                strategies = gap_complete_result['strategies']
                strategy_counts = {}
                for strategy in strategies.values():
                    strategy_counts[strategy] = strategy_counts.get(strategy, 0) + 1
                stats['strategy_counts'] = strategy_counts
            
            result.statistics = stats
            
            self.logger.info("Chromosome strategy statistics:")
            if 'strategies' in gap_complete_result:
                for chrom, strategy in gap_complete_result['strategies'].items():
                    self.logger.info(f"  {chrom}: {strategy}")
            
        except Exception as e:
            result.error = f"Stage 2 execution exception: {str(e)}"
        
        return result
    
    def _run_correct_refill(self, input_file: str) -> StageResult:
        result = StageResult(stage=Stage.CORRECT_REFILL, success=False)
        
        contigs_file = self.get_contigs_for_stage(3)
        
        if not contigs_file:
            self.logger.warning("Stage 3 has no contigs files, skipping")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {'skipped': True, 'reason': 'No contigs files'}
            return result
        
        if MODULES.get('GenomeGapFixerAdapter') is None:
            self.logger.warning("gap_fixer module not available, skipping error correction")
            return self._run_final_fill_only(input_file, result, contigs_file)
        
        genome_stats = analyze_genome_stats(input_file)
        chroms_with_gaps = genome_stats.get('has_gaps', [])
        
        if not chroms_with_gaps:
            self.logger.info("No chromosomes with GAPs, skipping error correction")
            self.logger.info("Proceeding to final refill...")
            return self._run_final_fill_only(input_file, result, contigs_file)
        
        output_dir = os.path.join(self.work_dir, "correct_refill")
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info(f"Starting error correction for {len(chroms_with_gaps)} chromosomes with GAPs")
        self.logger.info(f"Configured threads: {self.config.get('threads', 16)}")
        self.logger.info(f"Using contigs: {contigs_file}")
        
        fixed_chrom_files = []
        
        for chrom_id in chroms_with_gaps:
            self.logger.info(f"Processing chromosome: {chrom_id}")
            
            chrom_file = os.path.join(output_dir, f"{chrom_id}_input.fa")
            if not extract_single_chromosome(input_file, chrom_id, chrom_file):
                self.logger.warning(f"Cannot extract chromosome {chrom_id}, skipping")
                continue
            
            config_dict = {
                'query_fasta': chrom_file,
                'reference_fasta': contigs_file,
                'threads': self.config.get('threads', 16),
                'output_fasta': os.path.join(output_dir, f"{chrom_id}_fixed.fa"),
                'work_dir': os.path.join(output_dir, f"work_{chrom_id}"),
                'repair_mode': self.config.get('repair_mode', 'aggressive'),
                'verbose': self.config.get('verbose', True)
            }
            
            self.logger.debug(f"Chromosome {chrom_id} config: threads={config_dict['threads']}")
            
            try:
                fixer = MODULES['GenomeGapFixerAdapter'](config_dict)
                fixer_result = fixer.run()
                
                if fixer_result.get('status') == 'success':
                    fixed_file = fixer_result.get('output_files', {}).get('final_genome', '')
                    if fixed_file and os.path.exists(fixed_file):
                        fixed_chrom_files.append(fixed_file)
                        self.logger.info(f"  ✓ Chromosome {chrom_id} correction successful")
                    else:
                        if os.path.exists(chrom_file):
                            fixed_chrom_files.append(chrom_file)
                            self.logger.warning(f"  ⚠ Fixed file not found for {chrom_id}, using original")
                else:
                    if os.path.exists(chrom_file):
                        fixed_chrom_files.append(chrom_file)
                    self.logger.warning(f"  ⚠ Chromosome {chrom_id} correction failed: {fixer_result.get('error')}")
                    
            except Exception as e:
                if os.path.exists(chrom_file):
                    fixed_chrom_files.append(chrom_file)
                self.logger.warning(f"  ⚠ Chromosome {chrom_id} correction exception: {e}")
        
        if not fixed_chrom_files:
            self.logger.warning("No chromosomes successfully corrected, proceeding to final refill")
            corrected_genome = input_file
            correction_stats = {
                'message': 'No chromosomes corrected', 
                'skipped': True,
                'contigs_file': contigs_file
            }
        else:
            corrected_genome = os.path.join(output_dir, "corrected_genome.fasta")
            all_chrom_files = []
            
            all_chrom_files.extend(fixed_chrom_files)
            
            for chrom_id in genome_stats.get('no_gaps', []):
                chrom_file = os.path.join(output_dir, f"{chrom_id}_no_gap.fa")
                if extract_single_chromosome(input_file, chrom_id, chrom_file):
                    all_chrom_files.append(chrom_file)
            
            if merge_fasta_files(all_chrom_files, corrected_genome):
                correction_stats = {
                    'chromosomes_with_gaps': len(chroms_with_gaps),
                    'chromosomes_processed': len(fixed_chrom_files),
                    'chromosomes_without_gaps': len(genome_stats.get('no_gaps', [])),
                    'threads_used': self.config.get('threads', 16),
                    'contigs_file': contigs_file
                }
            else:
                self.logger.warning("Failed to merge corrected chromosomes, using original file for refill")
                corrected_genome = input_file
                correction_stats = {
                    'message': 'Merge failed, using original file', 
                    'skipped': True,
                    'contigs_file': contigs_file
                }
        
        self.logger.info("Starting final refill...")
        return self._run_final_fill_integrated(corrected_genome, result, contigs_file, correction_stats)
    
    def _run_final_fill_only(self, input_file: str, result: StageResult, contigs_file: Optional[str]) -> StageResult:
        if MODULES.get('run_assembly_gap_fill') is None:
            self.logger.warning("assembly_gap_fill module not available, skipping fill")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {'skipped': True, 'reason': 'Fill module not available'}
            return result
        
        return self._run_final_fill_integrated(input_file, result, contigs_file, {'correction_skipped': True})
    
    def _run_final_fill_integrated(self, input_file: str, result: StageResult, 
                                  contigs_file: Optional[str], correction_stats: Dict) -> StageResult:
        if not contigs_file:
            self.logger.warning("No contigs available for final refill, using input file")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {
                'correction_stats': correction_stats,
                'fill_status': 'no_contigs',
                'message': 'No contigs available for refill'
            }
            return result
        
        output_dir = os.path.join(self.work_dir, "correct_refill", "final_fill")
        os.makedirs(output_dir, exist_ok=True)
        
        config_dict = {
            'query_fasta': input_file,
            'output_prefix': os.path.join(output_dir, "final_filled"),
            'threads': self.config.get('threads', 16),
            'verbose': self.config.get('verbose', True),
            'use_tgs': self.config.get('use_tgs', False),
            'assembly_tools': None,
            'skip_assembly': True,
            'contigs_files': [contigs_file]
        }
        
        self.logger.info(f"Final fill using contigs file: {os.path.basename(contigs_file)}")
        self.logger.info("Fill mode: Contigs only (no re-assembly)")
        
        try:
            fill_result = MODULES['run_assembly_gap_fill'](config_dict)
            
            if fill_result.get('status') != 'success':
                self.logger.warning(f"Final refill failed: {fill_result.get('message')}, using input file")
                result.success = True
                result.output_files = {'output_genome': input_file}
                result.statistics = {
                    'correction_stats': correction_stats,
                    'fill_status': 'failed',
                    'fill_message': fill_result.get('message'),
                    'contigs_file': contigs_file,
                    'fill_mode': 'contigs_only'
                }
                return result
            
            filled_genome = fill_result.get('output_file', '')
            if not filled_genome or not os.path.exists(filled_genome):
                self.logger.warning("Final refill file not found, using input file")
                result.success = True
                result.output_files = {'output_genome': input_file}
                result.statistics = {
                    'correction_stats': correction_stats,
                    'fill_status': 'file_not_found',
                    'contigs_file': contigs_file,
                    'fill_mode': 'contigs_only'
                }
                return result
            
            renamed_genome = os.path.join(output_dir, "final_genome_renamed.fasta")
            if rename_sequences_in_fasta(filled_genome, renamed_genome):
                filled_genome = renamed_genome
                self.logger.info("✓ Sequences renamed to original chromosome names")
            
            result.success = True
            result.output_files = {
                'output_genome': filled_genome,
                'log_file': os.path.join(output_dir, "final_fill.log")
            }
            result.statistics = {
                'correction_stats': correction_stats,
                'fill_stats': fill_result.get('summary', {}),
                'contigs_file': contigs_file,
                'fill_mode': 'contigs_only'
            }
            
            self.logger.info(f"✓ Final refill completed (contigs only)")
            
        except Exception as e:
            self.logger.error(f"Final refill execution exception: {str(e)}, using input file")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {
                'correction_stats': correction_stats,
                'fill_status': 'exception',
                'fill_error': str(e),
                'contigs_file': contigs_file,
                'fill_mode': 'contigs_only'
            }
        
        return result
    
    def _run_telomere_recover(self, input_file: str) -> StageResult:
        result = StageResult(stage=Stage.TELOMERE_RECOVER, success=False)
        
        contigs_file = self.get_contigs_for_stage(4)
        
        if not contigs_file:
            self.logger.warning("Stage 4 requires contigs file, skipping")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {'skipped': True, 'reason': 'No contigs available'}
            return result
        
        if MODULES.get('TelomereRecoveryAdapter') is None:
            self.logger.warning("telomere_recovery_controller module not available, skipping stage 4")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {'skipped': True, 'reason': 'Module not available'}
            return result
        
        output_dir = os.path.join(self.work_dir, "telomere_recover")
        input_dir = os.path.join(output_dir, "input")
        intermediate_dir = os.path.join(output_dir, "intermediate")
        output_file_dir = os.path.join(output_dir, "output")
        
        for dir_path in [output_dir, input_dir, intermediate_dir, output_file_dir]:
            os.makedirs(dir_path, exist_ok=True)
        
        input_copy = os.path.join(input_dir, "complete_genome.fasta")
        shutil.copy2(input_file, input_copy)
        
        config_dict = {
            'query_fasta': input_copy,
            'contigs_file': contigs_file,
            'output_file': os.path.join(output_file_dir, "final_genome.fasta"),
            'threads': self.config.get('threads', 16),
            'verbose': self.config.get('verbose', True)
        }
        
        self.logger.info(f"Using contigs: {contigs_file}")
        
        try:
            recovery = MODULES['TelomereRecoveryAdapter'](config_dict)
            recovery_result = recovery.run()
            
            if recovery_result.get('status') != 'success':
                self.logger.warning(f"Telomere recovery failed: {recovery_result.get('error')}, using original file")
                result.success = True
                result.output_files = {'output_genome': input_file}
                result.statistics = {
                    'message': 'Telomere recovery failed, using original file', 
                    'skipped': True,
                    'contigs_file': contigs_file
                }
                return result
            
            final_genome = recovery_result.get('output_files', {}).get('final_genome', '')
            if not final_genome or not os.path.exists(final_genome):
                self.logger.warning("Telomere recovered file not found, using original file")
                result.success = True
                result.output_files = {'output_genome': input_file}
                result.statistics = {
                    'message': 'Output file not found, using original file', 
                    'skipped': True,
                    'contigs_file': contigs_file
                }
                return result
            
            renamed_genome = os.path.join(output_file_dir, "final_genome_renamed.fasta")
            if rename_sequences_in_fasta(final_genome, renamed_genome):
                final_genome = renamed_genome
                self.logger.info("✓ Sequences renamed to original chromosome names")
            
            result.success = True
            result.output_files = {
                'output_genome': final_genome,
                'log_file': os.path.join(output_dir, "telomere_recover.log")
            }
            
            recovery_stats = recovery_result.copy()
            recovery_stats['contigs_file'] = contigs_file
            result.statistics = recovery_stats
            
        except Exception as e:
            self.logger.error(f"Stage 4 execution exception: {str(e)}, using original file")
            result.success = True
            result.output_files = {'output_genome': input_file}
            result.statistics = {
                'message': f'Exception: {str(e)}', 
                'skipped': True,
                'contigs_file': contigs_file
            }
        
        return result
    
    def run_pipeline(self, stages_to_run: List[int] = None) -> Dict[str, Any]:
        start_time = time.time()
        
        if stages_to_run is None:
            stages_to_run = [1, 2, 3, 4]
        
        self.logger.info(f"Starting stages: {stages_to_run}")
        self.logger.info(f"Total threads configured: {self.config.get('threads', 16)}")
        
        assembly_tools = self.config.get('assembly_tools', None)
        if assembly_tools:
            self.logger.info(f"Assembly tools specified: {', '.join(assembly_tools)}")
        else:
            self.logger.info("No assembly tools specified, using all available tools by default")
        
        if self.external_contigs:
            self.logger.info(f"External contigs files: {len(self.external_contigs)}")
        
        is_valid, error_msg = self.validate_config()
        if not is_valid:
            return {"status": "error", "error": error_msg}
        
        for stage_num in stages_to_run:
            result = self.run_stage(stage_num)
            
            if not result.success and stage_num < max(stages_to_run):
                self.logger.warning(f"Stage {stage_num} failed, but continuing with remaining stages")
        
        total_time = time.time() - start_time
        report = self._generate_report(total_time)
        
        self.logger.info("\n" + "=" * 80)
        self.logger.info("🎉 Pipeline execution completed!")
        self.logger.info("=" * 80)
        
        return report
    
    def _generate_report(self, total_time: float) -> Dict[str, Any]:
        report = {
            "status": "success",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "total_time_seconds": round(total_time, 2),
            "work_dir": self.work_dir,
            "contigs_info": {
                "external_contigs": self.external_contigs,
                "final_contigs_file": self.final_contigs_file
            },
            "stage_results": {
                stage.value: result.to_dict() 
                for stage, result in self.stage_results.items()
            },
            "parameters": self.config
        }
        
        report_file = os.path.join(self.work_dir, "genomerepair_report.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        self.logger.info(f"Report saved: {report_file}")
        
        self._print_summary(report)
        
        return report
    
    def _print_summary(self, report: Dict[str, Any]):
        self.logger.info("\n📊 Execution Summary")
        self.logger.info("-" * 40)
        self.logger.info(f"Total runtime: {report['total_time_seconds']:.1f} seconds")
        
        contigs_info = report.get('contigs_info', {})
        external_count = len(contigs_info.get('external_contigs', []))
        final_contigs_file = contigs_info.get('final_contigs_file')
        
        self.logger.info(f"Contigs sources: {external_count} external files")
        if final_contigs_file:
            self.logger.info(f"Final contigs file: {os.path.basename(final_contigs_file)}")
        
        assembly_tools = report.get('parameters', {}).get('assembly_tools', None)
        if assembly_tools:
            self.logger.info(f"Assembly tools used: {', '.join(assembly_tools)} (Stage 1 only)")
        else:
            self.logger.info("Assembly tools used: All available tools (Stage 1 only)")
        
        success_count = 0
        skipped_count = 0
        for stage_name, result in report.get('stage_results', {}).items():
            if result.get('success'):
                success_count += 1
                if result.get('statistics', {}).get('skipped'):
                    skipped_count += 1
        self.logger.info(f"Successful stages: {success_count}/{len(report.get('stage_results', {}))}")
        if skipped_count > 0:
            self.logger.info(f"Skipped stages: {skipped_count}")
        
        for stage_name, result in report.get('stage_results', {}).items():
            status = "✓ Success" if result.get('success') else "✗ Failed"
            duration = result.get('duration', 0)
            
            if result.get('statistics', {}).get('skipped'):
                status = "⏭ Skipped"
                reason = result.get('statistics', {}).get('reason', result.get('statistics', {}).get('message', ''))
                if reason:
                    status += f" ({reason})"
            
            if stage_name == 'correct_refill' and result.get('success'):
                stats = result.get('statistics', {})
                if 'correction_stats' in stats:
                    correction = stats['correction_stats']
                    if not correction.get('skipped'):
                        threads_used = correction.get('threads_used', 'Not recorded')
                        chroms_processed = correction.get('chromosomes_processed', 0)
                        contigs_file = correction.get('contigs_file', 'Not specified')
                        
                        self.logger.info(f"  {stage_name}: {status} ({duration:.1f}s)")
                        self.logger.info(f"    - Error correction: {chroms_processed} chromosomes, {threads_used} threads")
                        if contigs_file and contigs_file != 'Not specified':
                            self.logger.info(f"    - Contigs used: {os.path.basename(contigs_file)}")
                    else:
                        self.logger.info(f"  {stage_name}: {status} ({duration:.1f}s) - Error correction skipped")
                        
                if 'fill_stats' in stats:
                    fill_mode = stats.get('fill_mode', 'contigs_only')
                    contigs_file = stats.get('contigs_file', 'Not specified')
                    self.logger.info(f"    - Final refill: Completed, using contigs file: {os.path.basename(contigs_file)}")
                    self.logger.info(f"    - Fill mode: {fill_mode}")
            else:
                stats = result.get('statistics', {})
                contigs_info_str = ""
                
                if 'contigs_file' in stats and stats['contigs_file']:
                    contigs_file = stats['contigs_file']
                    if isinstance(contigs_file, str):
                        contigs_info_str = f", contigs: {os.path.basename(contigs_file)}"
                
                assembly_info_str = ""
                if 'assembly_tools_used' in stats and stats['assembly_tools_used']:
                    assembly_info_str = f", assembly tools: {', '.join(stats['assembly_tools_used'])}"
                
                self.logger.info(f"  {stage_name}: {status} ({duration:.1f}s){contigs_info_str}{assembly_info_str}")
        
        final_stage = None
        for stage_name in ['telomere_recover', 'correct_refill', 'patch_repair', 'assemble_fill']:
            if stage_name in report.get('stage_results', {}):
                stage_result = report['stage_results'][stage_name]
                if stage_result.get('success') and not stage_result.get('statistics', {}).get('skipped'):
                    final_stage = stage_name
                    break
        
        if final_stage:
            final_result = report['stage_results'][final_stage]
            output_files = final_result.get('output_files', {})
            if 'output_genome' in output_files:
                self.logger.info(f"\n📁 Final output file: {output_files['output_genome']}")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Genome Repair Tools v2.5 with Simplified Contigs Management",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument('-q', '--query', dest='query_fasta', required=True,
                       help='Genome FASTA file to repair')
    
    parser.add_argument('--hifi', '--pacbio-hifi', action='append', default=[],
                       help='PacBio HiFi sequencing data file')
    parser.add_argument('--ont', '--ont-ul', action='append', default=[],
                       help='Oxford Nanopore sequencing data file')
    parser.add_argument('--clr', '--pacbio-clr', action='append', default=[],
                       help='PacBio CLR sequencing data file')
    parser.add_argument('-c', '--contigs', dest='contigs_file',
                       help='Reference contigs file (optional, will be merged with generated contigs)')
    
    parser.add_argument('-t', '--threads', type=int, default=16,
                       help='Number of threads (default: 16)')
    parser.add_argument('-o', '--output-dir', default='genomerepair_results',
                       help='Output directory (default: genomerepair_results)')
    parser.add_argument('--stages', default='1,2,3,4',
                       help='Stages to run, comma-separated (default: 1,2,3,4)')
    parser.add_argument('--tgs', action='store_true',
                       help='Enable TGS-GapCloser mode')
    parser.add_argument('--repair-mode', choices=['conservative', 'aggressive'],
                       default='aggressive', help='Repair mode (default: aggressive)')
    parser.add_argument('--quiet', action='store_true',
                       help='Quiet mode (less output)')
    parser.add_argument('--log-file', dest='log_file',
                       help='Log file path (optional)')
    parser.add_argument('--run-stage-1', action='store_true',
                       help='Run Stage 1 independently (AssembleFill)')
    parser.add_argument('--run-stage-2', action='store_true',
                       help='Run Stage 2 independently (PatchRepair)')
    parser.add_argument('--run-stage-3', action='store_true',
                       help='Run Stage 3 independently (CorrectRefill)')
    parser.add_argument('--run-stage-4', action='store_true',
                       help='Run Stage 4 independently (TelomereRecover)')
    
    parser.add_argument('--run-hifiasm', action='store_true',
                       help='Run hifiasm assembly tool (for HiFi data)')
    parser.add_argument('--run-verkko', action='store_true',
                       help='Run verkko assembly tool (for telomere-to-telomere assembly)')
    parser.add_argument('--run-nextdenovo', action='store_true',
                       help='Run nextDenovo assembly tool')
    parser.add_argument('--run-flye', action='store_true',
                       help='Run flye assembly tool (for long reads)')
    parser.add_argument('--run-shasta', action='store_true',
                       help='Run shasta assembly tool (for ONT data)')
    
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    if not os.path.exists(args.query_fasta):
        print(f"\n❌ Error: Input file does not exist {args.query_fasta}")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print("🧬 Genome Repair Tools v2.5 with Simplified Contigs Management")
    print(f"{'='*80}")
    
    print("\n📊 Input Data Statistics:")
    print(f"  Input genome: {args.query_fasta}")
    
    data_counts = {
        'HiFi': len(args.hifi),
        'ONT': len(args.ont), 
        'CLR': len(args.clr),
        'Contigs': 1 if args.contigs_file else 0
    }
    
    has_data = False
    for data_type, count in data_counts.items():
        if count > 0:
            print(f"  {data_type}: {count} file(s)")
            has_data = True
    
    assembly_tools = []
    if args.run_hifiasm:
        assembly_tools.append('hifiasm')
    if args.run_verkko:
        assembly_tools.append('verkko')
    if args.run_nextdenovo:
        assembly_tools.append('nextdenovo')
    if args.run_flye:
        assembly_tools.append('flye')
    if args.run_shasta:
        assembly_tools.append('shasta')
    
    if assembly_tools:
        print(f"  Assembly tools: {', '.join(assembly_tools)} (Stage 1 only)")
    else:
        print(f"  Assembly tools: All available (Stage 1 only)")
    
    if not has_data:
        print(f"\n❌ Error: No data files provided")
        print("Please provide at least one of the following:")
        print("  --hifi    HiFi sequencing data")
        print("  --ont     ONT sequencing data")
        print("  --clr     CLR sequencing data")
        print("  -c/--contigs   Reference contigs file (optional)")
        sys.exit(1)
    
    config = {
        'verbose': not args.quiet,
        'threads': args.threads,
        'use_tgs': args.tgs,
        'repair_mode': args.repair_mode,
        'work_dir': os.path.abspath(args.output_dir),
        'query_fasta': os.path.abspath(args.query_fasta),
    }
    
    if assembly_tools:
        config['assembly_tools'] = assembly_tools
    
    if args.hifi:
        config['hifi_files'] = [os.path.abspath(f) for f in args.hifi if os.path.exists(f)]
        if len(config['hifi_files']) != len(args.hifi):
            print(f"⚠  Warning: {len(args.hifi) - len(config['hifi_files'])} HiFi file(s) do not exist")
    
    if args.ont:
        config['ont_files'] = [os.path.abspath(f) for f in args.ont if os.path.exists(f)]
        if len(config['ont_files']) != len(args.ont):
            print(f"⚠  Warning: {len(args.ont) - len(config['ont_files'])} ONT file(s) do not exist")
    
    if args.clr:
        config['clr_files'] = [os.path.abspath(f) for f in args.clr if os.path.exists(f)]
        if len(config['clr_files']) != len(args.clr):
            print(f"⚠  Warning: {len(args.clr) - len(config['clr_files'])} CLR file(s) do not exist")
    
    if args.contigs_file:
        if os.path.exists(args.contigs_file):
            config['contigs_file'] = os.path.abspath(args.contigs_file)
        else:
            print(f"⚠  Warning: Contigs file does not exist {args.contigs_file}")
    
    try:
        stages_to_run = [int(stage.strip()) for stage in args.stages.split(',')]
    except:
        stages_to_run = [1, 2, 3, 4]
        print(f"⚠  Note: Stage parameter parsing failed, using default: {stages_to_run}")
    
    print(f"\n⚙️  Run Configuration:")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Threads: {args.threads}")
    print(f"  Stages to run: {stages_to_run}")
    print(f"  TGS mode: {'Enabled' if args.tgs else 'Disabled'}")
    print(f"  Repair mode: {args.repair_mode}")
    if assembly_tools:
        print(f"  Assembly tools: {', '.join(assembly_tools)} (Stage 1 only)")
    else:
        print(f"  Assembly tools: All available tools (Stage 1 only)")
    print(f"  Contigs management: Simplified (Stage 1 creates final contigs file)")
    print(f"  Stage 3 mode: Contigs only (no re-assembly)")
    print(f"  Sequence renaming: Enabled (complex names to simple ChrN format)")
    print(f"{'='*80}")
    
    independent_modes = [
        ('run_stage_1', args.run_stage_1, run_stage1_independently),
        ('run_stage_2', args.run_stage_2, run_stage2_independently),
        ('run_stage_3', args.run_stage_3, run_stage3_independently),
        ('run_stage_4', args.run_stage_4, run_stage4_independently),
    ]
    
    specified_stages = [(name.replace('_', ' ').title().replace('Run Stage', 'Stage'), func) 
                       for name, flag, func in independent_modes if flag]
    
    if specified_stages:
        print(f"\n🚀 Independent Execution Mode:")
        results = []
        
        for stage_name, stage_func in specified_stages:
            print(f"\n>>> Running {stage_name}...")
            result = stage_func(args)
            results.append(result)
            
            if result.get('status') != 'success':
                print(f"❌ {stage_name} failed: {result.get('error', 'Unknown error')}")
                sys.exit(1)
        
        if results:
            last_result = results[-1]
            print(f"\n{'✅'*3} Independent execution completed")
            if 'output_file' in last_result:
                print(f"  Final output file: {last_result['output_file']}")
            if 'final_contigs_file' in last_result:
                print(f"  Final contigs file: {last_result['final_contigs_file']}")
            if 'assembly_tools_used' in last_result and last_result['assembly_tools_used']:
                print(f"  Assembly tools: {', '.join(last_result['assembly_tools_used'])}")
            if 'fill_mode' in last_result and last_result['fill_mode']:
                print(f"  Fill mode: {last_result['fill_mode']}")
            sys.exit(0)
        
    else:
        try:
            controller = GenomeRepairController(config)
            result = controller.run_pipeline(stages_to_run)
            
            if result.get('status') == 'success':
                print(f"\n{'✅'*3} Pipeline execution completed!")
                print(f"  Output directory: {args.output_dir}")
                
                contigs_info = result.get('contigs_info', {})
                external_count = len(contigs_info.get('external_contigs', []))
                final_contigs_file = contigs_info.get('final_contigs_file')
                
                if external_count > 0:
                    print(f"  External contigs: {external_count} file(s)")
                if final_contigs_file:
                    print(f"  Final contigs file: {os.path.basename(final_contigs_file)}")
                
                if assembly_tools:
                    print(f"  Assembly tools used: {', '.join(assembly_tools)} (Stage 1 only)")
                else:
                    print(f"  Assembly tools used: All available tools (Stage 1 only)")
                
                if 'correct_refill' in result.get('stage_results', {}):
                    stage3_stats = result['stage_results']['correct_refill'].get('statistics', {})
                    if 'fill_mode' in stage3_stats:
                        print(f"  Stage 3 fill mode: {stage3_stats['fill_mode']}")
                
                final_stage = None
                for stage_name in ['telomere_recover', 'correct_refill', 'patch_repair', 'assemble_fill']:
                    if stage_name in result.get('stage_results', {}):
                        stage_result = result['stage_results'][stage_name]
                        if stage_result.get('success') and not stage_result.get('statistics', {}).get('skipped'):
                            final_stage = stage_name
                            break
                
                if final_stage:
                    final_output = result['stage_results'][final_stage].get('output_files', {}).get('output_genome')
                    if final_output:
                        print(f"  Final output: {os.path.basename(final_output)} (with original chromosome names)")
                
                sys.exit(0)
            else:
                print(f"\n❌ Failed: {result.get('error', 'Unknown error')}")
                sys.exit(1)
                
        except ImportError as e:
            print(f"\n❌ Required module unavailable: {e}")
            print("Please ensure gap_complete_controller_module_enhanced.py is in the same directory")
            sys.exit(1)
        except KeyboardInterrupt:
            print("\n\n⚠️ User interrupted operation")
            sys.exit(1)
        except Exception as e:
            print(f"\n❌ Execution failed: {e}")
            traceback.print_exc()
            sys.exit(1)

if __name__ == "__main__":
    main()
