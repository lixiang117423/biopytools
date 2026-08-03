"""Phytophthora CRN效应子鉴定模块|Phytophthora CRN Effector Identification Module"""

import os
import re
from datetime import datetime
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Set

from .config import PhytoEffectorConfig
from .utils import (
    build_conda_command,
    generate_software_versions, is_step_completed,
    merge_fasta_files, parse_domtblout_details,
    parse_fasta_to_dict, run_command, run_signalp_pipeline,
)


class CRNFinder:
    """Phytophthora CRN效应子鉴定器|Phytophthora CRN Effector Finder"""

    LFLAK_PATTERN = re.compile(r'LFLAK')

    def __init__(self, config: PhytoEffectorConfig, logger=None, signalp_results: Dict = None):
        self.config = config
        self.logger = logger
        self.signalp_results: Dict[str, Dict] = signalp_results if signalp_results is not None else {}
        self.source_map: Dict[str, str] = {}

    def find_effectors(self) -> bool:
        """运行CRN效应子鉴定流程|Run CRN effector identification pipeline"""
        start_time = datetime.now()
        self.logger.info("=" * 60)
        self.logger.info("开始CRN效应子鉴定|Starting CRN effector identification")
        self.logger.info("=" * 60)

        try:
            # 步骤0: 准备输入|Step 0: Prepare input
            self._prepare_input()

            # 步骤1: SignalP信号肽预测|Step 1: SignalP prediction (已由main.py统一运行|already run by main.py)
            if not self.signalp_results and not self.config.skip_signalp:
                self._run_signalp()
            elif not self.signalp_results:
                self.logger.info("跳过SignalP预测|Skipping SignalP prediction")

            # 步骤2: hmmsearch搜索|Step 2: hmmsearch
            hmm_details = self._run_hmmsearch()

            # 步骤3: LFLAK基序注释和输出|Step 3: LFLAK motif annotation and output
            self._validate_and_output(hmm_details)

            # 记录运行时间|Record runtime
            elapsed = (datetime.now() - start_time).total_seconds()
            self.logger.info(f"CRN效应子鉴定完成|CRN effector identification completed in {elapsed:.1f}秒|seconds")
            self.logger.info("=" * 60)

            self._generate_versions(start_time)
            return True

        except Exception as e:
            self.logger.error(f"CRN效应子鉴定失败|CRN effector identification failed: {e}")
            raise

    def _prepare_input(self):
        """准备输入文件|Prepare input files"""
        self.logger.info("步骤0: 准备输入文件|Step 0: Preparing input files")
        self.logger.info(f"输入文件数|Input file count: {len(self.config._input_files)}")

        if len(self.config._input_files) == 1:
            self.config._combined_fasta = self.config._input_files[0]
            stem = Path(self.config._input_files[0]).stem
            self.source_map = {sid: stem for sid, _ in self._prepare_gen()}
        else:
            total, self.source_map = merge_fasta_files(
                self.config._input_files, self.config._combined_fasta
            )
            self.logger.info(f"合并{total}条序列到|merged {total} sequences to {self.config._combined_fasta}")

    def _prepare_gen(self):
        """生成器辅助函数|Generator helper"""
        from .utils import parse_fasta
        yield from parse_fasta(self.config._combined_fasta)

    def _run_signalp(self):
        """运行SignalP预测|Run SignalP prediction"""
        step_dir = os.path.join(self.config.output_dir, '01_signalp')
        summary_file = os.path.join(step_dir, 'prediction_results.txt')
        self.signalp_results = run_signalp_pipeline(
            self.config, self.logger, step_dir, summary_file
        )

    def _run_hmmsearch(self) -> List[Dict]:
        """运行hmmsearch|Run hmmsearch"""
        step_dir = os.path.join(self.config.output_dir, '02_hmmsearch')
        domtblout = os.path.join(step_dir, 'crn.domtblout')

        if is_step_completed(domtblout):
            self.logger.info("跳过已完成步骤|Skipping completed step: hmmsearch CRN搜索|CRN hmmsearch")
            details = parse_domtblout_details(domtblout)
            self.logger.info(f"加载hmmsearch结果|Loaded hmmsearch results: {len(details)}条记录|records")
            return details

        self.logger.info("步骤2: hmmsearch搜索(CRN HMM)|Step 2: hmmsearch (CRN HMM, default parameters)")
        os.makedirs(step_dir, exist_ok=True)

        cmd = build_conda_command(self.config.hmmsearch_path, [
            '--domtblout', domtblout,
            '--cpu', str(self.config.threads),
            self.config.crn_hmm,
            self.config._combined_fasta,
        ])

        success, stdout, stderr = run_command(cmd, self.logger, "hmmsearch CRN HMM搜索|hmmsearch CRN HMM search")
        if not success:
            self.logger.error("hmmsearch运行失败|hmmsearch failed")
            return []

        details = parse_domtblout_details(domtblout)
        unique_targets = len(set(d['target'] for d in details))
        self.logger.info(f"hmmsearch命中|hmmsearch hits: {unique_targets}个唯一序列|unique sequences, {len(details)}条domain记录|domain records")
        return details

    def _validate_and_output(self, hmm_details: List[Dict]):
        """LFLAK基序注释并生成输出|Annotate LFLAK motif and generate output"""
        step_dir = os.path.join(self.config.output_dir, '03_candidates')
        candidates_file = os.path.join(step_dir, 'crn_candidates.tsv')

        self.logger.info("步骤3: LFLAK基序注释|Step 3: LFLAK motif annotation")
        os.makedirs(step_dir, exist_ok=True)

        seq_dict = parse_fasta_to_dict(self.config._combined_fasta)

        best_hits = {}
        for d in hmm_details:
            tid = d['target']
            if tid not in best_hits or d['domain_score'] > best_hits[tid]['domain_score']:
                best_hits[tid] = d

        rows = []
        for seq_id in sorted(best_hits.keys()):
            sequence = seq_dict.get(seq_id, '')
            if not sequence:
                continue

            detail = best_hits[seq_id]

            lflak_matches = list(self.LFLAK_PATTERN.finditer(sequence))
            has_lflak = len(lflak_matches) > 0
            lflak_positions = ';'.join(str(m.start() + 1) for m in lflak_matches) if lflak_matches else 'None'

            sp_info = self.signalp_results.get(seq_id, {})
            has_sp = sp_info.get('has_signal_peptide', None)
            sp_pos = sp_info.get('cs_position', 'N/A')

            rows.append({
                'Effector_Type': 'CRN',
                'Sample': self.source_map.get(seq_id, 'Unknown'),
                'Sequence_ID': seq_id,
                'Sequence': sequence,
                'SignalP': 'Yes' if has_sp else ('No' if has_sp is not None else 'N/A'),
                'SP_Cleavage_Site': sp_pos if has_sp else ('-' if has_sp is not None else 'N/A'),
                'HMM': 'Yes',
                'HMM_Evalue': f"{detail['domain_evalue']:.2e}",
                'HMM_Score': f"{detail['domain_score']:.1f}",
                'LFLAK': 'Yes' if has_lflak else 'No',
                'LFLAK_Position': lflak_positions,
            })

        df = pd.DataFrame(rows)
        df.to_csv(candidates_file, sep='\t', index=False)
        self.logger.info(f"候选CRN效应子已保存|CRN candidates saved: {candidates_file}")

        total = len(df)
        lflak_count = len(df[df['LFLAK'] == 'Yes'])
        sp_yes = len(df[df['SignalP'] == 'Yes'])
        self.logger.info(f"总候选数|Total candidates: {total}")
        self.logger.info(f"LFLAK=Yes: {lflak_count}, SignalP=Yes: {sp_yes}")

    def _generate_versions(self, start_time):
        """生成版本信息文件|Generate software versions file"""
        tools = {
            'signalp': self.config.signalp_path,
            'hmmsearch': self.config.hmmsearch_path,
        }
        params = {
            'mode': 'crn',
            'organism': self.config.organism,
            'skip_signalp': self.config.skip_signalp,
            'threads': self.config.threads,
            'crn_hmm': self.config.crn_hmm,
            'lflak_motif': 'LFLAK (annotation only)',
        }
        generate_software_versions(self.config.output_dir, tools, params, start_time)
