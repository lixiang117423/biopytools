"""Phytophthora通用效应子鉴定模块(基于HMM)|Phytophthora Generic Effector Finder (HMM-based)"""

import os
import re
from datetime import datetime
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Pattern, Set

from .config import PhytoEffectorConfig
from .utils import (
    PhytoEffectorLogger, build_conda_command,
    format_number, generate_software_versions, is_step_completed,
    merge_fasta_files, parse_domtblout_details, parse_domtblout_hits,
    parse_fasta_to_dict, parse_signalp_output, run_command,
)


# 各效应子类型的配置|Effector type configurations
EFFECTOR_TYPE_CONFIG = {
    'nlp': {
        'label': 'NLP',
        'hmm_attr': 'nlp_hmm',
        'default_hmm': 'paper_nlp.hmm',
        'motif': 'GHRHDWE',
        'motif_regex': re.compile(r'G.H.{0,2}D.W[ED]'),
        'motif_label': 'GHRHDWE',
        'output_dir': '04_nlp',
        'candidates_file': 'nlp_candidates.tsv',
        'score_threshold': None,
    },
    'protease': {
        'label': 'Protease',
        'hmm_attr': 'protease_hmm',
        'default_hmm': 'paper_protease.hmm',
        'motif': None,
        'motif_regex': None,
        'motif_label': None,
        'output_dir': '05_protease',
        'candidates_file': 'protease_candidates.tsv',
        'score_threshold': None,
    },
    'scp': {
        'label': 'SCP',
        'hmm_attr': 'scp_hmm',
        'default_hmm': 'paper_scp.hmm',
        'motif': None,
        'motif_regex': None,
        'motif_label': None,
        'output_dir': '06_scp',
        'candidates_file': 'scp_candidates.tsv',
        'score_threshold': 15.0,
    },
    'elicitin': {
        'label': 'elicitin',
        'hmm_attr': 'elicitin_hmm',
        'default_hmm': 'paper_elicitin.hmm',
        'motif': None,
        'motif_regex': None,
        'motif_label': None,
        'output_dir': '07_elicitin',
        'candidates_file': 'elicitin_candidates.tsv',
        'score_threshold': 20.0,
    },
    'yxsl': {
        'label': 'YxSL',
        'hmm_attr': 'yxsl_hmm',
        'default_hmm': 'paper_yxsl.hmm',
        'motif': None,
        'motif_regex': None,
        'motif_label': None,
        'output_dir': '08_yxsl',
        'candidates_file': 'yxsl_candidates.tsv',
        'score_threshold': 25.0,
    },
}


class GenericEffectorFinder:
    """通用HMM效应子鉴定器|Generic HMM-based Effector Finder"""

    def __init__(self, effector_type: str, config: PhytoEffectorConfig,
                 logger=None, signalp_results: Dict = None):
        if effector_type not in EFFECTOR_TYPE_CONFIG:
            raise ValueError(f"不支持的效应子类型|Unsupported effector type: {effector_type}")
        self.effector_type = effector_type
        self.type_config = EFFECTOR_TYPE_CONFIG[effector_type]
        self.config = config
        self.logger = logger
        self.signalp_results: Dict[str, Dict] = signalp_results if signalp_results is not None else {}
        self.source_map: Dict[str, str] = {}

    @property
    def hmm_file(self) -> str:
        """获取HMM文件路径|Get HMM file path"""
        return getattr(self.config, self.type_config['hmm_attr'], '')

    @property
    def label(self) -> str:
        return self.type_config['label']

    def find_effectors(self) -> bool:
        """运行效应子鉴定流程|Run effector identification pipeline"""
        start_time = datetime.now()
        label = self.label
        self.logger.info("=" * 60)
        self.logger.info(f"开始{label}效应子鉴定|Starting {label} effector identification")
        self.logger.info("=" * 60)

        try:
            self._prepare_input()
            self._validate_and_output()

            elapsed = (datetime.now() - start_time).total_seconds()
            self.logger.info(f"{label}效应子鉴定完成|{label} effector identification completed in {elapsed:.1f}秒|seconds")
            self.logger.info("=" * 60)

            self._generate_versions(start_time)
            return True

        except Exception as e:
            self.logger.error(f"{label}效应子鉴定失败|{label} effector identification failed: {e}")
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
        from .utils import parse_fasta
        yield from parse_fasta(self.config._combined_fasta)

    def _validate_and_output(self):
        """HMM搜索、基序注释并生成输出|HMM search, motif annotation and output"""
        label = self.label
        tc = self.type_config
        step_dir = os.path.join(self.config.output_dir, tc['output_dir'])
        domtblout = os.path.join(step_dir, f"{self.effector_type}.domtblout")
        candidates_file = os.path.join(step_dir, tc['candidates_file'])

        # 使用类型专属阈值或全局阈值|Use type-specific or global threshold
        score_threshold = tc.get('score_threshold') if tc.get('score_threshold') is not None else self.config.score_threshold

        # Step 1: hmmsearch
        self.logger.info(f"步骤1: hmmsearch搜索({label} HMM, -T {score_threshold})|Step 1: hmmsearch ({label} HMM, -T {score_threshold})")

        if is_step_completed(domtblout):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: hmmsearch {label}搜索|hmmsearch {label}")
            details = parse_domtblout_details(domtblout)
        else:
            os.makedirs(step_dir, exist_ok=True)

            cmd = build_conda_command(self.config.hmmsearch_path, [
                '--domtblout', domtblout,
                '-T', str(score_threshold),
                '--cpu', str(self.config.threads),
                self.hmm_file,
                self.config._combined_fasta,
            ])

            success, stdout, stderr = run_command(
                cmd, self.logger, f"hmmsearch {label} HMM搜索|hmmsearch {label} HMM search"
            )
            if not success:
                self.logger.error(f"hmmsearch运行失败|hmmsearch failed for {label}")
                return

            details = parse_domtblout_details(domtblout)

        # Best hit per target
        best_hits = {}
        for d in details:
            tid = d['target']
            if tid not in best_hits or d['domain_score'] > best_hits[tid]['domain_score']:
                best_hits[tid] = d

        self.logger.info(f"hmmsearch命中|hmmsearch hits: {len(best_hits)}个唯一序列|unique sequences")

        # Step 2: Motif annotation
        has_motif = tc['motif_regex'] is not None
        if has_motif:
            self.logger.info(f"步骤2: {tc['motif_label']}基序注释|Step 2: {tc['motif_label']} motif annotation")
        else:
            self.logger.info(f"步骤2: 无需基序注释|Step 2: No motif annotation needed")

        seq_dict = parse_fasta_to_dict(self.config._combined_fasta)
        rows = []

        for seq_id in sorted(best_hits.keys()):
            sequence = seq_dict.get(seq_id, '')
            if not sequence:
                continue

            detail = best_hits[seq_id]

            # SignalP annotation
            sp_info = self.signalp_results.get(seq_id, {})
            has_sp = sp_info.get('has_signal_peptide', None)
            sp_pos = sp_info.get('cs_position', 'N/A')

            row = {
                'Effector_Type': label,
                'Sample': self.source_map.get(seq_id, 'Unknown'),
                'Sequence_ID': seq_id,
                'Sequence': sequence,
                'SignalP': 'Yes' if has_sp else ('No' if has_sp is not None else 'N/A'),
                'SP_Cleavage_Site': sp_pos if has_sp else ('-' if has_sp is not None else 'N/A'),
                'HMM': 'Yes',
                'HMM_Evalue': f"{detail['domain_evalue']:.2e}",
                'HMM_Score': f"{detail['domain_score']:.1f}",
            }

            # Motif annotation
            if has_motif:
                motif_regex: Pattern = tc['motif_regex']
                motif_matches = list(motif_regex.finditer(sequence))
                has_m = len(motif_matches) > 0
                m_label = tc['motif_label']
                row[f'{m_label}'] = 'Yes' if has_m else 'No'
                row[f'{m_label}_Position'] = ';'.join(
                    str(m.start() + 1) for m in motif_matches
                ) if motif_matches else 'None'

            rows.append(row)

        if not rows:
            self.logger.warning(f"未找到{label}候选效应子|No {label} candidates found")
            os.makedirs(step_dir, exist_ok=True)
            pd.DataFrame(rows).to_csv(candidates_file, sep='\t', index=False)
            return

        df = pd.DataFrame(rows)
        df.to_csv(candidates_file, sep='\t', index=False)

        # Summary
        total = len(df)
        sp_yes = len(df[df['SignalP'] == 'Yes'])
        self.logger.info(f"总候选数|Total candidates: {total}")
        self.logger.info(f"含信号肽|With signal peptide: {sp_yes}")

        if has_motif:
            m_label = tc['motif_label']
            motif_yes = len(df[df[m_label] == 'Yes'])
            self.logger.info(f"{m_label}=Yes: {motif_yes}")

        self.logger.info(f"候选{label}效应子已保存|{label} candidates saved: {candidates_file}")

    def _generate_versions(self, start_time):
        """生成版本信息文件|Generate software versions file"""
        tools = {
            'signalp': self.config.signalp_path,
            'hmmsearch': self.config.hmmsearch_path,
        }
        type_threshold = self.type_config.get('score_threshold') if self.type_config.get('score_threshold') is not None else self.config.score_threshold
        params = {
            'mode': self.effector_type,
            'organism': self.config.organism,
            'score_threshold': type_threshold,
            'skip_signalp': self.config.skip_signalp,
            'threads': self.config.threads,
            f'{self.effector_type}_hmm': self.hmm_file,
        }
        if self.type_config['motif_label']:
            params['motif'] = self.type_config['motif_label']
        generate_software_versions(self.config.output_dir, tools, params, start_time)
