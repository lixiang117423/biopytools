"""Phytophthora RxLR效应子鉴定模块|Phytophthora RxLR Effector Identification Module"""

import os
from datetime import datetime
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .config import PhytoEffectorConfig
from .utils import (
    PhytoEffectorLogger, RxLRMotifScanner, build_conda_command,
    extract_sequences_by_id, format_number, generate_software_versions,
    has_tm_outside_sp, is_step_completed, merge_fasta_files,
    parse_blastp_tabular, parse_domtblout_details, parse_domtblout_hits,
    parse_fasta_to_dict, parse_signalp_output, parse_tmhmm_output,
    run_command, run_signalp3, save_signalp3_compatible, merge_signalp_results,
)


class RxLRFinder:
    """Phytophthora RxLR效应子鉴定器|Phytophthora RxLR Effector Finder"""

    def __init__(self, config: PhytoEffectorConfig, logger=None, signalp_results: Dict = None):
        self.config = config
        self.logger = logger
        self.motif_scanner = RxLRMotifScanner(window_start=20, window_end=120)
        self.signalp_results: Dict[str, Dict] = signalp_results if signalp_results is not None else {}
        self.source_map: Dict[str, str] = {}

    def find_effectors(self) -> bool:
        """运行RxLR效应子鉴定流程|Run RxLR effector identification pipeline"""
        start_time = datetime.now()
        self.logger.info("=" * 60)
        self.logger.info("开始RxLR效应子鉴定|Starting RxLR effector identification")
        self.logger.info("=" * 60)

        try:
            self._prepare_input()

            if not self.signalp_results and not self.config.skip_signalp:
                self._run_signalp()
            elif not self.signalp_results:
                self.logger.info("跳过SignalP预测|Skipping SignalP prediction")

            # 步骤2: hmmsearch搜索|Step 2: hmmsearch
            hmmsearch_hits = self._run_hmmsearch()

            # 步骤3: BLASTP搜索(可选)|Step 3: BLASTP search (optional)
            blastp_hits = self._run_blastp_search()

            # 步骤4: 正则表达式直接搜索(含信号肽但未命中HMM/BLASTP的序列)|Step 4: Regex search
            regex_hits = self._run_regex_search(hmmsearch_hits | blastp_hits)

            # 步骤5: WY结构域搜索(可选)|Step 5: WY domain search (optional)
            wy_hits = set()
            if self.config.use_wy_domain:
                wy_hits = self._run_wy_search()

            # 合并所有hits|Merge all hits
            all_hits = hmmsearch_hits | blastp_hits | regex_hits | wy_hits
            self.logger.info(f"合并后候选序列数(TMHMM过滤前)|Total merged candidates (pre-TMHMM): {len(all_hits)}")

            # 步骤6: TMHMM过滤|Step 6: TMHMM filtering
            filtered_hits = self._run_tmhmm_filter(all_hits)

            # 步骤7: 基序校验和输出|Step 7: Motif validation and output
            self._validate_and_output(filtered_hits, hmmsearch_hits, blastp_hits)

            elapsed = (datetime.now() - start_time).total_seconds()
            self.logger.info(f"RxLR效应子鉴定完成|RxLR effector identification completed in {elapsed:.1f}秒|seconds")
            self.logger.info("=" * 60)

            self._generate_versions(start_time)
            return True

        except Exception as e:
            self.logger.error(f"RxLR效应子鉴定失败|RxLR effector identification failed: {e}")
            raise

    def _prepare_input(self):
        """准备输入文件(合并多个FASTA)|Prepare input files (merge multiple FASTAs)"""
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

        if is_step_completed(summary_file):
            self.logger.info("跳过已完成步骤|Skipping completed step: SignalP预测|SignalP prediction")
            self.signalp_results = parse_signalp_output(step_dir, 'combined_input')
            self.logger.info(f"加载SignalP结果|Loaded SignalP results: {len(self.signalp_results)}条记录|records")
            return

        os.makedirs(step_dir, exist_ok=True)

        if self.config.signalp_version == 'both':
            self._run_signalp_both(step_dir, summary_file)
            return

        if self.config.signalp_version == '3':
            self.logger.info("步骤1: 运行SignalP 3.0信号肽预测|Step 1: Running SignalP 3.0 prediction")
            self.signalp_results = run_signalp3(
                self.config.signalp3_path, self.config._combined_fasta,
                self.logger, self.config.signalp3_sprob_threshold,
            )
            save_signalp3_compatible(self.signalp_results, summary_file)
            return

        self.logger.info("步骤1: 运行SignalP 6.0信号肽预测|Step 1: Running SignalP 6.0 prediction")

    def _run_signalp_both(self, step_dir, summary_file):
        """同时运行SignalP 3.0和6.0，取并集|Run both SP3 and SP6, take union"""
        sp6_results = {}
        sp3_results = {}

        # 运行SignalP 6.0
        self.logger.info("步骤1: 同时运行SignalP 3.0+6.0|Step 1: Running SignalP 3.0 + 6.0")
        try:
            sp6_cmd = build_conda_command(self.config.signalp_path, [
                '--fastafile', self.config._combined_fasta,
                '--output_dir', step_dir,
                '--format', 'txt',
                '--organism', self.config.organism,
                '--mode', 'fast',
                '--bsize', str(self.config.threads),
                '--write_procs', str(self.config.threads),
                '--torch_num_threads', str(self.config.threads),
            ])
            success, stdout, stderr = run_command(sp6_cmd, self.logger, "SignalP 6.0信号肽预测|SignalP 6.0 prediction")
            if success:
                sp6_results = parse_signalp_output(step_dir, 'combined_input')
                sp6_sp = sum(1 for v in sp6_results.values() if v['has_signal_peptide'])
                self.logger.info(f"SignalP 6.0: {len(sp6_results)}条预测|predictions, {sp6_sp}条SP+|with SP")
            else:
                self.logger.warning("SignalP 6.0运行失败|SignalP 6.0 failed")
        except Exception as e:
            self.logger.warning(f"SignalP 6.0异常|SignalP 6.0 error: {e}")

        # 运行SignalP 3.0
        if os.path.exists(self.config.signalp3_path):
            try:
                sp3_results = run_signalp3(
                    self.config.signalp3_path, self.config._combined_fasta,
                    self.logger, self.config.signalp3_sprob_threshold,
                )
            except Exception as e:
                self.logger.warning(f"SignalP 3.0异常|SignalP 3.0 error: {e}")
        else:
            self.logger.warning(f"SignalP 3.0未找到，跳过|SignalP 3.0 not found, skipping: {self.config.signalp3_path}")

        # 合并(并集)
        self.signalp_results = merge_signalp_results(sp6_results, sp3_results)
        sp_count = sum(1 for v in self.signalp_results.values() if v['has_signal_peptide'])
        self.logger.info(
            f"SignalP合并结果|Merged results: {len(self.signalp_results)}条预测|predictions, "
            f"{sp_count}条SP+|with SP (SP6={len(sp6_results)}, SP3={len(sp3_results)})"
        )

        # 保存合并结果用于断点续传
        save_signalp3_compatible(self.signalp_results, summary_file)

        cmd = build_conda_command(self.config.signalp_path, [
            '--fastafile', self.config._combined_fasta,
            '--output_dir', step_dir,
            '--format', 'txt',
            '--organism', self.config.organism,
            '--mode', 'fast',
            '--bsize', str(self.config.threads),
            '--write_procs', str(self.config.threads),
            '--torch_num_threads', str(self.config.threads),
        ])

        success, stdout, stderr = run_command(cmd, self.logger, "SignalP信号肽预测|SignalP signal peptide prediction")
        if not success:
            self.logger.warning("SignalP运行失败，继续但SignalP列将为空|SignalP failed, continuing but SP columns will be empty")
            return

        self.signalp_results = parse_signalp_output(step_dir, 'combined_input')
        sp_count = sum(1 for v in self.signalp_results.values() if v['has_signal_peptide'])
        self.logger.info(f"SignalP完成|SignalP completed: {len(self.signalp_results)}条预测|predictions, {sp_count}条含信号肽|with signal peptide")

    def _run_hmmsearch(self) -> Set[str]:
        """运行hmmsearch|Run hmmsearch"""
        step_dir = os.path.join(self.config.output_dir, '02_hmmsearch')
        domtblout = os.path.join(step_dir, 'rxlr.domtblout')

        if is_step_completed(domtblout):
            self.logger.info("跳过已完成步骤|Skipping completed step: hmmsearch搜索|hmmsearch")
            hits = parse_domtblout_hits(domtblout)
            self.logger.info(f"加载hmmsearch结果|Loaded hmmsearch results: {len(hits)}条命中|hits")
            return hits

        self.logger.info("步骤2: hmmsearch搜索(RxLR HMM)|Step 2: hmmsearch (RxLR HMM)")
        os.makedirs(step_dir, exist_ok=True)

        cmd_args = [
            '--domtblout', domtblout,
            '-T', str(self.config.score_threshold),
            '--cpu', str(self.config.threads),
            self.config.rxlr_hmm,
            self.config._combined_fasta,
        ]

        cmd = build_conda_command(self.config.hmmsearch_path, cmd_args)

        success, stdout, stderr = run_command(cmd, self.logger, "hmmsearch RxLR HMM搜索|hmmsearch RxLR HMM search")
        if not success:
            self.logger.error("hmmsearch运行失败|hmmsearch failed")
            return set()

        hits = parse_domtblout_hits(domtblout)
        self.logger.info(f"hmmsearch命中|hmmsearch hits: {len(hits)}")
        return hits

    def _run_blastp_search(self) -> Set[str]:
        """运行BLASTP搜索补充RxLR候选|Run BLASTP search to supplement RxLR candidates"""
        if not self.config.rxlr_blastp_queries:
            self.logger.info("步骤3: 跳过BLASTP搜索(未提供查询序列)|Step 3: Skipping BLASTP search (no queries provided)")
            return set()

        step_dir = os.path.join(self.config.output_dir, '03_blastp')
        db_prefix = os.path.join(step_dir, 'genome_db')
        blastp_out = os.path.join(step_dir, 'blastp_results.tsv')

        if is_step_completed(blastp_out):
            self.logger.info("跳过已完成步骤|Skipping completed step: BLASTP搜索|BLASTP search")
            hits = parse_blastp_tabular(blastp_out)
            self.logger.info(f"加载BLASTP结果|Loaded BLASTP results: {len(hits)}条命中|hits")
            return hits

        self.logger.info("步骤3: BLASTP搜索(RxLR同源补充)|Step 3: BLASTP search (RxLR homology)")
        os.makedirs(step_dir, exist_ok=True)

        # 构建BLAST数据库|Build BLAST database
        import shutil
        makeblastdb = self.config.blastp_path.replace('blastp', 'makeblastdb')
        if not os.path.exists(makeblastdb):
            makeblastdb_bin = os.path.join(os.path.dirname(self.config.blastp_path), 'makeblastdb')
            if os.path.exists(makeblastdb_bin):
                makeblastdb = makeblastdb_bin

        db_files = [f"{db_prefix}.phr", f"{db_prefix}.pin", f"{db_prefix}.psq"]
        if not all(os.path.exists(f) for f in db_files):
            cmd_makedb = build_conda_command(makeblastdb, [
                '-in', self.config._combined_fasta,
                '-dbtype', 'prot',
                '-out', db_prefix,
            ])
            success, _, _ = run_command(cmd_makedb, self.logger, "构建BLAST数据库|Building BLAST database")
            if not success:
                self.logger.error("构建BLAST数据库失败|Failed to build BLAST database")
                return set()

        # 运行BLASTP|Run BLASTP
        cmd_blastp = build_conda_command(self.config.blastp_path, [
            '-query', self.config.rxlr_blastp_queries,
            '-db', db_prefix,
            '-evalue', str(self.config.rxlr_blastp_evalue),
            '-outfmt', '6',
            '-out', blastp_out,
            '-num_threads', str(self.config.threads),
        ])

        success, stdout, stderr = run_command(cmd_blastp, self.logger, "BLASTP RxLR同源搜索|BLASTP RxLR homology search")
        if not success:
            self.logger.error("BLASTP运行失败|BLASTP failed")
            return set()

        # 解析BLASTP结果并检查RxLR/EER基序|Parse BLASTP and check RxLR/EER motifs
        all_blastp_hits = parse_blastp_tabular(blastp_out)
        self.logger.info(f"BLASTP原始命中|BLASTP raw hits: {len(all_blastp_hits)}")

        motif_positive_hits = set()
        seq_dict = parse_fasta_to_dict(self.config._combined_fasta)
        for seq_id in all_blastp_hits:
            sequence = seq_dict.get(seq_id, '')
            if not sequence:
                continue
            motif_result = self.motif_scanner.scan(seq_id, sequence)
            if motif_result['is_candidate']:
                motif_positive_hits.add(seq_id)

        self.logger.info(f"BLASTP命中中含RxLR/EER基序|BLASTP hits with RxLR/EER motif: {len(motif_positive_hits)}/{len(all_blastp_hits)}")
        return motif_positive_hits

    def _run_regex_search(self, existing_hits: Set[str]) -> Set[str]:
        """正则表达式搜索含SignalP的序列|Regex search on SP+ sequences"""
        self.logger.info("步骤4: 正则表达式搜索RxLR/EER基序(含SignalP)|Step 4: Regex search for RxLR/EER motifs (SP+)")

        seq_dict = parse_fasta_to_dict(self.config._combined_fasta)
        regex_hits = set()
        regex_total = 0

        for seq_id, sequence in seq_dict.items():
            if seq_id in existing_hits:
                continue

            sp_info = self.signalp_results.get(seq_id, {})
            if not sp_info.get('has_signal_peptide', False):
                continue

            motif_result = self.motif_scanner.scan(seq_id, sequence)
            if motif_result['is_candidate']:
                regex_hits.add(seq_id)
            regex_total += 1

        self.logger.info(f"SignalP+序列中正则匹配|Regex matches in SP+ sequences: {len(regex_hits)}/{regex_total}")
        return regex_hits

    def _run_wy_search(self) -> Set[str]:
        """运行WY结构域搜索|Run WY domain search"""
        step_dir = os.path.join(self.config.output_dir, '04_wy_domain')
        domtblout = os.path.join(step_dir, 'wy.domtblout')

        if is_step_completed(domtblout):
            self.logger.info("跳过已完成步骤|Skipping completed step: WY结构域搜索|WY domain search")
            hits = parse_domtblout_hits(domtblout)
            self.logger.info(f"加载WY搜索结果|Loaded WY search results: {len(hits)}条命中|hits")
            return hits

        self.logger.info("步骤5: hmmsearch WY结构域搜索|Step 5: hmmsearch WY domain search")
        os.makedirs(step_dir, exist_ok=True)

        cmd = build_conda_command(self.config.hmmsearch_path, [
            '--domtblout', domtblout,
            '-T', str(self.config.score_threshold),
            '--cpu', str(self.config.threads),
            self.config.rxlr_wy_hmm,
            self.config._combined_fasta,
        ])

        success, stdout, stderr = run_command(cmd, self.logger, "hmmsearch WY结构域搜索|hmmsearch WY domain search")
        if not success:
            self.logger.warning("WY搜索运行失败|WY domain search failed")
            return set()

        hits = parse_domtblout_hits(domtblout)
        self.logger.info(f"WY结构域命中|WY domain hits: {len(hits)}")
        return hits

    def _run_tmhmm_filter(self, all_hits: Set[str]) -> Set[str]:
        """TMHMM过滤: 移除信号肽区外有跨膜域的候选|TMHMM filter: remove candidates with TM outside SP"""
        step_dir = os.path.join(self.config.output_dir, '05_tmhmm')
        tmhmm_out = os.path.join(step_dir, 'tmhmm_results.txt')
        filtered_file = os.path.join(step_dir, 'tmhmm_filtered_ids.txt')

        self.logger.info("步骤6: TMHMM跨膜域过滤|Step 6: TMHMM transmembrane domain filtering")

        if is_step_completed(filtered_file):
            self.logger.info("跳过已完成步骤|Skipping completed step: TMHMM过滤|TMHMM filtering")
            with open(filtered_file, 'r') as f:
                filtered_ids = set(line.strip().split('/')[0] for line in f if line.strip())
            self.logger.info(f"加载TMHMM过滤结果|Loaded TMHMM filter results: {len(filtered_ids)}通过|passed")
            return filtered_ids

        os.makedirs(step_dir, exist_ok=True)

        # 提取候选序列|Extract candidate sequences
        candidates_fasta = os.path.join(step_dir, 'rxlr_candidates_all.faa')
        extract_sequences_by_id(self.config._combined_fasta, all_hits, candidates_fasta)

        # 运行TMHMM|Run TMHMM
        cmd = build_conda_command(self.config.tmhmm_path, [
            '<', candidates_fasta,
        ])
        # TMHMM reads from stdin
        import subprocess
        self.logger.info(f"执行|Executing: TMHMM跨膜域预测|TMHMM transmembrane prediction")
        self.logger.info(f"命令|Command: {self.config.tmhmm_path} < {candidates_fasta}")

        try:
            result = subprocess.run(
                self.config.tmhmm_path, shell=False, capture_output=True, text=True,
                check=True, timeout=86400,
                input=open(candidates_fasta, 'r').read()
            )
            with open(tmhmm_out, 'w') as f:
                f.write(result.stdout)
            success = True
        except Exception as e:
            self.logger.error(f"TMHMM运行失败|TMHMM failed: {e}")
            return all_hits

        # 解析TMHMM结果|Parse TMHMM results
        tm_helices = parse_tmhmm_output(tmhmm_out)
        self.logger.info(f"TMHMM预测含TM螺旋的蛋白|Proteins with TM helices: {len(tm_helices)}/{len(all_hits)}")

        # 过滤: 保留有SignalP且TM不在SP区外的候选|Filter: keep candidates with SP and no TM outside SP
        filtered_hits = set()
        filtered_out_no_sp = 0
        filtered_out_tm = 0
        passed_no_tm = 0

        for seq_id in all_hits:
            sp_info = self.signalp_results.get(seq_id, {})
            has_sp = sp_info.get('has_signal_peptide', False)
            cs_str = sp_info.get('cs_position', '-')

            # 必须有信号肽|Must have signal peptide
            if not has_sp:
                filtered_out_no_sp += 1
                continue

            # 解析切割位点|Parse cleavage site
            sp_cleavage = 0
            if cs_str and cs_str not in ('-', 'N/A'):
                try:
                    sp_cleavage = int(cs_str.split('-')[0]) if '-' in cs_str else int(cs_str)
                except ValueError:
                    sp_cleavage = 0

            # 检查TM是否在SP区外|Check if TM is outside SP region
            protein_tms = tm_helices.get(seq_id, [])
            if has_tm_outside_sp(protein_tms, sp_cleavage):
                filtered_out_tm += 1
                continue

            filtered_hits.add(seq_id)
            if not protein_tms:
                passed_no_tm += 1

        self.logger.info(f"TMHMM过滤结果|TMHMM filter results:")
        self.logger.info(f"  无信号肽被过滤|Filtered (no signal peptide): {filtered_out_no_sp}")
        self.logger.info(f"  SP区外有TM被过滤|Filtered (TM outside SP): {filtered_out_tm}")
        self.logger.info(f"  保留(无TM)|Passed (no TM): {passed_no_tm}")
        self.logger.info(f"  保留(TM仅在SP区内)|Passed (TM within SP only): {len(filtered_hits) - passed_no_tm}")
        self.logger.info(f"  总保留|Total passed: {len(filtered_hits)}/{len(all_hits)}")

        with open(filtered_file, 'w') as f:
            for sid in sorted(filtered_hits):
                f.write(f"{sid}\n")

        return filtered_hits

    def _validate_and_output(self, all_hits: Set[str], hmm_hits: Set[str], blastp_hits: Set[str]):
        """基序注释并生成输出|Annotate motifs and generate output"""
        step_dir = os.path.join(self.config.output_dir, '06_candidates')
        candidates_file = os.path.join(step_dir, 'rxlr_candidates.tsv')

        self.logger.info("步骤7: 基序注释|Step 7: Motif annotation")
        os.makedirs(step_dir, exist_ok=True)

        seq_dict = parse_fasta_to_dict(self.config._combined_fasta)

        # 加载TMHMM结果|Load TMHMM results
        tmhmm_file = os.path.join(self.config.output_dir, '05_tmhmm', 'tmhmm_results.txt')
        tm_helices = {}
        if os.path.exists(tmhmm_file):
            tm_helices = parse_tmhmm_output(tmhmm_file)

        rows = []

        for seq_id in sorted(all_hits):
            sequence = seq_dict.get(seq_id, '')
            if not sequence:
                continue

            motif_result = self.motif_scanner.scan(seq_id, sequence)

            sp_info = self.signalp_results.get(seq_id, {})
            has_sp = sp_info.get('has_signal_peptide', None)
            sp_pos = sp_info.get('cs_position', 'N/A')

            in_hmm = seq_id in hmm_hits
            in_blastp = seq_id in blastp_hits

            # TMHMM信息|TMHMM info
            protein_tms = tm_helices.get(seq_id, [])
            tm_count = len(protein_tms)
            tm_positions = ';'.join(f"{s}-{e}" for s, e in protein_tms) if protein_tms else 'None'

            # 判断TM是否在SP区外|Check if TM outside SP
            sp_cleavage = 0
            if has_sp and sp_pos and sp_pos not in ('-', 'N/A'):
                try:
                    sp_cleavage = int(sp_pos.split('-')[0]) if '-' in sp_pos else int(sp_pos)
                except ValueError:
                    sp_cleavage = 0
            tm_outside = has_tm_outside_sp(protein_tms, sp_cleavage) if protein_tms else False

            # 来源|Source
            sources = []
            if in_hmm: sources.append('HMM')
            if in_blastp: sources.append('BLASTP')
            if not in_hmm and not in_blastp: sources.append('Regex')

            rows.append({
                'Effector_Type': 'RxLR',
                'Sample': self.source_map.get(seq_id, 'Unknown'),
                'Sequence_ID': seq_id,
                'Sequence': sequence,
                'SignalP': 'Yes' if has_sp else ('No' if has_sp is not None else 'N/A'),
                'SP_Cleavage_Site': sp_pos if has_sp else ('-' if has_sp is not None else 'N/A'),
                'Source': '/'.join(sources),
                'HMM': 'Yes' if in_hmm else 'No',
                'BLASTP': 'Yes' if in_blastp else 'No',
                'RxLR': 'Yes' if motif_result['has_rxlr'] else 'No',
                'RxLR_Position': motif_result['rxlr_positions'],
                'EER': 'Yes' if motif_result['has_eer'] else 'No',
                'EER_Position': motif_result['eer_positions'],
                'TMHMM_TMs': tm_count,
                'TMHMM_Position': tm_positions,
            })

        df = pd.DataFrame(rows)
        df.to_csv(candidates_file, sep='\t', index=False)
        self.logger.info(f"候选RxLR效应子已保存|RxLR candidates saved: {candidates_file}")

        total = len(df)
        sp_yes = len(df[df['SignalP'] == 'Yes'])
        hmm_yes = len(df[df['HMM'] == 'Yes'])
        blastp_yes = len(df[df['BLASTP'] == 'Yes'])
        rxlr_yes = len(df[df['RxLR'] == 'Yes'])
        eer_yes = len(df[df['EER'] == 'Yes'])
        self.logger.info(f"总候选数|Total candidates: {total}")
        self.logger.info(f"SignalP=Yes: {sp_yes}, HMM=Yes: {hmm_yes}, BLASTP=Yes: {blastp_yes}, RxLR=Yes: {rxlr_yes}, EER=Yes: {eer_yes}")

    def _generate_versions(self, start_time):
        """生成版本信息文件|Generate software versions file"""
        tools = {
            'signalp': self.config.signalp_path,
            'hmmsearch': self.config.hmmsearch_path,
        }
        if self.config.rxlr_blastp_queries:
            tools['blastp'] = self.config.blastp_path
        tools['tmhmm'] = self.config.tmhmm_path

        params = {
            'mode': 'rxlr',
            'organism': self.config.organism,
            'evalue': self.config.evalue,
            'use_wy_domain': self.config.use_wy_domain,
            'skip_signalp': self.config.skip_signalp,
            'threads': self.config.threads,
            'rxlr_hmm': self.config.rxlr_hmm,
            'window': '21-120 (R/Q/GxLR or EER)',
            'tmhmm_filter': 'SP required, no TM outside SP',
        }
        if self.config.rxlr_blastp_queries:
            params['rxlr_blastp_queries'] = self.config.rxlr_blastp_queries
            params['rxlr_blastp_evalue'] = self.config.rxlr_blastp_evalue

        generate_software_versions(self.config.output_dir, tools, params, start_time)
