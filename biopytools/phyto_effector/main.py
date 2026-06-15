"""Phytophthora效应子鉴定主程序模块|Phytophthora Effector Identification Main Module"""

import argparse
import os
import sys
import time
from pathlib import Path

from .config import PhytoEffectorConfig
from .utils import (
    PhytoEffectorLogger, build_conda_command, generate_software_versions,
    is_step_completed, merge_candidate_files, parse_signalp_output,
    run_command, run_signalp3, save_signalp3_compatible, merge_signalp_results,
)
from .rxlr_finder import RxLRFinder
from .crn_finder import CRNFinder
from .generic_finder import GenericEffectorFinder, EFFECTOR_TYPE_CONFIG


# 通用效应子类型列表|Generic effector type list
GENERIC_TYPES = ['nlp', 'protease', 'scp', 'elicitin', 'yxsl']
# 所有效应子类型|All effector types
ALL_TYPES = ['rxlr', 'crn'] + GENERIC_TYPES

# 候选TSV路径映射(相对于效应子类型输出目录)|Candidate TSV path mapping
CANDIDATES_REL_PATH = {
    'rxlr': os.path.join('06_candidates', 'rxlr_candidates.tsv'),
    'crn': os.path.join('03_candidates', 'crn_candidates.tsv'),
}


def _get_candidates_rel_path(effector_type: str) -> str:
    """获取效应子类型的候选TSV相对路径|Get candidate TSV relative path for effector type"""
    if effector_type in CANDIDATES_REL_PATH:
        return CANDIDATES_REL_PATH[effector_type]
    tc = EFFECTOR_TYPE_CONFIG.get(effector_type, {})
    return os.path.join(tc.get('output_dir', ''), tc.get('candidates_file', ''))


FASTA_EXTENSIONS = {'.fa', '.fasta', '.faa', '.pep', '.protein', '.prot'}


def _determine_run_mode(input_path: str):
    """判断输入模式: 单文件或多样本|Determine run mode: single-file or multi-sample

    Returns:
        ('single', [file_path]) 或 ('multi', [file1, file2, ...])
    """
    p = Path(input_path)
    if p.is_file():
        return 'single', [str(p)]
    elif p.is_dir():
        files = sorted([
            str(f) for f in p.iterdir()
            if f.is_file() and f.suffix.lower() in FASTA_EXTENSIONS
        ])
        if not files:
            raise ValueError(f"目录中未找到FASTA文件|No FASTA files found in directory: {p}")
        if len(files) == 1:
            return 'single', files
        return 'multi', files
    else:
        raise ValueError(f"输入路径不存在|Input path not found: {input_path}")


def _check_duplicate_stems(files: list):
    """检测文件名stem冲突|Check for duplicate file stems"""
    stems = [Path(f).stem for f in files]
    seen = set()
    for s in stems:
        if s in seen:
            raise ValueError(
                f"目录中存在同名样本文件|Duplicate sample name detected: {s}.fa / {s}.fasta"
            )
        seen.add(s)


def _run_signalp(config, logger):
    """运行SignalP预测(共享步骤)|Run SignalP prediction (shared step)"""
    step_dir = os.path.join(config.output_dir, '01_signalp')
    summary_file = os.path.join(step_dir, 'prediction_results.txt')

    if is_step_completed(summary_file):
        logger.info("跳过已完成步骤|Skipping completed step: SignalP预测|SignalP prediction")
        return parse_signalp_output(step_dir)

    os.makedirs(step_dir, exist_ok=True)

    if config.signalp_version == 'both':
        return _run_signalp_both(config, logger, step_dir, summary_file)

    if config.signalp_version == '3':
        logger.info("步骤1: 运行SignalP 3.0信号肽预测|Step 1: Running SignalP 3.0 prediction")
        results = run_signalp3(
            config.signalp3_path, config._combined_fasta,
            logger, config.signalp3_sprob_threshold,
        )
        save_signalp3_compatible(results, summary_file)
        return results

    logger.info("步骤1: 运行SignalP 6.0信号肽预测|Step 1: Running SignalP 6.0 prediction")


def _run_signalp_both(config, logger, step_dir, summary_file):
    """同时运行SignalP 3.0和6.0，取并集|Run both SP3 and SP6, take union"""
    sp6_results = {}
    sp3_results = {}

    # SignalP 6.0
    try:
        sp6_cmd = build_conda_command(config.signalp_path, [
            '--fastafile', config._combined_fasta,
            '--output_dir', step_dir,
            '--format', 'txt',
            '--organism', config.organism,
            '--mode', config.signalp_mode,
            '--bsize', str(config.threads),
            '--write_procs', str(config.threads),
            '--torch_num_threads', str(config.threads),
        ])
        success, stdout, stderr = run_command(sp6_cmd, logger, "SignalP 6.0信号肽预测|SignalP 6.0 prediction")
        if success:
            sp6_results = parse_signalp_output(step_dir)
            sp6_sp = sum(1 for v in sp6_results.values() if v['has_signal_peptide'])
            logger.info(f"SignalP 6.0: {len(sp6_results)}条预测|predictions, {sp6_sp}条SP+|with SP")
        else:
            logger.warning("SignalP 6.0运行失败|SignalP 6.0 failed")
    except Exception as e:
        logger.warning(f"SignalP 6.0异常|SignalP 6.0 error: {e}")

    # SignalP 3.0
    if os.path.exists(config.signalp3_path):
        try:
            sp3_results = run_signalp3(
                config.signalp3_path, config._combined_fasta,
                logger, config.signalp3_sprob_threshold,
            )
        except Exception as e:
            logger.warning(f"SignalP 3.0异常|SignalP 3.0 error: {e}")
    else:
        logger.warning(f"SignalP 3.0未找到，跳过|SignalP 3.0 not found, skipping: {config.signalp3_path}")

    # 合并(并集)
    merged = merge_signalp_results(sp6_results, sp3_results)
    sp_count = sum(1 for v in merged.values() if v['has_signal_peptide'])
    logger.info(
        f"SignalP合并结果|Merged results: {len(merged)}条预测|predictions, "
        f"{sp_count}条SP+|with SP (SP6={len(sp6_results)}, SP3={len(sp3_results)})"
    )
    save_signalp3_compatible(merged, summary_file)
    return merged

    cmd = build_conda_command(config.signalp_path, [
        '--fastafile', config._combined_fasta,
        '--output_dir', step_dir,
        '--format', 'txt',
        '--organism', config.organism,
        '--mode', config.signalp_mode,
        '--bsize', str(config.threads),
        '--write_procs', str(config.threads),
        '--torch_num_threads', str(config.threads),
    ])

    success, stdout, stderr = run_command(cmd, logger, "SignalP信号肽预测|SignalP signal peptide prediction")
    if not success:
        logger.warning("SignalP运行失败，继续但SignalP列将为空|SignalP failed, continuing but SP columns will be empty")
        return {}

    signalp_results = parse_signalp_output(step_dir)
    sp_count = sum(1 for v in signalp_results.values() if v['has_signal_peptide'])
    logger.info(f"SignalP完成|SignalP completed: {len(signalp_results)}条预测|predictions, {sp_count}条含信号肽|with signal peptide")
    return signalp_results


def _make_rxlr_config(sample_fasta, sample_output_dir, args):
    """创建RxLR配置|Create RxLR config for a sample"""
    return PhytoEffectorConfig(
        mode='rxlr',
        input_path=sample_fasta,
        output_dir=os.path.join(sample_output_dir, 'rxlr'),
        skip_signalp=True,
        signalp_path=args.signalp_path,
        organism=args.organism,
        signalp_mode=args.signalp_mode,
        signalp_version=args.signalp_version,
        signalp3_path=args.signalp3_path,
        signalp3_sprob_threshold=args.signalp3_sprob_threshold,
        hmmsearch_path=args.hmmsearch_path,
        blastp_path=args.blastp_path,
        tmhmm_path=args.tmhmm_path,
        rxlr_hmm=args.rxlr_hmm,
        rxlr_blastp_queries=args.rxlr_blastp_queries,
        use_wy_domain=args.use_wy_domain,
        rxlr_wy_hmm=args.rxlr_wy_hmm,
        evalue=args.evalue,
        score_threshold=args.score_threshold,
        threads=args.threads,
    )


def _make_crn_config(sample_fasta, sample_output_dir, args):
    """创建CRN配置|Create CRN config for a sample"""
    return PhytoEffectorConfig(
        mode='crn',
        input_path=sample_fasta,
        output_dir=os.path.join(sample_output_dir, 'crn'),
        skip_signalp=True,
        signalp_path=args.signalp_path,
        organism=args.organism,
        signalp_version=args.signalp_version,
        signalp3_path=args.signalp3_path,
        hmmsearch_path=args.hmmsearch_path,
        crn_hmm=args.crn_hmm,
        threads=args.threads,
    )


def _make_type_config(etype, sample_fasta, sample_output_dir, args):
    """创建通用效应子类型配置|Create generic effector type config for a sample"""
    hmm_value = getattr(args, f'{etype}_hmm', None)
    return PhytoEffectorConfig(
        mode=etype,
        input_path=sample_fasta,
        output_dir=os.path.join(sample_output_dir, etype),
        skip_signalp=True,
        signalp_path=args.signalp_path,
        organism=args.organism,
        signalp_version=args.signalp_version,
        signalp3_path=args.signalp3_path,
        hmmsearch_path=args.hmmsearch_path,
        score_threshold=args.score_threshold,
        threads=args.threads,
        **{f'{etype}_hmm': hmm_value},
    )


def _run_single_sample(sample_fasta: str, sample_output_dir: str, args):
    """运行单个样本的完整效应子鉴定流程|Run full effector identification for a single sample

    Args:
        sample_fasta: 单个样本FASTA文件路径|Single sample FASTA file path
        sample_output_dir: 样本输出目录|Sample output directory
        args: 命令行参数|CLI arguments

    Returns:
        dict映射效应子类型到成功/失败|dict mapping effector type -> bool
    """
    sample_name = Path(sample_fasta).stem
    results = {}

    # 创建样品级logger|Create per-sample logger
    sample_log_dir = os.path.join(sample_output_dir, '99_logs')
    os.makedirs(sample_log_dir, exist_ok=True)
    sample_log_file = os.path.join(sample_log_dir, 'pipeline.log')
    sample_logger = PhytoEffectorLogger(log_file=sample_log_file).get_logger()

    sample_logger.info("=" * 60)
    sample_logger.info(f"开始样本效应子鉴定|Starting effector identification: {sample_name}")
    sample_logger.info("=" * 60)

    # SignalP (所有效应子类型共享)|SignalP (shared by all effector types)
    signalp_results = {}
    if not args.skip_signalp:
        rxlr_config = _make_rxlr_config(sample_fasta, sample_output_dir, args)
        rxlr_config.validate()
        signalp_results = _run_signalp(rxlr_config, sample_logger)

    # RxLR鉴定|RxLR identification
    try:
        rxlr_config = _make_rxlr_config(sample_fasta, sample_output_dir, args)
        rxlr_config.validate()
        rxlr_log = os.path.join(rxlr_config.output_dir, '99_logs', 'rxlr_pipeline.log')
        rxlr_logger = PhytoEffectorLogger(log_file=rxlr_log).get_logger()

        if signalp_results:
            finder = RxLRFinder(rxlr_config, rxlr_logger, signalp_results=signalp_results)
        else:
            finder = RxLRFinder(rxlr_config, rxlr_logger)
        results['rxlr'] = finder.find_effectors()
    except Exception as e:
        sample_logger.error(f"RxLR鉴定失败|RxLR identification failed: {e}")
        results['rxlr'] = False

    # CRN鉴定|CRN identification
    try:
        crn_config = _make_crn_config(sample_fasta, sample_output_dir, args)
        crn_config.validate()
        crn_log = os.path.join(crn_config.output_dir, '99_logs', 'crn_pipeline.log')
        crn_logger = PhytoEffectorLogger(log_file=crn_log).get_logger()
        finder = CRNFinder(crn_config, crn_logger, signalp_results=signalp_results)
        results['crn'] = finder.find_effectors()
    except Exception as e:
        sample_logger.error(f"CRN鉴定失败|CRN identification failed: {e}")
        results['crn'] = False

    # 其他效应子类型|Other effector types
    for etype in GENERIC_TYPES:
        try:
            type_config = _make_type_config(etype, sample_fasta, sample_output_dir, args)
            type_config.validate()
            type_log = os.path.join(type_config.output_dir, '99_logs', f'{etype}_pipeline.log')
            type_logger = PhytoEffectorLogger(log_file=type_log).get_logger()
            finder = GenericEffectorFinder(etype, type_config, type_logger,
                                          signalp_results=signalp_results)
            results[etype] = finder.find_effectors()
        except Exception as e:
            sample_logger.error(f"{etype}鉴定失败|{etype} identification failed: {e}")
            results[etype] = False

    sample_logger.info(f"样本效应子鉴定完成|Sample effector identification completed: {sample_name}")
    return results


def _generate_top_summary(output_dir, sample_results, fasta_files, args, start_time):
    """生成顶层汇总信息|Generate top-level summary"""
    os.makedirs(output_dir, exist_ok=True)

    tools = {
        'signalp': args.signalp_path,
        'hmmsearch': args.hmmsearch_path,
        'tmhmm': args.tmhmm_path,
    }
    if args.rxlr_blastp_queries:
        tools['blastp'] = args.blastp_path
    params = {
        'organism': args.organism,
        'evalue': args.evalue,
        'use_wy_domain': args.use_wy_domain,
        'skip_signalp': args.skip_signalp,
        'threads': args.threads,
        'samples': [Path(f).stem for f in fasta_files],
        'mode': 'multi-sample',
    }
    generate_software_versions(output_dir, tools, params, start_time)


def main():
    """效应子鉴定主函数(同时运行所有类型)|Main function (run all effector types together)"""
    parser = argparse.ArgumentParser(
        description="Phytophthora效应子鉴定工具(RxLR+CRN+NLP+...)|Phytophthora Effector Identification Tool",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入FASTA文件或目录|Input FASTA file or directory')
    parser.add_argument('-o', '--output-dir', default='./phyto_effector_output',
                        help='输出目录|Output directory (default: ./phyto_effector_output)')

    # SignalP参数|SignalP parameters
    parser.add_argument('--skip-signalp', action='store_true',
                        help='跳过SignalP预测|Skip SignalP prediction')
    parser.add_argument('--signalp-path', default='~/miniforge3/envs/signalp6/bin/signalp6',
                        help='SignalP程序路径|SignalP program path')
    parser.add_argument('--organism', default='eukarya', choices=['eukarya', 'other'],
                        help='生物类型|Organism type (default: eukarya)')
    parser.add_argument('--signalp-mode', default='slow-sequential',
                        choices=['fast', 'slow', 'slow-sequential'],
                        help='SignalP运行模式|SignalP run mode (default: slow-sequential)')

    # SignalP 3.0参数|SignalP 3.0 parameters
    parser.add_argument('--signalp-version', default='both', choices=['3', '6', 'both'],
                        help='SignalP版本|SignalP version: 3, 6, or both (default: both)')
    parser.add_argument('--signalp3-path',
                        default='~/miniforge3/envs/signalp_v.3.0b/bin/signalp',
                        help='SignalP 3.0程序路径|SignalP 3.0 program path')
    parser.add_argument('--signalp3-sprob-threshold', type=float, default=0.9,
                        help='SignalP 3.0 HMM Sprob阈值|SignalP 3.0 HMM Sprob threshold (default: 0.9)')

    # HMMER参数|HMMER parameters
    parser.add_argument('--hmmsearch-path', default='~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch',
                        help='hmmsearch程序路径|hmmsearch program path')

    # BLASTP参数|BLASTP parameters
    parser.add_argument('--blastp-path', default='~/miniforge3/envs/Blast_v.2.16.0/bin/blastp',
                        help='blastp程序路径|blastp program path')
    parser.add_argument('--rxlr-blastp-queries', default=None,
                        help='RxLR BLASTP查询序列FASTA(默认内置)|RxLR BLASTP query FASTA (default: bundled)')

    # TMHMM参数|TMHMM parameters
    parser.add_argument('--tmhmm-path', default='~/miniforge3/envs/tmmhmm_v.2.0c/bin/tmhmm',
                        help='tmhmm程序路径|tmhmm program path')

    # RxLR参数|RxLR parameters
    parser.add_argument('--rxlr-hmm', default=None,
                        help='RxLR HMM文件路径(默认内置)|RxLR HMM file path (default: bundled)')
    parser.add_argument('--use-wy-domain', action='store_true',
                        help='同时搜索WY结构域|Also search WY domain')
    parser.add_argument('--rxlr-wy-hmm', default=None,
                        help='WY HMM文件路径(默认内置)|WY HMM file path (default: bundled)')
    parser.add_argument('-e', '--evalue', type=float, default=1e-5,
                        help='E-value阈值(已弃用，使用--score-threshold)|E-value threshold (deprecated)')
    parser.add_argument('--score-threshold', type=float, default=0.0,
                        help='HMM score阈值|HMM score threshold (default: 0.0)')

    # CRN参数|CRN parameters
    parser.add_argument('--crn-hmm', default=None,
                        help='CRN HMM文件路径(默认内置)|CRN HMM file path (default: bundled)')

    # 其他效应子类型参数|Other effector type parameters
    for etype in GENERIC_TYPES:
        label = EFFECTOR_TYPE_CONFIG[etype]['label']
        parser.add_argument(f'--{etype}-hmm', default=None,
                            help=f'{label} HMM文件路径(默认内置)|{label} HMM file path (default: bundled)')

    # 通用参数|Common parameters
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads (default: 12)')

    args = parser.parse_args()
    start_time = time.time()

    try:
        run_mode, fasta_files = _determine_run_mode(args.input)

        if run_mode == 'single':
            # 单文件模式: 行为与之前完全一致|Single-file mode: identical to previous behavior
            sample_fasta = fasta_files[0]
            if Path(args.input).is_dir():
                sample_output_dir = args.output_dir
            else:
                sample_output_dir = args.output_dir

            os.makedirs(sample_output_dir, exist_ok=True)

            sample_results = _run_single_sample(sample_fasta, sample_output_dir, args)
            all_success = all(sample_results.values())

            # 版本信息|Version info
            elapsed = time.time() - start_time
            logger = PhytoEffectorLogger(
                log_file=os.path.join(sample_output_dir, 'rxlr', '99_logs', 'pipeline.log')
            ).get_logger()
            logger.info(f"{'='*60}")
            logger.info(f"效应子鉴定完成|Effector identification completed in {elapsed:.1f}秒|seconds")
            for etype in ALL_TYPES:
                type_outdirs = {
                    'rxlr': '06_candidates', 'crn': '03_candidates',
                    'nlp': '04_nlp', 'protease': '05_protease',
                    'scp': '06_scp', 'elicitin': '07_elicitin', 'yxsl': '08_yxsl',
                }
                out = os.path.join(sample_output_dir, etype, type_outdirs.get(etype, '06_candidates'))
                logger.info(f"{etype.upper()}结果|{etype.upper()} results: {out}/")
            logger.info(f"{'='*60}")

            sys.exit(0 if all_success else 1)

        else:
            # 多样本模式: 按样品独立运行|Multi-sample mode: run independently per sample
            _check_duplicate_stems(fasta_files)

            os.makedirs(args.output_dir, exist_ok=True)

            # 顶层logger|Top-level logger
            top_log_file = os.path.join(args.output_dir, '99_logs', 'pipeline.log')
            master_logger = PhytoEffectorLogger(log_file=top_log_file).get_logger()

            master_logger.info("=" * 60)
            master_logger.info("多样本效应子鉴定|Multi-sample effector identification")
            master_logger.info(f"样本数|Sample count: {len(fasta_files)}")
            master_logger.info("=" * 60)

            sample_results = {}
            failed_samples = []
            total_samples = len(fasta_files)

            for i, fasta_file in enumerate(fasta_files):
                sample_name = Path(fasta_file).stem
                sample_dir = os.path.join(args.output_dir, sample_name)

                master_logger.info(
                    f"处理样本|Processing sample [{i+1}/{total_samples}]: {sample_name}"
                )

                try:
                    results = _run_single_sample(fasta_file, sample_dir, args)
                    sample_results[sample_name] = results
                    master_logger.info(f"样本完成|Sample completed: {sample_name}")
                except Exception as e:
                    master_logger.error(f"样本处理失败|Sample processing failed: {sample_name}: {e}")
                    failed_samples.append(sample_name)

            # 合并各效应子类型的候选TSSV|Merge candidate TSVs across samples
            master_logger.info("合并各样本候选结果|Merging candidate results across samples")
            sample_dirs = [
                os.path.join(args.output_dir, Path(f).stem) for f in fasta_files
                if Path(f).stem not in failed_samples
            ]
            for etype in ALL_TYPES:
                merge_dir = os.path.join(args.output_dir, etype)
                merge_candidate_files(sample_dirs, etype, merge_dir, master_logger)

            # 汇总报告|Summary report
            elapsed = time.time() - start_time
            master_logger.info("=" * 60)
            master_logger.info(f"多样本效应子鉴定完成|Multi-sample effector identification completed in {elapsed:.1f}秒|seconds")
            master_logger.info(f"成功|Succeeded: {len(sample_results)}/{total_samples}")
            if failed_samples:
                master_logger.warning(f"失败|Failed samples: {', '.join(failed_samples)}")

            for sample_name, results in sample_results.items():
                status = ', '.join(
                    f"{k.upper()}={'OK' if v else 'FAIL'}" for k, v in results.items()
                )
                master_logger.info(f"  {sample_name}: {status}")
            master_logger.info("=" * 60)

            # 顶层版本信息|Top-level version info
            _generate_top_summary(args.output_dir, sample_results, fasta_files, args, start_time)

            sys.exit(0 if sample_results else 1)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)
