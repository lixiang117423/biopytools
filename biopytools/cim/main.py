"""
CIM分析主程序模块|CIM Analysis Main Module
"""

import argparse
import os
import sys
from datetime import datetime

from .config import CIMConfig
from .utils import CIMLogger, CommandRunner
from .converter import (
    parse_vcf_genotypes,
    filter_markers,
    filter_rf_markers,
    fix_geno_error,
    ld_prune_markers,
    save_filtered_vcf,
    build_tidy_files,
    build_mstmap_linkage_map,
    generate_r_cim_script,
    run_r_script
)


def save_pipeline_info(config: CIMConfig, output_file: str):
    """保存流程参数信息|Save pipeline parameters"""
    params = [
        f"# CIM Analysis Pipeline Info",
        f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"[Input]",
        f"VCF: {config.input_file}",
        f"Phenotype: {config.pheno_file}",
        f"",
        f"[Cross]",
        f"Type: {config.cross_type}",
        f"",
        f"[Map]",
        f"Mode: {config.map_mode}",
        f"",
        f"[Filtering]",
        f"MAF threshold: {config.maf}",
        f"Max missing rate: {config.max_missing}",
        f"",
        f"[LD Pruning]",
        f"Window: {config.ld_window}",
        f"Step: {config.ld_step}",
        f"R2 threshold: {config.ld_r2}",
        f"",
        f"[RF QC]",
        f"Max het rate: {config.max_het_rate}",
        f"Max mean RF: {config.max_mean_rf}",
        f"",
        f"[CIM]",
        f"Cofactors (n.marcovar): {config.n_marcovar}",
        f"Window: {config.window}",
        f"Method: {config.method}",
        f"Pseudomarker step: {config.step}",
        f"",
        f"[Permutation]",
        f"Permutations: {config.n_perm}",
        f"",
        f"[MSTmap]",
        f"P-value (initial, auto-tuned): {config.mstmap_pvalue}",
        f"Distance function: {config.mstmap_distfun}",
        f"MSTmap path: {config.mstmap_path}",
        f"",
        f"[Environment]",
        f"R env: {config.r_env}",
    ]
    with open(output_file, 'w') as f:
        f.write('\n'.join(params) + '\n')


def save_qc_stats(stats: dict, output_file: str):
    """保存QC统计信息|Save QC statistics"""
    lines = []
    for key, val in stats.items():
        lines.append(f"{key}\t{val}")
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="R/qtl复合区间作图(CIM)分析工具|R/qtl Composite Interval Mapping (CIM) Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  python -m biopytools.cim -i input.vcf.gz -p phe.txt -o output_dir
  python -m biopytools.cim -i input.vcf.gz -p phe.txt -o output_dir --map-mode estimate --n-perm 1000
  python -m biopytools.cim -i input.vcf.gz -p phe.txt -o output_dir --n-marcovar 5 --window 10 --method hk
        """
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input', required=True,
                        help='输入VCF文件|Input VCF file (supports .vcf and .vcf.gz)')
    parser.add_argument('-p', '--pheno', required=True,
                        help='表型文件|Phenotype file (TSV with sample, value columns)')
    parser.add_argument('-o', '--output', required=True,
                        help='输出目录|Output directory')

    # 群体类型|Cross type
    parser.add_argument('-t', '--type', dest='cross_type',
                        choices=['f2', 'bc'], default='f2',
                        help='群体类型|Cross type (f2/bc, default: f2)')

    # 遗传图谱|Genetic map
    parser.add_argument('--map-mode', choices=['physical', 'estimate', 'mstmap'],
                        default='mstmap',
                        help='cM位置来源|cM position source (physical/estimate/mstmap, default: mstmap)')

    # 标记过滤|Marker filtering
    parser.add_argument('--maf', type=float, default=0.05,
                        help='最小等位基因频率阈值|Minor allele frequency threshold (default: 0.05)')
    parser.add_argument('--missing', type=float, default=0.1,
                        help='最大缺失率|Maximum missing rate (default: 0.1)')

    # CIM参数|CIM parameters
    parser.add_argument('--n-marcovar', type=int, default=10,
                        help='协因子数量|Number of marker covariates (default: 10)')
    parser.add_argument('--window', type=float, default=10.0,
                        help='窗口大小(cM)|Window size in cM (default: 10)')
    parser.add_argument('--method', choices=['hk', 'em', 'imp'], default='hk',
                        help='扫描方法|Scanning method (default: hk)')
    parser.add_argument('--step', type=float, default=1.0,
                        help='伪标记步长(cM)|Pseudomarker step in cM (default: 1)')

    # 置换检验|Permutation test
    parser.add_argument('--n-perm', type=int, default=1000,
                        help='置换检验次数(0=跳过)|Number of permutations (0=skip, default: 1000)')

    # LD pruning参数|LD pruning parameters
    parser.add_argument('--ld-window', type=int, default=50,
                        help='LD计算窗口(SNP数)|LD window in SNP count (default: 50)')
    parser.add_argument('--ld-step', type=int, default=5,
                        help='LD计算步长(SNP数)|LD step in SNP count (default: 5)')
    parser.add_argument('--ld-r2', type=float, default=0.1,
                        help='LD r2阈值|LD r2 threshold (default: 0.1)')
    parser.add_argument('--skip-ld', action='store_true',
                        help='跳过LD降维步骤|Skip LD pruning')

    # 基因型纠错|Genotype error correction
    parser.add_argument('--fix-geno-error-size', type=int, default=10,
                        help='fixGenoError短片段阈值(0=跳过)|Minimum run length for fixGenoError (0=skip, default: 10)')

    # MSTmap参数|MSTmap parameters
    parser.add_argument('--mstmap-pvalue', type=float, default=1e-6,
                        help='MSTmap聚类p值起始值(自动调优)|MSTmap clustering p-value start value, auto-tuned (default: 1e-6)')
    parser.add_argument('--mstmap-distfun', choices=['kosambi', 'haldane'], default='kosambi',
                        help='MSTmap距离函数|MSTmap distance function (default: kosambi)')
    parser.add_argument('--mstmap-path', default='~/miniforge3/envs/Rqtl/bin/mstmap',
                        help='MSTmap二进制路径|MSTmap binary path (default: ~/miniforge3/envs/Rqtl/bin/mstmap)')

    # 环境参数|Environment parameters
    parser.add_argument('--r-env', default='Rqtl',
                        help='R conda环境名|R conda environment name (default: Rqtl)')
    parser.add_argument('--threads', type=int, default=1,
                        help='并行线程数|Number of parallel threads (default: 1)')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = CIMConfig(
            input_file=args.input,
            pheno_file=args.pheno,
            output_dir=args.output,
            cross_type=args.cross_type,
            map_mode=args.map_mode,
            maf=args.maf,
            max_missing=args.missing,
            n_marcovar=args.n_marcovar,
            window=args.window,
            method=args.method,
            step=args.step,
            n_perm=args.n_perm,
            ld_window=args.ld_window,
            ld_step=args.ld_step,
            ld_r2=args.ld_r2,
            skip_ld=args.skip_ld,
            fix_geno_error_size=args.fix_geno_error_size,
            r_env=args.r_env,
            threads=args.threads,
            mstmap_pvalue=args.mstmap_pvalue,
            mstmap_distfun=args.mstmap_distfun,
            mstmap_path=args.mstmap_path,
        )
        config.validate()

        # 创建日志|Create logger
        logger_manager = CIMLogger(log_file=config.log_file)
        logger = logger_manager.get_logger()
        cmd_runner = CommandRunner(logger)

        logger.info("=" * 60)
        logger.info("CIM Analysis Pipeline|R/qtl复合区间作图分析流程")
        logger.info("=" * 60)

        # 保存流程参数|Save pipeline info
        save_pipeline_info(config, os.path.join(config.info_dir, "pipeline_params.txt"))

        # Step 1: 解析VCF和表型|Parse VCF and phenotype
        genotype_matrix, marker_info, samples, pheno_values = parse_vcf_genotypes(
            config.input_file, config.pheno_file, logger
        )

        # Step 2: 标记过滤|Filter markers
        genotype_matrix, marker_info, filter_stats = filter_markers(
            genotype_matrix, marker_info, config.maf, config.max_missing, logger
        )
        save_qc_stats(filter_stats, os.path.join(config.qc_dir, "marker_filter_stats.txt"))

        if genotype_matrix.shape[0] == 0:
            logger.error("过滤后无标记保留，请检查MAF和缺失率参数|No markers left after filtering")
            sys.exit(1)

        # Step 3: LD降维|LD pruning
        if config.skip_ld:
            logger.info("跳过LD降维（--skip-ld）|Skipping LD pruning (--skip-ld)")
            prune_stats = {'markers_before_prune': genotype_matrix.shape[0],
                           'markers_after_prune': genotype_matrix.shape[0],
                           'pruning_method': 'skipped'}
        else:
            genotype_matrix, marker_info, prune_stats = ld_prune_markers(
                genotype_matrix, marker_info,
                config.qc_dir, config.ld_window, config.ld_step, config.ld_r2,
                cmd_runner, logger
            )
        save_qc_stats(prune_stats, os.path.join(config.qc_dir, "ld_prune_stats.txt"))

        # 保存LD降维后的VCF|Save VCF after LD pruning
        filtered_vcf_path = os.path.join(config.qc_dir, "filtered_markers.vcf.gz")
        if os.path.exists(filtered_vcf_path):
            logger.info("跳过已完成步骤|Skipping completed step: 过滤后VCF保存|filtered VCF save")
        else:
            save_filtered_vcf(config.input_file, marker_info, samples, filtered_vcf_path, cmd_runner, logger)

        # 对齐样本数|Align sample counts (filter_markers may remove all-NaN sample columns)
        if 'keep_sample_mask' in filter_stats:
            keep = filter_stats['keep_sample_mask']
            samples = [s for s, k in zip(samples, keep) if k]
            pheno_values = pheno_values[keep]

        # Step 3.5: 基因型纠错|Genotype error correction
        if config.fix_geno_error_size > 0:
            logger.info("-" * 60)
            logger.info("步骤3.5: 基因型纠错|Step 3.5: Genotype error correction")
            logger.info("-" * 60)
            genotype_matrix, marker_info, gec_stats = fix_geno_error(
                genotype_matrix, marker_info, config.fix_geno_error_size, logger
            )
            save_qc_stats(gec_stats, os.path.join(config.qc_dir, "geno_error_correction_stats.txt"))

        # Step 4: Pre-RF CIM分析（基线）|Pre-RF CIM analysis (baseline)
        logger.info("-" * 60)
        logger.info("步骤4: Pre-RF CIM分析 (基线)|Step 4: Pre-RF CIM analysis (baseline)")
        logger.info("-" * 60)

        pre_rf_tidy = build_tidy_files(
            genotype_matrix, marker_info, samples, pheno_values,
            config.pre_rf_tidy_dir, config.map_mode, logger
        )

        if config.map_mode == "mstmap":
            pre_rf_mstmap_map = os.path.join(config.pre_rf_mstmap_dir, 'mstmap_map.csv')
            if os.path.exists(pre_rf_mstmap_map):
                logger.info("跳过已完成步骤|Skipping completed step: MSTmap遗传图谱构建(Pre-RF)")
            else:
                logger.info("步骤4a: MSTmap遗传图谱构建|MSTmap linkage map construction")
                logger.info("-" * 60)
                build_mstmap_linkage_map(
                    genotype_matrix, marker_info, samples, pheno_values,
                    pre_rf_tidy, config.pre_rf_mstmap_dir, config,
                    cmd_runner, logger
                )

        pre_rf_done_marker = os.path.join(config.pre_rf_plots_dir, 'cim_lod_data_mstmap.tsv') \
            if config.map_mode == "mstmap" else os.path.join(config.pre_rf_plots_dir, 'cim_lod_data.tsv')
        if os.path.exists(pre_rf_done_marker):
            logger.info("跳过已完成步骤|Skipping completed step: Pre-RF CIM R脚本执行")
            pre_rf_success = True
        else:
            pre_rf_script = generate_r_cim_script(config, pre_rf_tidy, config.pre_rf_dir)
            logger.info(f"Pre-RF R脚本已生成|R script generated: {pre_rf_script}")

            pre_rf_success = run_r_script(pre_rf_script, config.r_env, cmd_runner, logger, timeout=None)

            if not pre_rf_success:
                logger.warning("Pre-RF CIM分析失败，继续RF过滤和Post-RF分析|"
                               "Pre-RF CIM failed, continuing with RF filtering and post-RF analysis")

        # Step 5: RF质控过滤|RF quality filtering
        logger.info("-" * 60)
        logger.info("步骤5: RF质控过滤|Step 5: RF quality filtering")
        logger.info("-" * 60)

        genotype_matrix, marker_info, rf_stats = filter_rf_markers(
            genotype_matrix, marker_info, config.max_het_rate, config.max_mean_rf, logger,
            qc_dir=config.qc_dir
        )
        rf_stats_file = os.path.join(config.qc_dir, "rf_filter_stats.txt")
        save_qc_stats(rf_stats, rf_stats_file)

        if genotype_matrix.shape[0] == 0:
            logger.error("RF过滤后无标记保留，请调整max_het_rate和max_mean_rf参数|"
                         "No markers left after RF filtering")
            sys.exit(1)

        # 校验样本数一致性|RF filtering不改变列数，若不一致则说明逻辑错误
        if genotype_matrix.shape[1] != len(samples):
            logger.error(f"样本数不匹配|Sample count mismatch: "
                         f"genotype={genotype_matrix.shape[1]}, phenotype={len(samples)}")
            sys.exit(1)

        # 保存RF过滤后的VCF|Save VCF after RF filtering
        rf_filtered_vcf_path = os.path.join(config.qc_dir, "rf_filtered_markers.vcf.gz")
        if os.path.exists(rf_filtered_vcf_path):
            logger.info("跳过已完成步骤|Skipping completed step: RF过滤后VCF保存|RF filtered VCF save")
        else:
            save_filtered_vcf(config.input_file, marker_info, samples, rf_filtered_vcf_path, cmd_runner, logger)

        # Step 6: Post-RF CIM分析（最终结果）|Post-RF CIM analysis (final results)
        logger.info("-" * 60)
        logger.info("步骤6: Post-RF CIM分析 (最终结果)|Step 6: Post-RF CIM analysis (final)")
        logger.info("-" * 60)

        post_rf_tidy = build_tidy_files(
            genotype_matrix, marker_info, samples, pheno_values,
            config.post_rf_tidy_dir, config.map_mode, logger
        )

        if config.map_mode == "mstmap":
            post_rf_mstmap_map = os.path.join(config.post_rf_mstmap_dir, 'mstmap_map.csv')
            if os.path.exists(post_rf_mstmap_map):
                logger.info("跳过已完成步骤|Skipping completed step: MSTmap遗传图谱构建(Post-RF)")
            else:
                logger.info("步骤6a: MSTmap遗传图谱构建|MSTmap linkage map construction")
                logger.info("-" * 60)
                build_mstmap_linkage_map(
                    genotype_matrix, marker_info, samples, pheno_values,
                    post_rf_tidy, config.post_rf_mstmap_dir, config,
                    cmd_runner, logger
                )

        post_rf_done_marker = os.path.join(config.post_rf_plots_dir, 'cim_lod_data_mstmap.tsv') \
            if config.map_mode == "mstmap" else os.path.join(config.post_rf_plots_dir, 'cim_lod_data.tsv')
        if os.path.exists(post_rf_done_marker):
            logger.info("跳过已完成步骤|Skipping completed step: Post-RF CIM R脚本执行")
            post_rf_success = True
        else:
            post_rf_script = generate_r_cim_script(config, post_rf_tidy, config.post_rf_dir)
            logger.info(f"Post-RF R脚本已生成|R script generated: {post_rf_script}")

            post_rf_success = run_r_script(post_rf_script, config.r_env, cmd_runner, logger, timeout=None)

        if post_rf_success:
            logger.info("=" * 60)
            logger.info("CIM分析完成!|CIM Analysis Complete!")
            logger.info(f"结果目录|Output directory: {config.output_dir}")
            logger.info(f"  - 过滤后VCF: {filtered_vcf_path}")
            logger.info(f"  - RF过滤后VCF: {rf_filtered_vcf_path}")
            logger.info(f"  - RF过滤统计: {os.path.join(config.qc_dir, 'rf_filter_stats.txt')}")
            if pre_rf_success:
                logger.info(f"  - Pre-RF结果目录: {config.pre_rf_dir}")
            if config.map_mode == "mstmap":
                logger.info(f"  - Post-RF物理位置CIM峰值: {os.path.join(config.post_rf_plots_dir, 'cim_peaks_physical.tsv')}")
                logger.info(f"  - Post-RF物理位置CIM LOD数据: {os.path.join(config.post_rf_plots_dir, 'cim_lod_data_physical.tsv')}")
                logger.info(f"  - Post-RF遗传图谱: {os.path.join(config.post_rf_mstmap_dir, 'linkage_map.csv')}")
                logger.info(f"  - Post-RF MSTmap CIM峰值: {os.path.join(config.post_rf_plots_dir, 'cim_peaks_mstmap.tsv')}")
                logger.info(f"  - Post-RF MSTmap CIM LOD数据: {os.path.join(config.post_rf_plots_dir, 'cim_lod_data_mstmap.tsv')}")
            else:
                logger.info(f"  - Post-RF QTL峰值表: {os.path.join(config.post_rf_plots_dir, 'cim_peaks.tsv')}")
                logger.info(f"  - Post-RF LOD数据: {os.path.join(config.post_rf_plots_dir, 'cim_lod_data.tsv')}")
                logger.info(f"  - Post-RF全基因组LOD图: {os.path.join(config.post_rf_dir, 'cim_genome_plot.pdf')}")
            logger.info("=" * 60)
            sys.exit(0)
        else:
            logger.error("Post-RF CIM分析失败|Post-RF CIM Analysis Failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("用户中断|User interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
