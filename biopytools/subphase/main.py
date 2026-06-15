"""SubPhaser亚基因组分离主程序|SubPhaser subgenome phasing main program"""

import glob
import os
import subprocess
import sys
import time
from pathlib import Path

from .config import SubPhaserConfig, SUBGENOME_LABELS
from .utils import (
    SubPhaserLogger,
    build_conda_command,
    generate_auto_config,
    generate_phased_config,
    assign_chr_names,
    build_subgenome_groups,
    find_chrom_subgenome_file,
    parse_chrom_subgenome,
    parse_fasta_info,
    write_phasing_result,
    write_renamed_fasta,
    validate_with_parents,
)


class SubPhaserRunner:
    """SubPhaser运行器|SubPhaser Runner"""

    def __init__(self, **kwargs):
        self.config = SubPhaserConfig(**kwargs)
        self.config.validate()

        log_dir = os.path.join(self.config.output_dir, "99_logs")
        Path(log_dir).mkdir(parents=True, exist_ok=True)
        log_file = os.path.join(log_dir, "subphaser.log")

        self.logger_manager = SubPhaserLogger(log_file=log_file)
        self.logger = self.logger_manager.get_logger()

        self.temp_config = None
        self.chroms = []

    def run(self):
        """运行SubPhaser|Run SubPhaser"""
        start_time = time.time()
        mode_label = "自动模式|Auto mode" if self.config.auto_mode else "配置模式|Config mode"
        self.logger.info("=" * 60)
        self.logger.info(f"SubPhaser亚基因组分离流程 ({mode_label})")
        self.logger.info(f"SubPhaser subgenome phasing pipeline ({mode_label})")
        self.logger.info("=" * 60)

        # 自动模式：生成临时配置文件
        if self.config.auto_mode:
            self._prepare_auto_config()

        # 运行SubPhaser
        cmd = self._build_command()
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                stdout=sys.stdout,
                stderr=sys.stderr,
                shell=False,
            )

            elapsed = time.time() - start_time
            if result.returncode == 0:
                self.logger.info(f"SubPhaser运行完成|SubPhaser completed ({elapsed:.1f}s)")
                self._record_versions()
                self._postprocess()
            else:
                self.logger.error(
                    f"SubPhaser运行失败|SubPhaser failed (exit code: {result.returncode})"
                )
                sys.exit(result.returncode)

        except FileNotFoundError:
            self.logger.error("subphaser命令未找到|subphaser command not found")
            self.logger.error(
                f"请确认conda环境 '{self.config.conda_env}' 已安装SubPhaser"
                f"|Please ensure SubPhaser is installed in conda env '{self.config.conda_env}'"
            )
            sys.exit(1)
        except KeyboardInterrupt:
            self.logger.warning("用户中断|User interrupted")
            sys.exit(130)
        finally:
            self._cleanup_temp()

    def _prepare_auto_config(self):
        """自动模式：解析基因组并生成临时配置|Auto mode: parse genome and generate temp config"""
        self.logger.info("自动模式：生成临时配置文件|Auto mode: generating temp config")

        self.chroms = parse_fasta_info(self.config.genomes[0])
        self.logger.info(
            f"检测到 {len(self.chroms)} 条染色体|"
            f"Detected {len(self.chroms)} chromosomes"
        )

        for name, length in self.chroms:
            self.logger.info(f"  {name}: {length:,} bp")

        temp_dir = os.path.join(self.config.output_dir, "00_auto_config")
        Path(temp_dir).mkdir(parents=True, exist_ok=True)

        self.temp_config = os.path.join(temp_dir, "auto_sg.config")
        filtered_chroms = [
            (name, length) for name, length in self.chroms
            if length >= self.config.min_chrom_size
        ] if self.config.min_chrom_size > 0 else self.chroms

        if len(filtered_chroms) < self.config.nsg:
            self.logger.error(
                f"过滤后仅剩 {len(filtered_chroms)} 条染色体 (< nsg={self.config.nsg})，"
                f"请降低 --min-chrom-size|"
                f"Only {len(filtered_chroms)} chromosomes left after filtering, "
                f"lower --min-chrom-size"
            )
            sys.exit(1)

        self.logger.info(
            f"过滤后保留 {len(filtered_chroms)} 条染色体 "
            f"(阈值: {self.config.min_chrom_size:,} bp)|"
            f"Kept {len(filtered_chroms)} chromosomes "
            f"(threshold: {self.config.min_chrom_size:,} bp)"
        )

        generate_auto_config(filtered_chroms, self.config.nsg, self.temp_config)
        self.config.temp_config = self.temp_config

        self.logger.info(f"临时配置已生成|Temp config generated: {self.temp_config}")

    def _build_command(self):
        """构建SubPhaser命令|Build SubPhaser command"""
        args = []

        # 必需参数|Required parameters
        args.extend(["-i"] + self.config.genomes)
        config_files = self.config.sg_cfgs if self.config.sg_cfgs else [self.temp_config]
        args.extend(["-c"] + config_files)

        # 输出|Output
        args.extend(["-o", self.config.output_dir])

        if self.config.prefix:
            args.extend(["-pre", self.config.prefix])

        # 资源|Resources
        if self.config.threads:
            args.extend(["-p", str(self.config.threads)])

        # 亚基因组数量|Number of subgenomes
        args.extend(["-nsg", str(self.config.nsg)])

        # K-mer参数|K-mer parameters
        if self.config.kmer_size != 15:
            args.extend(["-k", str(self.config.kmer_size)])
        if self.config.min_fold != 2.0:
            args.extend(["-f", str(self.config.min_fold)])
        if self.config.min_freq != 200:
            args.extend(["-q", str(self.config.min_freq)])

        # 聚类参数|Cluster parameters
        if self.config.max_pval != 0.05:
            args.extend(["-max_pval", str(self.config.max_pval)])
        if self.config.replicates != 1000:
            args.extend(["-replicates", str(self.config.replicates)])
        if self.config.test_method != "ttest_ind":
            args.extend(["-test_method", self.config.test_method])

        # 步骤控制|Step control
        if self.config.disable_ltr:
            args.append("-disable_ltr")
        if self.config.disable_circos:
            args.append("-disable_circos")
        if self.config.disable_blocks:
            args.append("-disable_blocks")
        if self.config.just_core:
            args.append("-just_core")

        # LTR参数|LTR parameters
        if self.config.ltr_detectors:
            args.extend(["-ltr_detectors"] + self.config.ltr_detectors)
        if self.config.mu != 13e-9:
            args.extend(["-mu", str(self.config.mu)])

        # Circos参数|Circos parameters
        if self.config.window_size != 1000000:
            args.extend(["-window_size", str(self.config.window_size)])
        if self.config.aligner != "minimap2":
            args.extend(["-aligner", self.config.aligner])

        # 高级参数|Advanced parameters
        if self.config.sg_assigned:
            args.extend(["-sg_assigned", self.config.sg_assigned])
        if self.config.target:
            args.extend(["-target", self.config.target])
        if self.config.labels:
            args.extend(["-labels"] + self.config.labels)
        if self.config.no_label:
            args.append("-no_label")
        if self.config.custom_features:
            args.extend(["-custom_features"] + self.config.custom_features)
        if self.config.figfmt != "pdf":
            args.extend(["-figfmt", self.config.figfmt])

        # 覆盖/清理|Overwrite/cleanup
        if self.config.overwrite:
            args.append("-overwrite")
        if self.config.cleanup:
            args.append("-cleanup")

        return build_conda_command(self.config.conda_env, "subphaser", args)

    def _postprocess(self):
        """后处理：解析结果并生成命名|Post-process: parse results and generate naming"""
        tsv_path = find_chrom_subgenome_file(
            self.config.output_dir, self.config.prefix
        )
        if not tsv_path:
            self.logger.warning(
                "未找到chrom-subgenome.tsv，跳过后处理|"
                "chrom-subgenome.tsv not found, skipping post-processing"
            )
            return

        assignments = parse_chrom_subgenome(tsv_path)
        self.logger.info(
            f"Phasing结果|Phasing result: {len(assignments)} chromosomes -> "
            f"{len(set(assignments.values()))} subgenomes"
        )

        if not self.chroms:
            self.chroms = parse_fasta_info(self.config.genomes[0])

        groups = build_subgenome_groups(assignments, self.chroms)
        for sg, chrom_list in sorted(groups.items()):
            names = [n for n, _ in chrom_list]
            self.logger.info(f"  {sg}: {', '.join(names)}")

        # 验证模式：用父本基因组确定A/B标签
        if self.config.parental_genomes:
            self.logger.info("验证模式：比对父本基因组|Validation: aligning to parental genomes")
            label_map = validate_with_parents(
                groups, self.config.parental_genomes,
                self.config.genomes[0], self.config.threads,
            )
            groups = _relabel_groups(groups, label_map)
            assignments = _relabel_assignments(assignments, label_map)
            self.logger.info(
                f"父本验证完成|Parental validation: {label_map}"
            )

        # 生成Chr1A/B/C命名
        name_map = assign_chr_names(groups, self.config.nsg)
        for original, new in sorted(name_map.items(), key=lambda x: x[1]):
            self.logger.info(f"  {original} -> {new}")

        # 输出重命名FASTA
        phased_dir = os.path.join(self.config.output_dir, "phased_genome")
        phased_fasta = write_renamed_fasta(
            self.config.genomes[0], name_map, phased_dir
        )
        self.logger.info(f"重命名FASTA已输出|Renamed FASTA: {phased_fasta}")

        # 输出phasing结果汇总
        result_tsv = os.path.join(self.config.output_dir, "phasing_result.tsv")
        write_phasing_result(groups, name_map, assignments, self.chroms, result_tsv)
        self.logger.info(f"Phasing结果汇总|Phasing result: {result_tsv}")

        # 生成最终配置文件
        config_path = os.path.join(self.config.output_dir, "phased_sg.config")
        generate_phased_config(groups, name_map, self.config.nsg, config_path)
        self.logger.info(f"最终配置文件|Final config: {config_path}")

    def _cleanup_temp(self):
        """清理临时文件|Clean up temporary files"""
        pass

    def _record_versions(self):
        """记录软件版本信息|Record software version info"""
        info_dir = os.path.join(self.config.output_dir, "00_pipeline_info")
        Path(info_dir).mkdir(parents=True, exist_ok=True)

        version_file = os.path.join(info_dir, "software_versions.yml")
        try:
            cmd = build_conda_command(
                self.config.conda_env, "subphaser", ["-v"]
            )
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            version = result.stdout.strip() if result.returncode == 0 else "unknown"
        except Exception:
            version = "unknown"

        mode = "auto" if self.config.auto_mode else "config"
        content = (
            f"pipeline:\n"
            f"  name: biopytools subphaser\n"
            f"  version: 1.0.0\n"
            f"  mode: {mode}\n"
            f"\n"
            f"tools:\n"
            f"  subphaser:\n"
            f"    version: '{version}'\n"
            f"    conda_env: {self.config.conda_env}\n"
            f"\n"
            f"parameters:\n"
            f"  genomes: {self.config.genomes}\n"
            f"  nsg: {self.config.nsg}\n"
            f"  auto_mode: {self.config.auto_mode}\n"
            f"  kmer_size: {self.config.kmer_size}\n"
            f"  min_fold: {self.config.min_fold}\n"
            f"  min_freq: {self.config.min_freq}\n"
            f"  threads: {self.config.threads}\n"
            f"  parental_genomes: {self.config.parental_genomes}\n"
        )
        with open(version_file, 'w') as f:
            f.write(content)
        self.logger.info(f"版本信息已记录|Version info saved: {version_file}")


def _relabel_groups(
    groups: dict, label_map: dict
) -> dict:
    """用父本验证结果重新标记亚基因组|Relabel subgenome groups with parental validation"""
    new_groups = {}
    for sg_key, chrom_list in groups.items():
        new_label = label_map.get(sg_key, sg_key)
        new_groups[new_label] = chrom_list
    return new_groups


def _relabel_assignments(
    assignments: dict, label_map: dict
) -> dict:
    """用父本验证结果重新标记分配|Relabel assignments with parental validation"""
    return {chrom: label_map.get(sg, sg) for chrom, sg in assignments.items()}


def main():
    """命令行入口|Command line entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="SubPhaser亚基因组分离工具|SubPhaser subgenome phasing tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # 必需参数|Required
    parser.add_argument("-i", "--genomes", nargs="+", required=True,
                        help="基因组FASTA文件|Genome FASTA files")
    parser.add_argument("--nsg", type=int, required=True,
                        help="亚基因组数量|Number of subgenomes")

    # 输入模式|Input mode
    parser.add_argument("-c", "--sg-cfgs", nargs="+", default=None,
                        help="亚基因组配置文件(可选)|Subgenome config files (optional)")
    parser.add_argument("--parental-genomes", nargs=2, default=None,
                        help="父本基因组(验证模式)|Parental genomes (validation mode)")

    # 输出|Output
    parser.add_argument("-o", "--output-dir", default="./subphaser_output",
                        help="输出目录|Output directory")
    parser.add_argument("--prefix", default=None,
                        help="输出前缀|Output prefix")

    # 资源|Resources
    parser.add_argument("-t", "--threads", type=int, default=24,
                        help="线程数|Number of threads")

    # 染色体过滤|Chromosome filtering
    parser.add_argument("--min-chrom-size", type=int, default=1_000_000,
                        help="最小染色体长度(bp)|Min chromosome size (bp)")

    # K-mer参数|K-mer
    parser.add_argument("-k", "--kmer-size", type=int, default=15,
                        help="K-mer大小|K-mer size")
    parser.add_argument("-f", "--min-fold", type=float, default=2.0,
                        help="最小倍数差异|Minimum fold difference")
    parser.add_argument("-q", "--min-freq", type=int, default=200,
                        help="最小k-mer频率|Minimum k-mer frequency")

    # 聚类|Cluster
    parser.add_argument("--max-pval", type=float, default=0.05,
                        help="最大P值|Maximum P-value")
    parser.add_argument("--replicates", type=int, default=1000,
                        help="Bootstrap重复次数|Bootstrap replicates")
    parser.add_argument("--test-method", default="ttest_ind",
                        choices=["ttest_ind", "kruskal", "wilcoxon", "mannwhitneyu"],
                        help="统计检验方法|Statistical test method")

    # 步骤控制|Step control
    parser.add_argument("--disable-ltr", action="store_true",
                        help="禁用LTR分析|Disable LTR analysis")
    parser.add_argument("--disable-circos", action="store_true",
                        help="禁用Circos图|Disable Circos plot")
    parser.add_argument("--disable-blocks", action="store_true",
                        help="禁用同源区块|Disable homologous blocks")
    parser.add_argument("--just-core", action="store_true",
                        help="仅运行核心phasing|Only run core phasing")

    # LTR参数|LTR
    parser.add_argument("--ltr-detectors", nargs="+", default=None,
                        choices=["ltr_finder", "ltr_harvest"],
                        help="LTR检测工具|LTR detection tools")
    parser.add_argument("--mu", type=float, default=13e-9,
                        help="替换率/年|Substitution rate per year")

    # Circos参数|Circos
    parser.add_argument("--window-size", type=int, default=1000000,
                        help="Circos窗口大小|Circos window size")
    parser.add_argument("--aligner", default="minimap2",
                        choices=["minimap2", "unimap"],
                        help="比对工具|Aligner")

    # 高级参数|Advanced
    parser.add_argument("--sg-assigned", default=None,
                        help="已知亚基因组分配文件|Pre-assigned subgenome file")
    parser.add_argument("--target", default=None,
                        help="目标染色体文件|Target chromosomes file")
    parser.add_argument("--labels", nargs="+", default=None,
                        help="基因组标签|Genome labels")
    parser.add_argument("--no-label", action="store_true",
                        help="不添加标签前缀|No label prefix")
    parser.add_argument("--custom-features", nargs="+", default=None,
                        help="自定义特征FASTA|Custom feature FASTA files")
    parser.add_argument("--figfmt", default="pdf", choices=["pdf", "png"],
                        help="图片格式|Figure format")

    # 其他|Other
    parser.add_argument("--overwrite", action="store_true",
                        help="覆盖已有结果|Overwrite existing results")
    parser.add_argument("--cleanup", action="store_true",
                        help="清理临时文件|Clean up temporary files")
    parser.add_argument("--conda-env", default="SubPhaser",
                        help="conda环境名称|Conda environment name")

    args = parser.parse_args()

    kwargs = vars(args)
    runner = SubPhaserRunner(**kwargs)
    runner.run()


if __name__ == "__main__":
    main()
