"""
选择性扫荡检测流水线模块|Selective Sweep Detection Pipeline Module
"""

import itertools
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from .utils import (CommandRunner, build_conda_command, check_dependencies,
                    generate_software_versions_yml)


def parse_pop_info(path: str) -> List[Tuple[str, str]]:
    """解析群体分组文件(样品ID<TAB>分组,无表头)|Parse population info file

    Args:
        path: pop_info文件路径|pop_info file path

    Returns:
        [(sample_id, group)] 保持文件顺序|preserve file order

    Raises:
        ValueError: 列数不足/空字段/重复样品,收集所有错误一次抛出
    """
    errors = []
    entries = []
    seen_samples = set()
    with open(path, encoding='utf-8') as fh:
        for lineno, line in enumerate(fh, start=1):
            # 每字段单独strip:样品ID尾随空白会被bcftools -S静默丢弃;
            # 分组名尾随空白会生成含空格文件名(命名违规)
            # |strip each field: trailing whitespace in a sample ID would be
            # silently dropped by bcftools -S; trailing whitespace in a group
            # name would create space-containing file names
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = [f.strip() for f in line.split('\t')]
            if len(fields) < 2:
                errors.append(
                    f"第{lineno}行列数不足(需要2列,tab分隔)|Line {lineno}: "
                    f"expected 2 tab-separated cols, got {len(fields)}")
                continue
            sample, group = fields[0], fields[1]
            if not sample or not group:
                errors.append(f"第{lineno}行存在空字段|Line {lineno}: empty sample or group")
                continue
            if sample in seen_samples:
                errors.append(f"重复样品|Duplicate sample: {sample}")
                continue
            seen_samples.add(sample)
            entries.append((sample, group))
    if errors:
        raise ValueError("\n".join(errors))
    return entries


def group_samples(entries: List[Tuple[str, str]]) -> Dict[str, List[str]]:
    """按群体分组样本(保持首次出现顺序)|Group samples by population (first-seen order)"""
    groups = {}
    for sample, group in entries:
        groups.setdefault(group, []).append(sample)
    return groups


def pairwise_combinations(pops: List[str]) -> List[Tuple[str, str]]:
    """所有两两组合,固定顺序|All pairwise combinations in fixed order"""
    return list(itertools.combinations(pops, 2))


class SweepPipeline:
    """选择性扫荡检测流水线|Selective Sweep Detection Pipeline"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.runner = CommandRunner(logger, Path(config.output_dir))
        self.pops: List[str] = []
        self.pop_samples: Dict[str, List[str]] = {}

    def _is_step_completed(self, output_file) -> bool:
        """断点续传:输出文件存在即跳过|Resume: skip if output exists"""
        return Path(output_file).exists()

    # ---- 1. 过滤|filter ----
    def run_filter(self) -> Path:
        """bcftools过滤:双等位SNP+MAF+缺失率|Filter: biallelic SNP + MAF + missing rate"""
        filtered_vcf = Path(self.config.filter_dir) / 'filtered.vcf.gz'
        if self._is_step_completed(filtered_vcf):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 过滤|filter")
            return filtered_vcf

        self.logger.info("开始步骤|Starting step: 过滤|filter")
        tmp_tagged = Path(self.config.tmp_dir) / 'filtered.tagged.vcf.gz'

        # 填充MAF/F_MISSING标签|fill MAF/F_MISSING tags
        fill_cmd = build_conda_command(
            self.config.bcftools_path,
            ['+fill-tags', self.config.input_vcf, '-Oz', '-o', str(tmp_tagged),
             '--', '-t', 'MAF,F_MISSING'])
        if not self.runner.run(fill_cmd, 'bcftools填充标签|bcftools fill-tags'):
            raise RuntimeError("bcftools fill-tags 失败|bcftools fill-tags failed")

        # 过滤:双等位SNP+MAF+缺失率|filter: biallelic SNP + MAF + missing
        expr = f'MAF>={self.config.min_maf} & F_MISSING<{self.config.max_missing}'
        view_cmd = build_conda_command(
            self.config.bcftools_path,
            ['view', '-m2', '-M2', '-v', 'snps', '-i', expr,
             str(tmp_tagged), '-Oz', '-o', str(filtered_vcf)])
        if not self.runner.run(view_cmd, 'bcftools过滤|bcftools view filter'):
            raise RuntimeError("bcftools view 过滤失败|bcftools view filter failed")

        # 建索引|index
        index_cmd = build_conda_command(self.config.bcftools_path,
                                        ['index', '-f', str(filtered_vcf)])
        if not self.runner.run(index_cmd, 'bcftools索引|bcftools index'):
            raise RuntimeError("bcftools index 失败|bcftools index failed")

        # 过滤结果为空则报错并提示降低阈值(spec§9)|error if filter emptied the VCF
        n_variants = self._count_variants(filtered_vcf)
        if n_variants == 0:
            raise RuntimeError(
                "过滤后VCF为空,请降低MAF/缺失率阈值或检查输入"
                "|Filtered VCF is empty; lower MAF/missing thresholds or check input")

        # 清理过滤中间文件(tmp下仅filtered.tagged.vcf.gz;删文件不删目录)
        # |clean filter intermediates (only the tagged VCF lives in tmp; files only, keep dir)
        try:
            for p in self.config.tmp_dir.iterdir():
                if p.is_file():
                    p.unlink()
            self.logger.info("已清理过滤中间文件|Filter intermediates cleaned")
        except OSError as e:
            self.logger.warning(
                f"清理过滤中间文件失败,不影响结果|Failed to clean filter intermediates: {e}")

        self.logger.info(f"过滤完成|Filtering completed: {filtered_vcf}")
        return filtered_vcf

    def _count_variants(self, vcf: Path) -> int:
        """统计VCF变异数(用于过滤后空检查)|Count variants in VCF (for empty check)

        失败返回-1(无法判断时放行,不阻塞主流程)|-1 on failure (pass through)
        """
        if not vcf.exists():
            self.logger.warning(f"变异统计输入不存在,放行|Variants-count input missing; pass through: {vcf}")
            return -1
        cmd = build_conda_command(self.config.bcftools_path, ['view', '-H', str(vcf)])
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, shell=False, capture_output=True,
                                    text=True, check=True)
            return sum(1 for _ in result.stdout.splitlines())
        except Exception as e:
            self.logger.warning(f"统计变异数失败|Failed to count variants: {e}")
            return -1

    # ---- 2. 拆分|split ----
    def _vcf_samples(self, vcf: Path):
        """获取VCF头中样本集合(用于pop_info样本校验)|Sample set from VCF header

        Returns:
            set 样本集合|set of samples; 失败返回None(校验放行,不阻塞)|None on
            failure (validation passes through, never blocks)
        """
        if not vcf.exists():
            self.logger.warning(f"样本校验输入缺失,放行|Sample-check input missing; pass through: {vcf}")
            return None
        cmd = build_conda_command(self.config.bcftools_path, ['query', '-l', str(vcf)])
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, shell=False, capture_output=True,
                                    text=True, check=True)
            return set(result.stdout.splitlines())
        except Exception as e:
            self.logger.warning(f"获取VCF样本列表失败,跳过样本校验|Failed to get VCF samples; skip validation: {e}")
            return None

    def split_pops(self, filtered_vcf: Path) -> Dict[str, List[str]]:
        """按群体拆分样本列表与子VCF|Split per-pop sample lists and sub-VCFs"""
        entries = parse_pop_info(self.config.pop_info)
        # pop_info中的样本需真实存在于VCF:表头行/拼写错误的样本会触发bcftools 255
        # |validate against VCF header: header rows or typo'd samples make bcftools -S fail
        vcf_samples = self._vcf_samples(filtered_vcf)
        valid_entries = []
        dropped_by_group = {}
        for sample, group in entries:
            if vcf_samples is not None and sample not in vcf_samples:
                self.logger.warning(
                    f"样本不在VCF中,已跳过|Sample not in VCF header, skipped: {sample} (group {group})")
                dropped_by_group[group] = dropped_by_group.get(group, 0) + 1
                continue
            valid_entries.append((sample, group))
        self.pop_samples = group_samples(valid_entries)
        self.pops = list(self.pop_samples.keys())
        if not self.pops:
            raise RuntimeError(
                "所有群体样本均不在VCF中,请检查pop_info与VCF样本名是否一致"
                "|No sample in pop_info exists in the VCF; check sample names against the VCF header")
        for group, n_dropped in dropped_by_group.items():
            if group not in self.pop_samples:
                self.logger.warning(
                    f"群体{group}的全部样本均不在VCF中,已跳过该群体"
                    f"|Pop {group} has no valid samples in VCF; pop skipped")

        for pop, samples in self.pop_samples.items():
            sample_file = Path(self.config.filter_dir) / f'{pop}.samples.txt'
            pop_vcf = Path(self.config.filter_dir) / f'{pop}.vcf.gz'
            n = len(samples)

            if not self._is_step_completed(pop_vcf):
                sample_file.write_text('\n'.join(samples) + '\n', encoding='utf-8')
                view_cmd = build_conda_command(
                    self.config.bcftools_path,
                    ['view', '-S', str(sample_file), str(filtered_vcf),
                     '-Oz', '-o', str(pop_vcf)])
                if not self.runner.run(view_cmd, f'提取群体子VCF|Extract pop VCF: {pop}'):
                    raise RuntimeError(f"提取群体子VCF失败|Failed to extract pop VCF: {pop}")
                index_cmd = build_conda_command(self.config.bcftools_path,
                                                ['index', '-f', str(pop_vcf)])
                if not self.runner.run(index_cmd, f'群体VCF索引|Index pop VCF: {pop}'):
                    raise RuntimeError(f"索引失败|Failed to index pop VCF: {pop}")
            else:
                self.logger.info(f"跳过已完成步骤|Skipping completed step: 群体VCF|pop VCF {pop}")

            self.logger.info(f"群体{pop}样本数|Pop {pop} sample count: {n}")
            if n < self.config.raisd_min_samples:
                self.logger.warning(
                    f"群体{pop}样本量{n}偏少(<{self.config.raisd_min_samples}),"
                    f"RAiSD μ统计噪音较大,composite_score中将排除其MU分量"
                    f"(--include-mu-low-n可强制加入)"
                    f"|Pop {pop} has few samples ({n}<{self.config.raisd_min_samples}); "
                    f"RAiSD mu noisy; MU component excluded from score "
                    f"(--include-mu-low-n to force)")
        return self.pop_samples

    # ---- 3. 每群体统计|per-pop stats ----
    def run_pop_stats(self, filtered_vcf: Path) -> None:
        """每群体计算 pi / Tajima's D / RAiSD μ|Per-pop pi, TajimaD, RAiSD mu"""
        for pop in self.pops:
            pop_vcf = Path(self.config.filter_dir) / f'{pop}.vcf.gz'

            # π|pi
            pi_out = Path(self.config.stats_dir) / f'{pop}.windowed.pi'
            if not self._is_step_completed(pi_out):
                pi_cmd = build_conda_command(
                    self.config.vcftools_path,
                    ['--gzvcf', str(pop_vcf),
                     '--window-pi', str(self.config.win),
                     '--window-pi-step', str(self.config.step),
                     '--out', str(Path(self.config.stats_dir) / pop)])
                if not self.runner.run(pi_cmd, f'vcftools窗口π|vcftools windowed pi: {pop}'):
                    raise RuntimeError(f"vcftools π计算失败|vcftools pi failed: {pop}")
            else:
                self.logger.info(f"跳过已完成步骤|Skipping completed step: π {pop}")

            # Tajima's D|Tajima's D
            tajd_out = Path(self.config.stats_dir) / f'{pop}.Tajima.D'
            if not self._is_step_completed(tajd_out):
                tajd_cmd = build_conda_command(
                    self.config.vcftools_path,
                    ['--gzvcf', str(pop_vcf),
                     '--TajimaD', str(self.config.win),
                     '--out', str(Path(self.config.stats_dir) / pop)])
                if not self.runner.run(tajd_cmd, f"vcftools窗口Tajima's D|vcftools Tajima's D: {pop}"):
                    raise RuntimeError(f"vcftools TajimaD失败|vcftools TajimaD failed: {pop}")
            else:
                self.logger.info(f"跳过已完成步骤|Skipping completed step: TajimaD {pop}")

            # RAiSD μ|cwd=stats_dir(报告文件无路径参数)|mu; cwd=stats_dir
            raisd_reports = list(Path(self.config.stats_dir).glob(f'RAiSD_Report.{pop}.*'))
            if raisd_reports:
                self.logger.info(f"跳过已完成步骤|Skipping completed step: RAiSD {pop}")
            else:
                raisd_cmd = build_conda_command(
                    self.config.raisd_path,
                    ['-n', pop, '-I', str(pop_vcf),
                     '-R', '-s', '-f', '-O', '-w', str(self.config.raisd_window)])
                if not self.runner.run(raisd_cmd, f'RAiSD μ统计|RAiSD mu statistic: {pop}',
                                       cwd=Path(self.config.stats_dir)):
                    raise RuntimeError(f"RAiSD运行失败|RAiSD failed: {pop}")
                raisd_reports = list(Path(self.config.stats_dir).glob(f'RAiSD_Report.{pop}.*'))
                if not raisd_reports:
                    raise RuntimeError(f"RAiSD未生成报告|RAiSD produced no reports for pop: {pop}")

            # SweeD CLR|cwd=stats_dir;SweeD不支持.gz输入,先解压到tmp,跑完清理
            sweed_reports = list(Path(self.config.stats_dir).glob(f'SweeD_Report.{pop}.*'))
            if sweed_reports:
                self.logger.info(f"跳过已完成步骤|Skipping completed step: SweeD {pop}")
            else:
                sweed_vcf = Path(self.config.tmp_dir) / f'{pop}.sweed.vcf'
                decompress_cmd = build_conda_command(
                    self.config.bcftools_path,
                    ['view', '-O', 'v', str(pop_vcf), '-o', str(sweed_vcf)])
                if not self.runner.run(decompress_cmd, f'SweeD输入解压|Decompress for SweeD: {pop}'):
                    raise RuntimeError(f"SweeD输入解压失败|Failed to decompress for SweeD: {pop}")
                sweed_args = ['-name', pop, '-input', str(sweed_vcf),
                              '-grid', str(self.config.sweed_grid),
                              '-sampleList', str(Path(self.config.filter_dir) / f'{pop}.samples.txt'),
                              '-reports', '-threads', str(self.config.threads)]
                if self.config.sweed_folded:
                    sweed_args.append('-folded')
                sweed_cmd = build_conda_command(self.config.sweed_path, sweed_args)
                if not self.runner.run(sweed_cmd, f'SweeD CLR统计|SweeD CLR statistic: {pop}',
                                       cwd=Path(self.config.stats_dir)):
                    raise RuntimeError(f"SweeD运行失败|SweeD failed: {pop}")
                try:
                    sweed_vcf.unlink()
                    self.logger.info(f"已清理SweeD解压文件|SweeD decompressed input cleaned: {sweed_vcf}")
                except OSError as e:
                    self.logger.warning(
                        f"清理SweeD解压文件失败,不影响结果|Failed to clean SweeD input: {e}")
                sweed_reports = list(Path(self.config.stats_dir).glob(f'SweeD_Report.{pop}.*'))
                if not sweed_reports:
                    raise RuntimeError(f"SweeD未生成报告|SweeD produced no reports for pop: {pop}")

    # ---- 4. Fst两两组合|pairwise Fst ----
    def run_fst(self, filtered_vcf: Path) -> None:
        """每两群体计算weir-fst(自动全部组合)|Pairwise weir-fst (all combinations)"""
        if len(self.pops) < 2:
            self.logger.info("群体数<2,跳过Fst|Fewer than 2 pops; skipping Fst")
            return
        for a, b in pairwise_combinations(self.pops):
            label = f'{a}_{b}'
            fst_out = Path(self.config.stats_dir) / f'{label}.windowed.weir.fst'
            if self._is_step_completed(fst_out):
                self.logger.info(f"跳过已完成步骤|Skipping completed step: Fst {label}")
                continue
            fst_cmd = build_conda_command(
                self.config.vcftools_path,
                ['--gzvcf', str(filtered_vcf),
                 '--weir-fst-pop', str(Path(self.config.filter_dir) / f'{a}.samples.txt'),
                 '--weir-fst-pop', str(Path(self.config.filter_dir) / f'{b}.samples.txt'),
                 '--fst-window-size', str(self.config.win),
                 '--fst-window-step', str(self.config.step),
                 '--out', str(Path(self.config.stats_dir) / label)])
            if not self.runner.run(fst_cmd, f'vcftools窗口Fst|vcftools windowed Fst: {a} vs {b}'):
                raise RuntimeError(f"vcftools Fst失败|vcftools Fst failed: {label}")

    # ---- 5. 元数据|metadata ----
    def write_pop_summary(self) -> None:
        """写群体样本量摘要(供merger低样本MU排除)|Write pop sample summary"""
        out_file = Path(self.config.info_dir) / 'pop_summary.tsv'
        lines = ['pop\tn_samples']
        for pop in self.pops:
            lines.append(f'{pop}\t{len(self.pop_samples[pop])}')
        out_file.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        self.logger.info(f"群体摘要|Pop summary -> {out_file}")

    # ---- 主入口|main entry ----
    def run(self) -> None:
        """运行全流程|Run the full pipeline"""
        self.logger.info("开始选择性扫荡统计计算|Starting selective sweep statistics pipeline")
        if not check_dependencies(self.config, self.logger):
            raise RuntimeError("依赖软件检查失败|Dependency check failed")

        filtered_vcf = self.run_filter()
        self.split_pops(filtered_vcf)
        self.run_pop_stats(filtered_vcf)
        self.run_fst(filtered_vcf)
        self.write_pop_summary()
        generate_software_versions_yml(self.config, self.logger)
        self.logger.info("流水线完成|Pipeline completed")
