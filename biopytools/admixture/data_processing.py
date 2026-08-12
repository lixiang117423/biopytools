"""
 ADMIXTURE数据处理模块|ADMIXTURE Data Processing Module

VCF预处理(去多等位/染色体重命名)与PLINK转换/质控/LD剪枝/染色体编号修复。
所有外部命令走 list + build_conda_command + shell=False(§13);
管道命令改Python实现,避免shell与conda-run-in-pipe(§13.2.1);
各步骤支持断点续传(§10.2)。
"""

import os
import re
import gzip
import shutil
from collections import Counter
import pandas as pd
from .utils import CommandRunner, build_conda_command


def _open_vcf(path: str):
    """打开VCF文件(自动处理.gz)|Open VCF file (auto-handle .gz)"""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def _iter_vcf_data_fields(path: str):
    """迭代VCF数据行(跳过#注释头)的字段列表|Iterate VCF data rows (skip # header) as field lists"""
    with _open_vcf(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if parts and parts[0]:
                yield parts


def _is_step_done(path: str, force: bool = False) -> bool:
    """断点续传:输出文件存在且非空则视为完成(§10.2)|Checkpoint: done if output exists and non-empty"""
    return (not force) and os.path.exists(path) and os.path.getsize(path) > 0


class VCFProcessor:
    """ VCF处理器|VCF Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def preprocess_vcf(self):
        """ 预处理VCF文件|Preprocess VCF file

        Returns:
            预处理后的VCF路径;skip_preprocessing时返回None|Preprocessed VCF path; None when skipped
        """
        if self.config.skip_preprocessing:
            self.logger.info("跳过VCF预处理|Skipping VCF preprocessing")
            return None

        vcf_file = self.config.vcf_file
        out_dir = self.config.preprocessing_dir

        self.logger.info(f"预处理VCF文件|Preprocessing VCF file: {vcf_file}")

        #  检查VCF文件格式|Check VCF file format
        self._check_vcf_format(vcf_file)

        #  移除多等位基因位点|Remove multi-allelic sites
        biallelic_vcf = os.path.join(out_dir, "biallelic.vcf.gz")
        self._remove_multiallelic(vcf_file, biallelic_vcf)

        #  检查染色体命名并重命名|Check chromosome naming and rename
        renamed_vcf = self._check_and_rename_chromosomes(biallelic_vcf)

        return renamed_vcf

    def _check_vcf_format(self, vcf_file: str):
        """ 检查VCF文件格式(读取前20行)|Check VCF file format (read first 20 lines)"""
        header_lines = []
        try:
            with _open_vcf(vcf_file) as f:
                for i, line in enumerate(f):
                    if i >= 20:
                        break
                    header_lines.append(line.rstrip('\n'))
        except Exception as e:
            raise FileNotFoundError(f"读取VCF失败|Failed to read VCF: {vcf_file}: {e}")

        self.logger.info(f"VCF文件头信息|VCF header info:\n" + "\n".join(header_lines))

    def _remove_multiallelic(self, input_vcf: str, output_vcf: str):
        """ 移除多等位基因位点|Remove multi-allelic sites (bcftools)"""
        # 断点续传|Checkpoint resume
        if _is_step_done(output_vcf, self.config.force):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 双等位过滤|biallelic filter ({output_vcf})")
            return

        args = ["view", "-m2", "-M2", "-v", "snps", input_vcf, "-Oz", "-o", output_vcf]
        cmd = build_conda_command(self.config.bcftools_path, args)
        self.cmd_runner.run(cmd, "移除多等位基因位点|Remove multi-allelic sites")

    def _check_and_rename_chromosomes(self, vcf_file: str):
        """ 检查染色体命名并重命名为数字|Check chromosome naming and rename to numbers"""
        #  获取染色体列表|Get chromosome list
        chr_info = self._get_chromosome_info(vcf_file)
        self._current_chr_info = chr_info  #  保存以便后续使用|Save for later use

        #  检查是否需要重命名|Check if renaming is needed
        need_rename = any(not chr_name.isdigit() for chr_name in chr_info.keys())

        if not need_rename:
            self.logger.info(" 染色体已经是数字格式，无需重命名|Chromosomes are already numeric, no renaming needed")
            return vcf_file

        #  生成重命名映射|Generate renaming mapping
        chr_mapping = self._generate_chromosome_mapping(chr_info)

        #  保存染色体对应关系表|Save chromosome mapping table
        self._save_chromosome_mapping(chr_mapping)

        #  执行重命名|Perform renaming
        renamed_vcf = self._rename_chromosomes(vcf_file, chr_mapping)

        #  验证重命名结果|Verify renaming results
        self._verify_renamed_chromosomes(renamed_vcf)

        return renamed_vcf

    def _get_chromosome_info(self, vcf_file: str):
        """ 获取染色体信息(CHROM->SNP数)|Get chromosome information (CHROM -> SNP count)"""
        counter = Counter()
        for parts in _iter_vcf_data_fields(vcf_file):
            counter[parts[0]] += 1

        chr_info = dict(counter)

        self.logger.info(f"染色体信息|Chromosome info:")
        for chr_name, count in chr_info.items():
            self.logger.info(f"{chr_name}: {count} SNPs")

        return chr_info

    def _generate_chromosome_mapping(self, chr_info: dict):
        """ 生成染色体重命名映射|Generate chromosome renaming mapping"""
        #  对染色体进行自然排序|Natural sorting of chromosomes
        sorted_chrs = self._natural_sort_chromosomes(list(chr_info.keys()))

        #  生成映射：原始名称 -> 数字|Generate mapping: original name -> number
        chr_mapping = {}
        for i, chr_name in enumerate(sorted_chrs, 1):
            chr_mapping[chr_name] = str(i)

        self.logger.info("染色体重命名映射|Chromosome renaming mapping:")
        for original, new in chr_mapping.items():
            self.logger.info(f"  {original} -> {new}")

        return chr_mapping

    def _natural_sort_chromosomes(self, chr_list):
        """ 自然排序染色体|Natural sort chromosomes"""
        def natural_key(chr_name):
            # 提取数字部分进行排序|Extract numeric part for sorting
            numbers = re.findall(r'\d+', chr_name)
            if numbers:
                return [int(num) for num in numbers]
            else:
                return [float('inf')]  # 非数字部分排在最后
        return sorted(chr_list, key=natural_key)

    def _save_chromosome_mapping(self, chr_mapping: dict):
        """ 保存染色体对应关系表|Save chromosome mapping table"""
        out_dir = self.config.preprocessing_dir
        mapping_file = os.path.join(out_dir, "chromosome_mapping.txt")

        # 获取染色体信息|Get chromosome info
        chr_info = getattr(self, '_current_chr_info', {})

        with open(mapping_file, 'w') as f:
            f.write("#  染色体重命名对应关系表|Chromosome Renaming Mapping Table\n")
            f.write("# 格式|Format: 原始染色体编号 -> 重命名后编号|Original -> Renamed\n")
            f.write("# 说明|Note: 用于ADMIXTURE分析的染色体重命名|Chromosome renaming for ADMIXTURE analysis\n")
            f.write("Original_Chr\tRenamed_Chr\tSNP_Count\n")
            for original, renamed in sorted(chr_mapping.items(), key=lambda x: int(x[1])):
                snp_count = chr_info.get(original, 0)
                f.write(f"{original}\t{renamed}\t{snp_count}\n")

        # bcftools需要的重命名文件|Renaming file for bcftools
        rename_file = os.path.join(out_dir, "chr_rename.txt")
        with open(rename_file, 'w') as f:
            for original, renamed in chr_mapping.items():
                f.write(f"{original}\t{renamed}\n")

        self.logger.info(f"染色体对应关系表已保存|Chromosome mapping table saved: {mapping_file}")
        self.logger.info(f"染色体重命名文件已保存|Chromosome rename file saved: {rename_file}")
        return mapping_file, rename_file

    def _rename_chromosomes(self, vcf_file: str, chr_mapping: dict):
        """ 重命名染色体|Rename chromosomes (bcftools annotate)"""
        out_dir = self.config.preprocessing_dir
        renamed_vcf = os.path.join(out_dir, "renamed_chromosomes.vcf.gz")
        rename_file = os.path.join(out_dir, "chr_rename.txt")

        # 断点续传|Checkpoint resume
        if _is_step_done(renamed_vcf, self.config.force):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 染色体重命名|chr rename ({renamed_vcf})")
            return renamed_vcf

        args = ["annotate", "--rename-chrs", rename_file, vcf_file, "-Oz", "-o", renamed_vcf]
        cmd = build_conda_command(self.config.bcftools_path, args)
        self.cmd_runner.run(cmd, "重命名染色体为数字|Rename chromosomes to numbers")

        self.logger.info(f"染色体重命名完成|Chromosome renaming completed: {renamed_vcf}")
        return renamed_vcf

    def _verify_renamed_chromosomes(self, renamed_vcf: str):
        """ 验证重命名结果|Verify renaming results"""
        self.logger.info("验证染色体重命名结果|Verifying chromosome renaming results")

        chr_counter = Counter()
        for parts in _iter_vcf_data_fields(renamed_vcf):
            chr_counter[parts[0]] += 1

        self.logger.info("重命名后染色体信息|Renamed chromosome info:")
        for chr_name in sorted(chr_counter, key=lambda x: (int(x) if x.isdigit() else float('inf'))):
            self.logger.info(f"  {chr_name}: {chr_counter[chr_name]} SNPs")

        # 检查是否所有染色体都是数字|Check if all chromosomes are numeric
        non_numeric = [c for c in chr_counter if not c.isdigit()]
        if non_numeric:
            self.logger.warning(f"警告：仍有非数字染色体|Warning: Still have non-numeric chromosomes: {non_numeric}")
        else:
            self.logger.info("所有染色体已成功重命名为数字|All chromosomes successfully renamed to numbers")


class PlinkProcessor:
    """ PLINK处理器|PLINK Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def convert_vcf_to_plink(self, vcf_file: str):
        """ VCF转换为PLINK格式|Convert VCF to PLINK format"""
        output_prefix = os.path.join(self.config.plink_dir, "raw_data")

        # 断点续传|Checkpoint resume
        if _is_step_done(f"{output_prefix}.bed", self.config.force):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: VCF转PLINK|VCF->PLINK ({output_prefix}.bed)")
            return output_prefix

        #  检测染色体数量|Detect chromosome count
        chr_count = self._detect_chromosome_count(vcf_file)

        # 构建参数列表(§13: list + build_conda_command)|Build arg list (§13)
        args = [
            "--vcf", vcf_file,
            "--make-bed", "--out", output_prefix,
            "--allow-extra-chr", "--double-id",
            "--threads", str(self.config.threads),
        ]

        # 如果染色体数量确定且小于95，添加autosome-num参数
        if chr_count and chr_count < 95:
            args += ["--autosome-num", str(chr_count)]
            self.logger.info(f"设置autosome-num参数为 {chr_count}|Set autosome-num parameter to {chr_count}")
        else:
            self.logger.info("未设置autosome-num参数|autosome-num parameter not set")

        cmd = build_conda_command(self.config.plink_path, args)
        self.cmd_runner.run(cmd, "VCF转换为PLINK格式|Convert VCF to PLINK format")
        return output_prefix

    def quality_control(self, input_prefix: str):
        """质量控制|Quality control"""
        output_prefix = os.path.join(self.config.plink_dir, self.config.base_name)

        # 断点续传|Checkpoint resume
        if _is_step_done(f"{output_prefix}.bed", self.config.force):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 质控|QC ({output_prefix}.bed)")
            return output_prefix

        #  记录质控前统计信息
        self._log_qc_stats(input_prefix, "质控前 (Before QC)")

        #  MAF过滤|MAF filtering
        maf_prefix = f"{output_prefix}_maf"
        cmd_maf = build_conda_command(self.config.plink_path, [
            "--bfile", input_prefix, "--maf", str(self.config.maf),
            "--make-bed", "--out", maf_prefix,
            "--allow-extra-chr", "--threads", str(self.config.threads),
        ])
        self.cmd_runner.run(cmd_maf, "MAF过滤|MAF filtering")

        #  HWE过滤|HWE filtering
        hwe_prefix = f"{output_prefix}_hwe"
        cmd_hwe = build_conda_command(self.config.plink_path, [
            "--bfile", maf_prefix, "--hwe", str(self.config.hwe_pvalue),
            "--make-bed", "--out", hwe_prefix,
            "--allow-extra-chr", "--threads", str(self.config.threads),
        ])
        self.cmd_runner.run(cmd_hwe, "HWE过滤|HWE filtering")

        #  缺失率过滤|Missing rate filtering
        cmd_missing = build_conda_command(self.config.plink_path, [
            "--bfile", hwe_prefix,
            "--geno", str(self.config.missing_rate), "--mind", str(self.config.missing_rate),
            "--make-bed", "--out", output_prefix,
            "--allow-extra-chr", "--threads", str(self.config.threads),
        ])
        self.cmd_runner.run(cmd_missing, "缺失率过滤|Missing rate filtering")

        #  记录质控后统计信息
        self._log_qc_stats(output_prefix, "质控后 (After QC)")

        #  清理中间文件
        if not self.config.keep_intermediate:
            self._cleanup_intermediate_files([maf_prefix, hwe_prefix])

        return output_prefix

    def ld_prune(self, input_prefix: str) -> str:
        """
        LD剪枝：去除高度连锁SNP，保留代表性SNP|LD pruning: remove linked SNPs, keep representative ones

        使用PLINK --indep-pairwise 去除高度连锁的SNP（参考ADMIXTURE手册§2.3）。
        Uses PLINK --indep-pairwise to remove highly linked SNPs (ref: ADMIXTURE manual §2.3).

        Args:
            input_prefix: 输入PLINK文件前缀（已质控）|Input PLINK prefix (already QC'd)

        Returns:
            剪枝后的PLINK文件前缀；未开启时原样返回|Pruned PLINK prefix; unchanged when disabled
        """
        # 未开启LD剪枝时直接返回|Return directly when LD pruning is disabled
        if not self.config.ld_prune:
            self.logger.info("已关闭LD剪枝，跳过|LD pruning disabled, skipping")
            return input_prefix

        output_prefix = f"{input_prefix}_ld"

        # 断点续传|Checkpoint resume
        if _is_step_done(f"{output_prefix}.bed", self.config.force):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: LD剪枝|LD pruning ({output_prefix}.bed)")
            return output_prefix

        self.logger.info("LD剪枝：去除高度连锁SNP|LD pruning: removing linked SNPs")

        # 记录剪枝前样本/SNP数|Log sample/SNP counts before pruning
        self._log_qc_stats(input_prefix, "LD剪枝前|Before LD pruning")

        # 步骤1：计算连锁，生成prune.in/prune.out|Step 1: compute LD, generate prune.in/prune.out
        prune_prefix = f"{output_prefix}_prune"
        cmd_prune = build_conda_command(self.config.plink_path, [
            "--bfile", input_prefix,
            "--keep-allele-order", "--real-ref-alleles",
            "--allow-extra-chr", "--double-id",
            "--threads", str(self.config.threads),
            "--indep-pairwise",
            self.config.ld_window, str(self.config.ld_step), str(self.config.ld_r2),
            "--out", prune_prefix,
        ])
        self.cmd_runner.run(cmd_prune, "LD剪枝计算|LD pruning computation")

        # 步骤2：提取非连锁SNP|Step 2: extract unlinked SNPs
        cmd_extract = build_conda_command(self.config.plink_path, [
            "--bfile", input_prefix,
            "--extract", f"{prune_prefix}.prune.in",
            "--keep-allele-order", "--real-ref-alleles",
            "--allow-extra-chr", "--double-id",
            "--threads", str(self.config.threads),
            "--make-bed", "--out", output_prefix,
        ])
        self.cmd_runner.run(cmd_extract, "提取非连锁SNP|Extracting unlinked SNPs")

        # 记录剪枝后样本/SNP数|Log sample/SNP counts after pruning
        self._log_qc_stats(output_prefix, "LD剪枝后|After LD pruning")

        # 清理中间文件|Clean intermediates
        if not self.config.keep_intermediate:
            self._cleanup_intermediate_files([prune_prefix])
            for ext in ['prune.in', 'prune.out']:
                prune_file = f"{prune_prefix}.{ext}"
                try:
                    if os.path.exists(prune_file):
                        os.remove(prune_file)
                except OSError as e:
                    self.logger.warning(f"清理失败|Failed to clean {prune_file}: {e}")

        return output_prefix

    def fix_chromosome_codes(self, plink_prefix: str):
        """ 修复染色体编号为整数格式|Fix chromosome codes to integer format"""
        self.logger.info("修复染色体编号为整数格式|Fixing chromosome codes to integer format")

        fixed_prefix = os.path.join(self.config.plink_dir, "admixture_chr_fixed")

        # 断点续传|Checkpoint resume
        if _is_step_done(f"{fixed_prefix}.bed", self.config.force):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 染色体编号修复|chr code fix ({fixed_prefix}.bed)")
            return fixed_prefix

        #  读取bim文件|Read bim file
        bim_file = f"{plink_prefix}.bim"
        if not os.path.exists(bim_file):
            raise FileNotFoundError(f"BIM文件不存在|BIM file not found: {bim_file}")

        bim_df = pd.read_csv(bim_file, sep='\t', header=None,
                            names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])

        # 统一转字符串处理|Unify to string handling
        unique_chrs = [str(chr_code) for chr_code in bim_df['CHR'].unique()]

        def chr_sort_key(chr_code):
            if chr_code.isdigit():
                return (0, int(chr_code))
            return (1, chr_code)

        sorted_chrs = sorted(unique_chrs, key=chr_sort_key)
        self.logger.info(f"当前染色体编号|Current chromosome codes: {sorted_chrs}")

        # 已全为数字则直接返回原前缀|All numeric -> return original prefix
        non_numeric_chrs = [c for c in unique_chrs if not c.isdigit()]
        if not non_numeric_chrs:
            self.logger.info("所有染色体编号已为数字格式|All chromosome codes are already numeric")
            return plink_prefix

        self.logger.warning(f"发现非数字染色体编号|Found non-numeric chromosome codes: {non_numeric_chrs}")

        # 获取染色体SNP计数信息|Get chromosome SNP count info
        chr_counts = {}
        for chr_code in unique_chrs:
            original_chr = None
            for orig_chr in bim_df['CHR'].unique():
                if str(orig_chr) == chr_code:
                    original_chr = orig_chr
                    break
            if original_chr is not None:
                chr_counts[chr_code] = (bim_df['CHR'] == original_chr).sum()
            else:
                chr_counts[chr_code] = 0

        # 创建染色体编号映射|Create chromosome code mapping
        chr_mapping = {}
        numeric_start = 1
        for chr_code in sorted_chrs:
            if chr_code.isdigit():
                chr_mapping[chr_code] = chr_code
                numeric_start = max(numeric_start, int(chr_code) + 1)
        for chr_code in sorted(non_numeric_chrs):
            chr_mapping[chr_code] = str(numeric_start)
            numeric_start += 1

        self.logger.info(f"染色体编号映射|Chromosome code mapping: {chr_mapping}")

        # 保存染色体对应关系表|Save chromosome mapping table
        self._save_chromosome_mapping_plink(chr_mapping, chr_counts)

        #  复制bed和fam文件(直接用shutil,系统命令不走conda)|Copy bed and fam files via shutil
        for ext in ['bed', 'fam']:
            src_file = f"{plink_prefix}.{ext}"
            dst_file = f"{fixed_prefix}.{ext}"
            try:
                shutil.copy(src_file, dst_file)
                self.logger.info(f"复制{ext.upper()}文件|Copy {ext.upper()} file: {dst_file}")
            except OSError as e:
                raise OSError(f"复制{ext}失败|Failed to copy {ext}: {src_file} -> {dst_file}: {e}")

        #  修改bim文件中的染色体编号|Modify chromosome codes in bim file
        mapping_dict = {}
        for orig_chr in bim_df['CHR'].unique():
            str_chr = str(orig_chr)
            if str_chr in chr_mapping:
                mapping_dict[orig_chr] = int(chr_mapping[str_chr])

        bim_df['CHR'] = bim_df['CHR'].map(mapping_dict)
        bim_df.to_csv(f"{fixed_prefix}.bim", sep='\t', header=False, index=False)

        #  验证修复结果|Verify fix results
        self._verify_chromosome_fix(fixed_prefix)

        self.logger.info(f"染色体编号修复完成|Chromosome code fix completed: {fixed_prefix}")
        return fixed_prefix

    def _save_chromosome_mapping_plink(self, chr_mapping: dict, chr_counts: dict):
        """ 保存染色体对应关系表（PLINK版本）| Save chromosome mapping table (PLINK version)"""
        mapping_file = os.path.join(self.config.plink_dir, "chromosome_mapping.txt")

        with open(mapping_file, 'w') as f:
            f.write("#   染色体重命名对应关系表|Chromosome Renaming Mapping Table\n")
            f.write("# 格式|Format: 原始染色体编号 -> 重命名后编号|Original -> Renamed\n")
            f.write("# 说明|Note: 用于ADMIXTURE分析的染色体重命名|Chromosome renaming for ADMIXTURE analysis\n")
            f.write("# 来源|Source: PLINK文件染色体编号修复|PLINK file chromosome code fix\n")
            f.write("Original_Chr\tRenamed_Chr\tSNP_Count\n")
            for original, renamed in sorted(chr_mapping.items(), key=lambda x: int(str(x[1]))):
                snp_count = chr_counts.get(original, 0)
                f.write(f"{original}\t{renamed}\t{snp_count}\n")

        rename_file = os.path.join(self.config.plink_dir, "chr_rename.txt")
        with open(rename_file, 'w') as f:
            for original, renamed in chr_mapping.items():
                f.write(f"{original}\t{renamed}\n")

        self.logger.info(f"染色体对应关系表已保存|Chromosome mapping table saved: {mapping_file}")
        self.logger.info(f"染色体重命名文件已保存|Chromosome rename file saved: {rename_file}")
        return mapping_file, rename_file

    def _verify_chromosome_fix(self, plink_prefix: str):
        """ 验证染色体编号修复结果|Verify chromosome code fix results"""
        bim_file = f"{plink_prefix}.bim"
        bim_df = pd.read_csv(bim_file, sep='\t', header=None,
                            names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])

        unique_chrs = sorted(bim_df['CHR'].unique())
        self.logger.info(f" 修复后染色体编号|Fixed chromosome codes: {unique_chrs}")

        non_numeric = [c for c in unique_chrs if not str(c).isdigit()]
        if non_numeric:
            self.logger.error(f"仍有非数字染色体编号|Still have non-numeric chromosome codes: {non_numeric}")
            raise ValueError("染色体编号修复失败|Chromosome code fix failed")
        else:
            self.logger.info("所有染色体编号已成功修复为数字格式|All chromosome codes successfully fixed to numeric format")

    def _detect_chromosome_count(self, vcf_file: str):
        """ 检测染色体数量|Detect chromosome count (Python实现,失败返回None)"""
        try:
            chrs = set()
            for parts in _iter_vcf_data_fields(vcf_file):
                chrs.add(parts[0])
            chr_count = len(chrs)
            self.logger.info(f"检测到 {chr_count} 个染色体|Detected {chr_count} chromosomes")
            return chr_count
        except Exception as e:
            self.logger.warning(f"无法检测染色体数量|Cannot detect chromosome count: {e}")
            return None

    def _log_qc_stats(self, prefix: str, stage: str):
        """ 记录质控统计信息|Log QC statistics"""
        try:
            # 统计样本数|Count samples
            fam_file = f"{prefix}.fam"
            if os.path.exists(fam_file):
                with open(fam_file, 'r') as f:
                    sample_count = len(f.readlines())
                self.logger.info(f"{stage}样本数|{stage} samples: {sample_count}")

            # 统计SNP数|Count SNPs
            bim_file = f"{prefix}.bim"
            if os.path.exists(bim_file):
                with open(bim_file, 'r') as f:
                    snp_count = len(f.readlines())
                self.logger.info(f"{stage}SNP数|{stage} SNPs: {snp_count}")

        except Exception as e:
            self.logger.warning(f"无法获取{stage}统计信息|Cannot get {stage} statistics: {e}")

    def _cleanup_intermediate_files(self, prefixes: list):
        """清理中间文件|Cleanup intermediate files"""
        try:
            for prefix in prefixes:
                for ext in ['bed', 'bim', 'fam', 'log', 'nosex']:
                    file_path = f"{prefix}.{ext}"
                    if os.path.exists(file_path):
                        os.remove(file_path)
            self.logger.info("清理中间文件完成|Intermediate files cleaned")
        except Exception as e:
            self.logger.warning(f"清理中间文件时出错|Error cleaning intermediate files: {e}")
