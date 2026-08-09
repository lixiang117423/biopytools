"""
VCF过滤器模块|VCF Filter Module
"""

from pathlib import Path

from .utils import build_conda_command


class VCFFilter:
    """VCF文件过滤器|VCF File Filter"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 定义输出文件名|Define output file names
        self.filtered_snp_file = Path(f"{self.config.base_name}.filtered.snp.vcf.gz")
        self.filtered_indel_file = Path(f"{self.config.base_name}.filtered.indel.vcf.gz")
        self.merged_file = Path(f"{self.config.base_name}.filtered.merged.vcf.gz")

        # 双等位位点过滤后的文件|Biallelic filtered SNP file
        self.biallelic_snp_file = Path(f"{self.config.base_name}.filtered.snp.biallelic.vcf.gz")

    def filter_snps(self, input_file: Path) -> bool:
        """过滤SNP|Filter SNPs"""
        self.logger.info("=" * 60)
        self.logger.info("步骤2a: 过滤SNP|Step 2a: Filtering SNPs")
        self.logger.info("=" * 60)

        # 构建过滤表达式(文献常用标准)|Build filter expression (literature standards)
        filters = [
            f"QUAL >= {self.config.snp_qual}",
            f"INFO/DP >= {self.config.snp_dp}",
            f"MQ >= {self.config.snp_mq}",
            f"QD >= {self.config.snp_qd}",
            f"FS <= {self.config.snp_fs}",
            f"SOR <= {self.config.snp_sor}",
            f"MQRankSum >= {self.config.snp_mqrs}",
            f"ReadPosRankSum >= {self.config.snp_rprs}",
            f"MAF >= {self.config.snp_maf}",
        ]
        filter_expr = " && ".join(filters)

        self.logger.info("SNP过滤条件|SNP filter criteria:")
        for f in filters:
            self.logger.info(f"  - {f}")

        # 执行过滤|Execute filtering
        filter_args = [
            "view",
            "-i", filter_expr,
            "--threads", str(self.config.threads),
            "-Oz", "-o", str(self.filtered_snp_file),
            str(input_file),
        ]
        filter_cmd = build_conda_command(self.config.bcftools_path, filter_args)
        if not self.cmd_runner.run(filter_cmd, "应用SNP过滤|Applying SNP filters"):
            return False

        # 创建索引|Create index
        index_cmd = build_conda_command(
            self.config.bcftools_path, ["index", "-f", str(self.filtered_snp_file)]
        )
        if not self.cmd_runner.run(index_cmd, "索引过滤后的SNP文件|Indexing filtered SNP file"):
            return False

        self.logger.info("SNP过滤完成|SNP filtering completed")
        return True

    def filter_indels(self, input_file: Path) -> bool:
        """过滤INDEL|Filter INDELs"""
        self.logger.info("=" * 60)
        self.logger.info("步骤2b: 过滤INDEL|Step 2b: Filtering INDELs")
        self.logger.info("=" * 60)

        # INDEL使用更宽松的标准|INDELs use more relaxed standards
        filters = [
            f"QUAL >= {self.config.indel_qual}",
            f"INFO/DP >= {self.config.indel_dp}",
            f"MQ >= {self.config.indel_mq}",
            f"QD >= {self.config.indel_qd}",
            f"FS <= {self.config.indel_fs}",
            f"SOR <= {self.config.indel_sor}",
            f"ReadPosRankSum >= {self.config.indel_rprs}",
        ]
        filter_expr = " && ".join(filters)

        self.logger.info("INDEL过滤条件|INDEL filter criteria:")
        for f in filters:
            self.logger.info(f"  - {f}")

        # 执行过滤|Execute filtering
        filter_args = [
            "view",
            "-i", filter_expr,
            "--threads", str(self.config.threads),
            "-Oz", "-o", str(self.filtered_indel_file),
            str(input_file),
        ]
        filter_cmd = build_conda_command(self.config.bcftools_path, filter_args)
        if not self.cmd_runner.run(filter_cmd, "应用INDEL过滤|Applying INDEL filters"):
            return False

        # 创建索引|Create index
        index_cmd = build_conda_command(
            self.config.bcftools_path, ["index", "-f", str(self.filtered_indel_file)]
        )
        if not self.cmd_runner.run(index_cmd, "索引过滤后的INDEL文件|Indexing filtered INDEL file"):
            return False

        self.logger.info("INDEL过滤完成|INDEL filtering completed")
        return True

    def filter_biallelic_snps(self) -> bool:
        """
        过滤双等位位点SNP|Filter biallelic SNPs

        只保留双等位位点且MAF大于阈值的SNP
        Only keep biallelic SNPs with MAF above threshold
        """
        if not self.config.snp_biallelic:
            self.logger.info("跳过双等位位点过滤|Skipping biallelic filtering (disabled)")
            # 不进行双等位过滤时,最终使用的SNP文件就是filtered_snp_file
            self.biallelic_snp_file = self.filtered_snp_file
            return True

        self.logger.info("=" * 60)
        self.logger.info("步骤2c: 过滤双等位位点SNP|Step 2c: Filtering biallelic SNPs")
        self.logger.info("=" * 60)

        # -m2 -M2: 只保留双等位位点|-m2 -M2: keep only biallelic sites
        maf_threshold = self.config.snp_maf
        self.logger.info("双等位位点过滤条件|Biallelic filter criteria:")
        self.logger.info("  - 双等位位点|Biallelic sites only (max 2 alleles)")
        self.logger.info(f"  - MAF > {maf_threshold}")

        # 执行双等位位点过滤|Execute biallelic filtering
        filter_args = [
            "view",
            "-m2", "-M2",                 # 只保留双等位位点|Only keep biallelic sites
            "-v", "snps",                 # 只保留SNP|Only keep SNPs
            "-i", f"MAF > {maf_threshold}",
            "--threads", str(self.config.threads),
            "-Oz", "-o", str(self.biallelic_snp_file),
            str(self.filtered_snp_file),
        ]
        filter_cmd = build_conda_command(self.config.bcftools_path, filter_args)
        if not self.cmd_runner.run(filter_cmd, "应用双等位位点过滤|Applying biallelic filters"):
            return False

        # 创建索引|Create index
        index_cmd = build_conda_command(
            self.config.bcftools_path, ["index", "-f", str(self.biallelic_snp_file)]
        )
        if not self.cmd_runner.run(index_cmd, "索引双等位位点SNP文件|Indexing biallelic SNP file"):
            return False

        self.logger.info("双等位位点SNP过滤完成|Biallelic SNP filtering completed")
        return True

    def merge_filtered_variants_with_biallelic(self) -> bool:
        """
        合并过滤后的变异(使用双等位位点过滤后的SNP)|
        Merge filtered variants (using biallelic-filtered SNPs)
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤3: 合并过滤后的变异|Step 3: Merging filtered variants")
        self.logger.info("=" * 60)

        # 根据variant_type决定合并策略|Determine merge strategy based on variant_type
        if self.config.variant_type == 'snp_only':
            # 只有SNP,直接复用SNP文件作为合并结果|Only SNPs, reuse SNP file
            self.logger.info("只有SNP变异,跳过合并|Only SNP variants, skipping merge")
            source_file = (self.biallelic_snp_file
                           if self.config.snp_biallelic else self.filtered_snp_file)
            import shutil
            try:
                shutil.copy(str(source_file), str(self.merged_file))
                self.logger.info("复制SNP文件作为合并结果|Copied SNP file as merged result")
            except Exception as e:
                self.logger.error(f"复制文件失败|Failed to copy file: {e}")
                return False

        elif self.config.variant_type == 'indel_only':
            # 只有INDEL,直接复用INDEL文件作为合并结果|Only INDELs, reuse INDEL file
            self.logger.info("只有INDEL变异,跳过合并|Only INDEL variants, skipping merge")
            import shutil
            try:
                shutil.copy(str(self.filtered_indel_file), str(self.merged_file))
                self.logger.info("复制INDEL文件作为合并结果|Copied INDEL file as merged result")
            except Exception as e:
                self.logger.error(f"复制文件失败|Failed to copy file: {e}")
                return False

        else:  # variant_type == 'both'
            # 正常合并SNP和INDEL|Normal merge of SNPs and INDELs
            snp_file_to_merge = (self.biallelic_snp_file
                                 if self.config.snp_biallelic else self.filtered_snp_file)
            if self.config.snp_biallelic:
                self.logger.info("使用双等位位点过滤后的SNP|Using biallelic-filtered SNPs")
            else:
                self.logger.info("使用常规过滤后的SNP|Using regular filtered SNPs")

            # 临时未排序合并文件名(在working_dir下)|Temp unsorted merge filename (under working_dir)
            unsorted_merged_filename = f"{self.config.base_name}.filtered.merged.unsorted.vcf.gz"

            # 步骤1: 合并到临时文件|Step 1: concat into temp file
            merge_args = [
                "concat",
                "-a", "-D",
                "--threads", str(self.config.threads),
                "-Oz", "-o", unsorted_merged_filename,
                str(snp_file_to_merge), str(self.filtered_indel_file),
            ]
            merge_cmd = build_conda_command(self.config.bcftools_path, merge_args)
            if not self.cmd_runner.run(merge_cmd, "合并SNP和INDEL|Merging SNPs and INDELs"):
                return False

            # 步骤2: 排序并输出到最终文件|Step 2: sort into final file
            sort_args = ["sort", "-Oz", "-o", str(self.merged_file), unsorted_merged_filename]
            sort_cmd = build_conda_command(self.config.bcftools_path, sort_args)
            if not self.cmd_runner.run(sort_cmd, "按位置排序|Sorting by position"):
                return False

            # 步骤3: 删除临时未排序文件|Step 3: remove temp unsorted file
            import os
            try:
                temp_file_full_path = self.config.output_path / unsorted_merged_filename
                self.logger.info(f"清理临时文件|Cleaning up temporary file: {temp_file_full_path}")
                os.remove(temp_file_full_path)
                self.logger.info("临时文件已删除|Temporary file deleted")
            except OSError as e:
                self.logger.warning(
                    f"删除临时文件失败(可能已不存在)|Failed to remove temp file (may not exist): {e}"
                )

        # 步骤4: 对最终文件创建索引|Step 4: index final file
        index_cmd = build_conda_command(
            self.config.bcftools_path, ["index", "-f", str(self.merged_file)]
        )
        if not self.cmd_runner.run(index_cmd, "索引合并文件|Indexing merged file"):
            return False

        self.logger.info("变异合并完成|Variant merging completed")
        return True
