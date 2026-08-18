"""
VCF数据处理模块|VCF Data Processing Module
"""

import os
import subprocess
from pathlib import Path
from typing import Dict
from .utils import CommandRunner
from ..common.conda_runner import build_conda_command, build_pipeline_command
from ..common.conda_runner import build_conda_command, build_pipeline_command

class VCFProcessor:
    """VCF文件处理器|VCF File Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def is_compressed_vcf(self, vcf_file: str) -> bool:
        """检查VCF文件是否压缩|Check if VCF file is compressed"""
        return vcf_file.endswith('.gz')
    
    def get_vcf_stats(self, vcf_file: str) -> Dict[str, int]:
        """获取VCF文件统计信息|Get VCF file statistics"""
        self.logger.info("获取VCF文件统计信息|Getting VCF file statistics")
        
        stats = {}
        
        # 获取样本数量|Get sample count
        if self.is_compressed_vcf(vcf_file):
            # 管道(方案B): bcftools解析到域环境二进制, wc为系统工具
            # |Pipeline (solution B): bcftools resolves to domain binary, wc stays bare
            cmd = build_pipeline_command([
                [self.config.bcftools_path, 'query', '-l', vcf_file],
                ['wc', '-l'],
            ])
        else:
            cmd = f"grep '^#CHROM' {vcf_file}|cut -f10-|tr '\\t' '\\n'|wc -l"
        
        try:
            self.logger.info(f"命令|Command: {cmd}")
            self.logger.info(f"命令|Command: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            stats['samples'] = int(result.stdout.strip())
        except:
            stats['samples'] = 0
        
        # 获取变异数量|Get variant count
        if self.is_compressed_vcf(vcf_file):
            # 管道(方案B): bcftools解析到域环境二进制, wc为系统工具
            # |Pipeline (solution B): bcftools resolves to domain binary, wc stays bare
            cmd = build_pipeline_command([
                [self.config.bcftools_path, 'view', '-H', vcf_file],
                ['wc', '-l'],
            ])
        else:
            cmd = f"grep -v '^#' {vcf_file}|wc -l"
        
        try:
            self.logger.info(f"命令|Command: {cmd}")
            self.logger.info(f"命令|Command: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            stats['variants'] = int(result.stdout.strip())
        except:
            stats['variants'] = 0
        
        self.logger.info(f"VCF统计|VCF statistics: {stats['samples']} 样本|samples, {stats['variants']} 变异|variants")
        
        return stats
    
    def vcf_to_plink(self):
        """VCF转换为PLINK格式|Convert VCF to PLINK format"""
        vcf_path = self.config.vcf_file
        output_prefix = self.config.output_path / self.config.base_name
        
        # 确保输出目录存在并获取绝对路径|Ensure output directory exists and get absolute path
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        output_prefix_abs = output_prefix.resolve()
        
        # 获取VCF统计信息|Get VCF statistics
        stats = self.get_vcf_stats(vcf_path)
        self.logger.info(f"开始处理包含 {stats['samples']} 个样本和 {stats['variants']} 个变异的VCF文件")
        
        # 检查目录是否可写|Check if directory is writable
        try:
            test_file = output_prefix.parent / "test_write.tmp"
            test_file.touch()
            test_file.unlink()
            self.logger.info(f"输出目录可写|Output directory is writable: {output_prefix.parent}")
        except Exception as e:
            self.logger.error(f"输出目录不可写|Output directory is not writable: {output_prefix.parent}, 错误: {e}")
            return False
        
        args = [
            "--vcf", vcf_path,
            "--make-bed", "--out", str(output_prefix_abs),
            "--allow-extra-chr", "--double-id",
        ]
        cmd = build_conda_command(self.config.plink_path, args)
        
        success = self.cmd_runner.run(cmd, "VCF转换为PLINK格式|Convert VCF to PLINK format")
        
        if success:
            self.config.plink_prefix = str(output_prefix_abs)
            self.logger.info(f"PLINK文件已生成|PLINK files generated: {output_prefix_abs}")
        
        return success

class QualityController:
    """质量控制器|Quality Controller"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def apply_quality_filters(self):
        """应用质量过滤|Apply quality filters"""
        # 检查是否跳过质控|Check if skipping QC
        if self.config.skip_qc:
            self.logger.info("跳过质量控制过滤 (默认不过滤)|Skipping quality control filtering (no filtering by default)")
            # 默认不过滤,但剔除零基因型样本(100%缺失)避免PLINK --pca崩溃
            # No filtering by default, but remove zero-genotype samples to avoid PLINK crash
            self.config.plink_prefix_qc = self._drop_zero_genotype_samples(self.config.plink_prefix)
            self.log_qc_stats()
            return True
        
        input_prefix = self.config.plink_prefix
        output_prefix = self.config.output_path / f"{self.config.base_name}_qc"
        
        self.logger.info("应用质量控制过滤|Applying quality control filters")
        self.logger.info(f"MAF阈值|MAF threshold: {self.config.maf}")
        self.logger.info(f"缺失率阈值|Missing rate threshold: {self.config.missing_rate}")
        self.logger.info(f"HWE p值阈值|HWE p-value threshold: {self.config.hwe_pvalue}")
        
        # Step 1: MAF过滤|MAF filtering
        maf_prefix = self.config.output_path / f"{self.config.base_name}_maf"
        maf_prefix.parent.mkdir(parents=True, exist_ok=True)
        maf_prefix_abs = maf_prefix.resolve()
        cmd_maf = build_conda_command(self.config.plink_path, [
            "--bfile", str(input_prefix),
            "--maf", str(self.config.maf), "--make-bed",
            "--out", str(maf_prefix_abs), "--allow-extra-chr",
        ])
        
        if not self.cmd_runner.run(cmd_maf, f"MAF过滤|MAF filtering (>={self.config.maf})"):
            return False
        
        # Step 2: HWE过滤|HWE filtering
        hwe_prefix = self.config.output_path / f"{self.config.base_name}_hwe"
        hwe_prefix.parent.mkdir(parents=True, exist_ok=True)
        hwe_prefix_abs = hwe_prefix.resolve()
        cmd_hwe = build_conda_command(self.config.plink_path, [
            "--bfile", str(maf_prefix_abs),
            "--hwe", str(self.config.hwe_pvalue), "--make-bed",
            "--out", str(hwe_prefix_abs), "--allow-extra-chr",
        ])
        
        if not self.cmd_runner.run(cmd_hwe, f"HWE过滤|HWE filtering (p>{self.config.hwe_pvalue})"):
            return False
        
        # Step 3: 缺失率过滤|Missing rate filtering
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        output_prefix_abs = output_prefix.resolve()
        cmd_missing = build_conda_command(self.config.plink_path, [
            "--bfile", str(hwe_prefix_abs),
            "--geno", str(self.config.missing_rate),
            "--mind", str(self.config.missing_rate),
            "--make-bed", "--out", str(output_prefix_abs), "--allow-extra-chr",
        ])
        
        success = self.cmd_runner.run(cmd_missing, f"缺失率过滤|Missing rate filtering (<{self.config.missing_rate})")
        
        if success:
            self.config.plink_prefix_qc = str(output_prefix_abs)
            self.log_qc_stats()

        return success

    def _drop_zero_genotype_samples(self, input_prefix: str) -> str:
        """
        剔除零基因型样本(F_MISS==1.0)以避免PLINK --pca崩溃|Remove zero-genotype samples

        默认不过滤时调用:只剔除100%缺失的样本(结构性无用),不剔除任何变异或部分缺失样本。
        Called under no-filter default: only removes 100%-missing samples, no variants or partial-missing samples.

        Args:
            input_prefix: PLINK文件前缀|PLINK file prefix

        Returns:
            处理后的PLINK前缀(若无零数据样本则原样返回)|
            Processed PLINK prefix (unchanged if no zero-genotype samples)
        """
        misscan_prefix = self.config.output_path / f"{self.config.base_name}_misscan"
        misscan_prefix.parent.mkdir(parents=True, exist_ok=True)
        misscan_prefix_abs = misscan_prefix.resolve()

        # 仅扫描缺失率,不修改数据|Scan missing rates only, no data modification
        cmd_scan = build_conda_command(self.config.plink_path, [
            "--bfile", str(input_prefix),
            "--missing", "--out", str(misscan_prefix_abs), "--allow-extra-chr",
        ])
        if not self.cmd_runner.run(cmd_scan, "扫描样本缺失率|Scan sample missing rates"):
            # 扫描失败则原样返回,交给后续PCA步骤处理|On scan failure return as-is
            self.logger.warning("缺失率扫描失败,保持原数据|Missing-rate scan failed, keeping original data")
            return input_prefix

        imiss_file = Path(f"{misscan_prefix_abs}.imiss")
        if not imiss_file.exists():
            self.logger.warning("缺失率文件未生成,保持原数据|.imiss not generated, keeping original data")
            return input_prefix

        # 解析.imiss收集F_MISS==1.0的样本(零基因型)|Parse .imiss for F_MISS==1.0 (zero-genotype)
        zero_samples = []
        try:
            with open(imiss_file) as f:
                f.readline()  # 跳过表头|Skip header
                for line in f:
                    parts = line.split()
                    if len(parts) < 6:
                        continue
                    try:
                        f_miss = float(parts[5])
                    except ValueError:
                        continue
                    if f_miss >= 1.0:
                        zero_samples.append((parts[0], parts[1]))  # (FID, IID)
        except Exception as e:
            self.logger.warning(f"解析缺失率文件失败|Failed to parse .imiss: {e}")
            return input_prefix

        # 清理扫描产物(诊断文件)|Clean up diagnostic scan files
        for suffix in ('.imiss', '.lmiss', '.log'):
            p = Path(f"{misscan_prefix_abs}{suffix}")
            if p.exists():
                try:
                    p.unlink()
                except OSError:
                    pass

        # 无零数据样本:真正零过滤|No zero-genotype samples: truly zero filtering
        if not zero_samples:
            self.logger.info("无零基因型样本,数据原样通过|No zero-genotype samples, data passed through")
            return input_prefix

        # 有零数据样本:--remove剔除|Remove zero-genotype samples
        # PLINK 以 cwd=resolve(outdir) 执行,须给绝对路径,否则相对路径找不到→降级跳过剔除
        # |PLINK runs with cwd=resolve(outdir); must pass absolute path or it won't be found
        remove_file = (self.config.output_path / f"{self.config.base_name}_zero_samples.txt").resolve()
        try:
            with open(remove_file, 'w') as f:
                for fid, iid in zero_samples:
                    f.write(f"{fid}\t{iid}\n")
        except Exception as e:
            self.logger.error(f"写入剔除清单失败|Failed to write remove list: {e}")
            return input_prefix

        safe_prefix = self.config.output_path / f"{self.config.base_name}_safe"
        safe_prefix.parent.mkdir(parents=True, exist_ok=True)
        safe_prefix_abs = safe_prefix.resolve()
        cmd_remove = build_conda_command(self.config.plink_path, [
            "--bfile", str(input_prefix),
            "--remove", str(remove_file), "--make-bed",
            "--out", str(safe_prefix_abs), "--allow-extra-chr",
        ])
        if not self.cmd_runner.run(
            cmd_remove,
            f"剔除零基因型样本|Remove zero-genotype samples ({len(zero_samples)} samples)"
        ):
            self.logger.warning(
                "剔除零基因型样本失败,保持原数据|Failed to remove zero samples, keeping original data"
            )
            return input_prefix

        self.logger.info(
            f"检测到 {len(zero_samples)} 个零基因型样本(100%缺失),已自动剔除以避免PLINK崩溃|"
            f"Detected {len(zero_samples)} zero-genotype (100% missing) samples, auto-removed to avoid PLINK crash"
        )
        return str(safe_prefix_abs)

    def log_qc_stats(self):
        """记录质控统计|Log QC statistics"""
        try:
            fam_file = Path(f"{self.config.plink_prefix_qc}.fam")
            bim_file = Path(f"{self.config.plink_prefix_qc}.bim")
            
            if fam_file.exists():
                sample_count = sum(1 for line in open(fam_file))
                if self.config.skip_qc:
                    self.logger.info(f"输入样本数|Input samples: {sample_count}")
                else:
                    self.logger.info(f"质控后样本数|Samples after QC: {sample_count}")
            
            if bim_file.exists():
                snp_count = sum(1 for line in open(bim_file))
                if self.config.skip_qc:
                    self.logger.info(f"输入SNP数|Input SNPs: {snp_count}")
                else:
                    self.logger.info(f"质控后SNP数|SNPs after QC: {snp_count}")
                
        except Exception as e:
            self.logger.warning(f"无法读取统计信息|Cannot read statistics: {e}")
