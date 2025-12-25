"""
VCF联合分型处理器模块 | VCF Joint Genotyping Processor Module
"""

import subprocess
from pathlib import Path
from typing import Dict, Any, List
from .utils import CommandRunner, VCFUtils

class VCFJointGenotypingProcessor:
    """VCF联合分型处理器 | VCF Joint Genotyping Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_paths = config.get_output_paths()
    
    def create_sample_map(self) -> bool:
        """📋 创建样本映射文件 | Create sample mapping file"""
        self.logger.info("📋 创建样本映射文件 | Creating sample mapping file")
        
        vcf_files = VCFUtils.get_vcf_files(self.config.vcf_input_dir)
        
        if not vcf_files:
            self.logger.error(f"❌ 在目录中没有找到*.vcf.gz文件 | No *.vcf.gz files found in directory: {self.config.vcf_input_dir}")
            return False
        
        try:
            with open(self.output_paths['sample_map_file'], 'w') as f:
                for vcf_file in vcf_files:
                    sample_name = vcf_file.stem.replace('.vcf', '')
                    f.write(f"{sample_name}\t{vcf_file.absolute()}\n")
            
            self.logger.info(f"✅ 样本映射文件已创建 | Sample mapping file created: {self.output_paths['sample_map_file']} ({len(vcf_files)} 个样本 | samples)")
            
            # 显示前几行内容 | Show first few lines
            with open(self.output_paths['sample_map_file'], 'r') as f:
                lines = f.readlines()[:5]
                self.logger.info("📋 样本映射文件内容预览 | Sample mapping file preview:")
                for line in lines:
                    self.logger.info(f"  📄 {line.strip()}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 创建样本映射文件失败 | Failed to create sample mapping file: {e}")
            return False
    
    def check_and_index_vcf_files(self) -> bool:
        """🔍 检查并创建VCF索引文件 | Check and create VCF index files"""
        self.logger.info("🔍 检查VCF文件索引 | Checking VCF file indexes")
        
        vcf_files = VCFUtils.get_vcf_files(self.config.vcf_input_dir)
        
        missing_indexes = []
        corrupted_files = []
        
        for vcf_file in vcf_files:
            index_file = Path(f"{vcf_file}.tbi")
            if not index_file.exists():
                # 先检查文件完整性
                if self._check_vcf_integrity(vcf_file):
                    missing_indexes.append(vcf_file)
                else:
                    corrupted_files.append(vcf_file)
        
        # 处理损坏的文件
        if corrupted_files:
            self.logger.warning(f"⚠️ 发现 {len(corrupted_files)} 个损坏的VCF文件 | Found {len(corrupted_files)} corrupted VCF files")
            for corrupted_file in corrupted_files:
                self.logger.warning(f"🚨 损坏文件 | Corrupted file: {corrupted_file.name}")
                
                # 尝试修复文件
                if self._try_repair_vcf(corrupted_file):
                    self.logger.info(f"🔧 文件修复成功 | File repaired successfully: {corrupted_file.name}")
                    missing_indexes.append(corrupted_file)
                else:
                    if self.config.skip_corrupted_files:
                        self.logger.warning(f"⚠️ 文件修复失败，将跳过此文件 | File repair failed, will skip this file: {corrupted_file.name}")
                        # 从样本映射中移除损坏的文件
                        self._remove_corrupted_file_from_map(corrupted_file)
                    else:
                        self.logger.error(f"❌ 文件修复失败，程序终止 | File repair failed, program terminated: {corrupted_file.name}")
                        return False
        
        if missing_indexes:
            self.logger.info(f"🔧 发现 {len(missing_indexes)} 个VCF文件缺少索引，开始创建索引 | Found {len(missing_indexes)} VCF files missing indexes, creating indexes")
            
            for vcf_file in missing_indexes:
                self.logger.info(f"📂 为文件创建索引 | Creating index for file: {vcf_file.name}")
                cmd = f"{self.config.tabix_path} -p vcf {vcf_file}"
                
                if not self.cmd_runner.run(cmd, f"创建VCF索引 | Create VCF index: {vcf_file.name}"):
                    self.logger.error(f"❌ 创建索引失败 | Failed to create index for: {vcf_file}")
                    
                    # 尝试其他方法创建索引
                    if self._try_alternative_indexing(vcf_file):
                        self.logger.info(f"✅ 使用替代方法创建索引成功 | Alternative indexing method succeeded: {vcf_file.name}")
                    else:
                        if self.config.skip_corrupted_files:
                            self.logger.warning(f"⚠️ 所有索引方法都失败，将跳过此文件 | All indexing methods failed, will skip this file: {vcf_file.name}")
                            self._remove_corrupted_file_from_map(vcf_file)
                            continue
                        else:
                            self.logger.error(f"❌ 所有索引方法都失败，程序终止 | All indexing methods failed, program terminated: {vcf_file.name}")
                            return False
                
                # 验证索引文件是否创建成功
                index_file = Path(f"{vcf_file}.tbi")
                if index_file.exists():
                    self.logger.info(f"✅ 索引创建成功 | Index created successfully: {index_file.name}")
                else:
                    if self.config.skip_corrupted_files:
                        self.logger.warning(f"⚠️ 索引文件未找到，将跳过此文件 | Index file not found, will skip this file: {index_file}")
                        self._remove_corrupted_file_from_map(vcf_file)
                        continue
                    else:
                        self.logger.error(f"❌ 索引文件未找到，程序终止 | Index file not found, program terminated: {index_file}")
                        return False
        else:
            self.logger.info("✅ 所有VCF文件都已有索引 | All VCF files already have indexes")
        
        # 最终验证：确保至少有一些有效的VCF文件
        final_vcf_count = len(VCFUtils.get_vcf_files(self.config.vcf_input_dir))
        if final_vcf_count == 0:
            self.logger.error("❌ 没有有效的VCF文件可以处理 | No valid VCF files available for processing")
            return False
        
        return True
    
    def _check_vcf_integrity(self, vcf_file: Path) -> bool:
        """🔍 检查VCF文件完整性 | Check VCF file integrity"""
        try:
            # 尝试读取文件头部
            cmd = f"zcat {vcf_file} | head -n 100 > /dev/null 2>&1"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                self.logger.debug(f"🔍 文件头部检查失败 | File header check failed: {vcf_file.name}")
                return False
            
            # 检查bgzip格式
            cmd = f"bgzip -t {vcf_file}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                self.logger.debug(f"🔍 bgzip格式检查失败 | bgzip format check failed: {vcf_file.name}")
                return False
            
            return True
            
        except Exception as e:
            self.logger.debug(f"🔍 文件完整性检查异常 | File integrity check exception: {vcf_file.name} - {e}")
            return False
    
    def _try_repair_vcf(self, vcf_file: Path) -> bool:
        """🔧 尝试修复VCF文件 | Try to repair VCF file"""
        self.logger.info(f"🔧 尝试修复VCF文件 | Trying to repair VCF file: {vcf_file.name}")
        
        backup_file = vcf_file.with_suffix(vcf_file.suffix + '.backup')
        repaired_file = vcf_file.with_suffix('.repaired.vcf.gz')
        
        try:
            # 备份原文件
            cmd = f"cp {vcf_file} {backup_file}"
            if not subprocess.run(cmd, shell=True, capture_output=True).returncode == 0:
                return False
            
            # 方法1: 重新压缩
            self.logger.info(f"🔧 方法1: 重新压缩文件 | Method 1: Re-compressing file")
            cmd = f"zcat {vcf_file} | bgzip -c > {repaired_file}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and repaired_file.exists():
                # 验证修复后的文件
                if self._check_vcf_integrity(repaired_file):
                    # 替换原文件
                    cmd = f"mv {repaired_file} {vcf_file}"
                    subprocess.run(cmd, shell=True, capture_output=True)
                    self.logger.info(f"✅ 重新压缩成功 | Re-compression successful")
                    return True
            
            # 方法2: 使用gunzip重新压缩
            if repaired_file.exists():
                repaired_file.unlink()
            
            self.logger.info(f"🔧 方法2: 使用gunzip重新压缩 | Method 2: Re-compressing with gunzip")
            temp_vcf = vcf_file.with_suffix('.temp.vcf')
            cmd = f"gunzip -c {vcf_file} > {temp_vcf} && bgzip -c {temp_vcf} > {repaired_file}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and repaired_file.exists():
                if self._check_vcf_integrity(repaired_file):
                    cmd = f"mv {repaired_file} {vcf_file}"
                    subprocess.run(cmd, shell=True, capture_output=True)
                    if temp_vcf.exists():
                        temp_vcf.unlink()
                    self.logger.info(f"✅ gunzip重新压缩成功 | gunzip re-compression successful")
                    return True
            
            # 清理临时文件
            if temp_vcf.exists():
                temp_vcf.unlink()
            if repaired_file.exists():
                repaired_file.unlink()
            
            return False
            
        except Exception as e:
            self.logger.error(f"❌ 文件修复过程中出现异常 | Exception during file repair: {e}")
            return False
    
    def _try_alternative_indexing(self, vcf_file: Path) -> bool:
        """🔧 尝试其他索引方法 | Try alternative indexing methods"""
        self.logger.info(f"🔧 尝试其他索引方法 | Trying alternative indexing methods: {vcf_file.name}")
        
        try:
            # 方法1: 使用bcftools创建索引
            self.logger.info(f"🔧 尝试使用bcftools创建索引 | Trying bcftools indexing")
            cmd = f"{self.config.bcftools_path} index -t {vcf_file}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                index_file = Path(f"{vcf_file}.tbi")
                if index_file.exists():
                    self.logger.info(f"✅ bcftools索引创建成功 | bcftools indexing successful")
                    return True
            
            # 方法2: 强制重建索引
            self.logger.info(f"🔧 尝试强制重建索引 | Trying force rebuild index")
            cmd = f"{self.config.tabix_path} -f -p vcf {vcf_file}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                index_file = Path(f"{vcf_file}.tbi")
                if index_file.exists():
                    self.logger.info(f"✅ 强制重建索引成功 | Force rebuild index successful")
                    return True
            
            return False
            
        except Exception as e:
            self.logger.error(f"❌ 其他索引方法异常 | Alternative indexing exception: {e}")
            return False
    
    def _remove_corrupted_file_from_map(self, corrupted_file: Path):
        """🗑️ 移除损坏的文件 | Remove corrupted file"""
        try:
            # 如果样本映射文件已存在，从中移除
            if self.output_paths['sample_map_file'].exists():
                # 读取现有映射文件
                with open(self.output_paths['sample_map_file'], 'r') as f:
                    lines = f.readlines()
                
                # 过滤掉损坏的文件
                filtered_lines = []
                for line in lines:
                    if corrupted_file.name not in line:
                        filtered_lines.append(line)
                
                # 重写映射文件
                with open(self.output_paths['sample_map_file'], 'w') as f:
                    f.writelines(filtered_lines)
                
                self.logger.info(f"🗑️ 已从样本映射中移除损坏文件 | Removed corrupted file from sample map: {corrupted_file.name}")
            
            # 可选：创建一个损坏文件的备份目录并移动文件
            corrupted_dir = Path(self.config.output_dir) / "corrupted_files"
            corrupted_dir.mkdir(exist_ok=True)
            
            backup_path = corrupted_dir / corrupted_file.name
            if not backup_path.exists():
                import shutil
                shutil.copy2(corrupted_file, backup_path)
                self.logger.info(f"📁 损坏文件已备份到 | Corrupted file backed up to: {backup_path}")
        
        except Exception as e:
            self.logger.error(f"❌ 处理损坏文件失败 | Failed to handle corrupted file: {e}")
    
    def joint_genotyping(self) -> bool:
        """🔗 GTX联合分型 | GTX joint genotyping"""
        if self.config.skip_joint:
            self.logger.info("⏭️ 跳过联合分型步骤 | Skipping joint genotyping step")
            return True
        
        self.logger.info("🔗 使用GTX进行联合分型 | Performing joint genotyping with GTX")
        
        if self.output_paths['merged_vcf'].exists():
            self.logger.info("⏭️ 合并VCF文件已存在，跳过联合分型 | Merged VCF file exists, skipping joint genotyping")
            file_size = VCFUtils.get_file_size(self.output_paths['merged_vcf'])
            self.logger.info(f"📏 已存在合并VCF大小 | Existing merged VCF size: {file_size}")
            return True
        
        cmd = (f"{self.config.gtx_path} joint "
               f"-r {self.config.reference_genome} "
               f"--sample-name-map {self.output_paths['sample_map_file']} "
               f"-o {self.output_paths['merged_vcf']}")
        
        success = self.cmd_runner.run(cmd, "🔗 GTX联合分型 | GTX joint genotyping")
        
        if success and self.output_paths['merged_vcf'].exists():
            file_size = VCFUtils.get_file_size(self.output_paths['merged_vcf'])
            self.logger.info(f"✅ 联合分型完成，输出文件大小 | Joint genotyping completed, output file size: {file_size}")
        
        return success
    
    def extract_variants(self) -> bool:
        """🧬 提取SNP和INDEL变异 | Extract SNP and INDEL variants"""
        if self.config.skip_extract:
            self.logger.info("⏭️ 跳过变异提取步骤 | Skipping variant extraction step")
            return True
        
        # 提取SNP | Extract SNPs
        self.logger.info("🧬 提取SNP变异 | Extracting SNP variants")
        if not self.output_paths['snp_vcf'].exists():
            cmd = (f"{self.config.gatk_path} --java-options \"-Xmx{self.config.memory}\" SelectVariants "
                   f"-V {self.output_paths['merged_vcf']} "
                   f"-select-type SNP "
                   f"-O {self.output_paths['snp_vcf']}")
            
            if not self.cmd_runner.run(cmd, "🧬 提取SNP变异 | Extract SNP variants"):
                return False
            
            if self.output_paths['snp_vcf'].exists():
                file_size = VCFUtils.get_file_size(self.output_paths['snp_vcf'])
                self.logger.info(f"✅ SNP提取完成，文件大小 | SNP extraction completed, file size: {file_size}")
        else:
            self.logger.info("⏭️ SNP VCF文件已存在，跳过SNP提取 | SNP VCF file exists, skipping SNP extraction")
            file_size = VCFUtils.get_file_size(self.output_paths['snp_vcf'])
            self.logger.info(f"📏 已存在SNP VCF大小 | Existing SNP VCF size: {file_size}")
        
        # 提取INDEL | Extract INDELs
        self.logger.info("🧬 提取INDEL变异 | Extracting INDEL variants")
        if not self.output_paths['indel_vcf'].exists():
            cmd = (f"{self.config.gatk_path} --java-options \"-Xmx{self.config.memory}\" SelectVariants "
                   f"-V {self.output_paths['merged_vcf']} "
                   f"-select-type INDEL "
                   f"-O {self.output_paths['indel_vcf']}")
            
            if not self.cmd_runner.run(cmd, "🧬 提取INDEL变异 | Extract INDEL variants"):
                return False
            
            if self.output_paths['indel_vcf'].exists():
                file_size = VCFUtils.get_file_size(self.output_paths['indel_vcf'])
                self.logger.info(f"✅ INDEL提取完成，文件大小 | INDEL extraction completed, file size: {file_size}")
        else:
            self.logger.info("⏭️ INDEL VCF文件已存在，跳过INDEL提取 | INDEL VCF file exists, skipping INDEL extraction")
            file_size = VCFUtils.get_file_size(self.output_paths['indel_vcf'])
            self.logger.info(f"📏 已存在INDEL VCF大小 | Existing INDEL VCF size: {file_size}")
        
        return True
    
    def filter_variants(self) -> bool:
        """🔍 过滤变异 | Filter variants"""
        if self.config.skip_filter:
            self.logger.info("⏭️ 跳过变异过滤步骤 | Skipping variant filtering step")
            return True
        
        # 过滤SNP | Filter SNPs
        self.logger.info("🔍 过滤SNP变异 | Filtering SNP variants")
        self.logger.info(f"🔧 SNP过滤参数 | SNP filtering parameters: "
                        f"--maf {self.config.snp_maf} "
                        f"--max-missing {self.config.snp_max_missing} "
                        f"--hwe {self.config.snp_hwe_pvalue} "
                        f"--min-meanDP {self.config.snp_min_mean_dp} "
                        f"--max-meanDP {self.config.snp_max_mean_dp}")
        
        if not self.output_paths['filtered_snp_vcf'].exists():
            cmd = (f"{self.config.vcftools_path} "
                   f"--gzvcf {self.output_paths['snp_vcf']} "
                   f"--maf {self.config.snp_maf} "
                   f"--max-missing {self.config.snp_max_missing} "
                   f"--hwe {self.config.snp_hwe_pvalue} "
                   f"--min-meanDP {self.config.snp_min_mean_dp} "
                   f"--max-meanDP {self.config.snp_max_mean_dp} "
                   f"--recode --recode-INFO-all "
                   f"--out {Path(self.config.output_dir) / self.config.filtered_snp_prefix}")
            
            if not self.cmd_runner.run(cmd, "🔍 过滤SNP变异 | Filter SNP variants"):
                return False
            
            if self.output_paths['filtered_snp_vcf'].exists():
                file_size = VCFUtils.get_file_size(self.output_paths['filtered_snp_vcf'])
                self.logger.info(f"✅ SNP过滤完成，文件大小 | SNP filtering completed, file size: {file_size}")
        else:
            self.logger.info("⏭️ 过滤后SNP VCF文件已存在，跳过SNP过滤 | Filtered SNP VCF file exists, skipping SNP filtering")
            file_size = VCFUtils.get_file_size(self.output_paths['filtered_snp_vcf'])
            self.logger.info(f"📏 已存在过滤后SNP VCF大小 | Existing filtered SNP VCF size: {file_size}")
        
        # 过滤INDEL | Filter INDELs  
        self.logger.info("🔍 过滤INDEL变异 | Filtering INDEL variants")
        self.logger.info(f"🔧 INDEL过滤参数 | INDEL filtering parameters: "
                        f"--maf {self.config.indel_maf} "
                        f"--max-missing {self.config.indel_max_missing} "
                        f"--hwe {self.config.indel_hwe_pvalue} "
                        f"--min-meanDP {self.config.indel_min_mean_dp} "
                        f"--max-meanDP {self.config.indel_max_mean_dp}")
        
        if not self.output_paths['filtered_indel_vcf'].exists():
            cmd = (f"{self.config.vcftools_path} "
                   f"--gzvcf {self.output_paths['indel_vcf']} "
                   f"--maf {self.config.indel_maf} "
                   f"--max-missing {self.config.indel_max_missing} "
                   f"--hwe {self.config.indel_hwe_pvalue} "
                   f"--min-meanDP {self.config.indel_min_mean_dp} "
                   f"--max-meanDP {self.config.indel_max_mean_dp} "
                   f"--recode --recode-INFO-all "
                   f"--out {Path(self.config.output_dir) / self.config.filtered_indel_prefix}")
            
            if not self.cmd_runner.run(cmd, "🔍 过滤INDEL变异 | Filter INDEL variants"):
                return False
            
            if self.output_paths['filtered_indel_vcf'].exists():
                file_size = VCFUtils.get_file_size(self.output_paths['filtered_indel_vcf'])
                self.logger.info(f"✅ INDEL过滤完成，文件大小 | INDEL filtering completed, file size: {file_size}")
        else:
            self.logger.info("⏭️ 过滤后INDEL VCF文件已存在，跳过INDEL过滤 | Filtered INDEL VCF file exists, skipping INDEL filtering")
            file_size = VCFUtils.get_file_size(self.output_paths['filtered_indel_vcf'])
            self.logger.info(f"📏 已存在过滤后INDEL VCF大小 | Existing filtered INDEL VCF size: {file_size}")
        
        return True
    
    def compress_and_index(self) -> bool:
        """📦 压缩和索引结果文件 | Compress and index result files"""
        if self.config.skip_compress:
            self.logger.info("⏭️ 跳过压缩和索引步骤 | Skipping compression and indexing step")
            return True
        
        # 压缩和索引SNP | Compress and index SNPs
        self.logger.info("📦 压缩和索引SNP结果 | Compressing and indexing SNP results")
        if self.output_paths['filtered_snp_vcf'].exists() and not self.output_paths['final_snp_vcf'].exists():
            cmd = f"{self.config.bgzip_path} {self.output_paths['filtered_snp_vcf']}"
            if not self.cmd_runner.run(cmd, "📦 压缩SNP VCF文件 | Compress SNP VCF file"):
                return False
        
        if self.output_paths['final_snp_vcf'].exists() and not Path(f"{self.output_paths['final_snp_vcf']}.tbi").exists():
            cmd = f"{self.config.tabix_path} -p vcf {self.output_paths['final_snp_vcf']}"
            if not self.cmd_runner.run(cmd, "🔗 创建SNP VCF索引 | Create SNP VCF index"):
                return False
        
        # 压缩和索引INDEL | Compress and index INDELs
        self.logger.info("📦 压缩和索引INDEL结果 | Compressing and indexing INDEL results")
        if self.output_paths['filtered_indel_vcf'].exists() and not self.output_paths['final_indel_vcf'].exists():
            cmd = f"{self.config.bgzip_path} {self.output_paths['filtered_indel_vcf']}"
            if not self.cmd_runner.run(cmd, "📦 压缩INDEL VCF文件 | Compress INDEL VCF file"):
                return False
        
        if self.output_paths['final_indel_vcf'].exists() and not Path(f"{self.output_paths['final_indel_vcf']}.tbi").exists():
            cmd = f"{self.config.tabix_path} -p vcf {self.output_paths['final_indel_vcf']}"
            if not self.cmd_runner.run(cmd, "🔗 创建INDEL VCF索引 | Create INDEL VCF index"):
                return False
        
        return True
    
    def generate_statistics(self) -> Dict[str, Any]:
        """📊 生成统计信息 | Generate statistics"""
        self.logger.info("📊 生成统计信息 | Generating statistics")
        
        stats = {}
        
        # 获取输入文件统计 | Get input file statistics
        vcf_count = VCFUtils.count_vcf_files(self.config.vcf_input_dir)
        stats['input_vcf_count'] = vcf_count
        
        # 统计变异数量 | Count variants
        for vcf_path, var_type in [
            (self.output_paths['final_snp_vcf'], "SNP"), 
            (self.output_paths['final_indel_vcf'], "INDEL")
        ]:
            if vcf_path.exists():
                count = VCFUtils.count_variants_in_vcf(vcf_path)
                stats[f"{var_type.lower()}_count"] = count
                self.logger.info(f"🔢 过滤后{var_type}总数 | Filtered {var_type} count: {count}")
                
                # 计算过滤效果 | Calculate filtering efficiency
                original_vcf = self.output_paths[f"{var_type.lower()}_vcf"]
                if original_vcf.exists():
                    original_count = VCFUtils.count_variants_in_vcf(original_vcf)
                    if original_count > 0:
                        filter_rate = (count / original_count) * 100
                        stats[f"{var_type.lower()}_filter_rate"] = filter_rate
                        self.logger.info(f"📊 {var_type}过滤保留率 | {var_type} filter retention rate: {filter_rate:.2f}% ({count}/{original_count})")
            else:
                stats[f"{var_type.lower()}_count"] = 0
        
        # 计算总变异数 | Calculate total variants
        total_variants = stats.get('snp_count', 0) + stats.get('indel_count', 0)
        stats['total_variants'] = total_variants
        
        if total_variants > 0:
            self.logger.info(f"🔢 过滤后变异总数 | Total filtered variants: {total_variants}")
        
        return stats
