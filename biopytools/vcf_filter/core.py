"""
VCF文件筛选核心模块 | VCF File Filtering Core Module
基于原始函数，集成所有筛选功能
"""

import os
import gzip
import subprocess
from pathlib import Path
from typing import Optional, List, Union
from .utils import PerformanceLogger, FastVCFStats, check_plink_availability

class VCFFilter:
    """VCF文件筛选类（基于原始函数逻辑）| VCF File Filter Class (Based on Original Function Logic)"""
    
    def __init__(self, vcf_file: str):
        """初始化VCF筛选器"""
        self.vcf_file = Path(vcf_file)
        if not self.vcf_file.exists():
            raise FileNotFoundError(f"VCF文件不存在: {vcf_file}")
    
    def filter_vcf(self, 
                   chr_name: Union[str, List[str]], 
                   start: Optional[int] = None, 
                   end: Optional[int] = None,
                   output_file: Optional[str] = None,
                   convert_format: bool = False,
                   plink_path: str = "plink",
                   allow_extra_chr: bool = True,
                   min_maf: Optional[float] = None,
                   max_missing: Optional[float] = None,
                   quality_threshold: Optional[float] = None,
                   min_depth: Optional[int] = None,
                   max_depth: Optional[int] = None,
                   keep_samples: Optional[List[str]] = None,
                   remove_samples: Optional[List[str]] = None,
                   biallelic_only: bool = False,
                   remove_indels: bool = False,
                   skip_validation: bool = True,
                   verbose: bool = False) -> str:
        """筛选VCF文件"""
        
        # 性能优化的日志器
        logger = PerformanceLogger(verbose)
        stats = FastVCFStats(logger)
        
        # 生成输出文件名
        if output_file is None:
            base_name = self.vcf_file.stem
            if isinstance(chr_name, list):
                chr_str = "_".join(chr_name)
            else:
                chr_str = str(chr_name)
            
            suffix = f"_{chr_str}"
            if start is not None and end is not None:
                suffix += f"_{start}_{end}"
            elif start is not None:
                suffix += f"_{start}_end"
            elif end is not None:
                suffix += f"start_{end}"
            
            output_file = f"{base_name}_filtered{suffix}.vcf"
        
        output_path = Path(output_file)
        
        # 性能优化：快速验证或跳过
        if not skip_validation:
            logger.info("快速验证输入文件...")
            variant_count = stats.get_variant_count(str(self.vcf_file), skip_count=False)
            if variant_count:
                logger.info(f"估算包含 {variant_count} 个变异位点")
        else:
            logger.info("跳过输入验证以提高性能")
        
        if convert_format:
            # 使用plink进行筛选和格式转换
            return self._filter_with_plink(
                chr_name, start, end, output_path, plink_path, 
                allow_extra_chr, min_maf, max_missing, quality_threshold, logger
            )
        else:
            # 使用Python进行筛选
            return self._filter_with_python(
                chr_name, start, end, output_path, quality_threshold,
                min_depth, max_depth, keep_samples, remove_samples,
                biallelic_only, remove_indels, logger
            )
    
    def _filter_with_plink(self, chr_name, start, end, output_path, plink_path, 
                          allow_extra_chr, min_maf, max_missing, quality_threshold, logger) -> str:
        """使用plink进行筛选"""
        
        logger.info("使用PLINK进行筛选和格式转换...")
        
        # 检查plink可用性
        if not check_plink_availability(plink_path, logger):
            logger.warning("PLINK不可用，切换到Python筛选")
            return self._filter_with_python(chr_name, start, end, output_path, 
                                          quality_threshold, None, None, None, None, 
                                          False, False, logger)
        
        # 构建plink命令
        cmd = [
            plink_path,
            "--vcf", str(self.vcf_file),
            "--recode", "vcf-iid",
            "--out", str(output_path.with_suffix(''))
        ]
        
        if allow_extra_chr:
            cmd.append("--allow-extra-chr")
        
        # 添加染色体筛选
        if isinstance(chr_name, list):
            chr_list = ",".join(str(c) for c in chr_name)
            cmd.extend(["--chr", chr_list])
        else:
            cmd.extend(["--chr", str(chr_name)])
        
        # 添加位置筛选
        if start is not None and end is not None:
            cmd.extend(["--from-bp", str(start), "--to-bp", str(end)])
        elif start is not None:
            cmd.extend(["--from-bp", str(start)])
        elif end is not None:
            cmd.extend(["--to-bp", str(end)])
        
        # 添加质量控制参数
        if min_maf is not None:
            cmd.extend(["--maf", str(min_maf)])
        
        if max_missing is not None:
            cmd.extend(["--geno", str(max_missing)])
        
        try:
            # 执行plink命令
            logger.info(f"执行PLINK命令: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("PLINK执行成功")
            
            # plink输出文件会自动添加.vcf后缀
            final_output = str(output_path.with_suffix('')) + ".vcf"
            
            # 检查输出文件是否存在且有内容
            if not os.path.exists(final_output):
                raise FileNotFoundError(f"PLINK输出文件未生成: {final_output}")
            
            # 检查文件是否有变异位点（除了头部）
            variant_count = self._count_variants_in_file(final_output)
            if variant_count == 0:
                logger.warning(f"PLINK筛选后没有变异位点，输出空VCF文件")
            else:
                logger.info(f"PLINK筛选完成，输出 {variant_count} 个变异位点")
            
            return final_output
            
        except subprocess.CalledProcessError as e:
            logger.error(f"PLINK执行失败，退出码: {e.returncode}")
            if e.stderr:
                logger.error(f"PLINK错误输出: {e.stderr}")
            if e.stdout:
                logger.info(f"PLINK标准输出: {e.stdout}")
            
            # 特殊处理"All variants excluded"错误
            if "All variants excluded" in str(e.stderr):
                logger.warning("PLINK报告所有变异都被排除，可能原因:")
                logger.warning(f"1. 染色体 '{chr_name}' 在VCF文件中不存在")
                logger.warning(f"2. 位置范围 {start}-{end} 内没有变异位点")
                logger.warning("3. 质量控制参数过于严格")
                logger.warning("尝试切换到Python筛选模式...")
                
                # 切换到Python筛选
                return self._filter_with_python(chr_name, start, end, output_path, 
                                              quality_threshold, None, None, None, None, 
                                              False, False, logger)
            else:
                raise
        except FileNotFoundError:
            raise FileNotFoundError(f"找不到plink可执行文件: {plink_path}")
    
    def _count_variants_in_file(self, vcf_file: str) -> int:
        """统计VCF文件中的变异位点数量"""
        count = 0
        try:
            opener = gzip.open if vcf_file.endswith('.gz') else open
            with opener(vcf_file, 'rt') as f:
                for line in f:
                    if not line.startswith('#') and line.strip():
                        count += 1
        except:
            pass
        return count
    
    def _filter_with_python(self, chr_name, start, end, output_path, 
                           quality_threshold, min_depth, max_depth, 
                           keep_samples, remove_samples, biallelic_only, 
                           remove_indels, logger) -> str:
        """使用Python进行VCF筛选"""
        
        logger.info("使用Python进行VCF筛选...")
        
        if isinstance(chr_name, str):
            chr_list = [chr_name]
        else:
            chr_list = chr_name
        
        filtered_count = 0
        total_count = 0
        
        opener = gzip.open if str(self.vcf_file).endswith('.gz') else open
        
        # 样本筛选相关变量
        sample_indices = None
        header_written = False
        
        # 首先检查是否能找到目标染色体的变异
        logger.info(f"检查染色体 {chr_list} 在位置范围 {start}-{end} 内的变异...")
        
        with opener(self.vcf_file, 'rt') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                # 处理注释行
                if line.startswith('##'):
                    outfile.write(line)
                    continue
                elif line.startswith('#CHROM'):
                    # 处理样本头部
                    if keep_samples or remove_samples:
                        fields = line.strip().split('\\t')
                        if len(fields) > 9:
                            samples = fields[9:]
                            sample_indices = self._get_sample_indices(samples, keep_samples, remove_samples)
                            # 重写头部
                            new_header = fields[:9] + [samples[i] for i in sample_indices]
                            outfile.write('\\t'.join(new_header) + '\\n')
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)
                    header_written = True
                    continue
                
                if not header_written:
                    continue
                
                total_count += 1
                fields = line.strip().split('\\t')
                
                if len(fields) < 8:
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                
                # 应用筛选条件
                if not self._should_keep_variant(chrom, pos, ref, alt, qual, chr_list, 
                                               start, end, quality_threshold, 
                                               biallelic_only, remove_indels):
                    continue
                
                # 处理样本筛选
                if sample_indices is not None and len(fields) > 9:
                    sample_fields = fields[9:]
                    filtered_samples = [sample_fields[i] for i in sample_indices if i < len(sample_fields)]
                    line = '\\t'.join(fields[:9] + filtered_samples) + '\\n'
                
                outfile.write(line)
                filtered_count += 1
                
                # 进度报告（降低频率以提高性能）
                if total_count % 100000 == 0:
                    logger.info(f"已处理 {total_count} 个变异位点，保留 {filtered_count} 个")
        
        # 结果报告
        if filtered_count == 0:
            logger.warning(f"警告: 筛选后没有找到任何变异位点!")
            logger.warning(f"可能的原因:")
            logger.warning(f"1. 染色体名称 '{chr_name}' 在VCF文件中不存在")
            logger.warning(f"2. 位置范围 {start}-{end} 内没有变异位点")
            logger.warning(f"3. 筛选条件过于严格")
            
            # 提供调试信息
            logger.info("\\n建议检查:")
            logger.info("1. 查看VCF文件中的染色体名称:")
            logger.info(f"   zcat {self.vcf_file} | grep -v '^##' | head -10")
            logger.info("2. 检查指定区域是否有变异:")
            # 使用变量来避免f-string中的转义问题
            chr_check_cmd = f"   zcat {self.vcf_file} | grep -v '^#' | awk '$1==\"{chr_name}\" && $2>={start} && $2<={end}' | head -5"
            logger.info(chr_check_cmd)
        else:
            logger.info(f"Python筛选完成: {total_count} -> {filtered_count} 个变异位点")
        
        return str(output_path)
    
    def _get_sample_indices(self, samples, keep_samples, remove_samples):
        """获取要保留的样本索引"""
        indices = []
        for i, sample in enumerate(samples):
            keep = True
            
            if keep_samples is not None:
                keep = sample in keep_samples
            
            if remove_samples is not None:
                if sample in remove_samples:
                    keep = False
            
            if keep:
                indices.append(i)
        
        return indices
    
    def _should_keep_variant(self, chrom, pos, ref, alt, qual, chr_list, 
                           start, end, quality_threshold, biallelic_only, remove_indels):
        """判断是否保留变异位点"""
        
        # 染色体筛选
        if chrom not in chr_list:
            return False
        
        # 位置筛选
        if start is not None and pos < start:
            return False
        if end is not None and pos > end:
            return False
        
        # 质量筛选
        if quality_threshold is not None:
            if qual != '.' and qual != 'PASS':
                try:
                    qual_value = float(qual)
                    if qual_value < quality_threshold:
                        return False
                except ValueError:
                    pass
        
        # 双等位基因位点筛选
        if biallelic_only:
            alt_alleles = alt.split(',')
            if len(alt_alleles) > 1:
                return False
        
        # 移除插入缺失
        if remove_indels:
            if len(ref) > 1 or any(len(a) > 1 for a in alt.split(',')):
                return False
        
        return True
