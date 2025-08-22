"""
🧬 覆盖度数据处理模块 | Depth Data Processing Module
"""

import os
import subprocess
from pathlib import Path
from typing import List, Dict
from .utils import CommandRunner, get_sample_name, format_region

class DepthProcessor:
    """🧬 覆盖度处理器 | Depth Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.temp_files = []
    
    def process_all_files(self) -> bool:
        """📊 处理所有输入文件 | Process all input files"""
        self.logger.info(f"🚀 开始处理 {len(self.config.input_files)} 个文件 | Starting to process {len(self.config.input_files)} files")
        
        # 创建输出文件并写入表头 | Create output file and write header
        with open(self.config.output_file, 'w') as f:
            f.write("Sample\tChromosome\tPosition\tDepth\n")
        
        # 如果启用窗口分析，创建窗口结果文件 | If window analysis is enabled, create window results file
        if self.config.enable_window_analysis:
            window_output_file = self._get_window_output_file()
            with open(window_output_file, 'w') as f:
                f.write("Sample\tChromosome\tWindow_Start\tWindow_End\tWindow_Center\tAverage_Depth\tData_Points\n")
        
        total_processed = 0
        
        for input_file in self.config.input_files:
            sample_name = get_sample_name(input_file)
            self.logger.info(f"📄 处理文件: {input_file} (样品: {sample_name}) | Processing file: {input_file} (sample: {sample_name})")
            
            if self._process_single_file(input_file, sample_name):
                total_processed += 1
            else:
                self.logger.error(f"❌ 文件处理失败: {input_file} | Failed to process file: {input_file}")
        
        # 清理临时文件 | Clean up temporary files
        self._cleanup_temp_files()
        
        self.logger.info(f"✅ 成功处理 {total_processed}/{len(self.config.input_files)} 个文件 | Successfully processed {total_processed}/{len(self.config.input_files)} files")
        
        return total_processed > 0
    
    def _process_single_file(self, input_file: str, sample_name: str) -> bool:
        """📄 处理单个文件 | Process single file"""
        try:
            # 检查并创建索引（如果需要区间筛选） | Check and create index (if region filtering is needed)
            region = format_region(self.config.chromosome, self.config.start_pos, self.config.end_pos)
            if region:
                if not self._check_and_create_index(input_file):
                    # 如果索引创建失败，移除区间参数继续分析 | If index creation fails, remove region parameter and continue
                    self.logger.warning(f"⚠️ 索引创建失败，将分析整个BAM文件（可能较慢） | Index creation failed, will analyze entire BAM file (may be slow)")
                    region = ""
            
            # 构建samtools depth命令 | Build samtools depth command
            cmd_parts = [
                self.config.samtools_path, "depth",
                f"-@ {self.config.threads}",
                f"-q {self.config.quality_threshold}",
                f"-Q {self.config.mapping_quality}",
                "-a"  # 输出所有位置包括0覆盖度 | Output all positions including zero coverage
            ]
            
            # 添加区间参数 | Add region parameter
            if region:
                cmd_parts.append(f"-r {region}")
            
            cmd_parts.append(input_file)
            
            cmd = " ".join(cmd_parts)
            
            # 创建临时文件 | Create temporary file
            temp_file = Path(self.config.output_file).parent / f"temp_{sample_name}.depth"
            self.temp_files.append(temp_file)
            
            # 执行samtools depth命令 | Execute samtools depth command
            self.logger.info(f"🔬 分析覆盖度: {sample_name} | Analyzing depth: {sample_name}")
            
            with open(temp_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    shell=True,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    cwd=self.cmd_runner.working_dir
                )
            
            if result.returncode != 0:
                self.logger.error(f"❌ samtools执行失败: {result.stderr}")
                return False
            
            # 处理输出并过滤 | Process output and filter
            return self._filter_and_append_data(temp_file, sample_name)
            
        except Exception as e:
            self.logger.error(f"❌ 处理文件时出错: {e}")
            return False
    
    def _filter_and_append_data(self, temp_file: Path, sample_name: str) -> bool:
        """🔍 过滤数据并追加到输出文件 | Filter data and append to output file"""
        try:
            filtered_count = 0
            window_data = {}  # 用于存储窗口数据 | Store window data
            
            with open(temp_file, 'r') as f_in, open(self.config.output_file, 'a') as f_out:
                for line in f_in:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        chrom, pos, depth = parts[0], parts[1], parts[2]
                        
                        # 位置过滤 | Position filtering
                        try:
                            pos_int = int(pos)
                            # 检查是否需要位置过滤 | Check if position filtering is needed
                            if self.config.start_pos is not None and self.config.end_pos is not None:
                                if not (self.config.start_pos <= pos_int <= self.config.end_pos):
                                    continue
                        except ValueError:
                            continue
                        
                        # 染色体过滤 | Chromosome filtering
                        if self.config.chromosome != 'all' and chrom != self.config.chromosome:
                            continue
                        
                        # 深度值验证 | Depth value validation
                        try:
                            depth_float = float(depth)
                        except ValueError:
                            continue
                        
                        # 写入输出文件 | Write to output file
                        f_out.write(f"{sample_name}\t{chrom}\t{pos}\t{depth}\n")
                        filtered_count += 1
                        
                        # 如果启用窗口分析，收集窗口数据 | If window analysis is enabled, collect window data
                        if self.config.enable_window_analysis:
                            self._collect_window_data(window_data, sample_name, chrom, pos_int, depth_float)
            
            # 如果启用窗口分析，处理窗口数据 | If window analysis is enabled, process window data
            if self.config.enable_window_analysis and window_data:
                self._write_window_data(window_data)
            
            self.logger.info(f"📊 {sample_name}: 过滤出 {filtered_count} 个位点 | {sample_name}: filtered {filtered_count} positions")
            
            if self.config.enable_window_analysis:
                window_count = sum(len(chroms) for chroms in window_data.values())
                self.logger.info(f"📊 {sample_name}: 生成 {window_count} 个窗口 | {sample_name}: generated {window_count} windows")
            
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 过滤数据时出错: {e}")
            return False
    
    def _cleanup_temp_files(self):
        """🧹 清理临时文件 | Clean up temporary files"""
        for temp_file in self.temp_files:
            try:
                if temp_file.exists():
                    temp_file.unlink()
            except Exception as e:
                self.logger.warning(f"⚠️ 清理临时文件失败: {temp_file}, 错误: {e}")

    def _check_and_create_index(self, bam_file: str) -> bool:
        """🔍 检查并创建BAM索引文件 | Check and create BAM index file"""
        # 检查可能的索引文件位置 | Check possible index file locations
        possible_index_files = [
            f"{bam_file}.bai",      # 标准位置 | Standard location
            f"{bam_file}.csi",      # CSI索引 | CSI index
            bam_file.replace('.bam', '.bai'),  # 同名但不同扩展名 | Same name different extension
        ]
        
        # 检查是否已存在索引 | Check if index already exists
        for index_file in possible_index_files:
            if os.path.exists(index_file):
                self.logger.info(f"✅ 找到索引文件: {index_file} | Found index file: {index_file}")
                return True
        
        # 如果没有索引，创建索引 | If no index exists, create one
        self.logger.info(f"📋 未找到索引文件，正在创建... | Index not found, creating...")
        
        # index_cmd = f"{self.config.samtools_path} index {bam_file}"
        index_cmd = f"{self.config.samtools_path} index -@ {self.config.threads} {bam_file}"
        
        self.logger.info(f"🔧 创建索引命令: {index_cmd} | Creating index command: {index_cmd}")
        
        try:
            result = subprocess.run(
                index_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.cmd_runner.working_dir
            )
            
            self.logger.info(f"✅ 索引创建成功: {bam_file}.bai | Index created successfully: {bam_file}.bai")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 索引创建失败: {e.stderr} | Index creation failed: {e.stderr}")
            self.logger.warning(f"⚠️ 尝试不使用索引进行分析... | Trying analysis without index...")
            return False  # 返回False但不终止，尝试不使用区间筛选 | Return False but don't terminate, try without region filtering

    def _get_window_output_file(self) -> str:
        """📊 获取窗口分析输出文件名 | Get window analysis output filename"""
        base_path = Path(self.config.output_file)
        window_file = base_path.with_name(f"{base_path.stem}_windows_{self.config.window_size}bp{base_path.suffix}")
        return str(window_file)

    def _collect_window_data(self, window_data: dict, sample_name: str, chrom: str, pos: int, depth: float):
        """📊 收集窗口数据 | Collect window data"""
        # 确保步长不为0，避免除零错误 | Ensure step size is not 0 to avoid division by zero
        step_size = self.config.window_step if self.config.window_step > 0 else self.config.window_size
        window_start = ((pos - 1) // step_size) * step_size + 1
        window_end = window_start + self.config.window_size - 1
        
        # 创建窗口键 | Create window key
        window_key = (sample_name, chrom, window_start, window_end)
        
        if window_key not in window_data:
            window_data[window_key] = {'depths': [], 'positions': []}
        
        window_data[window_key]['depths'].append(depth)
        window_data[window_key]['positions'].append(pos)

    def _write_window_data(self, window_data: dict):
        """📝 写入窗口分析数据 | Write window analysis data"""
        window_output_file = self._get_window_output_file()
        
        with open(window_output_file, 'a') as f:
            for (sample_name, chrom, window_start, window_end), data in window_data.items():
                if data['depths']:
                    avg_depth = sum(data['depths']) / len(data['depths'])
                    window_center = (window_start + window_end) / 2
                    data_points = len(data['depths'])
                    
                    f.write(f"{sample_name}\t{chrom}\t{window_start}\t{window_end}\t{window_center:.1f}\t{avg_depth:.3f}\t{data_points}\n")