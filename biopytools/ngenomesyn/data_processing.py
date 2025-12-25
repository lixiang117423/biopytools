"""
数据处理模块 | Data Processing Module
"""

import os
from typing import List, Dict, Optional
from Bio import SeqIO

class GenomeProcessor:
    """基因组处理器 | Genome Processor"""
    
    def __init__(self, logger):
        self.logger = logger
        self.processing_log = []  # 存储处理步骤
    
    def log_processing_step(self, step: str, description: str = "", file_path: str = ""):
        """记录数据处理步骤 | Log data processing step"""
        self.logger.info("📊" + "=" * 60)
        self.logger.info(f"📊 数据处理步骤 | Data processing step: {step}")
        if description:
            self.logger.info(f"📝 描述 | Description: {description}")
        if file_path:
            self.logger.info(f"📁 文件 | File: {file_path}")
        self.logger.info("📊" + "=" * 60)
        
        # 保存到处理历史
        self.processing_log.append({
            "step": step,
            "description": description,
            "file_path": file_path
        })
    
    def save_processing_script(self, output_dir: str):
        """保存数据处理脚本 | Save data processing script"""
        script_path = os.path.join(output_dir, "data_processing_steps.sh")
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 📊 数据处理步骤脚本 | Data processing steps script\n")
            f.write(f"# ⏰ 生成时间 | Generated at: $(date)\n\n")
            
            f.write("echo \"📊 数据处理步骤记录 | Data processing steps log\"\n")
            f.write("echo \"⚠️ 注意: 这些是Python内部处理步骤的记录 | Note: These are records of Python internal processing steps\"\n\n")
            
            for i, step_info in enumerate(self.processing_log, 1):
                f.write(f"# 📝 步骤 {i}: {step_info['step']}\n")
                f.write(f"echo \"📝 步骤 {i}: {step_info['description']}\"\n")
                if step_info['file_path']:
                    f.write(f"echo \"  📁 文件: {step_info['file_path']}\"\n")
                f.write("# 📝 这是Python内部处理，无法用shell命令重现\n\n")
            
            f.write("echo \"✅ 数据处理步骤记录完成 | Data processing steps log completed\"\n")
        
        # 设置可执行权限
        os.chmod(script_path, 0o755)
        
        self.logger.info(f"📝 数据处理脚本已保存 | Data processing script saved: {script_path}")
        return script_path
    
    def read_sample_map(self, sample_map_file: str) -> List[Dict[str, str]]:
        """读取sample map文件 | Read sample map file"""
        self.log_processing_step("read_sample_map", f"📂 读取样本映射文件", sample_map_file)
        
        samples = []
        try:
            with open(sample_map_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) != 2:
                        self.logger.warning(f"⚠️ 第{line_num}行格式错误，跳过 | Line {line_num} format error, skipping")
                        continue
                    
                    genome_file, genome_name = parts
                    
                    # 🔍 检查文件是否存在
                    if not os.path.exists(genome_file):
                        self.logger.error(f"❌ 基因组文件不存在 | Genome file not found: {genome_file}")
                        continue
                    
                    samples.append({
                        "order": len(samples) + 1,
                        "name": genome_name,
                        "file": os.path.abspath(genome_file)
                    })
                    
                    self.log_processing_step("sample_added", f"➕ 添加基因组: {genome_name}", genome_file)
            
            self.logger.info(f"✅ 成功读取 {len(samples)} 个基因组 | Successfully read {len(samples)} genomes")
            
            # 保存处理脚本
            try:
                output_dir = os.path.dirname(sample_map_file) if os.path.dirname(sample_map_file) else "./"
                self.save_processing_script(output_dir)
            except Exception as e:
                self.logger.warning(f"⚠️ 保存处理脚本失败 | Failed to save processing script: {e}")
            
            return samples
            
        except Exception as e:
            self.logger.error(f"💥 读取sample map文件失败 | Failed to read sample map file: {e}")
            raise
    
    def calculate_chromosome_lengths(self, genome_file: str) -> Dict[str, int]:
        """计算染色体长度 | Calculate chromosome lengths"""
        self.log_processing_step("calculate_chr_lengths", f"📏 计算染色体长度", genome_file)
        
        chr_lengths = {}
        try:
            for record in SeqIO.parse(genome_file, "fasta"):
                chr_lengths[record.id] = len(record.seq)
                self.logger.debug(f"🧬 染色体 {record.id}: {len(record.seq)} bp")
            
            self.logger.info(f"🧬 发现 {len(chr_lengths)} 条染色体 | Found {len(chr_lengths)} chromosomes")
            self.log_processing_step("chr_calculation_complete", f"🧬 完成染色体长度计算，共 {len(chr_lengths)} 条染色体")
            
            return chr_lengths
            
        except Exception as e:
            self.logger.error(f"💥 计算染色体长度失败 | Failed to calculate chromosome lengths: {e}")
            raise
    
    def filter_sequences_by_chromosome(self, input_fasta: str, output_fasta: str, selected_chromosomes: List[int]):
        """根据染色体索引过滤FASTA文件 | Filter FASTA file by chromosome indices"""
        self.log_processing_step("filter_sequences", f"🔧 过滤序列文件，选择染色体: {selected_chromosomes}", input_fasta)
        
        try:
            # 📖 读取所有记录
            records = list(SeqIO.parse(input_fasta, "fasta"))
            total_chromosomes = len(records)
            
            self.logger.info(f"📊 输入文件包含 {total_chromosomes} 条染色体 | Input file contains {total_chromosomes} chromosomes")
            
            # ✅ 验证索引范围
            valid_indices = []
            for idx in selected_chromosomes:
                if 1 <= idx <= total_chromosomes:
                    valid_indices.append(idx)
                else:
                    self.logger.warning(f"⚠️ 染色体索引超出范围，跳过 | Chromosome index out of range, skipping: {idx} (总数: {total_chromosomes})")
            
            if not valid_indices:
                raise ValueError("❌ 没有有效的染色体索引 | No valid chromosome indices")
            
            # 🎯 选择记录（转换为0-based索引）
            selected_records = []
            for idx in valid_indices:
                record = records[idx - 1]  # 转换为0-based索引
                selected_records.append(record)
                self.logger.debug(f"✅ 选择染色体 {idx}: {record.id} (长度: {len(record.seq)})")
                self.log_processing_step("chr_selected", f"🎯 选择染色体 {idx}: {record.id}", f"长度: {len(record.seq)} bp")
            
            # 💾 写入过滤后的文件
            SeqIO.write(selected_records, output_fasta, "fasta")
            
            self.logger.info(f"🎉 过滤完成 | Filtering completed:")
            self.logger.info(f"  📥 输入: {total_chromosomes} 条染色体")
            self.logger.info(f"  📤 输出: {len(selected_records)} 条染色体")
            self.logger.info(f"  💾 保存到: {output_fasta}")
            
            self.log_processing_step("filter_complete", f"🎉 过滤完成，从 {total_chromosomes} 条染色体中选择了 {len(selected_records)} 条", output_fasta)
            
            return output_fasta
            
        except Exception as e:
            self.logger.error(f"💥 过滤序列文件失败 | Failed to filter sequence file: {e}")
            raise
    
    def generate_len_file(self, genome_file: str, output_file: str, color_scheme: str = "paired", selected_chromosomes: Optional[List[int]] = None) -> str:
        """生成.len文件 | Generate .len file"""
        self.log_processing_step("generate_len_file", f"📏 生成基因组长度文件", genome_file)
        
        chr_lengths = self.calculate_chromosome_lengths(genome_file)
        
        # 如果指定了染色体，则过滤
        if selected_chromosomes:
            self.logger.debug(f"🧬 过滤染色体用于len文件 | Filtering chromosomes for len file: {selected_chromosomes}")
            chr_items = list(chr_lengths.items())
            filtered_chr = {}
            
            for idx in selected_chromosomes:
                if 1 <= idx <= len(chr_items):
                    chr_name, length = chr_items[idx - 1]  # 转换为0-based索引
                    filtered_chr[chr_name] = length
                    self.logger.debug(f"✅ 包含染色体 {idx}: {chr_name} (长度: {length})")
                else:
                    self.logger.warning(f"⚠️ 染色体索引超出范围，跳过 | Chromosome index out of range, skipping: {idx}")
            
            chr_lengths = filtered_chr
            self.logger.info(f"🧬 过滤后的染色体数量 | Filtered chromosome count: {len(chr_lengths)}")
            self.log_processing_step("chr_filtered", f"🧬 染色体过滤完成，保留 {len(chr_lengths)} 条染色体")
        
        # 生成颜色方案
        colors = self._generate_chromosome_colors(len(chr_lengths), color_scheme)
        self.log_processing_step("color_scheme", f"🎨 生成颜色方案: {color_scheme}，{len(colors)} 种颜色")
        
        with open(output_file, 'w') as f:
            for i, (chr_name, length) in enumerate(chr_lengths.items()):
                color = colors[i % len(colors)]
                f.write(f"{chr_name}\t1\t{length}\tfill=\"{color}\"\n")
        
        self.logger.info(f"✅ 生成.len文件 | Generated .len file: {output_file}")
        self.log_processing_step("len_file_complete", f"✅ 完成.len文件生成", output_file)
        
        return output_file
    
    def _generate_chromosome_colors(self, chr_count: int, scheme: str = "paired") -> List[str]:
        """生成染色体颜色 | Generate chromosome colors"""
        color_schemes = {
            "paired": ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
                      "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB"],
            "set1": ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF"],
            "blues": ["#08519C", "#3182BD", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7"],
            "greens": ["#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0"]
        }
        
        base_colors = color_schemes.get(scheme, color_schemes["paired"])
        
        # 如果染色体数量超过预定义颜色，则重复使用
        colors = []
        for i in range(chr_count):
            colors.append(base_colors[i % len(base_colors)])
        
        self.logger.debug(f"🎨 为 {chr_count} 条染色体生成颜色方案 | Generated color scheme for {chr_count} chromosomes")
        
        return colors
    
    def get_chromosome_info(self, genome_file: str, selected_chromosomes: Optional[List[int]] = None) -> Dict[str, Dict]:
        """获取染色体信息 | Get chromosome information"""
        self.log_processing_step("get_chr_info", f"📋 获取染色体信息", genome_file)
        
        chr_lengths = self.calculate_chromosome_lengths(genome_file)
        chr_info = {}
        
        chr_items = list(chr_lengths.items())
        indices_to_process = selected_chromosomes if selected_chromosomes else range(1, len(chr_items) + 1)
        
        for idx in indices_to_process:
            if 1 <= idx <= len(chr_items):
                chr_name, length = chr_items[idx - 1]
                chr_info[chr_name] = {
                    "index": idx,
                    "length": length,
                    "name": chr_name
                }
                self.logger.debug(f"📋 染色体信息 {idx}: {chr_name} ({length} bp)")
        
        self.logger.info(f"📋 获取 {len(chr_info)} 条染色体信息 | Retrieved {len(chr_info)} chromosome info")
        self.log_processing_step("chr_info_complete", f"📋 完成染色体信息获取，共 {len(chr_info)} 条染色体")
        
        return chr_info