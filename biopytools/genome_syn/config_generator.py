"""
配置文件生成模块 | Configuration File Generator Module
"""

import os
import yaml
import pandas as pd
from typing import List, Dict, Tuple
from .data_processing import GenomeProcessor

class ConfigGenerator:
    """配置文件生成器 | Configuration File Generator"""
    
    def __init__(self, logger):
        self.logger = logger
        self.processor = GenomeProcessor(logger)
        self.config_log = []  # 存储配置生成过程
    
    def log_config_step(self, step: str, description: str = ""):
        """记录配置生成步骤 | Log configuration generation step"""
        self.logger.info(f"📋 配置生成步骤 | Config step: {step} - {description}")
        self.config_log.append({
            "step": step,
            "description": description
        })
    
    def save_config_generation_script(self, output_dir: str):
        """保存配置生成脚本 | Save configuration generation script"""
        script_path = os.path.join(output_dir, "regenerate_config.sh")
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 📋 配置重新生成脚本 | Configuration regeneration script\n")
            f.write(f"# ⏰ 生成时间 | Generated at: $(date)\n\n")
            
            f.write("echo \"📋 重新生成配置文件 | Regenerating configuration files\"\n")
            f.write("echo \"⚠️ 注意: 这需要原始的sample map文件 | Note: This requires the original sample map file\"\n\n")
            
            for i, step_info in enumerate(self.config_log, 1):
                f.write(f"# 📝 步骤 {i}: {step_info['step']}\n")
                f.write(f"echo \"📝 步骤 {i}: {step_info['description']}\"\n")
                f.write("# 📝 这是Python内部处理步骤，无法用shell命令重现\n\n")
            
            f.write("echo \"✅ 完成配置重新生成 | Completed configuration regeneration\"\n")
        
        # 设置可执行权限
        os.chmod(script_path, 0o755)
        
        self.logger.info(f"📝 配置生成脚本已保存 | Config generation script saved: {script_path}")
        return script_path
    
    def generate_default_config(self, samples: List[Dict], output_dir: str, config: 'GenomeSynConfig' = None) -> Tuple[str, str]:
        """生成默认配置文件 | Generate default configuration files"""
        self.log_config_step("start", "开始生成默认配置文件")
        
        # 为每个基因组计算染色体信息
        self.log_config_step("chromosome_analysis", "分析基因组染色体信息")
        for sample in samples:
            if config and config.chromosomes:
                # 如果指定了染色体过滤，获取过滤后的信息
                self.log_config_step("chromosome_filter", f"🧬 过滤基因组 {sample['name']} 的染色体: {config.chromosomes}")
                chr_info = self.processor.get_chromosome_info(sample["file"], config.chromosomes)
                sample["chromosomes"] = {name: info["length"] for name, info in chr_info.items()}
                sample["chr_count"] = len(chr_info)
                sample["selected_chromosomes"] = config.chromosomes
                self.logger.info(f"🧬 基因组 {sample['name']} 过滤后包含 {len(chr_info)} 条染色体")
            else:
                # 计算所有染色体
                self.log_config_step("chromosome_all", f"📊 计算基因组 {sample['name']} 的所有染色体")
                chr_lengths = self.processor.calculate_chromosome_lengths(sample["file"])
                sample["chromosomes"] = chr_lengths
                sample["chr_count"] = len(chr_lengths)
                sample["selected_chromosomes"] = None
        
        # 生成链式比对关系
        self.log_config_step("generate_links", "🔗 生成基因组比对链接关系")
        links = self._generate_chain_links(samples, config)
        
        # 生成YAML配置 - 使用传入的配置参数而不是硬编码值
        self.log_config_step("yaml_config", "📝 构建YAML配置数据结构")
        config_data = {
            "genomes": samples,
            "alignment": {
                "mode": config.alignment_mode if config else "chain",
                "aligner": config.aligner if config else "minimap2",
                "threads": config.threads if config else 16,
                "min_length": config.min_length if config else 5000
            },
            "links": links,
            "visualization": {
                "canvas_width": config.canvas_width if config and config.canvas_width else self._calculate_canvas_width(len(samples)),
                "canvas_height": config.canvas_height if config and config.canvas_height else self._calculate_canvas_height(len(samples)),
                "output_formats": config.output_formats if config else ["svg", "png"]
            }
        }
        
        # 添加染色体过滤信息到配置
        if config and config.chromosomes:
            self.log_config_step("chromosome_config", f"🧬 添加染色体过滤配置: {config.chromosomes}")
            config_data["chromosome_filter"] = {
                "enabled": True,
                "chromosomes": config.chromosomes,
                "description": f"分析染色体: {', '.join(map(str, config.chromosomes))}"
            }
            self.logger.info(f"🧬 使用染色体过滤 | Using chromosome filter: {config.chromosomes}")
        else:
            config_data["chromosome_filter"] = {
                "enabled": False,
                "chromosomes": None,
                "description": "分析所有染色体"
            }
        
        # 记录使用的线程数
        if config:
            self.log_config_step("threads_config", f"⚡ 配置线程数: {config.threads}")
            self.logger.info(f"⚡ 使用配置线程数 | Using configured threads: {config.threads}")
        
        # 保存YAML文件
        self.log_config_step("save_yaml", "💾 保存YAML配置文件")
        yaml_file = os.path.join(output_dir, "genome_syn_config.yaml")
        with open(yaml_file, 'w', encoding='utf-8') as f:
            yaml.dump(config_data, f, default_flow_style=False, allow_unicode=True, indent=2)
        
        # 转换为Excel文件
        self.log_config_step("convert_excel", "📊 转换为Excel配置文件")
        excel_file = self._yaml_to_excel(config_data, output_dir)
        
        # 保存配置生成脚本
        self.save_config_generation_script(output_dir)
        
        self.logger.info(f"📄 YAML配置文件已生成 | YAML config file generated: {yaml_file}")
        self.logger.info(f"📊 Excel配置文件已生成 | Excel config file generated: {excel_file}")
        
        return yaml_file, excel_file
    
    def _generate_chain_links(self, samples: List[Dict], config: 'GenomeSynConfig' = None) -> List[Dict]:
        """生成链式链接关系 | Generate chain links"""
        self.log_config_step("chain_links", f"🔗 为 {len(samples)} 个基因组生成链式比对关系")
        
        links = []
        chr_suffix = config.get_chromosome_suffix() if config else ""
        
        for i in range(len(samples) - 1):
            link_data = {
                "from": samples[i]["order"],           # 保持原来的顺序
                "to": samples[i + 1]["order"],         # 保持原来的顺序
                "from_name": samples[i]["name"],       # 保持原来的顺序
                "to_name": samples[i + 1]["name"],     # 保持原来的顺序
                "aligner": config.aligner if config else "minimap2",
                "params": {
                    "preset": "asm5",
                    "min_length": config.min_length if config else 5000
                }
            }
            
            # 如果有染色体过滤，添加相关信息
            if config and config.chromosomes:
                link_data["chromosome_filter"] = config.chromosomes
                link_data["output_suffix"] = chr_suffix
                self.log_config_step("link_chromosome", f"🧬 为链接 {samples[i]['name']} vs {samples[i + 1]['name']} 添加染色体过滤")
            
            links.append(link_data)
        
        return links
    
    def _calculate_canvas_width(self, genome_count: int) -> int:
        """计算画布宽度 | Calculate canvas width"""
        base_width = 1200
        if genome_count <= 2:
            return base_width
        elif genome_count <= 5:
            return int(base_width * 1.2)
        else:
            return int(base_width * 1.5)
    
    def _calculate_canvas_height(self, genome_count: int) -> int:
        """计算画布高度 | Calculate canvas height"""
        base_height = 800
        height_per_genome = 150
        return base_height + (genome_count - 2) * height_per_genome
    
    def _yaml_to_excel(self, config_data: Dict, output_dir: str) -> str:
        """将YAML配置转换为Excel | Convert YAML config to Excel"""
        excel_file = os.path.join(output_dir, "genome_syn_config.xlsx")
        
        self.log_config_step("excel_conversion", "🔄 开始YAML到Excel的转换")
        
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            # 基因组信息表
            genomes_df = pd.DataFrame(config_data["genomes"])
            # 移除复杂的字典列
            if "chromosomes" in genomes_df.columns:
                genomes_df = genomes_df.drop("chromosomes", axis=1)
            genomes_df.to_excel(writer, sheet_name="Genomes", index=False)
            
            # 链接关系表
            links_df = pd.DataFrame(config_data["links"])
            # 展开params字典
            if not links_df.empty and "params" in links_df.columns:
                params_df = pd.json_normalize(links_df["params"])
                links_df = pd.concat([links_df.drop("params", axis=1), params_df], axis=1)
            links_df.to_excel(writer, sheet_name="Links", index=False)
            
            # 可视化设置表
            viz_df = pd.DataFrame([config_data["visualization"]])
            viz_df.to_excel(writer, sheet_name="Visualization", index=False)
            
            # 比对设置表
            align_df = pd.DataFrame([config_data["alignment"]])
            align_df.to_excel(writer, sheet_name="Alignment", index=False)
            
            # 染色体过滤设置表
            if "chromosome_filter" in config_data:
                chr_filter_df = pd.DataFrame([config_data["chromosome_filter"]])
                chr_filter_df.to_excel(writer, sheet_name="ChromosomeFilter", index=False)
        
        self.log_config_step("excel_complete", f"✅ Excel文件保存完成: {excel_file}")
        return excel_file
    
    def excel_to_yaml(self, excel_file: str, output_dir: str) -> str:
        """将Excel配置转换回YAML | Convert Excel config back to YAML"""
        self.log_config_step("excel_to_yaml", f"📊 读取Excel配置文件: {excel_file}")
        
        try:
            # 读取各个工作表
            excel_data = pd.read_excel(excel_file, sheet_name=None)
            
            self.log_config_step("excel_parse", "📊 解析Excel工作表数据")
            
            # 重构配置数据
            config_data = {
                "genomes": excel_data["Genomes"].to_dict("records"),
                "links": excel_data["Links"].to_dict("records"),
                "visualization": excel_data["Visualization"].iloc[0].to_dict(),
                "alignment": excel_data["Alignment"].iloc[0].to_dict()
            }
            
            # 读取染色体过滤设置（如果存在）
            if "ChromosomeFilter" in excel_data:
                chr_filter = excel_data["ChromosomeFilter"].iloc[0].to_dict()
                config_data["chromosome_filter"] = chr_filter
                self.log_config_step("chromosome_parse", "🧬 解析染色体过滤设置")
            
            # 保存为YAML
            yaml_file = os.path.join(output_dir, "genome_syn_config_from_excel.yaml")
            with open(yaml_file, 'w', encoding='utf-8') as f:
                yaml.dump(config_data, f, default_flow_style=False, allow_unicode=True, indent=2)
            
            self.log_config_step("yaml_save", f"📄 YAML配置文件已生成: {yaml_file}")
            self.logger.info(f"📄 YAML配置文件已生成 | YAML config file generated: {yaml_file}")
            return yaml_file
            
        except Exception as e:
            self.logger.error(f"💥 Excel转YAML失败 | Excel to YAML conversion failed: {e}")
            raise