"""
🧠 主分析模块 | Main Analysis Module
"""

import os
import sys
import yaml
import logging
import shutil
import argparse
from typing import Dict, List

class GenomeSynAnalyzer:
    """🧠 基因组共线性分析器 | Genome Synteny Analyzer"""
    
    def __init__(self, config):
        self.config = config
        self.logger = self._setup_logger()
        
        # 🔧 初始化各个组件
        from .environment import EnvironmentChecker
        from .data_processing import GenomeProcessor
        from .config_generator import ConfigGenerator
        from .visualizer import NGenomeSynVisualizer
        from .aligners import AlignerFactory
        
        self.env_checker = EnvironmentChecker(self.logger)
        self.processor = GenomeProcessor(self.logger)
        self.config_generator = ConfigGenerator(self.logger)
        self.visualizer = NGenomeSynVisualizer(self.logger)
        self.aligner_factory = AlignerFactory()
    
    def _setup_logger(self) -> logging.Logger:
        """🔧 设置日志记录器 | Setup logger"""
        logger = logging.getLogger("GenomeSyn")
        logger.setLevel(logging.INFO)
        
        if not logger.handlers:
            handler = logging.StreamHandler(sys.stdout)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        
        return logger
    
    def run_analysis(self):
        """🚀 运行完整分析流程 | Run complete analysis pipeline"""
        try:
            self.logger.info("🎯" + "=" * 60)
            self.logger.info("🧬 开始基因组共线性分析 | Starting genome synteny analysis")
            if self.config.chromosomes:
                self.logger.info(f"🧬 分析染色体 | Analyzing chromosomes: {self.config.chromosomes}")
            self.logger.info("🎯" + "=" * 60)
            
            # ✅ 验证配置
            self.config.validate()
            
            # 🔍 检查环境
            env_results = self.env_checker.check_environment(self.config.aligner)
            if not env_results["NGenomeSyn"]:
                raise RuntimeError("❌ NGenomeSyn未安装 | NGenomeSyn not installed")
            
            # 📋 处理配置文件生成模式
            if self.config.generate_config:
                return self._generate_config_only()
            
            # 🚀 运行完整分析
            return self._run_full_analysis()
            
        except Exception as e:
            self.logger.error(f"💥 分析失败 | Analysis failed: {e}")
            raise
    
    def _generate_config_only(self):
        """📋 仅生成配置文件模式 | Generate configuration only mode"""
        if not self.config.sample_map:
            raise ValueError("❌ 生成配置模式需要提供sample_map | Sample map required for config generation mode")
        
        self.logger.info("📋 配置文件生成模式 | Configuration generation mode")
        
        # 📂 读取样本信息
        samples = self.processor.read_sample_map(self.config.sample_map)
        
        # 🛠️ 生成配置文件 - 传递配置对象以获取正确的参数
        yaml_file, excel_file = self.config_generator.generate_default_config(
            samples, 
            self.config.output_dir,
            self.config  # 传递配置对象
        )
        
        self.logger.info("🎉" + "=" * 60)
        self.logger.info("✅ 配置文件已生成完成 | Configuration files generated")
        self.logger.info(f"📄 YAML配置文件 | YAML config: {yaml_file}")
        self.logger.info(f"📊 Excel配置文件 | Excel config: {excel_file}")
        self.logger.info("💡 请编辑Excel文件后使用 --config 参数运行分析 | Please edit Excel file and run with --config parameter")
        self.logger.info("🎉" + "=" * 60)
        
        return {"yaml_config": yaml_file, "excel_config": excel_file}
    
    def _run_full_analysis(self):
        """🚀 运行完整分析流程 | Run full analysis pipeline"""
        self.logger.info("🚀 开始完整分析流程 | Starting full analysis pipeline")
        
        # 📋 确定配置数据来源
        if self.config.config_file:
            # 📊 从配置文件读取
            if self.config.config_file.endswith('.xlsx'):
                yaml_file = self.config_generator.excel_to_yaml(self.config.config_file, self.config.output_dir)
            else:
                yaml_file = self.config.config_file
            
            with open(yaml_file, 'r', encoding='utf-8') as f:
                config_data = yaml.safe_load(f)
        else:
            # 📂 从sample_map生成 - 传递配置对象
            samples = self.processor.read_sample_map(self.config.sample_map)
            yaml_file, _ = self.config_generator.generate_default_config(
                samples, 
                self.config.output_dir,
                self.config  # 传递配置对象以获取正确的线程数等参数
            )
            with open(yaml_file, 'r', encoding='utf-8') as f:
                config_data = yaml.safe_load(f)
        
        # 🔗 步骤1: 执行比对
        self.logger.info("🔗" + "=" * 50)
        self.logger.info("🔗 步骤1: 执行基因组比对 | Step 1: Performing genome alignment")
        self.logger.info("🔗" + "=" * 50)
        
        self._run_alignments(config_data)
        
        # 🎨 步骤2: 生成可视化
        self.logger.info("\n" + "🎨" + "=" * 50) 
        self.logger.info("🎨 步骤2: 生成可视化图形 | Step 2: Generating visualization")
        self.logger.info("🎨" + "=" * 50)
        
        output_files = self._run_visualization(config_data)
        
        # 📝 步骤3: 生成完整的命令脚本
        self.logger.info("\n" + "📝" + "=" * 50)
        self.logger.info("📝 步骤3: 生成完整命令脚本 | Step 3: Generating complete command script")
        self.logger.info("📝" + "=" * 50)
        
        self._generate_complete_script(config_data)
        
        self.logger.info("🎉" + "=" * 60)
        self.logger.info("🎉 分析完成 | Analysis completed")
        for file_path in output_files:
            self.logger.info(f"📤 输出文件 | Output file: {file_path}")
        self.logger.info("🎉" + "=" * 60)
        
        return {"output_files": output_files}
    
    def _run_alignments(self, config_data: Dict):
        """🔗 执行比对 | Run alignments"""
        aligner_type = config_data["alignment"]["aligner"]
        threads = config_data["alignment"]["threads"]
        
        self.logger.info(f"🔧 创建比对器 | Creating aligner: {aligner_type} with {threads} threads")
        
        # 🔧 创建比对器时传递正确的线程数和输出目录
        aligner = self.aligner_factory.create_aligner(
            aligner_type, 
            self.logger, 
            threads,  # 使用配置中的线程数
            self.config.output_dir  # 添加输出目录参数
        )
        
        # 📋 创建基因组查找字典
        genomes_dict = {g["order"]: g for g in config_data["genomes"]}
        
        # 🔗 执行每个链接的比对
        for link in config_data["links"]:
            genome1 = genomes_dict[link["from"]]
            genome2 = genomes_dict[link["to"]]
            
            # 📄 生成文件名时考虑染色体过滤
            chr_suffix = self.config.get_chromosome_suffix()
            link_file = f"{link['from_name']}_vs_{link['to_name']}{chr_suffix}.link"
            paf_file = f"{link['from_name']}_vs_{link['to_name']}{chr_suffix}.paf"
            
            link_path = os.path.join(self.config.output_dir, link_file)
            paf_path = os.path.join(self.config.output_dir, paf_file)
            
            # 🔍 检查文件是否存在且非空
            skip_alignment = False
            if os.path.exists(link_path) and os.path.getsize(link_path) > 0:
                self.logger.info(f"⭐️ LINK文件已存在且有效，跳过比对 | Valid LINK file exists, skipping: {link_file}")
                skip_alignment = True
            elif os.path.exists(paf_path) and os.path.getsize(paf_path) > 0:
                self.logger.info(f"⭐️ PAF文件已存在且有效，跳过比对 | Valid PAF file exists, skipping: {paf_file}")
                skip_alignment = True
            elif os.path.exists(paf_path) and os.path.getsize(paf_path) == 0:
                self.logger.warning(f"⚠️ 发现空的PAF文件，将重新运行比对 | Found empty PAF file, will re-run alignment: {paf_file}")
                os.remove(paf_path)  # 删除空文件
            elif os.path.exists(link_path) and os.path.getsize(link_path) == 0:
                self.logger.warning(f"⚠️ 发现空的LINK文件，将重新运行比对 | Found empty LINK file, will re-run alignment: {link_file}")
                os.remove(link_path)  # 删除空文件
            
            if skip_alignment:
                continue
            
            # 🚀 执行比对，传递染色体过滤参数
            self.logger.info(f"🚀 开始比对 | Starting alignment: {genome1['name']} vs {genome2['name']}")
            result_file = aligner.align(
                genome1, 
                genome2, 
                selected_chromosomes=self.config.chromosomes,  # 传递染色体过滤参数
                **link.get("params", {})
            )
            
            self.logger.info(f"📄 比对文件已生成 | Alignment file generated: {result_file}")
            
            # ✅ 验证生成的文件不为空
            if os.path.exists(result_file) and os.path.getsize(result_file) == 0:
                self.logger.error(f"❌ 生成的比对文件为空 | Generated alignment file is empty: {result_file}")
                raise RuntimeError(f"❌ 比对失败，生成的文件为空 | Alignment failed, generated file is empty")
            
            self.logger.info(f"✅ 比对成功完成 | Alignment successfully completed: {os.path.basename(result_file)}")
    
    def _run_visualization(self, config_data: Dict) -> List[str]:
        """🎨 运行可视化 | Run visualization"""
        # 📋 生成NGenomeSyn配置文件，传递染色体过滤参数
        ngenomesyn_config = self.visualizer.generate_ngenomesyn_config(
            config_data, 
            self.config.output_dir,
            selected_chromosomes=self.config.chromosomes  # 传递染色体过滤参数
        )
        
        # 🎨 运行NGenomeSyn
        output_formats = config_data["visualization"]["output_formats"]
        output_files = self.visualizer.run_ngenomesyn(ngenomesyn_config, self.config.output_dir, output_formats)
        
        return output_files
    
    def _generate_complete_script(self, config_data: Dict):
        """📝 生成包含所有步骤的完整命令脚本 | Generate complete command script for all steps"""
        script_path = os.path.join(self.config.output_dir, "complete_genome_syn_pipeline.sh")
        
        # 获取染色体后缀用于文件命名
        chr_suffix = self.config.get_chromosome_suffix()
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 🧬 基因组共线性分析完整流程脚本 | Complete Genome Synteny Analysis Pipeline\n")
            f.write(f"# ⏰ 生成时间 | Generated at: $(date)\n")
            f.write(f"# 📁 输出目录 | Output directory: {self.config.output_dir}\n")
            if self.config.chromosomes:
                f.write(f"# 🧬 分析染色体 | Analyzed chromosomes: {', '.join(map(str, self.config.chromosomes))}\n")
            f.write("\n")
            
            f.write("set -e  # ⚠️ 遇到错误时停止执行 | Stop on error\n")
            f.write("set -u  # ⚠️ 使用未定义变量时报错 | Error on undefined variables\n\n")
            
            # 设置变量
            f.write("# 📝 设置变量 | Set variables\n")
            f.write(f"OUTPUT_DIR=\"{self.config.output_dir}\"\n")
            f.write(f"THREADS={self.config.threads}\n")
            f.write(f"MIN_LENGTH={self.config.min_length}\n")
            if self.config.chromosomes:
                f.write(f"CHROMOSOMES=\"{','.join(map(str, self.config.chromosomes))}\"\n")
            f.write("\n")
            
            f.write("echo \"🚀 开始基因组共线性分析流程 | Starting genome synteny analysis pipeline\"\n")
            f.write("echo \"📁 输出目录: $OUTPUT_DIR\"\n")
            if self.config.chromosomes:
                f.write("echo \"🧬 分析染色体: $CHROMOSOMES\"\n")
            f.write("\n")
            
            # 创建输出目录
            f.write("# 📁 创建输出目录 | Create output directory\n")
            f.write("mkdir -p \"$OUTPUT_DIR\"\n")
            f.write("cd \"$OUTPUT_DIR\"\n\n")
            
            # 第一步：生成基因组长度文件
            f.write("echo \"📋 步骤1: 生成基因组长度文件 | Step 1: Generate genome length files\"\n")
            for genome in config_data["genomes"]:
                len_filename = f"{genome['name']}{chr_suffix}.len"
                f.write(f"echo \"  📏 生成 {len_filename}\"\n")
                # 这里可以添加实际的生成命令，但由于这是Python函数调用，我们用注释表示
                f.write(f"# python -c \"from data_processing import GenomeProcessor; GenomeProcessor().generate_len_file('{genome['file']}', '{len_filename}', selected_chromosomes={self.config.chromosomes})\"\n")
            f.write("\n")
            
            # 第二步：执行基因组比对
            f.write("echo \"🔗 步骤2: 执行基因组比对 | Step 2: Perform genome alignment\"\n")
            for link in config_data["links"]:
                genome1_name = link['from_name']
                genome2_name = link['to_name']
                
                # 如果有染色体过滤，生成过滤文件
                if self.config.chromosomes:
                    f.write(f"echo \"  🧬 过滤染色体序列\"\n")
                    f.write(f"# 🧬 过滤 {genome1_name} 序列\n")
                    f.write(f"# 🧬 过滤 {genome2_name} 序列\n")
                
                # Minimap2比对命令
                paf_file = f"{genome1_name}_vs_{genome2_name}{chr_suffix}.paf"
                link_file = f"{genome1_name}_vs_{genome2_name}{chr_suffix}.link"
                
                f.write(f"echo \"  🔗 比对: {genome1_name} vs {genome2_name}\"\n")
                
                # 获取实际的基因组文件路径
                genome1_file = None
                genome2_file = None
                for genome in config_data["genomes"]:
                    if genome["name"] == genome1_name:
                        genome1_file = genome["file"]
                    if genome["name"] == genome2_name:
                        genome2_file = genome["file"]
                
                if genome1_file and genome2_file:
                    # 如果有染色体过滤，使用过滤后的文件
                    if self.config.chromosomes:
                        filtered_genome1 = f"{genome1_name}_filtered.fa"
                        filtered_genome2 = f"{genome2_name}_filtered.fa"
                        f.write(f"minimap2 -x asm5 -t $THREADS \"{filtered_genome2}\" \"{filtered_genome1}\" > \"{paf_file}\"\n")
                    else:
                        f.write(f"minimap2 -x asm5 -t $THREADS \"{genome2_file}\" \"{genome1_file}\" > \"{paf_file}\"\n")
                else:
                    f.write(f"# ⚠️ 注意: 请替换为实际的基因组文件路径\n")
                    f.write(f"# minimap2 -x asm5 -t $THREADS GENOME2_FILE GENOME1_FILE > \"{paf_file}\"\n")
                
                f.write("if [ $? -ne 0 ]; then\n")
                f.write("    echo \"❌ 错误: Minimap2比对失败\"\n")
                f.write("    exit 1\n")
                f.write("fi\n\n")
                
                # PAF转LINK格式
                f.write("echo \"  🔄 转换PAF到LINK格式\"\n")
                f.write("if command -v GetTwoGenomeSyn.pl >/dev/null 2>&1; then\n")
                f.write(f"    GetTwoGenomeSyn.pl Paf2Link \"{paf_file}\" $MIN_LENGTH \"{link_file}\"\n")
                f.write("    if [ $? -eq 0 ]; then\n")
                f.write(f"        echo \"    ✅ LINK文件生成成功: {link_file}\"\n")
                f.write("    else\n")
                f.write("        echo \"    ⚠️ 警告: PAF转换失败，将使用PAF文件\"\n")
                f.write("    fi\n")
                f.write("else\n")
                f.write("    echo \"    ⚠️ 警告: GetTwoGenomeSyn.pl不可用，将使用PAF文件\"\n")
                f.write("fi\n\n")
            
            # 第三步：生成NGenomeSyn配置文件
            config_filename = f"ngenomesyn{chr_suffix}.conf"
            f.write("echo \"📝 步骤3: 生成NGenomeSyn配置文件 | Step 3: Generate NGenomeSyn configuration\"\n")
            f.write(f"echo \"  📋 配置文件: {config_filename}\"\n")
            f.write("# 📋 配置文件内容将由Python脚本生成\n\n")
            
            # 第四步：运行NGenomeSyn可视化
            output_prefix = f"genome_synteny{chr_suffix}" if chr_suffix else "genome_synteny"
            f.write("echo \"🎨 步骤4: 生成可视化图形 | Step 4: Generate visualization\"\n")
            f.write("echo \"  🔍 检查NGenomeSyn是否可用\"\n")
            f.write("if ! command -v NGenomeSyn >/dev/null 2>&1; then\n")
            f.write("    echo \"❌ 错误: NGenomeSyn未安装或不在PATH中\"\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")
            
            f.write("echo \"  🚀 运行NGenomeSyn\"\n")
            f.write(f"NGenomeSyn -InConf \"{config_filename}\" -OutPut \"{output_prefix}\"\n")
            f.write("if [ $? -ne 0 ]; then\n")
            f.write("    echo \"⚠️ 警告: NGenomeSyn返回非零退出码，但可能已生成文件\"\n")
            f.write("fi\n\n")
            
            # 检查输出文件
            f.write("echo \"  🔍 检查生成的文件\"\n")
            f.write(f"if [ -f \"{output_prefix}.svg\" ] && [ -s \"{output_prefix}.svg\" ]; then\n")
            f.write(f"    echo \"    ✅ SVG文件生成成功: {output_prefix}.svg\"\n")
            f.write("else\n")
            f.write("    echo \"    ❌ 错误: SVG文件未生成或为空\"\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")
            
            # 第五步：转换为PNG（如果需要）
            if "png" in config_data["visualization"]["output_formats"]:
                f.write("echo \"🖼️ 步骤5: 转换为PNG格式 | Step 5: Convert to PNG format\"\n")
                f.write("echo \"  🔍 检查ImageMagick是否可用\"\n")
                f.write("if command -v convert >/dev/null 2>&1; then\n")
                f.write("    echo \"    🔄 转换SVG到PNG\"\n")
                f.write(f"    convert \"{output_prefix}.svg\" \"{output_prefix}.png\"\n")
                f.write(f"    if [ $? -eq 0 ] && [ -f \"{output_prefix}.png\" ] && [ -s \"{output_prefix}.png\" ]; then\n")
                f.write(f"        echo \"    ✅ PNG文件生成成功: {output_prefix}.png\"\n")
                f.write("    else\n")
                f.write("        echo \"    ⚠️ 警告: PNG转换失败\"\n")
                f.write("    fi\n")
                f.write("else\n")
                f.write("    echo \"    ⚠️ 警告: ImageMagick未安装，跳过PNG转换\"\n")
                f.write("fi\n\n")
            
            # 完成
            f.write("echo \"🎉 基因组共线性分析完成 | Genome synteny analysis completed\"\n")
            f.write("echo \"📂 输出文件列表 | Output files:\"\n")
            f.write(f"ls -la \"{output_prefix}\".*\n")
            
        # 设置可执行权限
        os.chmod(script_path, 0o755)
        
        self.logger.info(f"📝 完整流程脚本已生成 | Complete pipeline script generated: {script_path}")
        self.logger.info("🚀 您可以使用以下命令重新运行整个分析流程:")
        self.logger.info(f"bash {script_path}")
        self.logger.info("")
        self.logger.info("⚠️ 注意: 脚本中的某些步骤需要手动替换实际的输入文件路径")
        
        return script_path

def main():
    """💻 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="🧬 基因组共线性可视化工具 | Genome Synteny Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
🌟 使用示例 | Usage Examples:
  # 📋 生成配置文件 | Generate configuration file
  python run_genome_syn.py --sample-map genomes.tsv --generate-config --output-dir ./output
  
  # 🚀 运行分析 | Run analysis
  python run_genome_syn.py --config genomes_config.xlsx --output-dir ./output
  
  # ⚡ 一步完成 | One-step completion
  python run_genome_syn.py --sample-map genomes.tsv --output-dir ./output
  
  # 🧬 分析特定染色体 | Analyze specific chromosomes
  python run_genome_syn.py --sample-map genomes.tsv --output-dir ./output --chromosome "1,2,3"
  
  # 📊 分析染色体范围 | Analyze chromosome range
  python run_genome_syn.py --sample-map genomes.tsv --output-dir ./output --chromosome "1-5"
        """
    )
    
    # 📥 输入文件参数
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-s', '--sample-map', 
                           help='📂 样本映射文件 | Sample mapping file (tab-separated: genome_file\\tgenome_name)')
    input_group.add_argument('-c', '--config',
                           help='📋 配置文件 | Configuration file (.xlsx or .yaml)')
    
    # 📤 输出参数
    parser.add_argument('-o', '--output-dir', default='./genome_syn_output',
                       help='📁 输出目录 | Output directory (default: ./genome_syn_output)')
    
    # ⚙️ 模式参数
    parser.add_argument('--generate-config', action='store_true',
                       help='📋 仅生成配置文件 | Generate configuration file only')
    
    # 🔗 比对参数
    parser.add_argument('-a', '--aligner', default='minimap2', 
                       choices=['minimap2', 'mcscanx', 'syri', 'mummer'],
                       help='🔧 比对器类型 | Aligner type (default: minimap2)')
    parser.add_argument('--alignment-mode', default='chain',
                       choices=['chain', 'star', 'all_vs_all'],
                       help='🔗 比对模式 | Alignment mode (default: chain)')
    parser.add_argument('-t', '--threads', type=int, default=16,
                       help='⚡ 线程数 | Number of threads (default: 16)')
    parser.add_argument('--min-length', type=int, default=5000,
                       help='📏 最小比对长度 | Minimum alignment length (default: 5000)')
    
    # 🧬 染色体过滤参数
    parser.add_argument('--chromosome', type=str,
                       help='🧬 指定要分析的染色体 | Specify chromosomes to analyze (e.g., "1,2,3" or "1-5" or "1")')
    
    # 🎨 可视化参数
    parser.add_argument('--canvas-width', type=int,
                       help='📏 画布宽度 | Canvas width (auto-calculated if not specified)')
    parser.add_argument('--canvas-height', type=int,
                       help='📏 画布高度 | Canvas height (auto-calculated if not specified)')
    parser.add_argument('--output-formats', nargs='+', default=['svg', 'png'],
                       choices=['svg', 'png'],
                       help='📄 输出格式 | Output formats (default: svg png)')
    
    args = parser.parse_args()
    
    # ⚙️ 创建配置对象
    from .config import GenomeSynConfig
    config = GenomeSynConfig(
        sample_map=args.sample_map,
        config_file=args.config,
        output_dir=args.output_dir,
        aligner=args.aligner,
        alignment_mode=args.alignment_mode,
        threads=args.threads,
        min_length=args.min_length,
        canvas_width=args.canvas_width,
        canvas_height=args.canvas_height,
        output_formats=args.output_formats,
        generate_config=args.generate_config,
        _chromosome_str=args.chromosome  # 传递染色体字符串
    )
    
    # 🚀 运行分析
    analyzer = GenomeSynAnalyzer(config)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()