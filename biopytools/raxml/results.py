"""
🌳 RAxML分析结果处理模块 | RAxML Analysis Results Processing Module
"""

import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List

class ResultsManager:
    """结果管理器 | Results Manager"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def collect_output_files(self) -> Dict[str, str]:
        """收集输出文件 | Collect output files"""
        self.logger.info("📂 收集分析结果文件 | Collecting analysis result files")
        
        output_files = {}
        
        # RAxML输出文件类型 | RAxML output file types
        file_types = {
            'best_tree': f"RAxML_bestTree.{self.config.output_name}",
            'bipartitions': f"RAxML_bipartitions.{self.config.output_name}",
            'bipartitions_labels': f"RAxML_bipartitionsBranchLabels.{self.config.output_name}",
            'bootstrap': f"RAxML_bootstrap.{self.config.output_name}",
            'info': f"RAxML_info.{self.config.output_name}",
            'log': f"RAxML_log.{self.config.output_name}",
            'parsimony_tree': f"RAxML_parsimonyTree.{self.config.output_name}",
            'result': f"RAxML_result.{self.config.output_name}",
            'reduced_alignment': f"{self.config.sequence_file}.reduced"
        }
        
        for file_type, filename in file_types.items():
            file_path = self.config.output_path / filename
            if file_path.exists():
                output_files[file_type] = str(file_path)
                size = file_path.stat().st_size
                self.logger.info(f"  ✅ {file_type}: {filename} ({size} bytes)")
            else:
                self.logger.debug(f"  ➖ {file_type}: {filename} (未生成 | not generated)")
        
        return output_files
    
    def parse_raxml_info(self) -> Dict[str, str]:
        """解析RAxML信息文件 | Parse RAxML info file"""
        info_file = self.config.output_path / f"RAxML_info.{self.config.output_name}"
        
        if not info_file.exists():
            self.logger.warning("⚠️ RAxML信息文件未找到 | RAxML info file not found")
            return {}
        
        info_data = {}
        
        try:
            with open(info_file, 'r') as f:
                content = f.read()
                
                # 提取关键信息 | Extract key information
                lines = content.split('\n')
                for line in lines:
                    if 'Final GAMMA  likelihood' in line:
                        info_data['final_likelihood'] = line.strip()
                    elif 'alpha' in line and 'GAMMA' in line:
                        info_data['gamma_alpha'] = line.strip()
                    elif 'Time for' in line:
                        info_data['runtime'] = line.strip()
                    elif 'Overall execution time' in line:
                        info_data['total_time'] = line.strip()
        
        except Exception as e:
            self.logger.warning(f"⚠️ 解析信息文件失败 | Failed to parse info file: {e}")
        
        return info_data
    
    def generate_summary_report(self):
        """生成总结报告 | Generate summary report"""
        report_file = self.config.output_path / f"{self.config.output_name}_summary.txt"
        
        self.logger.info("📝 生成分析总结报告 | Generating analysis summary report")
        
        # 收集输出文件信息 | Collect output file information
        output_files = self.collect_output_files()
        info_data = self.parse_raxml_info()
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("🌳 RAxML系统发育分析总结报告 | RAxML Phylogenetic Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 基本信息 | Basic information
            f.write("📋 基本信息 | Basic Information:\n")
            f.write(f"  分析时间 | Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"  输出名称 | Output name: {self.config.output_name}\n")
            f.write(f"  工作目录 | Working directory: {self.config.output_dir}\n")
            f.write(f"  RAxML版本 | RAxML version: {self.config.raxml_path}\n\n")
            
            # 输入参数 | Input parameters
            f.write("🔧 分析参数 | Analysis Parameters:\n")
            f.write(f"  序列文件 | Sequence file: {self.config.sequence_file}\n")
            f.write(f"  替换模型 | Substitution model: {self.config.model}\n")
            f.write(f"  算法类型 | Algorithm type: {self.config.algorithm}\n")
            f.write(f"  线程数 | Threads: {self.config.threads}\n")
            f.write(f"  运行次数 | Number of runs: {self.config.runs}\n")
            
            if self.config.outgroup:
                f.write(f"  外群 | Outgroup: {self.config.outgroup}\n")
            
            f.write(f"  速率异质性类别 | Rate heterogeneity categories: {self.config.categories}\n")
            f.write(f"  似然精度 | Likelihood epsilon: {self.config.likelihood_epsilon}\n\n")
            
            # 输出文件 | Output files
            f.write("📂 输出文件 | Output Files:\n")
            
            important_files = {
                'best_tree': '最佳树 | Best tree',
                'bipartitions': '支持值标注树 | Tree with support values',
                'bootstrap': 'Bootstrap树集合 | Bootstrap tree collection',
                'info': '分析信息 | Analysis information',
                'result': '结果树 | Result tree'
            }
            
            for file_key, description in important_files.items():
                if file_key in output_files:
                    filename = Path(output_files[file_key]).name
                    f.write(f"  - {description}: {filename}\n")
            
            # 分析结果 | Analysis results
            if info_data:
                f.write(f"\n🎯 分析结果 | Analysis Results:\n")
                for key, value in info_data.items():
                    f.write(f"  - {key}: {value}\n")
            
            f.write(f"\n🔍 使用说明 | Usage Notes:\n")
            f.write("  1. 最佳树文件可用于进一步的系统发育分析\n")
            f.write("     Best tree file can be used for further phylogenetic analysis\n")
            f.write("  2. 如有bootstrap分析，支持值已标注在相应树文件中\n")
            f.write("     If bootstrap analysis performed, support values are annotated in corresponding tree files\n")
            f.write("  3. 建议使用FigTree、iTOL等软件进行树的可视化\n")
            f.write("     Recommend using FigTree, iTOL and other software for tree visualization\n\n")
        
        self.logger.info(f"📄 总结报告已生成 | Summary report generated: {report_file}")
        return str(report_file)
    
    def validate_results(self) -> bool:
        """验证结果文件 | Validate result files"""
        self.logger.info("✅ 验证输出结果 | Validating output results")
        
        success = True
        
        # 检查主要输出文件 | Check main output files
        essential_files = [
            f"RAxML_info.{self.config.output_name}",
            f"RAxML_result.{self.config.output_name}"
        ]
        
        for filename in essential_files:
            file_path = self.config.output_path / filename
            if not file_path.exists():
                self.logger.error(f"❌ 关键文件未生成 | Essential file not generated: {filename}")
                success = False
            else:
                self.logger.info(f"✅ 关键文件存在 | Essential file exists: {filename}")
        
        # 检查最佳树文件 | Check best tree file
        best_tree_file = self.config.output_path / f"RAxML_bestTree.{self.config.output_name}"
        if best_tree_file.exists():
            try:
                with open(best_tree_file, 'r') as f:
                    tree_content = f.read().strip()
                    if tree_content and tree_content.endswith(';'):
                        self.logger.info("✅ 最佳树文件格式正确 | Best tree file format is correct")
                    else:
                        self.logger.warning("⚠️ 最佳树文件可能格式不正确 | Best tree file may have incorrect format")
            except Exception as e:
                self.logger.error(f"❌ 无法验证最佳树文件 | Cannot validate best tree file: {e}")
                success = False
        
        return success
