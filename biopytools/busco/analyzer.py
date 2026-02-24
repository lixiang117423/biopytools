"""
BUSCO分析核心模块|BUSCO Analysis Core Module
"""

import json
import os
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
from .utils import CommandRunner

class BUSCORunner:
    """BUSCO运行器|BUSCO Runner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_busco_analysis(self, input_file: str, sample_name: str) -> Tuple[bool, Optional[Dict[str, Any]]]:
        """运行BUSCO分析|Run BUSCO analysis"""
        output_name = f"{sample_name}_busco"

        self.logger.info(f"开始BUSCO分析|Starting BUSCO analysis: {sample_name}")
        self.logger.info(f"分析模式|Analysis mode: {self.config.mode}")
        self.logger.info(f"数据库|Database: {self.config.lineage}")

        # 构建BUSCO命令|Build BUSCO command
        cmd = self.build_busco_command(input_file, output_name)

        # 执行BUSCO分析|Execute BUSCO analysis
        success, stdout, stderr = self.cmd_runner.run(
            cmd, f"BUSCO分析|BUSCO analysis: {sample_name}"
        )

        if success:
            # 解析结果|Parse results
            result_data = self.parse_busco_results(output_name, sample_name)
            if result_data:
                self.logger.info(f"BUSCO分析完成|BUSCO analysis completed: {sample_name}")
                return True, result_data
            else:
                self.logger.error(f"无法解析BUSCO结果|Cannot parse BUSCO results: {sample_name}")
                return False, None
        else:
            self.logger.error(f"BUSCO分析失败|BUSCO analysis failed: {sample_name}")
            return False, None
    
    def build_busco_command(self, input_file: str, output_name: str) -> str:
        """构建BUSCO命令|Build BUSCO command"""
        cmd_parts = [
            self.config.busco_path,
            f"-i {input_file}",
            f"-l {self.config.lineage}",
            f"-o {output_name}",
            f"-m {self.config.mode}",
            f"-c {self.config.threads}"
        ]
        
        # 添加高级参数|Add advanced parameters
        if self.config.force:
            cmd_parts.append("-f")
        
        if self.config.augustus:
            cmd_parts.append("--augustus")
        
        if self.config.augustus_parameters:
            cmd_parts.append(f'--augustus_parameters "{self.config.augustus_parameters}"')
        
        if self.config.augustus_species:
            cmd_parts.append(f"--augustus_species {self.config.augustus_species}")
        
        if self.config.auto_lineage:
            cmd_parts.append("--auto-lineage")
        
        if self.config.auto_lineage_euk:
            cmd_parts.append("--auto-lineage-euk")
        
        if self.config.auto_lineage_prok:
            cmd_parts.append("--auto-lineage-prok")
        
        if self.config.contig_break != 10:
            cmd_parts.append(f"--contig_break {self.config.contig_break}")
        
        if self.config.datasets_version != 'odb12':
            cmd_parts.append(f"--datasets_version {self.config.datasets_version}")
        
        if self.config.download_path:
            cmd_parts.append(f"--download_path {self.config.download_path}")
        
        if self.config.evalue != 1e-3:
            cmd_parts.append(f"-e {self.config.evalue}")
        
        if self.config.limit != 3:
            cmd_parts.append(f"--limit {self.config.limit}")
        
        if self.config.long:
            cmd_parts.append("--long")
        
        if self.config.metaeuk:
            cmd_parts.append("--metaeuk")
        
        if self.config.metaeuk_parameters:
            cmd_parts.append(f'--metaeuk_parameters "{self.config.metaeuk_parameters}"')
        
        if self.config.metaeuk_rerun_parameters:
            cmd_parts.append(f'--metaeuk_rerun_parameters "{self.config.metaeuk_rerun_parameters}"')
        
        if self.config.miniprot:
            cmd_parts.append("--miniprot")
        
        if self.config.skip_bbtools:
            cmd_parts.append("--skip_bbtools")
        
        if self.config.offline:
            cmd_parts.append("--offline")
        
        if self.config.restart:
            cmd_parts.append("-r")
        
        if self.config.quiet:
            cmd_parts.append("-q")
        
        if self.config.scaffold_composition:
            cmd_parts.append("--scaffold_composition")
        
        if self.config.tar:
            cmd_parts.append("--tar")
        
        return " ".join(cmd_parts)
    
    def parse_busco_results(self, output_name: str, sample_name: str) -> Optional[Dict[str, Any]]:
        """解析BUSCO结果|Parse BUSCO results"""
        # 查找JSON结果文件|Find JSON result file
        output_dir = self.config.output_path / output_name
        
        # 可能的JSON文件位置|Possible JSON file locations
        json_patterns = [
            output_dir / f"short_summary.specific.*.{output_name}.json",
            output_dir / f"short_summary.generic.*.{output_name}.json",
            output_dir / "short_summary.*.json"
        ]
        
        json_file = None
        for pattern in json_patterns:
            import glob
            matches = glob.glob(str(pattern))
            if matches:
                json_file = matches[0]
                break
        
        if not json_file or not os.path.exists(json_file):
            self.logger.error(f" 未找到BUSCO结果JSON文件|BUSCO result JSON file not found: {output_name}")
            return None
        
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # 提取关键结果信息|Extract key result information
            results = data.get('results', {})
            lineage_dataset = data.get('lineage_dataset', {})
            
            return {
                'sample_name': sample_name,
                'lineage_name': lineage_dataset.get('name', 'Unknown'),
                'complete_percentage': results.get('Complete percentage', 0),
                'complete_buscos': results.get('Complete BUSCOs', 0),
                'single_copy_percentage': results.get('Single copy percentage', 0),
                'single_copy_buscos': results.get('Single copy BUSCOs', 0),
                'multi_copy_percentage': results.get('Multi copy percentage', 0),
                'multi_copy_buscos': results.get('Multi copy BUSCOs', 0),
                'fragmented_percentage': results.get('Fragmented percentage', 0),
                'fragmented_buscos': results.get('Fragmented BUSCOs', 0),
                'missing_percentage': results.get('Missing percentage', 0),
                'missing_buscos': results.get('Missing BUSCOs', 0),
                'n_markers': results.get('n_markers', 0),
                'domain': results.get('domain', 'Unknown'),
                'status': '成功|Success'
            }
            
        except Exception as e:
            self.logger.error(f" 解析JSON文件失败|Failed to parse JSON file: {e}")
            return None
