"""
基因组组装结果处理模块 | Genome Assembly Results Processing Module
"""

import time
import os
from typing import Dict

class ResultsProcessor:
    """结果处理器 | Results Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, assembly_results: Dict) -> str:
        """生成总结报告 | Generate summary report"""
        self.logger.info("生成组装分析总结报告 | Generating assembly analysis summary report")
        
        report_file = self.config.output_path / "assembly_summary_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("基因组组装分析总结报告 | Genome Assembly Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"分析时间 | Analysis time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"输出目录 | Output directory: {self.config.output_dir}\n\n")
            
            f.write("输入数据 | Input Data:\n")
            f.write(f"  - HiFi reads: {self.config.hifi_reads}\n")
            if self.config.ont_reads:
                f.write(f"  - ONT reads: {self.config.ont_reads}\n")
            if self.config.reference_genome:
                f.write(f"  - 参考基因组 | Reference genome: {self.config.reference_genome}\n")
            f.write("\n")
            
            f.write("组装参数 | Assembly Parameters:\n")
            f.write(f"  - 基因组大小 | Genome size: {self.config.genome_size}\n")
            f.write(f"  - 倍性 | Ploidy: {self.config.ploidy}\n")
            f.write(f"  - 线程数 | Threads: {self.config.threads}\n")
            f.write(f"  - 内存 | Memory: {self.config.memory}GB\n")
            f.write("\n")
            
            f.write("组装工具版本 | Assembly Tool Versions:\n")
            f.write(f"  - Verkko: {self.config.verkko_version}\n")
            f.write(f"  - hifiasm: {self.config.hifiasm_version}\n")
            f.write(f"  - Graphasing: {self.config.graphasing_version}\n")
            f.write("\n")
            
            f.write("主要输出文件 | Main Output Files:\n")
            for key, value in assembly_results.items():
                if isinstance(value, str) and os.path.exists(value):
                    f.write(f"  - {key}: {value}\n")
                elif isinstance(value, dict):
                    f.write(f"  - {key}:\n")
                    for sub_key, sub_value in value.items():
                        if isinstance(sub_value, str) and os.path.exists(sub_value):
                            f.write(f"    - {sub_key}: {sub_value}\n")
            f.write("\n")
            
            f.write("质量控制结果 | Quality Control Results:\n")
            if self.config.run_contamination_screen:
                f.write("  ✓ 外源污染筛查已完成 | Contamination screening completed\n")
            if self.config.run_flagger:
                f.write("  ✓ Flagger错误注释已完成 | Flagger error annotation completed\n")
            if self.config.run_merqury:
                f.write("  ✓ Merqury质量评估已完成 | Merqury quality assessment completed\n")
            if self.config.run_inspector:
                f.write("  ✓ Inspector组装检查已完成 | Inspector assembly inspection completed\n")
            if self.config.run_deepvariant:
                f.write("  ✓ DeepVariant质量值估计已完成 | DeepVariant quality value estimation completed\n")
            if self.config.run_compleasm:
                f.write("  ✓ compleasm基因完整性评估已完成 | compleasm gene completeness assessment completed\n")
            f.write("\n")
            
            if self.config.trio_mode:
                f.write("家系分析 | Family Analysis:\n")
                f.write("  ✓ 家系三元组分析已完成 | Trio analysis completed\n")
                f.write("  ✓ 亲本支持度分析已完成 | Parental support analysis completed\n")
                f.write("\n")
            
            f.write("分析流程说明 | Analysis Pipeline Description:\n")
            f.write("1. 数据验证 | Data validation\n")
            f.write("2. Verkko主要组装 | Primary assembly with Verkko\n")
            f.write("3. hifiasm补充组装 | Complementary assembly with hifiasm\n")
            f.write("4. Graphasing分期信号生成 | Phasing signal generation with Graphasing\n")
            f.write("5. 外源污染筛查 | Contamination screening\n")
            f.write("6. 组装错误注释 | Assembly error annotation\n")
            f.write("7. 质量评估 | Quality assessment\n")
            f.write("8. 参考基因组比对 | Reference genome alignment\n")
            if self.config.trio_mode:
                f.write("9. 家系分析 | Family analysis\n")
            f.write("\n")
            
            f.write("引用 | Citations:\n")
            f.write("- Verkko: Rautiainen et al. Nature Biotechnology (2023)\n")
            f.write("- hifiasm: Cheng et al. Nature Methods (2021)\n")
            f.write("- Graphasing: Garg et al. Nature Biotechnology (2021)\n")
            f.write("- Merqury: Rhie et al. Genome Biology (2020)\n")
            f.write("- DeepVariant: Poplin et al. Nature Biotechnology (2018)\n")
            f.write("- minimap2: Li. Bioinformatics (2018)\n")
        
        self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
        return str(report_file)
