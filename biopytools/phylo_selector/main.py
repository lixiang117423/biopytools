"""
系统发育树样品选择主模块|Phylogenetic Tree Sample Selector Main Module

基于均匀间隔和PCA去重的智能选择算法|Intelligent selection based on uniform interval and PCA deduplication
"""

import sys
import os
from typing import List, Dict, Optional
import pandas as pd

from .config import PhyloSelectorConfig
from .calculator import PhyloSelectorCalculator
from .utils import PhyloSelectorLogger
from .parser import GroupTableParser, PCAFileParser


class PhyloSelectorRunner:
    """系统发育树样品选择运行器|Phylogenetic Tree Sample Selector Runner"""

    def __init__(self,
                 hierarchy_file: str,
                 pca_file: str,
                 output_prefix: str,
                 n_samples: int = 150,
                 group_file: Optional[str] = None,
                 newick_file: Optional[str] = None,
                 min_samples_per_group: int = 1,
                 hierarchy_level: int = 10,
                 pca_dedup_threshold: float = 0.0001,
                 generate_report: bool = True,
                 generate_csv: bool = True,
                 generate_visualization: bool = True):
        """初始化运行器|Initialize runner

        Args:
            hierarchy_file: 层级关系文件（必需）|Hierarchy file (required)
            pca_file: PCA坐标文件（必需）|PCA coordinates file (required)
            output_prefix: 输出文件前缀|Output file prefix
            n_samples: 选择样品总数|Total number of samples to select
            group_file: 分组文件（可选）|Group file (optional)
            newick_file: Newick树文件（可选，仅用于兼容性）|Newick tree file (optional, for compatibility only)
            min_samples_per_group: 每组最小样品数|Minimum samples per group
            hierarchy_level: 层级深度|Hierarchy level
            pca_dedup_threshold: PCA去重阈值|PCA dedup threshold
            generate_report: 是否生成报告|Whether to generate report
            generate_csv: 是否生成CSV|Whether to generate CSV
            generate_visualization: 是否生成可视化|Whether to generate visualization
        """
        self.config = PhyloSelectorConfig(
            hierarchy_file=hierarchy_file,
            pca_file=pca_file,
            output_prefix=output_prefix,
            newick_file=newick_file,
            group_file=group_file,
            n_samples=n_samples,
            min_samples_per_group=min_samples_per_group,
            hierarchy_level=hierarchy_level,
            pca_dedup_threshold=pca_dedup_threshold,
            generate_report=generate_report,
            generate_csv=generate_csv,
            generate_visualization=generate_visualization
        )

        logger_manager = PhyloSelectorLogger()
        self.logger = logger_manager.get_logger()
        self.calculator = PhyloSelectorCalculator(self.config, self.logger)

    def run(self) -> bool:
        """运行选择流程|Run selection process

        Returns:
            bool: 是否成功|Success status
        """
        try:
            # 验证配置|Validate configuration
            self.config.validate()

            self.logger.info("="*70)
            self.logger.info("系统发育树样品选择|Phylogenetic Tree Sample Selector")
            self.logger.info("基于均匀间隔和PCA去重的智能选择|Intelligent Selection Based on Uniform Interval and PCA Dedup")
            self.logger.info("="*70)

            # 步骤1: 读取层级文件|Step 1: Read hierarchy file
            self.logger.info("\n[步骤1|Step 1] 读取层级关系文件|Reading hierarchy file...")
            hierarchy_df = pd.read_excel(self.config.hierarchy_file)
            self.logger.info(f"  成功加载|Successfully loaded: {hierarchy_df.shape[0]} 行|rows x {hierarchy_df.shape[1]} 列|columns")
            self.logger.info(f"  层级列|Hierarchy columns: {list(hierarchy_df.columns[:5])}... (parent_1 to parent_{hierarchy_df.shape[1]-1})")

            # 步骤2: 读取PCA文件|Step 2: Read PCA file
            self.logger.info("\n[步骤2|Step 2] 读取PCA坐标文件|Reading PCA coordinates file...")
            pca_parser = PCAFileParser(self.logger)
            sample_to_pca = pca_parser.parse_pca_file(self.config.pca_file)
            self.logger.info(f"  已加载|Loaded PCA data for {len(sample_to_pca)} 个样品|samples")

            # 步骤3: 读取分组文件（可选）|Step 3: Read group file (optional)
            sample_to_group = None
            if self.config.group_file:
                self.logger.info("\n[步骤3|Step 3] 读取分组文件|Reading group file...")
                group_parser = GroupTableParser(self.logger)
                sample_to_group = group_parser.parse_group_file(self.config.group_file)
                self.logger.info(f"  已加载|Loaded group data for {len(sample_to_group)} 个样品|samples")
            else:
                self.logger.info("\n[步骤3|Step 3] 无分组文件，跳过|No group file, skip")

            # 步骤4: 构建样品列表|Step 4: Build sample list
            self.logger.info("\n[步骤4|Step 4] 构建样品列表|Building sample list...")
            all_samples = []
            for label in hierarchy_df['label']:
                if label in sample_to_pca:
                    all_samples.append({'name': label})

            self.logger.info(f"  有效样品数|Valid samples: {len(all_samples)}")

            if len(all_samples) < self.config.n_samples:
                self.logger.warning(f"  有效样品数({len(all_samples)})少于目标数量({self.config.n_samples})|Valid samples less than target")
                self.logger.warning(f"  将选择所有可用样品|Will select all available samples")

            # 步骤5: 执行选择|Step 5: Execute selection
            self.logger.info("\n[步骤5|Step 5] 执行智能选择|Executing intelligent selection...")
            selected_samples = self.calculator.hierarchy_based_selection(
                all_samples,
                hierarchy_df,
                sample_to_pca,
                sample_to_group
            )

            # 步骤6: 保存结果|Step 6: Save results
            self.logger.info("\n[步骤6|Step 6] 保存结果|Saving results...")
            self._save_results(selected_samples)

            self.logger.info("\n"+"="*70)
            self.logger.info("样品选择完成|Sample selection completed!")
            self.logger.info("="*70)
            self.logger.info(f"\n最终结果|Final results:")
            self.logger.info(f"  选择样品数|Selected samples: {len(selected_samples)}")
            self.logger.info(f"  输出文件|Output files:")
            self.logger.info(f"    - {self.config.output_prefix}.txt")
            self.logger.info(f"    - {self.config.output_prefix}_report.txt")
            if self.config.generate_csv:
                self.logger.info(f"    - {self.config.output_prefix}.csv")
            if self.config.generate_visualization:
                self.logger.info(f"    - {self.config.output_prefix}_visualization.txt")

            return True

        except Exception as e:
            self.logger.error(f"错误|Error: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _save_results(self, selected_samples: List[Dict]):
        """保存结果|Save results

        Args:
            selected_samples: 选中的样品|Selected samples
        """
        # 保存样品列表|Save sample list
        output_file = self.config.output_prefix + '.txt'
        with open(output_file, 'w') as f:
            f.write("# 系统发育树样品选择结果|Phylogenetic Tree Sample Selection Results\n")
            f.write(f"# 策略|Strategy: 基于层级关系和PCA的智能选择|Hierarchy-based and PCA intelligent selection\n")
            f.write(f"# 层级深度|Hierarchy level: parent_{self.config.hierarchy_level}\n")
            f.write(f"# 总样品数|Total samples: {len(selected_samples)}\n")
            f.write(f"#\n")
            for sample in selected_samples:
                f.write(f"{sample['name']}\n")
        self.logger.info(f"  样品列表|Sample list: {output_file}")

        # 保存详细报告|Save detailed report
        if self.config.generate_report:
            report_file = self.config.output_prefix + '_report.txt'
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("="*70 + "\n")
                f.write("样品选择详细报告|Detailed Sample Selection Report\n")
                f.write("="*70 + "\n\n")
                f.write(f"选择策略|Selection Strategy:\n")
                f.write(f"  - 基于层级关系|Based on: Hierarchy relationships\n")
                f.write(f"  - PCA去重|PCA dedup: Yes (threshold={self.config.pca_dedup_threshold})\n")
                f.write(f"  - 层级深度|Hierarchy level: parent_{self.config.hierarchy_level}\n")
                f.write(f"\n")
                f.write(f"选择结果|Selection Results:\n")
                f.write(f"  总样品数|Total samples: {len(selected_samples)}\n")
                f.write(f"\n")
                f.write(f"样品列表|Sample List:\n")
                for i, sample in enumerate(selected_samples, 1):
                    f.write(f"  {i:3d}. {sample['name']}\n")
            self.logger.info(f"  详细报告|Detailed report: {report_file}")

        # 保存CSV文件|Save CSV file
        if self.config.generate_csv:
            csv_file = self.config.output_prefix + '.csv'
            with open(csv_file, 'w', encoding='utf-8') as f:
                f.write("Index,Sample\n")
                for i, sample in enumerate(selected_samples, 1):
                    f.write(f"{i},{sample['name']}\n")
            self.logger.info(f"  CSV文件|CSV file: {csv_file}")

        # 保存可视化|Save visualization
        if self.config.generate_visualization:
            viz_file = self.config.output_prefix + '_visualization.txt'
            with open(viz_file, 'w', encoding='utf-8') as f:
                f.write("样品索引|Sample Index\n")
                for i, sample in enumerate(selected_samples, 1):
                    f.write(f"{i}\n")
            self.logger.info(f"  可视化|Visualization: {viz_file}")


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='系统发育树样品选择工具（基于层级关系和PCA）|Phylogenetic Tree Sample Selector (Hierarchy-based and PCA)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 基本用法（必需参数）
  %(prog)s -x hierarchy.xlsx -p pca.txt -o selected_samples

  # 指定选择数量
  %(prog)s -x hierarchy.xlsx -p pca.txt -o selected_samples -n 100

  # 包含分组信息
  %(prog)s -x hierarchy.xlsx -p pca.txt -g groups.txt -o selected_samples -n 150

  # 自定义层级深度和PCA阈值
  %(prog)s -x hierarchy.xlsx -p pca.txt -o selected_samples --hierarchy-level 15 --pca-threshold 0.001
        '''
    )

    # 必需参数|Required parameters
    parser.add_argument('-x', '--hierarchy-file', required=True,
                        help='层级关系文件（必需）|Hierarchy file (required)')
    parser.add_argument('-p', '--pca-file', required=True,
                        help='PCA坐标文件（必需）|PCA coordinates file (required)')
    parser.add_argument('-o', '--output-prefix', required=True,
                        help='输出文件前缀|Output file prefix')

    # 可选参数|Optional parameters
    parser.add_argument('-n', '--n-samples',
                        type=int,
                        default=150,
                        help='选择样品总数|Total number of samples to select (default: 150)')

    parser.add_argument('-g', '--group-file',
                        help='样品分组表文件|Sample group table file')

    parser.add_argument('--hierarchy-level',
                        type=int,
                        default=10,
                        help='层级深度（用于分组）|Hierarchy level for grouping (default: 10)')

    parser.add_argument('--pca-threshold',
                        type=float,
                        default=0.0001,
                        help='PCA去重阈值|PCA dedup threshold (default: 0.0001)')

    parser.add_argument('--min-samples-per-group',
                        type=int,
                        default=1,
                        help='每组最小样品数|Minimum samples per group (default: 1)')

    parser.add_argument('--no-report',
                        action='store_true',
                        help='不生成详细报告|Do not generate detailed report')

    parser.add_argument('--no-csv',
                        action='store_true',
                        help='不生成CSV文件|Do not generate CSV file')

    parser.add_argument('--no-visualization',
                        action='store_true',
                        help='不生成可视化|Do not generate visualization')

    args = parser.parse_args()

    try:
        # 创建运行器|Create runner
        runner = PhyloSelectorRunner(
            hierarchy_file=args.hierarchy_file,
            pca_file=args.pca_file,
            output_prefix=args.output_prefix,
            n_samples=args.n_samples,
            group_file=args.group_file,
            hierarchy_level=args.hierarchy_level,
            pca_dedup_threshold=args.pca_threshold,
            min_samples_per_group=args.min_samples_per_group,
            generate_report=not args.no_report,
            generate_csv=not args.no_csv,
            generate_visualization=not args.no_visualization
        )

        # 运行|Run
        success = runner.run()

        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
