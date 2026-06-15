"""
VCF2PCA分析主程序模块|VCF2PCA Analysis Main Module
"""

import argparse
import sys
from .config import VCF2PCAConfig
from .utils import VCF2PCALogger
from .v2p_analyzer import V2PAnalyzer


class VCF2PCARunner:
    """VCF2PCA运行器|VCF2PCA Runner - 支持多backend|Support multiple backends"""

    def __init__(self, **kwargs):
        """
        初始化VCF2PCA运行器|Initialize VCF2PCA runner

        Args:
            **kwargs: 配置参数
        """
        # 初始化配置|Initialize configuration
        self.config = VCF2PCAConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = VCF2PCALogger(
            self.config.output_path / 'vcf2pca.log'
        )
        self.logger = self.logger_manager.get_logger()

        # 根据backend选择分析器|Select analyzer based on backend
        if self.config.backend == 'v2p':
            from .v2p_analyzer import V2PAnalyzer
            self.analyzer = V2PAnalyzer(self.config, self.logger)
        elif self.config.backend == 'plink':
            from .pca_analysis import PCAAnalyzer
            from .data_processing import VCFProcessor, QualityController
            from .visualization import PCAVisualizer
            from .results import SummaryGenerator
            from .utils import CommandRunner

            # PLINK backend需要更多组件|PLINK backend needs more components
            self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
            self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
            self.quality_controller = QualityController(
                self.config, self.logger, self.cmd_runner
            )
            self.pca_analyzer = PCAAnalyzer(self.config, self.logger, self.cmd_runner)
            self.visualizer = PCAVisualizer(self.config, self.logger)
            self.summary_generator = SummaryGenerator(self.config, self.logger)
            self.analyzer = None  # PLINK backend直接调用各组件

    def run_analysis(self):
        """运行PCA分析流程|Run PCA analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始VCF2PCA分析流程|Starting VCF2PCA analysis pipeline")
            self.logger.info("=" * 80)
            self.logger.info(f"Backend|后端: {self.config.backend}")
            self.logger.info(f"VCF文件|VCF file: {self.config.vcf_file}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"主成分数量|Number of PCs: {self.config.components}")

            if self.config.backend == 'v2p':
                self._run_v2p_backend()
            elif self.config.backend == 'plink':
                self._run_plink_backend()

            self.logger.info("=" * 80)
            self.logger.info("VCF2PCA分析完成！|VCF2PCA analysis completed!")
            self.logger.info("=" * 80)
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")

            return True

        except Exception as e:
            self.logger.error(f"分析流程失败|Analysis pipeline failed: {e}")
            sys.exit(1)

    def _run_v2p_backend(self):
        """运行VCF2PCACluster后端|Run VCF2PCACluster backend"""
        self.logger.info("使用VCF2PCACluster后端|Using VCF2PCACluster backend")

        if not self.analyzer.run_analysis():
            raise RuntimeError("VCF2PCACluster分析失败|VCF2PCACluster analysis failed")

    def _run_plink_backend(self):
        """运行PLINK后端|Run PLINK backend"""
        self.logger.info("使用PLINK后端|Using PLINK backend")

        # 检查可视化库|Check visualization libraries
        if self.config.plot:
            try:
                import matplotlib.pyplot as plt
                import seaborn as sns
                self.logger.info("可视化库可用|Visualization libraries available")
            except ImportError:
                self.logger.warning(
                    "可视化库不可用，将跳过图表生成|"
                    "Visualization libraries unavailable, skipping plots"
                )
                self.config.plot = False

        # 步骤1: VCF转PLINK|Step 1: VCF to PLINK conversion
        self.logger.info("=" * 80)
        self.logger.info("步骤1: VCF转PLINK格式|Step 1: VCF to PLINK conversion")
        self.logger.info("=" * 80)

        if not self.vcf_processor.vcf_to_plink():
            raise RuntimeError("VCF转换失败|VCF conversion failed")

        # 步骤2: 质量控制|Step 2: Quality control
        self.logger.info("=" * 80)
        self.logger.info("步骤2: 质量控制|Step 2: Quality control")
        self.logger.info("=" * 80)

        if not self.quality_controller.apply_quality_filters():
            raise RuntimeError("质量控制失败|Quality control failed")

        # 步骤3: PCA分析|Step 3: PCA analysis
        self.logger.info("=" * 80)
        self.logger.info("步骤3: PCA分析|Step 3: PCA analysis")
        self.logger.info("=" * 80)

        if not self.pca_analyzer.run_pca_analysis():
            raise RuntimeError("PCA分析失败|PCA analysis failed")

        # 步骤4: 可视化|Step 4: Visualization
        if self.config.plot:
            self.logger.info("=" * 80)
            self.logger.info("步骤4: 生成可视化图表|Step 4: Generate visualization plots")
            self.logger.info("=" * 80)

            self.visualizer.create_plots()

        # 生成总结报告|Generate summary report
        self.summary_generator.generate_summary_report()


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='VCF2PCA分析脚本|VCF2PCA Analysis Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 使用VCF2PCACluster后端（默认）|Use VCF2PCACluster backend (default)
  %(prog)s -i variants.vcf -o pca_results

  # 使用PLINK后端|Use PLINK backend
  %(prog)s -i variants.vcf -o pca_results --backend plink

  # 启用聚类分析|Enable clustering
  %(prog)s -i variants.vcf -o pca_results --cluster --cluster-method kmeans --cluster-k 3
        """
    )

    # 必需参数|Required arguments
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-i', '--input', dest='vcf_file',
        help='输入VCF文件路径|Input VCF file path'
    )
    input_group.add_argument('-v', '--vcf', dest='vcf_file', help=argparse.SUPPRESS)

    # 可选参数|Optional arguments
    parser.add_argument(
        '-o', '--output', dest='output_dir', default='./pca_output',
        help='输出目录|Output directory'
    )
    parser.add_argument(
        '-s', '--sample-info', dest='sample_info_file',
        help='样本信息文件|Sample information file'
    )

    # Backend选择|Backend selection
    parser.add_argument(
        '-b', '--backend', dest='backend', choices=['v2p', 'plink'], default='v2p',
        help='分析后端|Analysis backend: v2p (VCF2PCACluster, default) or plink'
    )

    # PCA参数|PCA parameters
    parser.add_argument(
        '-c', '--components', dest='components', type=int, default=10,
        help='主成分数量|Number of principal components'
    )

    # 质控参数（仅PLINK后端|PLINK backend only）
    parser.add_argument(
        '--maf', dest='maf', type=float, default=0.05,
        help='最小等位基因频率阈值|Minor allele frequency threshold (PLINK backend only)'
    )
    parser.add_argument(
        '--missing', dest='missing_rate', type=float, default=0.1,
        help='最大缺失率阈值|Maximum missing rate threshold (PLINK backend only)'
    )
    parser.add_argument(
        '--hwe', dest='hwe_pvalue', type=float, default=1e-6,
        help='Hardy-Weinberg平衡p值阈值|Hardy-Weinberg equilibrium p-value (PLINK backend only)'
    )
    parser.add_argument(
        '--skip-qc', dest='skip_qc', action='store_true',
        help='跳过质量控制过滤|Skip quality control filtering (PLINK backend only)'
    )

    # 聚类参数（仅V2P后端|V2P backend only）
    parser.add_argument(
        '--cluster', dest='cluster', action='store_true',
        help='启用聚类分析|Enable clustering analysis (V2P backend only)'
    )
    parser.add_argument(
        '--cluster-method', dest='cluster_method',
        choices=['kmeans', 'dbscan', 'em'], default='kmeans',
        help='聚类方法|Clustering method: kmeans, dbscan, em (V2P backend only)'
    )
    parser.add_argument(
        '--cluster-k', dest='cluster_k', type=int, default=3,
        help='K-means聚类数|Number of clusters for K-means (V2P backend only)'
    )

    # 可视化参数|Visualization parameters
    parser.add_argument(
        '-P', '--plot', dest='plot', action='store_true',
        help='生成PCA可视化图表|Generate PCA visualization plots'
    )

    # 线程数|Threads
    parser.add_argument(
        '-t', '--threads', dest='threads', type=int, default=12,
        help='线程数|Number of threads'
    )

    # 工具路径|Tool paths
    parser.add_argument(
        '--vcf2pca-path', dest='vcf2pca_path',
        help='VCF2PCACluster路径|VCF2PCACluster path'
    )
    parser.add_argument(
        '--plink-path', dest='plink_path', default='plink',
        help='PLINK软件路径|PLINK software path'
    )
    parser.add_argument(
        '--bcftools-path', dest='bcftools_path', default='bcftools',
        help='BCFtools软件路径|BCFtools software path'
    )

    args = parser.parse_args()

    # 构建配置参数字典|Build configuration kwargs
    kwargs = {
        'vcf_file': args.vcf_file,
        'output_dir': args.output_dir,
        'backend': args.backend,
        'components': args.components,
        'maf': args.maf,
        'missing_rate': args.missing_rate,
        'hwe_pvalue': args.hwe_pvalue,
        'skip_qc': args.skip_qc,
        'cluster': args.cluster,
        'cluster_method': args.cluster_method,
        'cluster_k': args.cluster_k,
        'plot': args.plot,
        'threads': args.threads,
        'sample_info_file': args.sample_info_file
    }

    # 添加工具路径参数|Add tool path parameters
    if args.vcf2pca_path:
        kwargs['vcf2pca_path'] = args.vcf2pca_path
    if args.plink_path:
        kwargs['plink_path'] = args.plink_path
    if args.bcftools_path:
        kwargs['bcftools_path'] = args.bcftools_path

    # 创建运行器并运行|Create runner and run
    runner = VCF2PCARunner(**kwargs)
    runner.run_analysis()


if __name__ == "__main__":
    main()
