"""
NGenomeSyn Main Module
NGenomeSyn主模块
"""

import os
import sys
import time
import argparse
import logging
from .config import NGenomeSynConfig
from .environment import EnvironmentChecker
from .data_processing import GenomeProcessor
from .aligners import AlignerFactory, SyriProcessor
from .visualizer import NGenomeSynVisualizer


def setup_logger(output_dir: str) -> logging.Logger:
    """设置日志 | Setup logging"""
    logger = logging.getLogger("NGenomeSyn")
    logger.setLevel(logging.INFO)

    # 清除已有的handlers
    logger.handlers.clear()

    # 创建formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # 添加文件handler
    log_file = os.path.join(output_dir, "ngenomesyn_pipeline.log")
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # 添加控制台handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='NGenomeSynteny Analysis Pipeline - NGenomeSyn基因组共线性分析流程',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 基本用法
  %(prog)s -s samples.txt -o output_dir

  # 使用MUMmer比对器
  %(prog)s -s samples.txt -o output_dir --aligner mummer

  # 分析特定染色体
  %(prog)s -s samples.txt -o output_dir --chromosome "1,2,3"

  # 自定义线程数和最小长度
  %(prog)s -s samples.txt -o output_dir --threads 32 --min-length 10000
        '''
    )

    # 输入参数（互斥）
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-s', '--sample-map',
                            help='[FILE] 样本映射文件 | Sample mapping file (genome_file\\tgenome_name)')
    input_group.add_argument('-c', '--config',
                            help='[FILE] 配置文件 | Configuration file')

    # 输出参数
    parser.add_argument('-o', '--output', required=True,
                       help='[DIR] 输出目录 | Output directory')

    # 比对参数
    parser.add_argument('-a', '--aligner', default='minimap2',
                       choices=['minimap2', 'mummer'],
                       help='[STR] 比对器类型 (默认: minimap2) | Aligner type (default: minimap2)')
    parser.add_argument('-t', '--threads', type=int, default=16,
                       help='[INT] 线程数 (默认: 16) | Number of threads (default: 16)')
    parser.add_argument('--min-length', type=int, default=5000,
                       help='[INT] 最小比对长度 (默认: 5000) | Minimum alignment length (default: 5000)')

    # Minimap2参数
    parser.add_argument('--minimap-preset', default='asm5',
                       help='[STR] Minimap2预设模式 (默认: asm5) | Minimap2 preset (default: asm5)')

    # MUMmer参数
    parser.add_argument('--mummer-match-type', default='mumreference',
                       choices=['mum', 'mumreference', 'maxmatch'],
                       help='[STR] MUMmer匹配类型 (默认: mumreference) | MUMmer match type (default: mumreference)')
    parser.add_argument('--mummer-min-match', type=int, default=20,
                       help='[INT] MUMmer最小匹配长度 (默认: 20) | MUMmer min match (default: 20)')

    # 染色体过滤
    parser.add_argument('--chromosome', type=str,
                       help='[STR] 指定染色体 (如: "1,2,3" 或 "1-5") | Specify chromosomes')

    # 可视化参数
    parser.add_argument('--output-formats', nargs='+', default=['svg', 'png'],
                       choices=['svg', 'png'],
                       help='[STR] 输出格式 (默认: svg png) | Output formats (default: svg png)')
    parser.add_argument('--ngenomesyn-bin',
                       help='[FILE] NGenomeSyn二进制文件路径 | NGenomeSyn binary path')

    # SyRI参数
    parser.add_argument('--use-syri', action='store_true',
                       help='[FLAG] 使用SyRI进行结构变异分析 | Use SyRI for structural variation analysis')
    parser.add_argument('--syri-bin',
                       help='[FILE] SyRI二进制文件路径 | SyRI binary path')

    args = parser.parse_args()

    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)

    # 设置日志
    logger = setup_logger(args.output)

    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("NGenomeSynteny Analysis Pipeline")
        logger.info("=" * 60)

        # 初始化配置
        config = NGenomeSynConfig(
            sample_map=args.sample_map,
            config_file=args.config,
            output_dir=args.output,
            aligner=args.aligner,
            threads=args.threads,
            min_length=args.min_length,
            output_formats=args.output_formats,
            ngenomesyn_bin=args.ngenomesyn_bin,
            _chromosome_str=args.chromosome,
            minimap_preset=args.minimap_preset,
            mummer_match_type=args.mummer_match_type,
            mummer_min_match=args.mummer_min_match,
            use_syri=args.use_syri,
            syri_bin=args.syri_bin
        )

        config.validate()

        # 初始化组件
        env_checker = EnvironmentChecker(logger)
        processor = GenomeProcessor(logger)
        visualizer = NGenomeSynVisualizer(logger)

        # 检查环境
        logger.info("检查环境")
        env_results = env_checker.check_environment(config.aligner)
        if not env_results["NGenomeSyn"]:
            logger.error("NGenomeSyn未安装或不在PATH中")
            return 1
        if not env_results.get(config.aligner, False):
            logger.error(f"{config.aligner} 未安装或不在PATH中")
            return 1

        # 检查SyRI（如果使用）
        if config.use_syri:
            if not env_results.get("syri", False):
                logger.error("SyRI未安装或不在PATH中，但已启用--use-syri选项")
                return 1
            logger.info("✅ SyRI已启用")

        # 读取样本信息
        if config.sample_map:
            logger.info("读取样本映射文件")
            samples = processor.read_sample_map(config.sample_map)

            if len(samples) < 2:
                logger.error("至少需要2个基因组进行共线性分析")
                return 1

            logger.info(f"成功读取 {len(samples)} 个基因组")

            # 生成比对链（chain模式）
            links = []
            for i in range(len(samples) - 1):
                links.append({
                    "from": i + 1,
                    "to": i + 2,
                    "from_name": samples[i]["name"],
                    "to_name": samples[i + 1]["name"]
                })

            # 构建配置数据
            config_data = {
                "genomes": samples,
                "links": links,
                "alignment": {
                    "aligner": config.aligner,
                    "threads": config.threads,
                    "min_length": config.min_length,
                    "minimap_preset": config.minimap_preset,
                    "mummer_match_type": config.mummer_match_type,
                    "mummer_min_match": config.mummer_min_match
                },
                "visualization": {
                    "output_formats": config.output_formats
                }
            }
        else:
            logger.error("暂不支持从配置文件读取模式")
            return 1

        # 步骤1: 执行比对
        logger.info("")
        logger.info("步骤1: 执行基因组比对")
        logger.info("=" * 60)

        # 如果使用SyRI，初始化SyRI处理器
        if config.use_syri:
            syri_processor = SyriProcessor(logger, config.syri_bin, config.output_dir)

        aligner = AlignerFactory.create_aligner(
            config.aligner,
            logger,
            config.threads,
            config.output_dir
        )

        for link in config_data["links"]:
            genome1 = next(g for g in config_data["genomes"] if g["name"] == link["from_name"])
            genome2 = next(g for g in config_data["genomes"] if g["name"] == link["to_name"])

            logger.info(f"比对: {genome1['name']} vs {genome2['name']}")

            if config.use_syri:
                # 使用SyRI流程
                if config.aligner == "minimap2":
                    # 对于Minimap2 + SyRI，需要生成SAM格式
                    logger.info("🧬 使用Minimap2 + SyRI流程")
                    sam_file = aligner.align_for_syri(
                        genome1,
                        genome2,
                        selected_chromosomes=config.chromosomes,
                        preset=config.minimap_preset
                    )
                    # 运行SyRI分析
                    link_file = syri_processor.run_syri_from_minimap2(
                        genome1,
                        genome2,
                        sam_file,
                        min_length=config.min_length,
                        selected_chromosomes=config.chromosomes
                    )
                    logger.info(f"✅ SyRI分析完成: {link_file}")

                elif config.aligner == "mummer":
                    # 对于MUMmer + SyRI，使用coords文件
                    logger.info("🧬 使用MUMmer + SyRI流程")
                    coords_file = aligner.align(
                        genome1,
                        genome2,
                        selected_chromosomes=config.chromosomes,
                        match_type=config.mummer_match_type,
                        min_match=config.mummer_min_match
                    )
                    # 运行SyRI分析
                    link_file = syri_processor.run_syri_from_mummer_coords(
                        genome1,
                        genome2,
                        coords_file,
                        min_length=config.min_length,
                        selected_chromosomes=config.chromosomes
                    )
                    logger.info(f"✅ SyRI分析完成: {link_file}")
            else:
                # 标准流程（不使用SyRI）
                align_params = {}
                if config.aligner == "minimap2":
                    align_params["preset"] = config.minimap_preset
                    align_params["min_length"] = config.min_length
                elif config.aligner == "mummer":
                    align_params["match_type"] = config.mummer_match_type
                    align_params["min_match"] = config.mummer_min_match

                link_file = aligner.align(
                    genome1,
                    genome2,
                    selected_chromosomes=config.chromosomes,
                    **align_params
                )

                logger.info(f"比对完成: {link_file}")

        # 步骤2: 生成可视化
        logger.info("")
        logger.info("步骤2: 生成可视化")
        logger.info("=" * 60)

        # 生成NGenomeSyn配置
        ngenomesyn_config = visualizer.generate_ngenomesyn_config(
            config_data,
            config.output_dir,
            selected_chromosomes=config.chromosomes
        )

        # 运行NGenomeSyn
        output_files = visualizer.run_ngenomesyn(
            ngenomesyn_config,
            config.output_dir,
            config.output_formats
        )

        # 完成
        elapsed_time = time.time() - start_time

        logger.info("")
        logger.info("=" * 60)
        logger.info("分析总结")
        logger.info("=" * 60)
        logger.info(f"输入基因组: {len(samples)}")
        logger.info(f"比对数量: {len(links)}")
        logger.info(f"输出文件:")
        for file_path in output_files:
            logger.info(f"  - {file_path}")
        logger.info(f"运行时间: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        logger.info("=" * 60)
        logger.info("分析完成!")
        logger.info("=" * 60)

        return 0

    except KeyboardInterrupt:
        logger.warning("用户中断操作")
        return 130
    except ValueError as e:
        logger.error(f"配置错误: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Pipeline失败: {str(e)}", exc_info=True)
        return 1


class NGenomeSynAnalyzer:
    """NGenomeSyn分析器类 | NGenomeSyn Analyzer Class"""

    def __init__(self, **kwargs):
        """
        初始化分析器 | Initialize analyzer

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        self.config = NGenomeSynConfig(**kwargs)
        self.config.validate()

        # 初始化日志
        self.logger = setup_logger(self.config.output_dir)

        # 初始化组件
        self.env_checker = EnvironmentChecker(self.logger)
        self.processor = GenomeProcessor(self.logger)
        self.visualizer = NGenomeSynVisualizer(self.logger)

    def run(self):
        """
        运行分析 | Run analysis

        Returns:
            是否成功 | Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("NGenomeSynteny Analysis Pipeline")
            self.logger.info("=" * 60)

            # 检查环境
            env_results = self.env_checker.check_environment(self.config.aligner)
            if not env_results["NGenomeSyn"]:
                raise RuntimeError("NGenomeSyn未安装")
            if not env_results.get(self.config.aligner, False):
                raise RuntimeError(f"{self.config.aligner} 未安装")

            # 读取样本信息
            samples = self.processor.read_sample_map(self.config.sample_map)

            if len(samples) < 2:
                raise ValueError("至少需要2个基因组")

            # 生成比对链
            links = []
            for i in range(len(samples) - 1):
                links.append({
                    "from": i + 1,
                    "to": i + 2,
                    "from_name": samples[i]["name"],
                    "to_name": samples[i + 1]["name"]
                })

            # 构建配置数据
            config_data = {
                "genomes": samples,
                "links": links,
                "alignment": {
                    "aligner": self.config.aligner,
                    "threads": self.config.threads,
                    "min_length": self.config.min_length,
                    "minimap_preset": self.config.minimap_preset,
                    "mummer_match_type": self.config.mummer_match_type,
                    "mummer_min_match": self.config.mummer_min_match
                },
                "visualization": {
                    "output_formats": self.config.output_formats
                }
            }

            # 执行比对
            aligner = AlignerFactory.create_aligner(
                self.config.aligner,
                self.logger,
                self.config.threads,
                self.config.output_dir
            )

            for link in config_data["links"]:
                genome1 = next(g for g in config_data["genomes"] if g["name"] == link["from_name"])
                genome2 = next(g for g in config_data["genomes"] if g["name"] == link["to_name"])

                align_params = {}
                if self.config.aligner == "minimap2":
                    align_params["preset"] = self.config.minimap_preset
                    align_params["min_length"] = self.config.min_length
                elif self.config.aligner == "mummer":
                    align_params["match_type"] = self.config.mummer_match_type
                    align_params["min_match"] = self.config.mummer_min_match

                aligner.align(
                    genome1,
                    genome2,
                    selected_chromosomes=self.config.chromosomes,
                    **align_params
                )

            # 生成可视化
            ngenomesyn_config = self.visualizer.generate_ngenomesyn_config(
                config_data,
                self.config.output_dir,
                selected_chromosomes=self.config.chromosomes
            )

            output_files = self.visualizer.run_ngenomesyn(
                ngenomesyn_config,
                self.config.output_dir,
                self.config.output_formats
            )

            elapsed_time = time.time() - start_time

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("分析总结")
            self.logger.info("=" * 60)
            self.logger.info(f"运行时间: {elapsed_time:.2f} seconds")
            self.logger.info("分析完成!")

            return True

        except Exception as e:
            self.logger.error(f"分析失败: {str(e)}", exc_info=True)
            return False


if __name__ == "__main__":
    sys.exit(main())
