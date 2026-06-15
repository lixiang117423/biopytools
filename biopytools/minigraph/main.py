"""
Minigraph主程序模块|Minigraph Main Module
"""

import argparse
import sys
import os
from .config import (MinigraphBuildConfig, MinigraphCallConfig,
                     MinigraphBubbleConfig, MinigraphMapConfig)
from .utils import MinigraphLogger, CommandRunner, check_dependencies
from .graph_builder import (MinigraphGraphBuilder, MinigraphSVCaller,
                             MinigraphBubbleExtractor, MinigraphMapper)


class MinigraphRunner:
    """Minigraph运行器|Minigraph Runner"""

    def __init__(self, command: str, **kwargs):
        """初始化|Initialize"""
        self.command = command
        self.config = None
        self.logger = None
        self.cmd_runner = None

        # 根据命令类型初始化配置|Initialize config based on command type
        if command == 'build':
            self.config = MinigraphBuildConfig(**kwargs)
            self.config.validate()
            self._setup_logging(self.config.output_gfa)

        elif command == 'call':
            self.config = MinigraphCallConfig(**kwargs)
            self.config.validate()
            self._setup_logging(self.config.output_dir)

        elif command == 'bubble':
            self.config = MinigraphBubbleConfig(**kwargs)
            self.config.validate()
            self._setup_logging(self.config.output_bed)

        elif command == 'map':
            self.config = MinigraphMapConfig(**kwargs)
            self.config.validate()
            self._setup_logging(self.config.output_gaf)

        # 初始化命令执行器|Initialize command runner
        if self.logger:
            self.cmd_runner = CommandRunner(self.logger, os.path.dirname(self._get_output_path()))

        # 检查依赖|Check dependencies
        if self.logger:
            if not check_dependencies(self.config, self.logger):
                self.logger.error("依赖检查失败|Dependency check failed")
                sys.exit(1)

    def _setup_logging(self, output_path: str):
        """设置日志|Setup logging"""
        output_dir = os.path.dirname(output_path)
        if not output_dir:
            output_dir = '.'
        self.logger_manager = MinigraphLogger(
            output_dir=output_dir,
            log_file='minigraph_pipeline.log',
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

    def _get_output_path(self) -> str:
        """获取输出路径|Get output path"""
        if self.command == 'build':
            return self.config.output_gfa
        elif self.command == 'call':
            return self.config.output_dir
        elif self.command == 'bubble':
            return self.config.output_bed
        elif self.command == 'map':
            return self.config.output_gaf
        return '.'

    def run(self):
        """运行分析|Run analysis"""
        try:
            if self.command == 'build':
                builder = MinigraphGraphBuilder(self.config, self.logger, self.cmd_runner)
                success = builder.build_graph()

            elif self.command == 'call':
                caller = MinigraphSVCaller(self.config, self.logger, self.cmd_runner)
                success = caller.call_samples()

            elif self.command == 'bubble':
                extractor = MinigraphBubbleExtractor(self.config, self.logger, self.cmd_runner)
                success = extractor.extract_bubbles()

            elif self.command == 'map':
                mapper = MinigraphMapper(self.config, self.logger, self.cmd_runner)
                success = mapper.map_sequences()

            else:
                self.logger.error(f"未知命令|Unknown command: {self.command}")
                return False

            return success

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return False


def main():
    """主函数|Main function"""
    # 创建顶层解析器|Create top-level parser
    parser = argparse.ArgumentParser(
        description='Minigraph泛基因组图构建和分析工具|Minigraph Pangenome Graph Construction and Analysis Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='子命令|Sub-commands')

    # ========== build 子命令 ==========
    parser_build = subparsers.add_parser('build', help='构建泛基因组图|Build pangenome graph')

    parser_build.add_argument('--ref', required=True,
                             help='参考基因组FASTA文件|Reference genome FASTA file')
    parser_build.add_argument('--samples', required=True, nargs='+',
                             help='样本基因组FASTA文件列表|Sample genome FASTA file list')
    parser_build.add_argument('-o', '--output-gfa',
                             default='./pangenome.gfa',
                             help='输出GFA文件路径|Output GFA file path')
    parser_build.add_argument('--preset',
                             default='ggs',
                             choices=['g', 'gs', 'ggs'],
                             help='图构建预设|Graph building preset')
    parser_build.add_argument('--min-identity',
                             type=float, default=0.9,
                             help='最小序列相似度|Minimum sequence identity')
    parser_build.add_argument('--min-aln-len',
                             type=int, default=100000,
                             help='最小比对长度|Minimum alignment length')
    parser_build.add_argument('--max-gap',
                             type=int, default=1000000,
                             help='最大gap大小|Maximum gap size')
    parser_build.add_argument('-t', '--threads',
                             type=int, default=16,
                             help='线程数|Number of threads')
    parser_build.add_argument('--batch-size',
                             type=int,
                             help='批处理大小(MB)|Batch size (MB)')
    parser_build.add_argument('--minigraph-path',
                             default='minigraph',
                             help='minigraph工具路径|minigraph tool path')
    parser_build.add_argument('--gfatools-path',
                             default='gfatools',
                             help='gfatools工具路径|gfatools tool path')
    parser_build.add_argument('--keep-intermediate',
                             action='store_true',
                             help='保留中间文件|Keep intermediate files')
    parser_build.add_argument('--append-mode',
                             action='store_true',
                             help='追加模式|Append mode')

    # ========== call 子命令 ==========
    parser_call = subparsers.add_parser('call', help='调用SV|Call structural variants')

    parser_call.add_argument('--graph-gfa', required=True,
                            help='泛基因组图GFA文件|Pangenome graph GFA file')
    parser_call.add_argument('--samples', required=True, nargs='+',
                            help='样本基因组FASTA文件列表|Sample genome FASTA file list')
    parser_call.add_argument('-o', '--output-dir',
                            default='./minigraph_call',
                            help='输出目录|Output directory')
    parser_call.add_argument('--preset',
                            default='asm',
                            choices=['asm'],
                            help='SV调用预设|SV calling preset')
    parser_call.add_argument('-t', '--threads',
                            type=int, default=16,
                            help='线程数|Number of threads')
    parser_call.add_argument('--minigraph-path',
                            default='minigraph',
                            help='minigraph工具路径|minigraph tool path')

    # ========== bubble 子命令 ==========
    parser_bubble = subparsers.add_parser('bubble', help='提取SV bubbles|Extract SV bubbles')

    parser_bubble.add_argument('--graph-gfa', required=True,
                              help='泛基因组图GFA文件|Pangenome graph GFA file')
    parser_bubble.add_argument('-o', '--output-bed',
                              default='./sv_bubbles.bed',
                              help='输出BED文件路径|Output BED file path')
    parser_bubble.add_argument('--gfatools-path',
                              default='gfatools',
                              help='gfatools工具路径|gfatools tool path')

    # ========== map 子命令 ==========
    parser_map = subparsers.add_parser('map', help='序列映射|Map sequences to graph')

    parser_map.add_argument('--graph-gfa', required=True,
                           help='泛基因组图GFA文件|Pangenome graph GFA file')
    parser_map.add_argument('--queries', required=True, nargs='+',
                           help='查询序列FASTA文件列表|Query sequence FASTA file list')
    parser_map.add_argument('-o', '--output-gaf',
                           default='./mapping.gaf',
                           help='输出GAF文件路径|Output GAF file path')
    parser_map.add_argument('--preset',
                           default='lr',
                           choices=['sr', 'lr', 'map-pb', 'map-ont', 'asm'],
                           help='映射预设|Mapping preset')
    parser_map.add_argument('--max-intron-len',
                           type=int,
                           help='最大内含子长度|Maximum intron length')
    parser_map.add_argument('-t', '--threads',
                           type=int, default=16,
                           help='线程数|Number of threads')
    parser_map.add_argument('--batch-size',
                           type=int,
                           help='批处理大小(MB)|Batch size (MB)')
    parser_map.add_argument('--minigraph-path',
                           default='minigraph',
                           help='minigraph工具路径|minigraph tool path')

    # 解析参数|Parse arguments
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        # 根据命令类型创建运行器|Create runner based on command type
        if args.command == 'build':
            runner = MinigraphRunner(
                command='build',
                ref_fasta=args.ref,
                sample_fastas=args.samples,
                output_gfa=args.output_gfa,
                preset=args.preset,
                min_identity=args.min_identity,
                min_aln_len=args.min_aln_len,
                max_gap=args.max_gap,
                threads=args.threads,
                batch_size=args.batch_size,
                minigraph_path=args.minigraph_path,
                gfatools_path=args.gfatools_path,
                keep_intermediate=args.keep_intermediate,
                append_mode=args.append_mode
            )

        elif args.command == 'call':
            runner = MinigraphRunner(
                command='call',
                graph_gfa=args.graph_gfa,
                sample_fastas=args.samples,
                output_dir=args.output_dir,
                preset=args.preset,
                threads=args.threads,
                minigraph_path=args.minigraph_path
            )

        elif args.command == 'bubble':
            runner = MinigraphRunner(
                command='bubble',
                graph_gfa=args.graph_gfa,
                output_bed=args.output_bed,
                gfatools_path=args.gfatools_path
            )

        elif args.command == 'map':
            runner = MinigraphRunner(
                command='map',
                graph_gfa=args.graph_gfa,
                query_fastas=args.queries,
                output_gaf=args.output_gaf,
                preset=args.preset,
                max_intron_len=args.max_intron_len,
                threads=args.threads,
                batch_size=args.batch_size,
                minigraph_path=args.minigraph_path
            )

        else:
            print(f"未知命令|Unknown command: {args.command}")
            sys.exit(1)

        # 运行分析|Run analysis
        success = runner.run()

        if success:
            sys.exit(0)
        else:
            sys.exit(1)

    except Exception as e:
        print(f"错误|Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
