"""
PanMAN主程序模块|PanMAN Main Program Module
"""

import argparse
import sys
from pathlib import Path

from .config import PanMANConfig
from .utils import PanMANLogger, PanMANCommandRunner
from .panman_builder import PanMANBuilder
from .panman_extractor import PanMANExtractor
from .pangraph_builder import PanGraphBuilder


class PanMANBuildRunner:
    """PanMAN构建运行器|PanMAN Build Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PanMANConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = PanMANLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = PanMANCommandRunner(
            self.logger,
            self.config.conda_env,
            self.config.backend,
            self.config.threads,
            self.config.conda_base,
            self.config.sif_image,
            self.config.singularity_path
        )

        # 初始化构建器|Initialize builder
        self.builder = PanMANBuilder(self.config, self.logger, self.cmd_runner)

    def run_analysis(self):
        """运行构建分析|Run build analysis"""
        try:
            # 检查Conda环境|Check Conda environment
            if self.config.backend == "conda":
                if not self.cmd_runner.check_conda_env():
                    self.logger.error(
                        f"Conda环境不存在|Conda environment does not exist: {self.config.conda_env}\n"
                        f"请使用以下命令创建环境|"
                        f"Please create environment using: conda create -n {self.config.conda_env} panman"
                    )
                    return False

            # 执行构建|Execute build
            output_file = self.builder.build()

            if output_file:
                self.logger.info("分析完成|Analysis completed successfully")
                return True
            else:
                self.logger.error("分析失败|Analysis failed")
                return False

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return False


class PanMANExtractRunner:
    """PanMAN提取运行器|PanMAN Extract Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PanMANConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = PanMANLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = PanMANCommandRunner(
            self.logger,
            self.config.conda_env,
            self.config.backend,
            self.config.threads,
            self.config.conda_base,
            self.config.sif_image,
            self.config.singularity_path
        )

        # 初始化提取器|Initialize extractor
        self.extractor = PanMANExtractor(self.config, self.logger, self.cmd_runner)

    def run_analysis(self):
        """运行提取分析|Run extraction analysis"""
        try:
            # 检查Conda环境|Check Conda environment
            if self.config.backend == "conda":
                if not self.cmd_runner.check_conda_env():
                    self.logger.error(
                        f"Conda环境不存在|Conda environment does not exist: {self.config.conda_env}\n"
                        f"请使用以下命令创建环境|"
                        f"Please create environment using: conda create -n {self.config.conda_env} panman"
                    )
                    return False

            # 执行提取|Execute extraction
            results = self.extractor.extract()

            if results:
                self.logger.info("分析完成|Analysis completed successfully")
                return True
            else:
                self.logger.error("分析失败|Analysis failed")
                return False

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return False


class PanGraphGenRunner:
    """PanGraph生成运行器|PanGraph Generation Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PanMANConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = PanMANLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化PanGraph构建器|Initialize PanGraph builder
        self.pangraph_builder = PanGraphBuilder(
            self.logger,
            pangraph_path=self.config.pangraph_path,
            pangraph_sif=self.config.pangraph_sif,
            singularity_path=self.config.singularity_path
        )

    def run_analysis(self):
        """运行PanGraph生成分析|Run PanGraph generation analysis"""
        try:
            # 验证FASTA文件|Validate FASTA file
            if not self.pangraph_builder.validate_fasta(self.config.fasta_file):
                return False

            # 执行PanGraph构建|Execute PanGraph build
            success, outputs = self.pangraph_builder.build_from_fasta(
                fasta_file=self.config.fasta_file,
                output_prefix=self.config.output_prefix,
                output_dir=self.config.output_dir,
                threads=self.config.threads
            )

            if success:
                self.logger.info("PanGraph生成成功|PanGraph generation completed successfully")
                self.logger.info(f"输出文件|Output files:")
                self.logger.info(f"  JSON|JSON: {outputs.get('json', 'N/A')}")
                self.logger.info(f"  Newick|Newick: {outputs.get('newick', 'N/A')}")
                self.logger.info(f"  日志|Log: {outputs.get('log', 'N/A')}")
                return True
            else:
                self.logger.error("PanGraph生成失败|PanGraph generation failed")
                return False

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return False


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='PanMAN泛基因组分析工具 (构建模式和提取模式)|PanMAN Pangenome Analysis Tool (Build and Extract Modes)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 模式选择|Mode selection
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--build', action='store_true',
                            help='构建模式|Build mode: build PanMAN from alignments')
    mode_group.add_argument('--extract', action='store_true',
                            help='提取模式|Extract mode: extract data from PanMAN')
    mode_group.add_argument('--generate-pangraph', action='store_true',
                            help='PanGraph生成模式|Generate PanGraph mode: generate PanGraph JSON from FASTA')

    # PanGraph生成模式参数|PanGraph generation mode parameters
    parser.add_argument('-i', '--fasta',
                       help='输入FASTA文件路径 (PanGraph生成模式必需)|Input FASTA file path (required for generate-pangraph mode)')

    # 构建模式参数|Build mode parameters
    parser.add_argument('-P', '--pangraph',
                       help='PanGraph JSON文件路径|PanGraph JSON file path')
    parser.add_argument('-G', '--gfa',
                       help='GFA文件路径|GFA file path')
    parser.add_argument('-M', '--msa',
                       help='MSA文件路径 (FASTA格式)|MSA file path (FASTA format)')
    parser.add_argument('-N', '--newick',
                       help='Newick树文件路径 (构建模式必需)|Newick tree file path (required for build mode)')

    # 提取模式参数|Extract mode parameters
    parser.add_argument('-I', '--input-panman',
                       help='PanMAN文件路径 (提取模式必需)|PanMAN file path (required for extract mode)')

    # 通用参数|Common parameters
    parser.add_argument('-o', '--output-prefix', default='output',
                       help='输出文件前缀|Output file prefix')
    parser.add_argument('--output-dir', default='./panman_output',
                       help='输出目录|Output directory')
    parser.add_argument('-r', '--reference',
                       help='参考序列名称 (VCF提取需要)|Reference sequence name (required for VCF extraction)')

    # 提取选项|Extraction options
    parser.add_argument('--summary', action='store_true',
                       help='提取摘要统计|Extract summary statistics')
    parser.add_argument('--extract-fasta', action='store_true',
                       help='提取FASTA序列|Extract FASTA sequences')
    parser.add_argument('--extract-msa', action='store_true',
                       help='提取MSA比对|Extract MSA alignment')
    parser.add_argument('--vcf', action='store_true',
                       help='提取VCF变异|Extract VCF variants')
    parser.add_argument('--extract-gfa', action='store_true',
                       help='提取GFA格式|Extract GFA format')
    parser.add_argument('--extract-newick', action='store_true',
                       help='提取Newick树|Extract Newick tree')
    parser.add_argument('--extended-newick', action='store_true',
                       help='提取扩展Newick格式|Extract extended Newick format')
    parser.add_argument('--maf', action='store_true',
                       help='提取MAF格式|Extract MAF format')
    parser.add_argument('--aa', action='store_true',
                       help='提取氨基酸翻译|Extract amino acid translations')
    parser.add_argument('--subnet', action='store_true',
                       help='提取子网络|Extract subnet')
    parser.add_argument('--annotate', action='store_true',
                       help='注释节点|Annotate nodes')
    parser.add_argument('--reroot', action='store_true',
                       help='重新扎根树|Reroot tree')
    parser.add_argument('--create-network', action='store_true',
                       help='创建网络|Create network')
    parser.add_argument('--print-mutations', action='store_true',
                       help='打印突变信息|Print mutations')

    # 高级选项|Advanced options
    parser.add_argument('--backend', choices=['conda', 'docker', 'singularity'], default='conda',
                       help='后端选择|Backend selection (default: conda)')
    parser.add_argument('--conda-env', default='panman_v.0.1.4',
                       help='Conda环境名称|Conda environment name (default: panman_v.0.1.4)')
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')
    parser.add_argument('--pangraph-path',
                       help='PanGraph可执行文件路径|PanGraph executable path')
    parser.add_argument('--pangraph-sif',
                       help='PanGraph Singularity SIF镜像路径|PanGraph Singularity SIF image path')
    parser.add_argument('--sif-image',
                       help='Singularity SIF镜像路径|Singularity SIF image path')
    parser.add_argument('--singularity-path',
                       help='Singularity可执行文件路径|Singularity executable path')
    parser.add_argument('--acr', choices=['fitch', 'mppa'], default='fitch',
                       help='ACR方法|ACR method (fitch/mppa)')
    parser.add_argument('--tree-id',
                       help='树ID (VCF提取需要)|Tree ID (required for VCF extraction)')
    parser.add_argument('--input-file',
                       help='输入文件路径 (用于subnet/annotate/create-network)|Input file path (for subnet/annotate/create-network)')
    parser.add_argument('--range-index',
                       help='范围查询index参数|Range query index parameter')
    parser.add_argument('--range-start', type=int,
                       help='范围查询起始坐标|Range query start coordinate')
    parser.add_argument('--range-end', type=int,
                       help='范围查询结束坐标|Range query end coordinate')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 确定模式|Determine mode
        if args.build:
            mode = "build"
        elif args.generate_pangraph:
            mode = "generate_pangraph"
        else:
            mode = "extract"

        # 准备配置参数|Prepare configuration parameters
        config_params = {
            'mode': mode,
            'output_prefix': args.output_prefix,
            'output_dir': args.output_dir,
            'backend': args.backend,
            'conda_env': args.conda_env,
            'threads': args.threads,
            'sif_image': args.sif_image,
            'singularity_path': args.singularity_path,
            'pangraph_sif': args.pangraph_sif
        }

        # 添加PanGraph生成模式参数|Add pangraph generation mode parameters
        if mode == "generate_pangraph":
            if args.fasta:
                config_params['fasta_file'] = args.fasta
            if args.pangraph_path:
                config_params['pangraph_path'] = args.pangraph_path

        # 添加构建模式参数|Add build mode parameters
        if mode == "build":
            if args.pangraph:
                config_params['pangraph_file'] = args.pangraph
            if args.gfa:
                config_params['gfa_file'] = args.gfa
            if args.msa:
                config_params['msa_file'] = args.msa
            if args.newick:
                config_params['newick_file'] = args.newick

        # 添加提取模式参数|Add extract mode parameters
        if mode == "extract":
            if args.input_panman:
                config_params['panman_file'] = args.input_panman
            if args.reference:
                config_params['reference'] = args.reference

            config_params['extract_summary'] = args.summary
            config_params['extract_fasta'] = args.extract_fasta
            config_params['extract_msa'] = args.extract_msa
            config_params['extract_vcf'] = args.vcf
            config_params['extract_gfa'] = args.extract_gfa
            config_params['extract_newick'] = args.extract_newick
            config_params['extract_extended_newick'] = args.extended_newick
            config_params['extract_maf'] = args.maf
            config_params['extract_aa'] = args.aa
            config_params['extract_subnet'] = args.subnet
            config_params['extract_annotate'] = args.annotate
            config_params['extract_reroot'] = args.reroot
            config_params['extract_create_network'] = args.create_network
            config_params['extract_print_mutations'] = args.print_mutations

            # 高级参数|Advanced parameters
            if args.tree_id:
                config_params['tree_id'] = args.tree_id
            if args.input_file:
                config_params['input_file'] = args.input_file
            config_params['acr_method'] = args.acr
            if args.range_index:
                config_params['range_query_index'] = args.range_index
            if args.range_start is not None:
                config_params['range_start'] = args.range_start
            if args.range_end is not None:
                config_params['range_end'] = args.range_end

        # 根据模式创建运行器|Create runner based on mode
        if mode == "build":
            runner = PanMANBuildRunner(**config_params)
        elif mode == "generate_pangraph":
            runner = PanGraphGenRunner(**config_params)
        else:
            runner = PanMANExtractRunner(**config_params)

        # 运行分析|Run analysis
        success = runner.run_analysis()

        if success:
            sys.exit(0)
        else:
            sys.exit(1)

    except ValueError as e:
        print(f"配置错误|Configuration Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"程序错误|Program Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
