"""
InterProScan注释主程序模块|InterProScan Annotation Main Module
"""

import argparse
import sys
import os
from pathlib import Path
from .config import InterProScanConfig
from .utils import InterProScanLogger, CommandRunner, count_sequences, format_number
from .parser import InterProScanParser
from .formatter import InterProScanFormatter


class InterProScanAnnotator:
    """InterProScan注释主类|Main InterProScan Annotator Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = InterProScanConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        output_dir = Path(os.path.dirname(self.config.output_prefix))
        self.logger_manager = InterProScanLogger(output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器，设置Python环境变量|Initialize command runner with Python env vars
        env_vars = {}
        if hasattr(self.config, 'python_path') and self.config.python_path:
            # 设置PATH环境变量，让兼容的Python优先|Set PATH to prioritize compatible Python
            python_dir = os.path.dirname(self.config.python_path)
            env_vars['PATH'] = f"{python_dir}:{os.environ.get('PATH', '')}"
            env_vars['PYTHON_PATH'] = self.config.python_path
            self.logger.info(f"使用Python解释器|Using Python interpreter: {self.config.python_path}")

        self.cmd_runner = CommandRunner(self.logger, working_dir=output_dir, env_vars=env_vars)

    def build_command(self) -> str:
        """构建InterProScan命令|Build InterProScan command"""
        cmd_parts = [self.config.interproscan_path]

        # 输入文件|Input file
        cmd_parts.extend(['-i', self.config.input_file])

        # 输出格式（支持多个格式）|Output format (support multiple formats)
        formats = [f.strip().upper() for f in self.config.output_format.split(',')]
        for fmt in formats:
            cmd_parts.extend(['-f', fmt])

        # 输出前缀|Output prefix
        cmd_parts.extend(['-b', self.config.output_prefix])

        # 线程数|Thread count
        cmd_parts.extend(['-cpu', str(self.config.threads)])

        # 禁用预计算查找服务|Disable precalc lookup service
        if self.config.disable_precalc:
            cmd_parts.append('-dp')

        # GO术语|GO terms
        if self.config.goterms:
            cmd_parts.append('--goterms')

        # Pathway信息|Pathway information
        if self.config.pathways:
            cmd_parts.append('--pathways')

        # 指定应用|Specify applications
        if self.config.applications:
            cmd_parts.extend(['-appl', self.config.applications])

        # 临时目录|Temporary directory
        if self.config.temp_dir:
            cmd_parts.extend(['-T', self.config.temp_dir])

        return ' '.join(cmd_parts)

    def run_analysis(self):
        """运行InterProScan分析|Run InterProScan analysis"""
        self.logger.info("开始InterProScan注释流程|Starting InterProScan annotation pipeline")

        # 统计序列数量|Count sequences
        seq_count = count_sequences(self.config.input_file)
        self.logger.info(f"输入序列数量|Input sequence count: {format_number(seq_count)}")

        # 构建命令|Build command
        cmd = self.build_command()
        self.logger.info(f"InterProScan命令|InterProScan command: {cmd}")

        # 执行命令|Execute command
        success = self.cmd_runner.run(
            cmd,
            description="InterProScan蛋白质注释|InterProScan protein annotation"
        )

        if success:
            self.logger.info(f"InterProScan注释完成|InterProScan annotation completed")

            # 显示所有输出文件|Show all output files
            output_files = self.config.get_output_files()
            for fmt, filepath in output_files.items():
                if Path(filepath).exists():
                    self.logger.info(f"输出文件 ({fmt})|Output file ({fmt}): {filepath}")

            # 解析和格式化结果|Parse and format results
            if self.config.generate_report:
                self._process_results()

            return True
        else:
            self.logger.error("InterProScan注释失败|InterProScan annotation failed")
            return False

    def _process_results(self):
        """处理和格式化结果|Process and format results"""
        self.logger.info("开始结果解析和格式化|Starting result parsing and formatting")

        # 获取输出文件|Get output files
        output_files = self.config.get_output_files()

        # 优先使用TSV解析（更可靠，包含完整GO/Pathway信息）|Prefer TSV parsing (more reliable, contains complete GO/Pathway info)
        tsv_file = output_files.get('TSV')
        xml_file = output_files.get('XML')
        json_file = output_files.get('JSON')

        # 创建解析器|Create parser
        parser = InterProScanParser(
            self.logger,
            parse_goterms=self.config.goterms,
            parse_pathways=self.config.pathways,
            go_database_path=self.config.go_database_path
        )

        # 解析结果（优先TSV）|Parse results (prefer TSV)
        parsed = False
        if tsv_file and Path(tsv_file).exists():
            self.logger.info("使用TSV格式解析|Parsing using TSV format")
            parser.parse_tsv(tsv_file)
            parsed = True
        elif xml_file and Path(xml_file).exists():
            self.logger.info("使用XML格式解析|Parsing using XML format")
            parser.parse_xml(xml_file)
            parsed = True
        elif json_file and Path(json_file).exists():
            self.logger.info("使用JSON格式解析|Parsing using JSON format")
            parser.parse_json(json_file)
            parsed = True

        if not parsed or not parser.matches:
            self.logger.warning("未找到可解析的输出文件或解析失败|No parsable output file found or parsing failed")
            return

        # 生成格式化报告|Generate formatted report
        formatter = InterProScanFormatter(self.logger)
        formatter.generate_report(
            parser=parser,
            output_prefix=self.config.output_prefix,
            format_type=self.config.report_format,
            include_summary=self.config.include_summary
        )

        self.logger.info("结果解析和格式化完成|Result parsing and formatting completed")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='InterProScan蛋白质功能注释自动化脚本|InterProScan Protein Function Annotation Automation Script',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入蛋白质FASTA文件路径|Input protein FASTA file path')
    parser.add_argument('-o', '--output-prefix', required=True,
                       help='输出文件前缀(不含扩展名)|Output file prefix (without extension)')

    # 可选参数|Optional arguments
    parser.add_argument('--interproscan-path',
                       default='~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh',
                       help='InterProScan软件路径|InterProScan software path')
    parser.add_argument('-f', '--format',
                       default='TSV,XML',
                       help='输出格式，支持多个格式逗号分隔 (例如: TSV,XML,JSON)|Output format, support multiple formats comma-separated (e.g., TSV,XML,JSON)')
    parser.add_argument('-t', '--threads',
                       type=int, default=24,
                       help='线程数|Number of threads')

    # 处理选项|Processing options
    parser.add_argument('--disable-precalc', action='store_true', default=True,
                       help='禁用预计算查找服务(默认启用，解决网络问题)|'
                            'Disable precalc lookup service (enabled by default, solves network issues)')
    parser.add_argument('--enable-precalc', action='store_true',
                       help='启用预计算查找服务(需要网络连接)|'
                            'Enable precalc lookup service (requires network connection)')
    parser.add_argument('--goterms', action='store_true',
                       help='获取GO术语(默认不获取)|Get GO terms (disabled by default)')
    parser.add_argument('--no-goterms', action='store_true', default=True,
                       help='不获取GO术语(默认)|Do not get GO terms (default)')
    parser.add_argument('--pathways', action='store_true',
                       help='获取Pathway信息(默认不获取)|Get pathway information (disabled by default)')
    parser.add_argument('--no-pathways', action='store_true', default=True,
                       help='不获取Pathway信息(默认)|Do not get pathway information (default)')

    # 应用选项|Application options
    parser.add_argument('-appl', '--applications',
                       help='指定运行的应用，逗号分隔(例如: Pfam,SMART,Gene3D)|'
                            'Specify applications to run, comma-separated (e.g., Pfam,SMART,Gene3D)')

    # 临时目录|Temporary directory
    parser.add_argument('--temp-dir',
                       help='临时目录路径|Temporary directory path')

    # 结果整理选项|Result formatting options
    parser.add_argument('--no-report', action='store_true',
                       help='不生成整理后的报告|Do not generate formatted report')
    parser.add_argument('--report-format',
                       default='both',
                       choices=['excel', 'csv', 'both'],
                       help='报告格式 (默认: both，同时生成Excel和CSV)|Report format (default: both, generate both Excel and CSV)')
    parser.add_argument('--no-summary', action='store_true',
                       help='不包含汇总统计表|Do not include summary statistics')

    # GO数据库选项|GO database options
    parser.add_argument('--go-database',
                       help='GO术语JSON数据库文件路径（用于填充GO名称和ontology）|GO term JSON database file path (for filling GO names and ontologies)')

    # Python解释器选项|Python interpreter options
    parser.add_argument('--python-path',
                       help='Python解释器路径（用于指定兼容的Python 3.8-3.11版本）|Python interpreter path (for specifying compatible Python 3.8-3.11)')

    args = parser.parse_args()

    # 处理互斥选项|Handle mutually exclusive options
    disable_precalc = args.enable_precalc is False  # 如果启用precalc，则不禁用
    goterms = args.no_goterms is False  # 如果不获取goterms，则不获取
    pathways = args.no_pathways is False  # 如果不获取pathways，则不获取
    generate_report = not args.no_report  # 如果不生成报告，则不生成
    include_summary = not args.no_summary  # 如果不包含汇总，则不包含

    # 创建注释器并运行|Create annotator and run
    annotator = InterProScanAnnotator(
        input_file=args.input,
        output_prefix=args.output_prefix,
        interproscan_path=args.interproscan_path,
        output_format=args.format,
        threads=args.threads,
        disable_precalc=disable_precalc,
        goterms=goterms,
        pathways=pathways,
        applications=args.applications,
        temp_dir=args.temp_dir,
        generate_report=generate_report,
        report_format=args.report_format,
        include_summary=include_summary,
        go_database_path=args.go_database,
        python_path=args.python_path
    )

    success = annotator.run_analysis()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
