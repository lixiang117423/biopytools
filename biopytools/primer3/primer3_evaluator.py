"""
Primer3引物设计评估器|Primer3 Primer Design Evaluator
"""

import os
import sys
import tempfile
from pathlib import Path
from .config import Primer3Config
from .utils import Primer3Logger, CommandRunner, build_conda_command
from .parser import FastaParser, Primer3InputGenerator, Primer3OutputParser, ResultsFormatter


class Primer3Evaluator:
    """Primer3引物设计评估器|Primer3 Primer Design Evaluator"""

    def __init__(self, **kwargs):
        """
        初始化评估器|Initialize evaluator

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = Primer3Config(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = Primer3Logger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        self.logger.info(f"Primer3引物设计评估器初始化完成|Primer3 evaluator initialized")
        self.logger.info(f"输入文件|Input file: {self.config.input_fasta}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

    def run_design(self):
        """
        运行引物设计流程|Run primer design pipeline

        Returns:
            bool: 是否成功|Success or not
        """
        try:
            self.logger.info("开始Primer3引物设计流程|Starting Primer3 primer design pipeline")

            # 步骤1: 解析FASTA文件|Step 1: Parse FASTA file
            self.logger.info("步骤1: 解析FASTA文件|Step 1: Parsing FASTA file")
            sequences = self._parse_fasta()
            self.logger.info(f"成功解析|Successfully parsed {len(sequences)} 条序列|sequences")

            # 步骤2: 生成Primer3输入文件|Step 2: Generate Primer3 input file
            self.logger.info("步骤2: 生成Primer3输入格式|Step 2: Generating Primer3 input format")
            primer3_input = self._generate_primer3_input(sequences)

            # 步骤3: 运行Primer3|Step 3: Run Primer3
            self.logger.info("步骤3: 运行Primer3设计引物|Step 3: Running Primer3 to design primers")
            primer3_output = self._run_primer3(primer3_input)

            # 步骤4: 解析Primer3输出|Step 4: Parse Primer3 output
            self.logger.info("步骤4: 解析Primer3输出结果|Step 4: Parsing Primer3 output")
            parsed_results = self._parse_primer3_output(primer3_output)

            # 步骤5: 格式化并保存结果|Step 5: Format and save results
            self.logger.info("步骤5: 格式化并保存结果|Step 5: Formatting and saving results")
            self._format_and_save_results(parsed_results, sequences)

            self.logger.info("Primer3引物设计流程完成|Primer3 primer design pipeline completed")
            return True

        except Exception as e:
            self.logger.error(f"引物设计流程失败|Primer design pipeline failed: {e}")
            return False

    def _parse_fasta(self):
        """解析FASTA文件|Parse FASTA file"""
        parser = FastaParser(self.config.input_fasta)
        return parser.parse()

    def _generate_primer3_input(self, sequences):
        """生成Primer3输入格式|Generate Primer3 input format"""
        settings = self.config.get_primer3_settings()
        generator = Primer3InputGenerator(settings)
        return generator.generate_batch_input(sequences, self.config)

    def _run_primer3(self, primer3_input: str) -> str:
        """
        运行Primer3|Run Primer3

        Args:
            primer3_input: Primer3输入格式字符串|Primer3 input format string

        Returns:
            Primer3输出字符串|Primer3 output string
        """
        # 创建临时输入文件|Create temporary input file
        tmp_root = os.path.join(str(self.config.output_path), 'tmp')
        os.makedirs(tmp_root, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, dir=tmp_root) as tmp_input:
            tmp_input.write(primer3_input)
            tmp_input_path = tmp_input.name

        try:
            # 构建Primer3命令|Build Primer3 command
            cmd = build_conda_command(
                self.config.primer3_core_path,
                [tmp_input_path]
            )

            # 运行Primer3|Run Primer3
            result = self.cmd_runner.run(cmd, "Primer3引物设计|Primer3 primer design")

            if not result:
                raise RuntimeError("Primer3运行失败|Primer3 execution failed")

            # 读取输出文件（Primer3会输出到stdout）|Read output (Primer3 outputs to stdout)
            # 需要重新运行以捕获输出|Need to re-run to capture output
            import subprocess
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                cwd=self.config.output_path
            )

            if result.returncode != 0:
                raise RuntimeError(f"Primer3运行失败|Primer3 execution failed: {result.stderr}")

            return result.stdout

        finally:
            # 清理临时文件|Clean up temporary file
            import os
            if os.path.exists(tmp_input_path):
                os.unlink(tmp_input_path)

    def _parse_primer3_output(self, primer3_output: str):
        """解析Primer3输出|Parse Primer3 output"""
        parser = Primer3OutputParser()
        return parser.parse(primer3_output)

    def _format_and_save_results(self, parsed_results, sequences):
        """格式化并保存结果|Format and save results"""
        # 将sequences转换为字典以便查找|Convert sequences to dict for lookup
        seq_dict = {seq_id: sequence for seq_id, sequence in sequences}

        # 转换为DataFrame|Convert to DataFrame
        formatter = ResultsFormatter()
        df = formatter.to_dataframe(parsed_results, seq_dict, self.config.output_header_lang)

        # 保存结果|Save results
        output_file = self.config.output_path / f"primers_result.{self.config.output_format}"
        formatter.save_table(df, str(output_file), self.config.output_format)

        self.logger.info(f"结果已保存到|Results saved to: {output_file}")
        self.logger.info(f"共设计|Total primers designed: {len(df)} 对引物|primer pairs")

        return df
