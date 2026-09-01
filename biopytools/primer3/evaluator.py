"""
Primer3引物设计评估器|Primer3 Primer Design Evaluator
"""

import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from ..common.conda_runner import CommandRunner, build_conda_command
from . import __version__
from .config import Primer3Config
from .parser import FastaParser, Primer3InputGenerator, Primer3OutputParser, ResultsFormatter
from .utils import Primer3Logger, format_number


class Primer3Evaluator:
    """Primer3引物设计评估器|Primer3 Primer Design Evaluator"""

    def __init__(self, **kwargs):
        """
        初始化评估器|Initialize evaluator

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        self.config = Primer3Config(**kwargs)
        self.config.validate()

        # 输出目录结构(§12): 01_primer_design/ 00_pipeline_info/ 99_logs/ tmp/
        # |Output layout (§12): design files, metadata, logs, tmp
        self.output_path = self.config.output_path
        self.design_dir = self.output_path / '01_primer_design'
        self.info_dir = self.output_path / '00_pipeline_info'
        self.tmp_dir = self.output_path / 'tmp'
        for directory in (self.design_dir, self.info_dir, self.tmp_dir):
            directory.mkdir(parents=True, exist_ok=True)

        # 断点续传 done-marker(§10.2): primer3 原始输出
        # |Checkpoint done-marker (§10.2): raw primer3 output
        self.primer3_output_file = self.design_dir / 'primer3_output.txt'

        self.logger_manager = Primer3Logger(self.output_path / '99_logs')
        self.logger = self.logger_manager.get_logger()

        self.cmd_runner = CommandRunner(self.logger)

        self.logger.info("Primer3引物设计评估器初始化完成|Primer3 evaluator initialized")
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

            # 步骤1: 记录软件版本(§12.5)|Step 1: record software versions (§12.5)
            self._record_software_versions()

            # 步骤2: 解析FASTA文件|Step 2: Parse FASTA file
            self.logger.info("步骤2: 解析FASTA文件|Step 2: Parsing FASTA file")
            sequences = self._parse_fasta()
            self.logger.info(
                f"成功解析|Successfully parsed {format_number(len(sequences))} 条序列|sequences"
            )

            # 步骤3: 运行Primer3(断点续传; 达阈值自动并行)
            # |Step 3: Run Primer3 (resume; auto-parallel at threshold)
            self._run_primer3(sequences)

            # 步骤4: 解析Primer3输出|Step 4: Parse Primer3 output
            self.logger.info("步骤4: 解析Primer3输出结果|Step 4: Parsing Primer3 output")
            parsed_results = self._parse_primer3_output()

            # 步骤5: 格式化并保存结果|Step 5: Format and save results
            self.logger.info("步骤5: 格式化并保存结果|Step 5: Formatting and saving results")
            self._format_and_save_results(parsed_results, sequences)

            self.logger.info("Primer3引物设计流程完成|Primer3 primer design pipeline completed")
            return True

        except Exception as e:
            self.logger.error(f"引物设计流程失败|Primer design pipeline failed: {e}")
            return False

        finally:
            self._cleanup_tmp()

    def _parse_fasta(self):
        """解析FASTA文件|Parse FASTA file"""
        parser = FastaParser(self.config.input_fasta)
        return parser.parse()

    def _generate_primer3_input(self, sequences):
        """生成Primer3输入格式|Generate Primer3 input format"""
        settings = self.config.get_primer3_settings()
        generator = Primer3InputGenerator(settings)
        return generator.generate_batch_input(sequences, self.config)

    def _run_primer3(self, sequences):
        """
        运行Primer3(断点续传; 达阈值自动并行)|Run Primer3 (resume; auto-parallel at threshold)

        primer3_core 本身单线程, 并行=把序列均分给多个 primer3_core 进程。
        序列数 < parallel_threshold 或 threads=1 时单次执行, stdout 直写
        01_primer_design/primer3_output.txt(兼断点续传 done-marker);
        达到阈值则按 threads 切块并行, 完成后按块顺序合并, 输出与单进程
        完全一致; 任一块失败则整体失败且不落合并产物。
        |primer3_core is single-threaded; parallelism splits sequences across
        multiple primer3_core processes. Below threshold or threads=1: one
        run whose stdout lands on the resume done-marker; at threshold:
        chunked concurrent runs merged in chunk order (byte-identical to the
        single-process output); any failed chunk fails the whole run.

        Args:
            sequences: 序列列表[(seq_id, seq), ...]|Sequence list
        """
        if self._is_step_completed(self.primer3_output_file):
            self.logger.info(
                f"跳过已完成步骤|Skipping completed step: primer3 ({self.primer3_output_file})"
            )
            return

        n_seqs = len(sequences)
        use_parallel = (
            self.config.threads > 1 and n_seqs >= self.config.parallel_threshold
        )

        if not use_parallel:
            primer3_input = self._generate_primer3_input(sequences)
            self._run_primer3_single(primer3_input)
            return

        n_chunks = min(self.config.threads, n_seqs)
        self.logger.info(
            f"序列数 {n_seqs} 达到阈值 {self.config.parallel_threshold}, "
            f"启动 {n_chunks} 个primer3进程并行|"
            f"{n_seqs} sequences >= threshold {self.config.parallel_threshold}, "
            f"running {n_chunks} primer3 processes in parallel"
        )
        chunks = self._split_sequences(sequences, n_chunks)
        self._run_primer3_parallel(chunks)

    @staticmethod
    def _split_sequences(sequences, n_chunks):
        """按全局顺序均分为至多n块(每块非空)|Split into <=n_chunks ordered non-empty chunks"""
        n_seqs = len(sequences)
        chunk_size = (n_seqs + n_chunks - 1) // n_chunks
        return [sequences[i:i + chunk_size] for i in range(0, n_seqs, chunk_size)]

    def _run_primer3_single(self, primer3_input: str):
        """
        单进程执行primer3(stdout直写done-marker)|Single-process run (stdout to done-marker)

        Args:
            primer3_input: Primer3输入格式字符串|Primer3 input format string
        """
        # 临时输入文件写入 output_dir/tmp(§12.4, 防超算 /tmp 爆满)
        # |Temp input under output_dir/tmp (§12.4, avoid filling HPC /tmp)
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.txt', delete=False, dir=self.tmp_dir
        ) as tmp_input:
            tmp_input.write(primer3_input)
            tmp_input_path = tmp_input.name

        try:
            cmd = build_conda_command(self.config.primer3_core_path, [tmp_input_path])
            success, _, stderr = self.cmd_runner.run(
                cmd,
                description="Primer3引物设计|Primer3 primer design",
                output_file=str(self.primer3_output_file),
            )

            if not success:
                # 失败残留的半截输出会骗过断点续传, 必须删除
                # |A partial output would fake resume success; delete it
                self.primer3_output_file.unlink(missing_ok=True)
                raise RuntimeError(f"Primer3运行失败|Primer3 execution failed: {stderr}")
        finally:
            Path(tmp_input_path).unlink(missing_ok=True)

    def _run_primer3_parallel(self, chunks):
        """
        并行执行多个primer3_core并按块序合并|Run chunks concurrently, merge in order

        Args:
            chunks: 序列块列表|List of sequence chunks
        """
        n_chunks = len(chunks)
        chunk_outputs = [
            self.tmp_dir / f"primer3_chunk_{i:03d}.txt" for i in range(n_chunks)
        ]

        def _run_chunk(chunk_index: int):
            """执行单块并返回(块号, 是否成功, 错误)|Run one chunk, return (index, success, stderr)"""
            chunk_input = self._generate_primer3_input(chunks[chunk_index])
            with tempfile.NamedTemporaryFile(
                mode='w', suffix='.txt', delete=False, dir=self.tmp_dir,
                prefix=f"chunk_{chunk_index:03d}_",
            ) as tmp_input:
                tmp_input.write(chunk_input)
                tmp_input_path = tmp_input.name
            try:
                cmd = build_conda_command(self.config.primer3_core_path, [tmp_input_path])
                success, _, stderr = self.cmd_runner.run(
                    cmd,
                    description=(
                        f"Primer3引物设计(块 {chunk_index + 1}/{n_chunks})|"
                        f"Primer3 primer design (chunk {chunk_index + 1}/{n_chunks})"
                    ),
                    output_file=str(chunk_outputs[chunk_index]),
                )
                return chunk_index, success, stderr
            finally:
                Path(tmp_input_path).unlink(missing_ok=True)

        try:
            with ThreadPoolExecutor(max_workers=n_chunks) as pool:
                futures = [pool.submit(_run_chunk, i) for i in range(n_chunks)]
                failures = []
                for future in as_completed(futures):
                    chunk_index, success, stderr = future.result()
                    if not success:
                        failures.append((chunk_index, stderr))

            if failures:
                # 任一块失败即整体失败; 合并产物不落盘, 断点续传只认最终文件
                # |Any failed chunk fails the run; no merged marker is written
                chunk_index, stderr = sorted(failures)[0]
                raise RuntimeError(
                    f"Primer3运行失败(块 {chunk_index + 1}/{n_chunks})|"
                    f"Primer3 execution failed (chunk {chunk_index + 1}/{n_chunks}): {stderr}"
                )

            # 按块顺序合并, 输出与单进程执行完全一致|Merge in chunk order
            with open(self.primer3_output_file, 'w', encoding='utf-8') as out_fh:
                for chunk_output in chunk_outputs:
                    out_fh.write(chunk_output.read_text(encoding='utf-8'))
            self.logger.info(
                f"并行执行完成, 已合并到|Parallel run merged into: {self.primer3_output_file}"
            )
        finally:
            for chunk_output in chunk_outputs:
                chunk_output.unlink(missing_ok=True)

    def _parse_primer3_output(self):
        """解析Primer3输出文件|Parse Primer3 output file"""
        parser = Primer3OutputParser()
        return parser.parse(self.primer3_output_file.read_text(encoding='utf-8'))

    def _format_and_save_results(self, parsed_results, sequences):
        """格式化并保存结果|Format and save results"""
        seq_dict = {seq_id: sequence for seq_id, sequence in sequences}

        formatter = ResultsFormatter()
        df = formatter.to_dataframe(parsed_results, seq_dict, self.config.output_header_lang)

        output_file = self.design_dir / f"primers_result.{self.config.output_format}"
        formatter.save_table(df, str(output_file), self.config.output_format)

        self.logger.info(f"结果已保存到|Results saved to: {output_file}")
        self.logger.info(f"共设计|Total primers designed: {format_number(len(df))} 对引物|primer pairs")

        return df

    def _probe_version(self) -> str:
        """探测 primer3 版本(失败降级为 unknown)|Probe primer3 version (degrade to unknown)"""
        cmd = build_conda_command(self.config.primer3_core_path, ['--about'])
        success, stdout, _ = self.cmd_runner.run(
            cmd, description="探测primer3版本|Probing primer3 version"
        )
        if success and stdout.strip():
            return stdout.strip().split('\n')[0]
        self.logger.warning("无法获取primer3版本, 记录为unknown|Cannot probe primer3 version, use unknown")
        return 'unknown'

    def _record_software_versions(self):
        """生成00_pipeline_info/software_versions.yml(§12.5)|Generate software_versions.yml (§12.5)"""
        version = self._probe_version()
        yml_file = self.info_dir / 'software_versions.yml'
        lines = [
            'software_versions:',
            f'  primer3: "{version}"',
            'module:',
            f'  primer3_biopytools: "{__version__}"',
            'params:',
            f'  method: {self.config.method}',
            f'  primer_num_return: {self.config.primer_num_return}',
            f'  output_format: {self.config.output_format}',
        ]
        yml_file.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        self.logger.info(f"版本信息已记录|Software versions recorded: {yml_file}")

    @staticmethod
    def _is_step_completed(output_file: Path) -> bool:
        """断点续传判断(§10.2)|Checkpoint resume check (§10.2)"""
        return output_file.exists()

    def _cleanup_tmp(self):
        """清理 output_dir/tmp(§12.4)|Clean up output_dir/tmp (§12.4)"""
        if not self.tmp_dir.exists():
            return
        for tmp_file in self.tmp_dir.iterdir():
            if tmp_file.is_file():
                tmp_file.unlink(missing_ok=True)
