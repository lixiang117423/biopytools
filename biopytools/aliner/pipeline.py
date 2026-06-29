"""a-liner pipeline 核心流程|aliner pipeline core flow"""

import os
import shlex
import subprocess
from datetime import datetime
from pathlib import Path

from .config import AlinerConfig
from .utils import (AlinerLogger, build_conda_command, extract_fasta_lengths,
                    check_dependencies, get_tool_version, get_aliner_version, format_number)


class AlinerPipeline:
    """a-liner pipeline 主类|aliner pipeline main class"""

    def __init__(self, config: AlinerConfig):
        """初始化|Initialize"""
        self.config = config
        config.validate()
        self.log_file = os.path.join(config.output_dir, '99_logs', 'aliner_pipeline.log')
        self.logger_manager = AlinerLogger(self.log_file)
        self.logger = self.logger_manager.get_logger()
        self.logger.info(f"输出目录|Output directory: {config.output_dir}")
        self.paf_path = os.path.join(config.output_dir, '01_alignment', f"{config.out_prefix}.paf")
        self.seq_config_path = os.path.join(config.output_dir, '01_alignment', 'sequence_config.txt')
        self.pdf_path = os.path.join(config.output_dir, '02_aliner', f"{config.out_prefix}.pdf")
        self.start_time = datetime.now()   # 记录pipeline开始时间|record pipeline start time

    def run(self):
        """运行完整流程|Run full pipeline"""
        check_dependencies(self.config, self.logger)

        # Step 1-2: 序列长度与校验|lengths + validation
        self.logger.info("获取序列长度|Fetching sequence lengths")
        ref_lengths = extract_fasta_lengths(self.config.ref_fasta, self.config.samtools_path, self.logger)
        query_lengths = extract_fasta_lengths(self.config.query_fasta, self.config.samtools_path, self.logger)
        self._validate_seqs(ref_lengths, query_lengths)

        # Step 3: minimap2（断点续传）|minimap2 with checkpoint
        self.run_minimap2()

        # Step 4: 生成sequence_config|generate sequence_config
        self.generate_sequence_config(ref_lengths, query_lengths)

        # Step 5: 调a-liner|run a-liner
        self.run_aliner()

        # Step 6: 版本元数据|version metadata
        self.generate_software_versions()

        self.logger.info(f"完成|Done. PDF: {self.pdf_path}")

    def _validate_seqs(self, ref_lengths, query_lengths):
        """校验序列存在性与区段越界|Validate seqs exist and regions in bounds"""
        for spec in self.config.ref_specs:
            self._check_one(spec, ref_lengths, 'ref')
        for spec in self.config.query_specs:
            self._check_one(spec, query_lengths, 'query')

    def _check_one(self, spec, lengths, side):
        """校验单条序列规格|Validate one seq spec"""
        seq_id, start, end = spec
        if seq_id not in lengths:
            raise ValueError(f"{side} FASTA中找不到序列|seq not in {side} FASTA: '{seq_id}'。"
                             f"可用|Available: {list(lengths.keys())}")
        if start is not None and end > lengths[seq_id]:
            raise ValueError(f"{side} 区段越界|{side} region out of bounds: "
                             f"'{seq_id}:{start}-{end}' > 全长|full length {lengths[seq_id]}")

    def run_minimap2(self):
        """运行minimap2（断点续传）|Run minimap2 with checkpoint"""
        if os.path.exists(self.paf_path) and os.path.getsize(self.paf_path) > 0:
            self.logger.info(f"跳过已完成步骤|Skipping completed step: minimap2 ({self.paf_path})")
            return
        args = ['-c', '-x', self.config.preset, '-t', str(self.config.threads),
                self.config.ref_fasta, self.config.query_fasta]
        # 关键：传递完整路径，禁止提取命令名（§13.6）|full path, no basename
        cmd = build_conda_command(self.config.minimap2_path, args)
        # 必须记录完整命令到INFO级别（§2.2.1）|log full command at INFO
        self.logger.info(f"执行|Executing: minimap2 比对|minimap2 alignment")
        self.logger.info(f"命令|Command: {' '.join(cmd)} > {self.paf_path}")
        with open(self.paf_path, 'w') as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            self.logger.error(f"minimap2失败|minimap2 failed: {result.stderr}")
            raise RuntimeError(f"minimap2失败|minimap2 failed (returncode={result.returncode})")
        self.logger.info(f"minimap2完成|minimap2 done. PAF: {self.paf_path}")

    def generate_sequence_config(self, ref_lengths, query_lengths):
        """生成sequence_config.txt（配对交错）|Generate sequence_config (paired interleave)"""
        self.logger.info("生成sequence_config|Generating sequence_config")
        rows = []
        n = 0
        for ref_spec, query_spec in zip(self.config.ref_specs, self.config.query_specs):
            for side, spec, lengths in (('ref', ref_spec, ref_lengths),
                                        ('query', query_spec, query_lengths)):
                seq_id, start, end = spec
                full = lengths[seq_id]
                s = start if start is not None else 1
                e = end if end is not None else full
                name = f"{side}-{seq_id}"
                rows.append((n, seq_id, s, e, '+', name))
                n += 1
        with open(self.seq_config_path, 'w') as f:
            f.write("n\tID\tstart\tend\tstrand\tname\n")
            for row in rows:
                f.write("\t".join(str(x) for x in row) + "\n")
        self.logger.info(f"sequence_config已写入|Written: {self.seq_config_path}（{len(rows)} 行|rows）")

    def run_aliner(self):
        """运行a-liner出图|Run a-liner to produce PDF"""
        out_prefix = os.path.splitext(self.pdf_path)[0]   # a-liner自动加.pdf|a-liner appends .pdf
        args = ['-i', self.seq_config_path,
                '--minimap2', self.paf_path,
                '--min_identity', str(self.config.min_identity),
                '--min_alignment_len', str(self.config.min_alignment_len),
                '--colormap', str(self.config.colormap),
                '--figure_size', str(self.config.figure_size[0]), str(self.config.figure_size[1]),
                '--out', out_prefix]
        if self.config.extra_args:
            args.extend(shlex.split(self.config.extra_args))
        # a-liner固定conda环境，显式构造（不走自动检测）|a-liner fixed env, explicit
        cmd = ['conda', 'run', '-n', self.config.aliner_env, '--no-capture-output', 'a-liner'] + args
        # 必须记录完整命令到INFO级别（§2.2.1）|log full command at INFO
        self.logger.info(f"执行|Executing: a-liner 可视化|a-liner visualization")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.logger.error(f"a-liner失败|a-liner failed: {result.stderr}")
            raise RuntimeError(f"a-liner失败|a-liner failed (returncode={result.returncode})")
        self.logger.info(f"a-liner完成|a-liner done. PDF: {self.pdf_path}")

    def generate_software_versions(self):
        """生成software_versions.yml|Generate software_versions.yml"""
        # yaml在方法内导入（测试环境可能无yaml）|import inside method
        import yaml
        end = datetime.now()
        start = self.start_time
        info = {
            'pipeline': {'name': 'biopytools aliner', 'version': '1.0.0'},
            'tools': {
                'minimap2': {'version': get_tool_version(self.config.minimap2_path),
                             'path': self.config.minimap2_path},
                'samtools': {'version': get_tool_version(self.config.samtools_path),
                             'path': self.config.samtools_path},
                'a-liner': {'version': get_aliner_version(self.config.aliner_env),
                            'path': f"conda run -n {self.config.aliner_env}"},
            },
            'parameters': {
                'preset': self.config.preset,
                'min_identity': self.config.min_identity,
                'min_alignment_len': self.config.min_alignment_len,
                'colormap': self.config.colormap,
                'figure_size': list(self.config.figure_size),
                'threads': self.config.threads,
                'ref_fasta': self.config.ref_fasta,
                'query_fasta': self.config.query_fasta,
                'ref_seqs': self.config.ref_seqs,
                'query_seqs': self.config.query_seqs,
            },
            'execution': {
                'start_time': start.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': int((end - start).total_seconds()),
            },
        }
        out = Path(self.config.output_dir, '00_pipeline_info', 'software_versions.yml')
        with open(out, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
        # 参数快照|params snapshot
        params_out = Path(self.config.output_dir, '00_pipeline_info', 'pipeline_params.yaml')
        with open(params_out, 'w') as f:
            yaml.dump(info['parameters'], f, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"版本元数据已写入|Metadata written: {out}")
