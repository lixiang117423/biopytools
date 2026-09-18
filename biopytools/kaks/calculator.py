"""
Ka/Ks Calculator计算引擎模块|Ka/Ks Calculator Calculation Engine Module
功能: 调用KaKs_Calculator2.0进行Ka/Ks计算|Ka/Ks calculation engine using KaKs_Calculator2.0

重要说明|Important Notes:
- KaKs_Calculator2.0 需要AXT格式输入文件，不是FASTA格式|Requires AXT format input, not FASTA
- AXT格式每条记录4行: 序列名称 / 序列1 / 序列2 / 空行
  AXT format: 4 lines per record - sequence_name / sequence_1 / sequence_2 / blank_line
- 方法名称必须使用官方名称，如GMYN而不是gamma-MYN|Method names must use official names
"""

import os
import subprocess
import tempfile
from typing import Dict, List, Tuple
from .config import KaKsConfig
from .aligner import CodonAligner, PairAlignmentResult
from .utils import KaKsLogger, build_conda_command


class KaKsCalculator:
    """Ka/Ks计算引擎|Ka/Ks calculation engine"""

    def __init__(self, logger: KaKsLogger, kaks_path: str = "KaKs_Calculator",
                 check_installation: bool = True):
        """
        初始化计算引擎|Initialize calculation engine

        Args:
            logger: 日志器实例|Logger instance
            kaks_path: KaKs_Calculator可执行文件路径|Path to KaKs_Calculator executable
            check_installation: 是否检查安装(测试时关闭)|Whether to check installation (off in tests)
        """
        self.logger = logger
        self.kaks_path = kaks_path
        self.config = KaKsConfig()
        self.aligner = CodonAligner()
        if check_installation:
            self._check_kaks_installation()

    def _check_kaks_installation(self):
        """检查KaKs_Calculator是否可用|Check if KaKs_Calculator is available"""
        try:
            self.logger.info("检查KaKs_Calculator安装|Checking KaKs_Calculator installation")
            cmd = build_conda_command(self.kaks_path, ["-h"])
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            combined_output = (result.stdout + result.stderr).lower()
            if "kaks_calculator" not in combined_output and "usage" not in combined_output:
                raise FileNotFoundError(f"KaKs_Calculator not found or unexpected output")
            self.logger.success("KaKs_Calculator2.0 可用|KaKs_Calculator2.0 is available")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.logger.error("未找到KaKs_Calculator2.0|KaKs_Calculator2.0 not found")
            self.logger.info("请安装KaKs_Calculator2.0或使用--kaks-path指定路径|Please install KaKs_Calculator2.0 or specify path with --kaks-path")
            raise RuntimeError("KaKs_Calculator2.0 not available")

    def prepare_input_file(self, seq1_dict: Dict[str, str], seq2_dict: Dict[str, str],
                          pairs: List[Tuple[str, str, str]], temp_dir: str,
                          align: bool = True) -> Tuple[str, List[PairAlignmentResult]]:
        """
        准备KaKs_Calculator输入文件(AXT格式)|Prepare input file for KaKs_Calculator (AXT format)

        Args:
            seq1_dict: 第一个FASTA的序列字典|Sequences from first FASTA
            seq2_dict: 第二个FASTA的序列字典|Sequences from second FASTA
            pairs: 配对列表|Pairs list
            temp_dir: 临时目录|Temporary directory
            align: True=翻译→逐对蛋白比对→回译(含gap,标准pal2nal式流程);
                False=严格直通模式,长度不等直接报错,永不静默截断
                True=translate, pairwise-align and back-translate (pal2nal-style);
                False=strict passthrough, unequal lengths raise, never silently truncate

        Returns:
            (AXT文件路径, 全部配对比对记录)|AXT file path and alignment records for all pairs
        """
        input_file = os.path.join(temp_dir, self.config.output_files['temp_axt'])

        try:
            mode_desc = "比对模式|alignment mode" if align else "严格直通模式|strict passthrough mode"
            self.logger.info(
                f"准备AXT格式输入文件|Preparing AXT format input file with {len(pairs)} pairs ({mode_desc})"
            )

            if align:
                return self._prepare_aligned(seq1_dict, seq2_dict, pairs, input_file)
            return self._prepare_strict(seq1_dict, seq2_dict, pairs, input_file)

        except Exception as e:
            self.logger.error(f"准备AXT输入文件失败|Failed to prepare AXT input file: {e}")
            raise

    def _prepare_aligned(self, seq1_dict: Dict[str, str], seq2_dict: Dict[str, str],
                         pairs: List[Tuple[str, str, str]], input_file: str
                         ) -> Tuple[str, List[PairAlignmentResult]]:
        """比对模式:全部配对过密码子比对|Alignment mode: all pairs pass through codon alignment"""
        records: List[PairAlignmentResult] = []
        skipped = 0
        with open(input_file, 'w') as f:
            for seq1_id, seq2_id, pair_name in pairs:
                record = self.aligner.align_pair(
                    seq1_id, seq2_id, pair_name, seq1_dict[seq1_id], seq2_dict[seq2_id]
                )
                records.append(record)
                if record.skip_reason:
                    skipped += 1
                    continue
                f.write(f"{pair_name}\n{record.codon_aln1}\n{record.codon_aln2}\n\n")

        written = len(pairs) - skipped
        if skipped > 0:
            self.logger.warning(
                f"跳过无法比对的配对|Skipped unalignable pairs: {skipped}/{len(pairs)} "
                f"(原因见kaks_detailed的skip_reason列|see skip_reason column in kaks_detailed)"
            )
        self.logger.info(
            f"比对完成|Alignment completed: {written}/{len(pairs)} 对写入AXT|pairs written to AXT"
        )
        self.logger.success(f"AXT格式输入文件准备完成|AXT format input file prepared: {input_file}")
        return input_file, records

    def _prepare_strict(self, seq1_dict: Dict[str, str], seq2_dict: Dict[str, str],
                        pairs: List[Tuple[str, str, str]], input_file: str
                        ) -> Tuple[str, List[PairAlignmentResult]]:
        """严格模式:长度不等收集后一次性报错|Strict mode: collect unequal-length errors, raise once"""
        records: List[PairAlignmentResult] = []
        unequal: List[str] = []
        with open(input_file, 'w') as f:
            for seq1_id, seq2_id, pair_name in pairs:
                seq1 = seq1_dict[seq1_id]
                seq2 = seq2_dict[seq2_id]
                if len(seq1) != len(seq2):
                    unequal.append(
                        f"{pair_name}({seq1_id}={len(seq1)}bp, {seq2_id}={len(seq2)}bp)"
                    )
                    continue
                records.append(PairAlignmentResult(
                    seq1_id=seq1_id, seq2_id=seq2_id, pair_name=pair_name,
                    codon_aln1=seq1, codon_aln2=seq2,
                    aln_codons=len(seq1) // 3,
                    cds1_len=len(seq1), cds2_len=len(seq2),
                ))
                f.write(f"{pair_name}\n{seq1}\n{seq2}\n\n")

        if unequal:
            # 一次性抛出全部错误(§六)|Raise all collected errors at once
            sample = "; ".join(unequal[:5])
            raise ValueError(
                f"序列长度不等|Unequal sequence lengths: {len(unequal)}/{len(pairs)} 个配对|pairs "
                f"(严格模式要求等长,请用默认比对模式|strict mode requires equal lengths, "
                f"use default alignment mode). 前|first {min(5, len(unequal))}: {sample}"
            )

        self.logger.success(
            f"AXT格式输入文件准备完成(等长直通)|AXT input prepared (equal-length passthrough): {input_file}"
        )
        return input_file, records

    def _verify_axt_file(self, axt_file: str) -> bool:
        """
        验证AXT文件格式|Verify AXT file format

        Args:
            axt_file: AXT文件路径|AXT file path

        Returns:
            是否有效|Whether valid
        """
        try:
            with open(axt_file, 'r') as f:
                lines = f.readlines()

            if len(lines) < 4:
                self.logger.warning("AXT文件行数不足|AXT file has insufficient lines")
                return False

            first_line = lines[0].strip()
            if not first_line:
                self.logger.warning("AXT头部格式不正确|Invalid AXT header format: empty first line")
                return False

            self.logger.debug(f"AXT文件验证通过|AXT file validation passed: {len(lines)} lines")
            return True

        except Exception as e:
            self.logger.warning(f"AXT文件验证失败|AXT file validation failed: {e}")
            return False

    def run_calculation(self, input_file: str, method: str, temp_dir: str) -> str:
        """
        运行Ka/Ks计算|Run Ka/Ks calculation

        Args:
            input_file: 输入文件路径|Input file path
            method: 计算方法|Calculation method
            temp_dir: 临时目录|Temporary directory

        Returns:
            输出文件路径|Output file path
        """
        output_file = os.path.join(temp_dir, "kaks_output.txt")

        try:
            if not self._verify_axt_file(input_file):
                raise ValueError(f"无效的AXT文件格式|Invalid AXT file format: {input_file}")

            self.logger.info(f"运行Ka/Ks计算|Running Ka/Ks calculation")
            self.logger.info(f"使用方法|Using method: {self.config.get_method_description(method)}")

            args = ["-i", input_file, "-o", output_file, "-m", method]
            cmd = build_conda_command(self.kaks_path, args)
            self.logger.info(f"命令|Command: {' '.join(cmd)}")

            result = subprocess.run(cmd, capture_output=True, text=True)

            self.logger.debug(f"KaKs_Calculator返回码|Return code: {result.returncode}")
            if result.stdout:
                self.logger.debug(f"KaKs_Calculator标准输出|Stdout: {result.stdout}")
            if result.stderr:
                self.logger.debug(f"KaKs_Calculator错误输出|Stderr: {result.stderr}")

            if result.returncode != 0:
                error_msg = f"KaKs_Calculator执行失败|KaKs_Calculator failed (code {result.returncode})"
                if result.stderr:
                    error_msg += f": {result.stderr}"
                self.logger.error(error_msg)
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)

            if not os.path.exists(output_file):
                raise FileNotFoundError(f"输出文件未生成|Output file not created: {output_file}")

            with open(output_file, 'r') as f:
                lines = f.readlines()
                if len(lines) <= 1:
                    raise ValueError("输出文件为空或无有效数据|Output file is empty or has no valid data")

            self.logger.success(f"计算完成|Calculation completed: {output_file}")

            return output_file

        except Exception as e:
            self.logger.error(f"Ka/Ks计算失败|Ka/Ks calculation failed: {e}")
            raise
