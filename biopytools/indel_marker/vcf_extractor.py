"""VCF INDEL提取模块|VCF INDEL Extractor (bcftools)

用bcftools view -v indels 过滤INDEL，管道到bcftools query 导出
CHROM/POS/END/REF/ALT + 每样本GT 矩阵，再在Python侧按strlen精确过滤大小。
"""

import shlex
import subprocess
from typing import List

from .genotype_analyzer import IndelRecord


class VCFExtractor:
    """用bcftools提取INDEL与GT矩阵|Extract INDELs + GT matrix via bcftools"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd = cmd_runner

    def get_sample_names(self) -> List[str]:
        """bcftools query -l 取样本名|List sample names via bcftools query -l"""
        from .utils import build_conda_command
        cmd = build_conda_command(
            self.config.bcftools_path,
            ['query', '-l', self.config.vcf_file],
        )
        out = self.cmd.run_capture(cmd, description='bcftools query -l | list VCF samples')
        if out is None:
            raise RuntimeError("无法读取VCF样本名|Failed to read VCF sample names")
        return [ln.strip() for ln in out.splitlines() if ln.strip()]

    @staticmethod
    def _strip_conda_prefix(wrapped_cmd: List[str]) -> List[str]:
        """
        从conda run包装中提取实际命令|Extract actual cmd from conda run wrapper

        build_conda_command 在conda环境下返回:
          ['conda','run','-n',ENV,'--no-capture-output',cmd,...args]
        管道中禁止 "conda run | conda run" (CLAUDE.md 13.2.1)，
        故提取实际命令后用普通shell管道连接。
        """
        if len(wrapped_cmd) >= 2 and wrapped_cmd[0] == 'conda' and wrapped_cmd[1] == 'run':
            idx = 4  # 跳过 'conda','run','-n',ENV
            if len(wrapped_cmd) > 4 and wrapped_cmd[4] == '--no-capture-output':
                idx = 5
            return wrapped_cmd[idx:]
        return wrapped_cmd

    def extract(self, out_tsv: str, samples: List[str] = None) -> List[IndelRecord]:
        """
        提取INDEL到out_tsv并解析|Extract INDELs to out_tsv and parse

        - bcftools view -v indels 过滤INDEL，管道到query导出GT矩阵
        - 大小/类型过滤在 _parse_matrix 中按strlen(REF)/strlen(ALT)精确完成
          （避免bcftools表达式跨版本差异）
        - 输出列：CHROM POS END REF ALT QUAL + 每样本GT
          (QUAL位于ALT与GT之间，供min_quality过滤)
        - samples可由调用方预先获取后传入，避免重复调用get_sample_names
          Caller may pre-fetch samples to avoid a redundant get_sample_names call

        Args:
            out_tsv: 输出矩阵TSV路径|Output matrix TSV path
            samples: 预取的样本名列表，None则内部调用get_sample_names|
                     Pre-fetched sample names; None => call get_sample_names here
        """
        from .utils import build_conda_command

        bcftools = self.config.bcftools_path
        vcf = self.config.vcf_file

        # view -v indels 取全部INDEL，Python侧按strlen过滤大小
        view_cmd = build_conda_command(bcftools, ['view', '-v', 'indels', vcf])
        # raw string: bcftools解析\t\n为制表/换行|raw str so bcftools interprets \t\n
        # QUAL位于ALT与GT之间，供min_quality过滤|QUAL between ALT and GT for min_quality
        fmt = r'%CHROM\t%POS\t%INFO/END\t%REF\t%ALT\t%QUAL[\t%GT]\n'
        query_cmd = build_conda_command(bcftools, ['query', '-f', fmt])

        # 提取实际bcftools命令，避免"conda run | conda run"
        view_actual = self._strip_conda_prefix(view_cmd)
        query_actual = self._strip_conda_prefix(query_cmd)
        # shlex正确引用含\t[]等特殊字符的format串|quote format string with special chars
        pipeline = (
            f"{shlex.join(view_actual)} | {shlex.join(query_actual)} "
            f"> {shlex.quote(out_tsv)}"
        )

        self.logger.info("执行|Executing: bcftools提取INDEL矩阵|extract INDEL matrix")
        self.logger.info(f"命令|Command: {pipeline}")
        try:
            subprocess.run(
                pipeline, shell=True, check=True,
                capture_output=True, text=True,
            )
        except subprocess.CalledProcessError as e:
            self.logger.error(f"bcftools提取失败|bcftools extract failed: {e.stderr}")
            raise

        if samples is None:
            samples = self.get_sample_names()
        return self._parse_matrix(out_tsv, samples)

    def _parse_matrix(self, matrix_tsv: str, samples: List[str]) -> List[IndelRecord]:
        """
        解析GT矩阵，按min_quality与min/max_indel_size过滤|Parse matrix, filter by quality+size

        矩阵列顺序|Column order: CHROM POS END REF ALT QUAL [sample GT...]
        - QUAL位于cols[5]，GT从cols[6]开始|QUAL at cols[5], GT columns start at cols[6]
        - QUAL为'.'或非数值时保留(bcftools有时省略QUAL，不过度过滤)
          Keep QUAL='.' or non-numeric (bcftools may omit QUAL; don't over-filter)
        - QUAL为有效数值且 < config.min_quality 时丢弃

        genotypes 存储原始GT字符串(如'1/1', './.')，不在此处parse_gt。
        原因：Task 4 的 compute_group_stats 会在统计时调用 parse_gt(raw)，
        若此处预解析为分类('hom_alt')会导致双重解析破坏(parse_gt('hom_alt')→'missing')。
        """
        records = []
        with open(matrix_tsv, encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line or line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) < 5:
                    continue
                chrom = cols[0]
                pos = int(cols[1])
                # INFO/END 可能缺失('.')（常见于 insertion 或部分 caller），用 POS+len(REF)-1 兜底
                # INFO/END may be missing ('.') for insertions/some callers; fall back to POS+len(REF)-1
                try:
                    end = int(cols[2])
                except (ValueError, TypeError):
                    end = pos + len(cols[3]) - 1
                ref = cols[3]
                alt = cols[4]
                # 跳过多等位(ALT含,)与 * allele：不利于PCR标记，且len/GT解析会错
                # skip multi-allelic (ALT has,) and * allele: bad for PCR markers, breaks len/GT parsing
                if ',' in alt or '*' in alt:
                    continue
                # QUAL位于cols[5]；缺失时保留|QUAL at cols[5]; keep when missing
                qual_str = cols[5] if len(cols) > 5 else '.'
                try:
                    qual = float(qual_str)
                    if qual < self.config.min_quality:
                        continue
                except (ValueError, TypeError):
                    pass  # QUAL='.'或非数值，保留|QUAL='.' or non-numeric: keep
                gt_cols = cols[6:]

                # 判定INDEL类型与大小(按strlen)|determine type & size by strlen
                if len(alt) > len(ref):
                    indel_type, indel_size = 'insertion', len(alt) - len(ref)
                elif len(alt) < len(ref):
                    indel_type, indel_size = 'deletion', len(ref) - len(alt)
                else:
                    continue  # SNP/complex，跳过|skip non-indel

                # 大小过滤|size filter
                if not (self.config.min_indel_size <= indel_size <= self.config.max_indel_size):
                    continue

                # 存储原始GT字符串|store RAW gt strings (Task 4 contract)
                n = min(len(samples), len(gt_cols))
                genotypes = {samples[i]: gt_cols[i] for i in range(n)}
                records.append(IndelRecord(
                    chrom=chrom, pos=pos, end=end, ref=ref, alt=alt,
                    indel_type=indel_type, indel_size=indel_size,
                    genotypes=genotypes,
                ))

        self.logger.info(f"解析得到|Parsed {len(records)} 个合格INDEL|qualified INDELs")
        return records
