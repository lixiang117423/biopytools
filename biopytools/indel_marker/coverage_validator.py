"""覆盖度验证模块（samtools depth + deletion骤降）|Coverage Validator (samtools depth + deletion drop)

对每个候选 INDEL 区间、每个样品用 samtools depth 计算平均覆盖度，
再做两层验证：
  1. 质控层：所有样品平均覆盖度 >= min_depth；
  2. deletion 骤降层：对长度达到 deletion_size_for_coverage_check 的 deletion，
     present 组（携带该 deletion）覆盖度应明显低于 absent 组
     （present_mean / absent_mean < deletion_depth_ratio）。
"""

from typing import Dict, List

from .genotype_analyzer import Candidate
from .samplesheet import SampleInfo
from .utils import build_conda_command


def candidate_id(c: Candidate) -> str:
    """
    生成候选ID|Generate candidate id

    格式|Format: ``{chrom}:{pos}-{end}:{DEL|INS}:{R|S}_spec``
    （R/S = present_group 首字母|R/S = present_group initial）
    被 Task 7/9 复用|Reused by Task 7/9
    """
    t = 'DEL' if c.indel.indel_type == 'deletion' else 'INS'
    g = 'R' if c.present_group == 'resistant' else 'S'
    return f"{c.indel.chrom}:{c.indel.pos}-{c.indel.end}:{t}:{g}_spec"


class CoverageValidator:
    """区间覆盖度计算与骤降验证|Region coverage & drop validation"""

    def __init__(self, config, logger, cmd_runner, samples: Dict[str, SampleInfo]):
        self.config = config
        self.logger = logger
        self.cmd = cmd_runner
        self.samples = samples

    def _region(self, c: Candidate) -> str:
        """samtools depth -r 用的区间字符串|Region string for samtools depth -r"""
        return f"{c.indel.chrom}:{c.indel.pos}-{c.indel.end}"

    def compute(self, candidates: List[Candidate]) -> Dict[str, Dict[str, float]]:
        """
        对每个候选区间、每个样品算平均覆盖度|Mean depth per candidate per sample

        Args:
            candidates: 候选列表|Candidate list

        Returns:
            key1=candidate_id, key2=sample_name, val=mean_depth
        """
        coverage: Dict[str, Dict[str, float]] = {}
        for c in candidates:
            region = self._region(c)
            cid = candidate_id(c)
            coverage[cid] = {}
            for name, info in self.samples.items():
                cmd = build_conda_command(self.config.samtools_path, [
                    'depth', '-r', region,
                    '-Q', str(self.config.min_mapq),
                    '-q', str(self.config.min_baseq),
                    info.bam_path])
                # 完整命令日志由 CommandRunner.run_capture 内部记录|Full command logged inside run_capture
                out = self.cmd.run_capture(cmd, description=f"depth {region} {name}")
                coverage[cid][name] = self._mean_depth(out)
        return coverage

    @staticmethod
    def _mean_depth(stdout) -> float:
        """
        samtools depth输出求平均|Mean of samtools depth output

        兼容2列(pos depth)与3列(chrom pos depth)输出格式，取每行最后一列作为深度。
        Robust to both 2-col (pos depth) and 3-col (chrom pos depth) output;
        takes the last column as depth.
        """
        if not stdout:
            return 0.0
        depths = []
        for ln in stdout.splitlines():
            parts = ln.split('\t')
            if len(parts) < 2:
                continue
            try:
                depths.append(float(parts[-1]))
            except ValueError:
                continue
        return sum(depths) / len(depths) if depths else 0.0

    def validate(self, candidates: List[Candidate], coverage: Dict[str, Dict[str, float]]) -> None:
        """
        质控过滤 + deletion骤降验证（就地更新candidate）|QC + drop validation (in-place)

        为每个 candidate 设置：r_depth_mean / s_depth_mean / depths /
        passes_coverage_qc / deletion_depth_ratio / passes_deletion_drop。
        insertion 或短 deletion 的骤降字段置为 'NA'。
        """
        for c in candidates:
            cid = candidate_id(c)
            cov = coverage.get(cid, {})

            # 按组聚合平均|group means
            r_samples = [n for n, s in self.samples.items() if s.group == 'resistant']
            s_samples = [n for n, s in self.samples.items() if s.group == 'susceptible']
            r_depths = [cov.get(n, 0.0) for n in r_samples]
            s_depths = [cov.get(n, 0.0) for n in s_samples]
            c.r_depth_mean = sum(r_depths) / len(r_depths) if r_depths else 0.0
            c.s_depth_mean = sum(s_depths) / len(s_depths) if s_depths else 0.0
            c.depths = cov

            # 质控层：所有样品>=min_depth|QC: all samples >= min_depth
            all_depths = list(cov.values())
            c.passes_coverage_qc = bool(all_depths) and all(
                d >= self.config.min_depth for d in all_depths)

            # deletion骤降层（仅>=阈值的deletion）|drop (deletion only, >= size threshold)
            if (c.indel.indel_type == 'deletion'
                    and c.indel.indel_size >= self.config.deletion_size_for_coverage_check):
                # present 组携带该 deletion，覆盖度应骤降|present group carries the deletion → lower depth
                present_depths = r_depths if c.present_group == 'resistant' else s_depths
                absent_depths = s_depths if c.present_group == 'resistant' else r_depths
                pmean = sum(present_depths) / len(present_depths) if present_depths else 0.0
                amean = sum(absent_depths) / len(absent_depths) if absent_depths else 0.0
                c.deletion_depth_ratio = pmean / amean if amean > 0 else 1.0
                c.passes_deletion_drop = c.deletion_depth_ratio < self.config.deletion_depth_ratio
            else:
                c.deletion_depth_ratio = 'NA'
                c.passes_deletion_drop = 'NA'
