"""侧翼序列提取模块（samtools faidx）|Flank sequence extractor"""

from typing import Dict, List
from .genotype_analyzer import Candidate
from .coverage_validator import candidate_id
from .utils import build_conda_command


class SequenceExtractor:
    """用samtools faidx提取区间+侧翼序列|Extract region+flank via samtools faidx"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd = cmd_runner

    def extract_flank(self, candidates: List[Candidate]) -> Dict[str, str]:
        """提取每个候选的区间+两侧flank_length序列|Extract region + flank each side"""
        flank = self.config.flank_length
        seqs: Dict[str, str] = {}
        for c in candidates:
            start = max(1, c.indel.pos - flank)
            end = c.indel.end + flank
            region = f"{c.indel.chrom}:{start}-{end}"
            cmd = build_conda_command(self.config.samtools_path,
                                      ['faidx', self.config.genome_fasta, region])
            out = self.cmd.run_capture(cmd, description=f"faidx {region}")
            seq = self._parse_fasta(out)
            if not seq:
                self.logger.warning(f"侧翼序列为空|empty flank: {region}")
            seqs[candidate_id(c)] = seq.upper()
        return seqs

    @staticmethod
    def _parse_fasta(stdout: str) -> str:
        """单条FASTA解析为纯序列|Parse single FASTA to sequence"""
        if not stdout:
            return ""
        seq_lines = [ln.strip() for ln in stdout.splitlines() if not ln.startswith('>')]
        return ''.join(seq_lines)
