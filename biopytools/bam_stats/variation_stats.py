"""
变异与特异性统计分析模块|Variation & Specificity Statistics Analysis Module

包含: Mismatch Rate, GC Content Bias, Soft-clipping Rate
"""

import statistics
from pathlib import Path
from typing import Dict, Any
from .utils import CommandRunner


class VariationStatsAnalyzer:
    """变异与特异性统计分析器|Variation & Specificity Statistics Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def get_variation_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取变异与特异性统计|Get variation & specificity statistics"""
        self.logger.info(
            f"分析变异与特异性统计|Analyzing variation & specificity stats: "
            f"{Path(bam_file).name}"
        )

        stats = {}

        mismatch = self._get_mismatch_rate(bam_file)
        stats.update(mismatch)

        gc_bias = self._get_gc_bias(bam_file)
        stats.update(gc_bias)

        softclip = self._get_softclip_rate(bam_file)
        stats.update(softclip)

        return stats

    def _get_mismatch_rate(self, bam_file: str) -> Dict[str, Any]:
        """获取错配率|Get mismatch rate via samtools stats"""
        cmd = f"{self.config.samtools_path} stats {bam_file}"
        success, output = self.cmd_runner.run(
            cmd, f"samtools stats for {Path(bam_file).name}"
        )
        if not success:
            return {'mismatch_rate': None}

        mismatch_rate = None
        for line in output.strip().split('\n'):
            if line.startswith('SN\tbases mapped (cigar):\t'):
                bases_mapped = int(line.split('\t')[2])
            elif line.startswith('SN\terror rate:\t'):
                mismatch_rate = float(line.split('\t')[2])

        if mismatch_rate is not None:
            return {
                'mismatch_rate': round(mismatch_rate * 100, 4),
                'mismatch_rate_raw': mismatch_rate
            }

        # 回退方案: 通过NM标签计算|Fallback: calculate via NM tag
        return self._get_mismatch_rate_from_nm(bam_file)

    def _get_mismatch_rate_from_nm(self, bam_file: str) -> Dict[str, Any]:
        """通过NM标签计算错配率|Calculate mismatch rate from NM tag"""
        awk_nm = (
            '{nm=0; for(i=12;i<=NF;i++) '
            'if($i ~ /^NM:/){split($i,a,":"); nm=a[2]} '
            'seq_len=length($10); if(seq_len>0) print nm, seq_len}'
        )
        cmd = (
            f"{self.config.samtools_path} view -F 4 -d 0 {bam_file}"
            f"| head -50000"
            f"| awk '{awk_nm}'"
        )
        success, output = self.cmd_runner.run(cmd, "Calculate mismatch rate via NM tag")
        if not success:
            return {'mismatch_rate': None}

        total_nm = 0
        total_bases = 0
        for line in output.strip().split('\n'):
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    total_nm += int(parts[0])
                    total_bases += int(parts[1])
                except ValueError:
                    continue

        if total_bases == 0:
            return {'mismatch_rate': None}

        rate = total_nm / total_bases * 100
        return {
            'mismatch_rate': round(rate, 4),
            'mismatch_rate_raw': rate / 100
        }

    def _get_gc_bias(self, bam_file: str) -> Dict[str, Any]:
        """获取GC含量偏倚|Get GC content bias"""
        cmd = (
            f"{self.config.samtools_path} view -F 4 {bam_file}"
            f"| cut -f10 | head -10000"
        )
        success, output = self.cmd_runner.run(
            cmd, "GC content bias analysis"
        )
        if not success:
            return {'gc_bias': None}

        gc_contents = []
        for line in output.strip().split('\n'):
            seq = line.strip()
            if seq and seq != '*':
                gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
                gc_contents.append(gc)

        if not gc_contents:
            return {'gc_bias': None}

        avg_gc = statistics.mean(gc_contents)
        gc_std = statistics.stdev(gc_contents) if len(gc_contents) > 1 else 0

        # 参考基因组GC对比|Compare with reference genome GC
        ref_gc = None
        if self.config.reference_file:
            ref_gc = self._get_reference_gc(self.config.reference_file)

        bias_info = {
            'average_gc_content': round(avg_gc, 2),
            'gc_content_std': round(gc_std, 2),
            'gc_bias': round(gc_std / avg_gc * 100, 2) if avg_gc > 0 else None,
            'sampled_reads_for_gc': len(gc_contents)
        }
        if ref_gc is not None:
            bias_info['reference_gc'] = ref_gc
            bias_info['gc_deviation'] = round(avg_gc - ref_gc, 2)

        return bias_info

    def _get_reference_gc(self, ref_file: str) -> float:
        """计算参考基因组GC含量|Calculate reference genome GC content"""
        awk_gc = (
            '{gc=gsub(/[GCgc]/,""); at=gsub(/[ATat]/,""); '
            'if(gc+at>0) print gc/(gc+at)*100}'
        )
        cmd = (
            f"grep -v '^>' {ref_file} | tr -d '\\n' | head -c 10000000"
            f"| awk '{awk_gc}'"
        )
        success, output = self.cmd_runner.run(
            cmd, "Calculate reference genome GC content"
        )
        if success and output.strip():
            try:
                return round(float(output.strip()), 2)
            except ValueError:
                return None
        return None

    def _get_softclip_rate(self, bam_file: str) -> Dict[str, Any]:
        """获取软剪切率|Get soft-clipping rate"""
        awk_script = (
            '{cigar=$6; total_len=0; clip_len=0; '
            'while(match(cigar,/[0-9]+[MIDNSHP=X]/)){'
            '  op=substr(cigar,RSTART+RLENGTH-1,1); '
            '  n=substr(cigar,RSTART,RLENGTH-1)+0; '
            '  if(op=="M"||op=="="||op=="X") total_len+=n; '
            '  if(op=="S") clip_len+=n; '
            '  cigar=substr(cigar,RSTART+RLENGTH);'
            '} '
            'if(total_len>0) print clip_len, total_len;}'
        )
        cmd = (
            f"{self.config.samtools_path} view -F 4 -d 0 {bam_file}"
            f"| awk '{awk_script}'"
        )
        success, output = self.cmd_runner.run(
            cmd, "Calculate soft-clipping rate"
        )
        if not success:
            return {'softclip_rate': None}

        total_clip = 0
        total_aligned = 0
        reads_with_clip = 0
        total_reads = 0

        for line in output.strip().split('\n'):
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    clip = int(parts[0])
                    aligned = int(parts[1])
                    total_clip += clip
                    total_aligned += aligned
                    total_reads += 1
                    if clip > 0:
                        reads_with_clip += 1
                except ValueError:
                    continue

        if total_aligned == 0 or total_reads == 0:
            return {'softclip_rate': None}

        return {
            'softclip_rate': round(total_clip / (total_clip + total_aligned) * 100, 4),
            'softclip_bases': total_clip,
            'total_bases_analyzed': total_clip + total_aligned,
            'reads_with_softclip': reads_with_clip,
            'reads_with_softclip_rate': round(reads_with_clip / total_reads * 100, 2)
        }
