"""
annorefine SSP(小分泌蛋白)回收通道辅助|SSP (Small Secreted Protein) helpers
- translate_hit_protein: 命中 CDS → 蛋白(复用 evidence 密码子表)|hit CDS -> protein
- run_signalp_batch: 批量 signalp6(复用 signalp 子模块, 一次跑完)|batch signalp6
- parse_signalp_sp_set: 解析 prediction_results.txt → 有信号肽的蛋白 ID 集合|-> SP set

设计依据|Design: 09.二次调试/ssp_lane_design.md
默认关闭(enable_ssp=False), 开启后为短分泌候选放宽长度/同源阈值,
用 信号肽+真实ORF+表达+同源 四重门控压噪声。
|Default off; when on, relaxes len/homology thresholds for short secreted
candidates, gated by signal-peptide + real-ORF + expression + homology.
"""

import os
from typing import Optional, Set

from .evidence import _CODON_TABLE, _revcomp


def translate_hit_protein(hit, genome: dict) -> str:
    """
    命中 CDS 拼接 + 翻译为蛋白|splice hit CDS + translate to protein
    复用 evidence 密码子表; 与 has_complete_orf 同款拼接(仅对完整 ORF 调用)
    |Reuses evidence codon table; same splicing as has_complete_orf (complete ORF only)
    """
    chrom = hit.chrom
    if chrom not in genome:
        return ''
    seq = genome[chrom]
    exons = sorted(hit.cds_exons, key=lambda x: x[0])
    spliced = ''.join(seq[s - 1:e] for s, e, _ in exons).upper()
    if hit.strand == '-':
        spliced = _revcomp(spliced)
    if len(spliced) % 3 != 0:
        # 截到 3 倍数(非完整 ORF 兜底, 正常 SSP 候选已过 has_complete_orf)
        # |trim to codon boundary (fallback; SSP candidates already passed has_complete_orf)
        spliced = spliced[:len(spliced) - (len(spliced) % 3)]
    return ''.join(_CODON_TABLE.get(spliced[i:i + 3], 'X')
                   for i in range(0, len(spliced), 3))


def run_signalp_batch(fasta: str, out_dir: str, config, logger) -> Optional[str]:
    """
    批量跑 signalp6(复用 signalp 子模块 SignalPPredictor, 一次跑完全候选)
    |Run signalp6 in batch via SignalPPredictor (one call for all candidates)
    断点续传: prediction_results.txt 存在则跳过(§10.2)|checkpoint: skip if output exists

    Returns:
        prediction_results.txt 路径; 失败/不可用返回 None
        |path to prediction_results.txt; None on failure/unavailable
    """
    pred = os.path.join(out_dir, 'prediction_results.txt')
    if os.path.exists(pred):
        logger.info(f"复用已有 signalp 结果|reuse existing signalp output: {pred}")
        return pred
    if not os.path.exists(fasta):
        logger.warning(f"SSP 候选 FASTA 不存在|SSP candidate FASTA missing: {fasta}")
        return None
    try:
        from ..signalp.main import SignalPPredictor
    except ImportError as e:
        logger.error(f"无法导入 signalp 子模块|cannot import signalp module: {e}")
        return None
    logger.info(f"批量信号肽预测(SSP 候选)|batch signalp: "
                f"{os.path.basename(fasta)} -> {out_dir}")
    try:
        ok = SignalPPredictor(
            fasta_file=fasta, output_dir=out_dir, organism='eukarya',
            signalp_path=config.signalp_path).run()
    except Exception as e:
        logger.error(f"signalp 运行失败|signalp failed: {e}")
        return None
    if ok and os.path.exists(pred):
        logger.info(f"signalp 完成|signalp done: {pred}")
        return pred
    logger.error("signalp 未产出 prediction_results.txt|signalp produced no output")
    return None


def parse_signalp_sp_set(prediction_results_txt: str,
                         secreted_types=('SP', 'SPI')) -> Set[str]:
    """
    解析 prediction_results.txt → 有经典信号肽的蛋白 ID 集合
    |Parse to set of IDs with a classical signal peptide (Sec/SPI)
    第2列 Prediction ∈ secreted_types 视为分泌(eukarya 模式下基本只有 SP)
    |col-2 Prediction in secreted_types => secreted
    """
    sp_set: Set[str] = set()
    if not prediction_results_txt or not os.path.exists(prediction_results_txt):
        return sp_set
    with open(prediction_results_txt) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 2 and cols[1] in secreted_types:
                sp_set.add(cols[0])
    return sp_set
