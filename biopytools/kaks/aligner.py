"""
密码子比对器模块|Codon Aligner Module
功能: 翻译CDS→逐对蛋白全局比对→回译为等长密码子比对|Translate CDS, align protein pairs, back-translate to codon alignments

背景|Background: KaKs_Calculator要求已对齐的密码子序列;原先不等长直接截断会导致
含indel的配对全部错位。本模块实现标准pal2nal式流程(needle同参BLOSUM62仿射gap),
蛋白gap回译为---三连,KaKs_Calculator自动剔除gap密码子。
"""

import warnings
from dataclasses import dataclass
from typing import Optional, Tuple

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio import BiopythonWarning


@dataclass
class PairAlignmentResult:
    """单对序列比对结果|Alignment result for one sequence pair"""

    seq1_id: str
    seq2_id: str
    pair_name: str
    codon_aln1: str = ""          # 等长密码子比对,含---|Equal-length codon alignment with ---
    codon_aln2: str = ""
    skip_reason: str = ""         # 空=成功|Empty means success
    protein_identity: float = 0.0 # 双非gap列中相同残基占比|Fraction of matching residues in non-gap overlap
    aln_codons: int = 0           # 比对列数(密码子)|Alignment columns in codons
    gap_codons: int = 0           # 任一侧含gap的列数|Columns with a gap on either side
    cds1_len: int = 0
    cds2_len: int = 0


def translate_cds(cds: str) -> Tuple[Optional[str], Optional[str], str]:
    """
    翻译CDS为蛋白并剥离末端终止子|Translate CDS to protein and strip terminal stop

    Returns:
        (protein, cds_stripped, reason): 成功时reason为空|reason empty on success
    """
    cds = cds.strip().upper()
    if not cds:
        return None, None, "empty_sequence"
    if len(cds) % 3 != 0:
        return None, None, "length_not_multiple_of_3"

    with warnings.catch_warnings():
        # NNN等模糊密码子译为X并保留原CDS(无损)|Ambiguous codons become X, CDS kept losslessly
        warnings.simplefilter("ignore", BiopythonWarning)
        protein = str(Seq(cds).translate(table=1))

    if protein.endswith("*"):
        protein = protein[:-1]
        cds = cds[:-3]
    if "*" in protein:
        return None, None, "internal_stop"

    return protein, cds, ""


def back_translate(aligned_prot: str, cds_stripped: str) -> str:
    """
    将蛋白比对回译为密码子比对|Back-translate a protein alignment to codons

    蛋白gap(-)扩展为---,其余残基取对应密码子|Protein gaps expand to ---, others take their codon
    """
    # 比对列含gap,只数非gap残基|Alignment columns include gaps; count non-gap residues only
    residue_count = sum(1 for residue in aligned_prot if residue != "-")
    if residue_count != len(cds_stripped) // 3:
        raise ValueError(
            f"非gap残基数与密码子数不符|Non-gap residues != codon count: "
            f"{residue_count} vs {len(cds_stripped) // 3}"
        )
    codons = []
    pos = 0
    for residue in aligned_prot:
        if residue == "-":
            codons.append("---")
        else:
            codons.append(cds_stripped[pos:pos + 3])
            pos += 3
    return "".join(codons)


class CodonAligner:
    """逐对密码子比对器|Pairwise codon aligner"""

    def __init__(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        # EMBOSS needle经典同参|Classic needle-equivalent parameters
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -11.0
        aligner.extend_gap_score = -1.0
        self._aligner = aligner

    def align_pair(self, seq1_id: str, seq2_id: str, pair_name: str,
                   cds1: str, cds2: str) -> PairAlignmentResult:
        """
        端到端比对一对CDS|Align one CDS pair end to end

        失败(内部终止子等)不抛异常,以skip_reason标记|Failures are flagged via skip_reason, not raised
        """
        result = PairAlignmentResult(
            seq1_id=seq1_id, seq2_id=seq2_id, pair_name=pair_name,
            cds1_len=len(cds1), cds2_len=len(cds2),
        )

        prot1, cds1_stripped, reason1 = translate_cds(cds1)
        prot2, cds2_stripped, reason2 = translate_cds(cds2)
        if reason1:
            result.skip_reason = reason1
            return result
        if reason2:
            result.skip_reason = reason2
            return result

        alignment = self._aligner.align(prot1, prot2)[0]
        aln_prot1, aln_prot2 = alignment[0], alignment[1]

        result.codon_aln1 = back_translate(aln_prot1, cds1_stripped)
        result.codon_aln2 = back_translate(aln_prot2, cds2_stripped)

        n_columns = len(aln_prot1)
        matches = 0
        overlap = 0
        gaps = 0
        for residue1, residue2 in zip(aln_prot1, aln_prot2):
            if residue1 == "-" or residue2 == "-":
                gaps += 1
            else:
                overlap += 1
                if residue1 == residue2:
                    matches += 1

        result.aln_codons = n_columns
        result.gap_codons = gaps
        result.protein_identity = matches / overlap if overlap > 0 else 0.0
        return result
