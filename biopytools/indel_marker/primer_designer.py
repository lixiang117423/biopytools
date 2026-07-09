"""引物设计模块（primer3-py）|Primer Designer via primer3-py

封装primer3-py的design_primers调用，为一个候选INDEL侧翼序列设计PCR引物对|
Wraps primer3-py design_primers to design a PCR primer pair for one candidate's
flank sequence.

注意|Note:
    primer3-py在本机可能未安装，且repo内存在同名biopytools/primer3/包会覆盖包名。
    本模块对import失败做优雅降级：self.primer3设为None，design()直接返回fail结构。
    |primer3-py may not be installed locally and the repo's biopytools/primer3/
    package shadows the name. On import failure this module degrades gracefully:
    self.primer3 is set to None and design() returns a fail-structured dict.
"""

import sys
import traceback


class PrimerDesigner:
    """用primer3-py设计PCR引物|Design PCR primers via primer3-py"""

    def __init__(self, config, logger):
        """绑定配置与日志，尝试加载primer3|Bind config/logger, attempt to load primer3

        Args:
            config: IndelMarkerConfig（使用flank_length/primer_product_min/max）|
                    IndelMarkerConfig (uses flank_length/primer_product_min/max)
            logger: 日志器|logger instance
        """
        self.config = config
        self.logger = logger
        self.primer3 = None
        try:
            import primer3
            self.primer3 = primer3
        except ImportError as e:
            self._log_import_failure(e)

    def _log_import_failure(self, err: ImportError):
        """记录primer3导入失败的诊断信息，区分未安装与被本地包覆盖|
        Log diagnostic for primer3 import failure; distinguish not-installed
        from shadowed-by-local-biopytools/primer3/."""
        # 导入失败后Python会清理sys.modules，无法直接据此判断；查traceback源文件|
        # Python clears sys.modules on failure, so inspect traceback source files.
        shadowed = False
        try:
            for frame in traceback.extract_tb(err.__traceback__):
                fpath = (frame.filename or '').replace('\\', '/').lower()
                if 'biopytools/primer3' in fpath:
                    shadowed = True
                    break
        except Exception:
            shadowed = False

        if shadowed:
            self.logger.error(
                "primer3导入失败：本地biopytools/primer3/包覆盖了primer3-py包名，"
                "请在不包含该子包的目录运行，或安装primer3-py并在包外调用|"
                f"primer3 import failed: shadowed by local biopytools/primer3/ package: "
                f"{type(err).__name__}: {err}"
            )
        else:
            self.logger.error(
                "primer3-py未安装，引物设计将被跳过（请安装: pip install primer3-py）|"
                f"primer3-py not installed, primer design will be skipped "
                f"({type(err).__name__}: {err}). Please install: pip install primer3-py"
            )
        self.primer3 = None

    def design(self, candidate_id: str, flank_seq: str) -> dict:
        """为一个候选设计引物|Design primers for one candidate

        Args:
            candidate_id: 候选标识（透传为SEQUENCE_ID）|candidate id (passed as SEQUENCE_ID)
            flank_seq: 侧翼序列字符串|flank sequence string

        Returns:
            dict: {left_primer, right_primer, product_size, tm_left, tm_right, primer_status}
                  设计失败时 primer_status='fail'，引物字段为''|
                  on failure primer_status='fail' and primer fields are ''
        """
        empty = {
            'left_primer': '',
            'right_primer': '',
            'product_size': '',
            'tm_left': '',
            'tm_right': '',
            'primer_status': 'fail',
        }

        # primer3不可用直接降级|Degrade when primer3 unavailable
        if self.primer3 is None:
            return empty

        # 序列过短无法容纳双侧引物+最小产物|Sequence too short for both primers + min product
        min_len = 2 * 18 + self.config.primer_product_min
        if not flank_seq or len(flank_seq) < min_len:
            self.logger.warning(
                f"侧翼序列过短，无法设计引物|flank too short to design primers: "
                f"{candidate_id} (len={len(flank_seq) if flank_seq else 0}, need>={min_len})"
            )
            return empty

        seq_args = {
            'SEQUENCE_ID': candidate_id,
            'SEQUENCE_TEMPLATE': flank_seq,
        }
        global_args = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_PRODUCT_SIZE_RANGE': [
                [self.config.primer_product_min, self.config.primer_product_max]],
            'PRIMER_NUM_RETURN': 1,
        }

        try:
            result = self.primer3.bindings.design_primers(seq_args, global_args)
        except Exception as e:
            self.logger.warning(
                f"primer3调用异常|primer3 design_primers error {candidate_id}: "
                f"{type(e).__name__}: {e}"
            )
            return empty

        n_returned = result.get('PRIMER_PAIR_NUM_RETURNED', 0)
        if not n_returned:
            self.logger.info(
                f"未找到合适引物|no suitable primer found: {candidate_id}"
            )
            return empty

        return {
            'left_primer': result.get('PRIMER_LEFT_0_SEQUENCE', ''),
            'right_primer': result.get('PRIMER_RIGHT_0_SEQUENCE', ''),
            'product_size': int(result.get('PRIMER_PAIR_0_PRODUCT_SIZE', 0)),
            'tm_left': round(result.get('PRIMER_LEFT_0_TM', 0.0), 2),
            'tm_right': round(result.get('PRIMER_RIGHT_0_TM', 0.0), 2),
            'primer_status': 'ok',
        }
