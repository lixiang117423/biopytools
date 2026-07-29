"""
PredGPI GPI锚定蛋白预测主模块|PredGPI GPI-anchor prediction main module
"""

import argparse
import os
import sys
import time

from .config import PredGPIConfig
from .utils import (
    PredGPILogger,
    setup_predgpi_env,
    classify_gpi,
    write_tsv,
    format_number,
    generate_software_version_yml,
)


class PredGPIPredictor:
    """PredGPI预测主类|PredGPI predictor main class"""

    # 非标准氨基酸替换映射|Non-standard amino acid substitution map
    _AA_REPLACE = str.maketrans({"U": "C", "Z": "A", "B": "A", "X": "A"})

    def __init__(self, **kwargs):
        # 初始化配置(含validate)|Init config (with validate)
        self.config = PredGPIConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Init logging
        self.logger_manager = PredGPILogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 设置predgpi环境并加载模型|Setup predgpi env and load models
        self._setup_and_load_models()

    def _setup_and_load_models(self):
        """设置predgpi环境并加载HMM/SVM模型|Setup env and load HMM/SVM models"""
        cfg = self.config
        setup_predgpi_env(cfg.predgpi_home)

        # 延迟导入predgpilib模块(sys.path已就绪)|Deferred import (sys.path ready)
        from predgpilib.hmm import HMM_IO
        from predgpilib.svm import SVMLike

        # 选择模型文件|Select model files
        gpidat = os.path.join(cfg.predgpi_home, "GPIDAT")
        if cfg.conservative:
            hmm_file = os.path.join(gpidat, "PHMM.TOT.ss.mod_CSDGN")
        else:
            hmm_file = os.path.join(gpidat, "PHMM.TOT.ss.mod")

        self.logger.info(f"加载模型|Loading models from: {gpidat}")
        self.hmmmod = HMM_IO.get_hmm(hmm_file)
        self.svm = SVMLike.getSVMLight(os.path.join(gpidat, "MOD"))
        self.logger.info("模型加载完成|Models loaded")

    def _predict_sequence(self, seq_id: str, sequence: str) -> dict:
        """
        预测单条序列|Predict single sequence

        Args:
            seq_id: 序列ID|sequence ID
            sequence: 氨基酸序列|amino acid sequence

        Returns:
            预测结果字典|prediction result dict
        """
        # 非标准氨基酸替换|Non-standard AA substitution
        seq_clean = sequence.translate(self._AA_REPLACE)

        if len(seq_clean) <= 40:
            # 序列太短,标记为非GPI|Sequence too short, mark as non-GPI
            return self._build_record(seq_id, sequence, failed=True)

        # 运行预测(延迟导入predGpipe)|Run prediction (lazy import)
        try:
            from predgpi import predGpipe
            lprob, cut, svmout, fitFPR = predGpipe(seq_clean, self.svm, self.hmmmod)
        except Exception:
            # predgpi内部HMM/SVM对某些序列可能抛异常(如numpy数组比较),标记为预测失败
            # |Internal HMM/SVM error on certain sequences — mark as prediction failed
            self.logger.warning(f"预测异常,标记为非GPI|Prediction failed, marking as non-GPI: {seq_id}")
            return self._build_record(seq_id, sequence, failed=True)

        # 强制转换为Python原生类型,避免numpy数组truth value歧义|Coerce to Python scalars
        return self._build_record(
            seq_id, sequence,
            lprob=float(lprob), cut=int(cut),
            svmout=float(svmout), fitFPR=float(fitFPR),
        )

    def _build_record(self, seq_id: str, sequence: str,
                      lprob: float = 0.0, cut: int = 0,
                      svmout: float = 0.0, fitFPR: float = 1.0,
                      failed: bool = False) -> dict:
        """
        根据预测原始值构建TSV记录(纯逻辑,无predgpi依赖)|Build TSV record (pure logic)

        Args:
            seq_id: 序列ID|sequence ID
            sequence: 原始氨基酸序列(用于计算长度/切割位点)|original AA sequence
            lprob: HMM log概率(归一化)|HMM log probability
            cut: HMM预测的GPI锚长度|predicted GPI-anchor length
            svmout: SVM原始预测值|raw SVM score
            fitFPR: 估计的假阳性率|estimated FPR
            failed: 是否预测失败(标记非GPI)|whether prediction failed

        Returns:
            预测结果字典|prediction result dict
        """
        gpi_anchored = fitFPR <= 0.01
        cleavage_site = len(sequence) - cut if gpi_anchored else "-"
        classification, probability = classify_gpi(fitFPR)

        return {
            "ID": seq_id,
            "Length": len(sequence),
            "GPI_Anchored": gpi_anchored,
            "Cleavage_Site": str(cleavage_site),
            "FPR": float(fitFPR),
            "HMM_LogProb": float(lprob),
            "SVM_Score": float(svmout),
            "Probability": probability,
            "Classification": classification,
        }

    def run(self) -> bool:
        """运行预测|Run prediction"""
        cfg = self.config
        tsv_file = os.path.join(cfg.output_dir, f"{cfg.output_prefix}.predgpi.tsv")
        start_time = time.time()

        self.logger.info("=" * 60)
        self.logger.info("PredGPI GPI锚定蛋白预测|PredGPI GPI-anchor prediction")
        self.logger.info("=" * 60)
        self.logger.info(f"输入文件|Input file: {cfg.input_file}")
        self.logger.info(f"输出目录|Output dir: {cfg.output_dir}")
        self.logger.info(f"输出前缀|Prefix: {cfg.output_prefix}")
        self.logger.info(f"predgpi目录|predgpi home: {cfg.predgpi_home}")
        self.logger.info(f"保守模型|Conservative model: {cfg.conservative}")

        # 断点续传:主输出已存在则跳过(§10.2)|Checkpoint resume: skip if TSV exists
        if os.path.exists(tsv_file):
            self.logger.info(f"跳过(已完成)|Skip (already done): {os.path.basename(tsv_file)}")
            return True

        # 读取FASTA|Read FASTA
        from predgpi import readFasta
        sequences = readFasta(cfg.input_file)
        n_total = len(sequences)
        self.logger.info(f"读取序列|Read sequences: {format_number(n_total)}")

        # 逐序列预测|Predict per sequence
        records = []
        n_gpi = 0
        n_failed = 0
        for seq_id, seq in sequences.items():
            rec = self._predict_sequence(seq_id, seq)
            records.append(rec)
            if rec["GPI_Anchored"]:
                n_gpi += 1
            # FPR=1是predgpi对短序列/预测失败的约定标记|FPR=1 marks short/failed (tool convention)
            if rec["FPR"] == 1.0:
                n_failed += 1

        # 写TSV|Write TSV
        write_tsv(records, tsv_file)

        # 统计|Stats
        elapsed = time.time() - start_time
        self.logger.info(f"预测完成|TSV written: {tsv_file}")
        self.logger.info(f"蛋白总数|Total proteins: {format_number(n_total)}")
        self.logger.info(f"GPI锚定蛋白|GPI-anchored: {n_gpi}")
        self.logger.info(f"非GPI蛋白|Non-GPI: {n_total - n_gpi}")
        if n_failed:
            self.logger.warning(f"预测失败序列|Failed sequences: {n_failed} (已标记非GPI|marked non-GPI)")
        self.logger.info(f"耗时|Elapsed: {elapsed:.2f}s")
        self.logger.info("=" * 60)
        self.logger.info("PredGPI预测完成|PredGPI prediction completed")
        self.logger.info("=" * 60)

        # 记录软件版本信息到00_pipeline_info(§12.5)|Record software versions
        try:
            generate_software_version_yml(
                cfg,
                execution_info={
                    "runtime_seconds": round(elapsed, 2),
                },
            )
        except Exception as e:
            self.logger.warning(f"软件版本信息写入失败|Failed to write software versions: {e}")
        return True


def main():
    """命令行入口|CLI entry"""
    parser = argparse.ArgumentParser(
        description="PredGPI GPI锚定蛋白预测|PredGPI GPI-anchor prediction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Examples: biopytools predgpi -i proteins.fa -o output_dir/",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="[FILE] 输入蛋白质FASTA|Input protein FASTA")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="[DIR] 输出目录|Output directory")
    parser.add_argument("--predgpi-home",
                        default="~/software/predgpi",
                        help="predgpi安装目录|predgpi install directory")
    parser.add_argument("--conservative", action="store_true",
                        help="使用保守omega模型|Use conservative omega model")
    parser.add_argument("--prefix", default=None,
                        help="[STR] 输出前缀(默认输入文件名)|Output prefix (default: input filename)")
    args = parser.parse_args()

    try:
        predictor = PredGPIPredictor(
            input_file=args.input,
            output_dir=args.output_dir,
            predgpi_home=args.predgpi_home,
            conservative=args.conservative,
            output_prefix=args.prefix,
        )
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        success = predictor.run()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception as e:
        print(f"运行失败|Run failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
