"""
MEME Parser主程序|MEME Parser Main Module
"""

import pandas as pd
import argparse
from pathlib import Path
from .config import MemeParserConfig
from .parser import MemeParser
from .runner import MemeRunner
from .sequence_extractor import MemeSequenceExtractor
from .utils import MemeParserLogger


class MemeParserPipeline:
    """MEME Parser主流程|MEME Parser Pipeline"""

    def __init__(self, config: MemeParserConfig, logger):
        self.config = config
        self.logger = logger
        self.parser = MemeParser(logger, config)
        self.runner = MemeRunner(logger, config)
        self.extractor = MemeSequenceExtractor(logger, config)

    def run(self):
        """运行完整的解析流程|Run complete parsing pipeline"""
        self.logger.info("="*60)
        self.logger.info("MEME Parser流程开始|MEME Parser pipeline starting")
        self.logger.info("="*60)

        # 1. 运行MEME（如果不是解析模式）|Run MEME (if not parse-only mode)
        if not self.config.parse_only:
            # 检查MEME输出是否已存在|Check if MEME output already exists
            xml_file = Path(self.config.get_meme_xml_path())

            if xml_file.exists():
                self.logger.info(f"检测到已存在的MEME输出|Found existing MEME output: {xml_file}")
                self.logger.info("跳过MEME运行，直接解析|Skipping MEME run, parsing existing output")
            else:
                self.logger.info("未检测到MEME输出，开始运行|No MEME output found, starting run")
                success = self.runner.run()
                if not success:
                    self.logger.error("MEME运行失败，终止流程|MEME run failed, aborting pipeline")
                    return None
        else:
            self.logger.info("解析模式，跳过MEME运行|Parse-only mode, skipping MEME run")

        # 2. 解析MEME输出|Parse MEME output
        df = self.parser.parse()

        # 3. 提取并添加motif序列（如果需要）|Extract and add motif sequences (if requested)
        if self.config.extract_motif_seqs:
            df = self._extract_and_add_motif_sequences(df)

        # 4. 保存结果|Save results
        self._save_results(df)

        self.logger.info("="*60)
        self.logger.info("MEME Parser流程完成|MEME Parser pipeline completed")
        self.logger.info(f"共解析|Total parsed: {len(df)}条motif记录")
        self.logger.info("="*60)

        return df

    def _save_results(self, df: pd.DataFrame):
        """
        保存结果到文件|Save results to files

        Args:
            df: 结果DataFrame|Result DataFrame
        """
        if df is None or len(df) == 0:
            self.logger.warning("没有数据需要保存|No data to save")
            return

        self.logger.info("保存结果|Saving results")

        # 保存TSV|Save TSV
        if self.config.output_tsv:
            tsv_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.tsv"
            df.to_csv(tsv_file, index=False, sep='\t', encoding='utf-8')
            self.logger.info(f"TSV已保存|TSV saved: {tsv_file}")

        # 保存CSV|Save CSV
        if self.config.output_csv:
            csv_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.csv"
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            self.logger.info(f"CSV已保存|CSV saved: {csv_file}")

        # 保存Excel|Save Excel
        if self.config.output_excel:
            try:
                excel_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.xlsx"
                df.to_excel(excel_file, index=False, engine='openpyxl')
                self.logger.info(f"Excel已保存|Excel saved: {excel_file}")
            except ImportError:
                self.logger.warning("openpyxl未安装，跳过Excel输出|openpyxl not installed, skipping Excel output")
            except Exception as e:
                self.logger.error(f"Excel保存失败|Excel save failed: {e}")

    def _extract_and_add_motif_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        提取并添加motif序列到DataFrame|Extract and add motif sequences to DataFrame

        Args:
            df: motif信息DataFrame|Motif information DataFrame

        Returns:
            pd.DataFrame: 包含motif序列的DataFrame|DataFrame with motif sequences
        """
        if df is None or len(df) == 0:
            self.logger.warning("没有motif信息，跳过序列提取|No motif information, skipping sequence extraction")
            return df

        self.logger.info("提取motif序列并添加到表格|Extracting motif sequences and adding to table")

        # 确定输入FASTA文件|Determine input FASTA file
        if self.config.parse_only:
            # 解析模式：优先使用指定的input_fasta
            if self.config.input_fasta:
                input_fasta = self.config.input_fasta
            else:
                # 尝试自动查找FASTA文件
                input_dir = Path(self.config.input_file).parent
                possible_fasta = [
                    input_dir / "*.fa",
                    input_dir / "*.fasta",
                    input_dir / "*.faa",
                    input_dir / "*.pep"
                ]

                import glob
                input_fasta = None
                for pattern in possible_fasta:
                    matches = glob.glob(str(pattern))
                    if matches:
                        input_fasta = matches[0]
                        break

                if input_fasta is None:
                    self.logger.warning("无法找到输入FASTA文件，跳过序列提取|Cannot find input FASTA file, skipping sequence extraction")
                    self.logger.warning("请使用--input-fasta参数指定|Please use --input-fasta parameter to specify")
                    return df
        else:
            # 运行模式：input_file就是FASTA文件
            input_fasta = self.config.input_file

        # 提取motif序列|Extract motif sequences
        motif_seqs = self.extractor.extract_motif_sequences(df, input_fasta)

        if not motif_seqs:
            return df

        # motif_seqs已经使用(seq_id, motif_id)作为tuple key
        # motif_seqs already uses (seq_id, motif_id) as tuple key
        seq_dict = motif_seqs

        # 添加motif_sequence列到DataFrame|Add motif_sequence column to DataFrame
        df['motif_sequence'] = df.apply(
            lambda row: seq_dict.get((row['sequence_id'], row['motif_id']), ''),
            axis=1
        )

        # 重新排列列顺序，将motif_sequence放在最后|Reorder columns, put motif_sequence at the end
        cols = list(df.columns)
        if 'motif_sequence' in cols:
            cols.remove('motif_sequence')
            cols.append('motif_sequence')
            df = df[cols]

        # 保存FASTA文件（额外保存）|Save FASTA file (additional save)
        output_file = f"{self.config.output_prefix}_motifs.fasta"
        self.extractor.write_motif_fasta(motif_seqs, output_file)

        self.logger.info(f"已添加motif_sequence列到表格|Added motif_sequence column to table")

        return df


def run_meme_parser(config: MemeParserConfig, logger):
    """
    运行MEME Parser的便捷函数|Convenience function to run MEME Parser

    Args:
        config: 配置对象|Configuration object
        logger: 日志对象|Logger object

    Returns:
        pd.DataFrame: 结果DataFrame|Result DataFrame
    """
    pipeline = MemeParserPipeline(config, logger)
    return pipeline.run()


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='MEME Parser工具：运行MEME并解析结果|MEME Parser Tool: Run MEME and parse results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input-file', required=True,
                       help='输入FASTA文件或MEME输出文件(xml/txt)|Input FASTA file or MEME output file (xml/txt)')
    parser.add_argument('--parse-only', action='store_true',
                       help='仅解析模式，不运行MEME|Parse-only mode, do not run MEME')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-prefix', default='meme_results',
                       help='输出文件前缀|Output file prefix')
    parser.add_argument('--output-dir', default='.',
                       help='输出目录|Output directory')
    parser.add_argument('--no-tsv', action='store_false', dest='output_tsv',
                       help='不输出TSV文件|Do not output TSV file')
    parser.add_argument('--no-csv', action='store_false', dest='output_csv',
                       help='不输出CSV文件|Do not output CSV file')
    parser.add_argument('--no-excel', action='store_false', dest='output_excel',
                       help='不输出Excel文件|Do not output Excel file')

    # MEME软件路径|MEME software path
    parser.add_argument('--meme-path',
                       default='~/miniforge3/envs/meme_v.5.5.9/bin/meme',
                       help='MEME可执行文件路径|MEME executable path')

    # MEME参数|MEME parameters
    parser.add_argument('-protein', action='store_true', default=True,
                       help='输入序列为蛋白质|Input sequences are protein (default: True)')
    parser.add_argument('-dna', action='store_true',
                       help='输入序列为DNA|Input sequences are DNA')
    parser.add_argument('-mod', default='zoops',
                       choices=['zoops', 'anr', 'oor'],
                       help='Motif分布模式|Motif distribution mode')
    parser.add_argument('-nmotifs', type=int, default=10,
                       help='Motif数量|Number of motifs')
    parser.add_argument('-minw', type=int, default=6,
                       help='最小motif宽度|Minimum motif width')
    parser.add_argument('-maxw', type=int, default=50,
                       help='最大motif宽度|Maximum motif width')
    parser.add_argument('-objfun', default='classic',
                       choices=['classic', 'de', 'ce', 'cd'],
                       help='目标函数|Objective function')
    parser.add_argument('-markov-order', type=int, default=0,
                       help='Markov链阶数|Markov chain order')
    parser.add_argument('--no-extract-motifs', action='store_false', dest='extract_motif_seqs',
                       help='不提取motif序列|Do not extract motif sequences')
    parser.add_argument('--input-fasta',
                       help='原始FASTA文件路径（解析模式时用于提取motif序列）|Original FASTA file path (for motif extraction in parse-only mode)')

    args = parser.parse_args()

    # 初始化日志|Initialize logging
    logger_manager = MemeParserLogger()
    logger = logger_manager.get_logger()

    # 创建配置|Create configuration
    config = MemeParserConfig(
        input_file=args.input_file,
        output_prefix=args.output_prefix,
        output_dir=args.output_dir,
        parse_only=args.parse_only,
        meme_path=args.meme_path,
        protein=args.protein,
        dna=args.dna,
        mod=args.mod,
        nmotifs=args.nmotifs,
        minw=args.minw,
        maxw=args.maxw,
        objfun=args.objfun,
        markov_order=args.markov_order,
        output_tsv=args.output_tsv,
        output_csv=args.output_csv,
        output_excel=args.output_excel,
        extract_motif_seqs=args.extract_motif_seqs,
        input_fasta=args.input_fasta
    )

    # 验证配置|Validate configuration
    try:
        config.validate()
    except ValueError as e:
        logger.error(f"配置错误|Configuration error:\n{e}")
        return 1

    # 运行分析|Run analysis
    try:
        run_meme_parser(config, logger)
        return 0
    except Exception as e:
        logger.error(f"运行失败|Run failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
