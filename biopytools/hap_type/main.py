"""
单倍型分析主程序模块|Haplotype Analysis Main Module
"""

from pathlib import Path

from .config import HapTypeConfig
from .utils import HapTypeLogger
from .builder import HaplotypeBuilder


def _process_single_region(builder, config, logger, chrom, start, end, name, output_files):
    """处理单个基因组区间|Process a single genomic region"""
    logger.info(f"--- 处理区间|Processing region: {chrom}:{start}-{end} ({name}) ---")

    result = builder.build_from_vcf(
        vcf_file=str(config.vcf_path),
        chrom=chrom,
        start=start,
        end=end,
        hetero_remove=config.hetero_remove,
        na_drop=config.na_drop,
        hap_prefix=config.hap_prefix,
        pad=config.pad,
    )

    if not result:
        logger.warning(f"区间 {name} 无变异位点，跳过|Region {name} has no variants, skipping")
        return None

    builder.export_table(result["result"], output_files['hap_result'])
    builder.export_xlsx(result["result"], output_files['hap_result_xlsx'])
    builder.export_table(result["summary"], output_files['hap_summary'])
    builder.export_xlsx(result["summary"], output_files['hap_summary_xlsx'])

    # 样本-单倍型映射表: chr, start, end, name, sample, haplotype
    header = ["Chr", "Start", "End", "Name", "Sample", "Haplotype"]
    sample_hap_rows = [header]
    for sample_name, hap_id in result["sample_hap"]:
        sample_hap_rows.append([chrom, start, end, name, sample_name, hap_id])
    builder.export_table(sample_hap_rows, output_files['sample_hap'])
    builder.export_xlsx(sample_hap_rows, output_files['sample_hap_xlsx'])

    logger.info(f"变异位点数|Variants: {result['n_variants']}")
    logger.info(f"单倍型数|Haplotypes: {result['n_haplotypes']}")
    logger.info(f"样本数|Samples: {result['n_samples']}")
    logger.info(f"hapResult: {output_files['hap_result']}")
    logger.info(f"hapSummary: {output_files['hap_summary']}")
    logger.info(f"sampleHap: {output_files['sample_hap']}")
    return result


def main(args) -> int:
    """主函数入口|Main function entry point"""
    config = HapTypeConfig(
        vcf_file=args.vcf,
        region=args.region,
        output=args.output,
        hetero_remove=args.hetero_remove,
        na_drop=args.na_drop,
        hap_prefix=args.hap_prefix,
        pad=args.pad,
    )

    try:
        config.validate()
    except ValueError as e:
        print(f"配置错误|Configuration error:\n{e}")
        return 1

    builder = HaplotypeBuilder

    if config.is_bed_mode:
        return _run_bed_mode(config, builder)
    else:
        return _run_single_mode(config, builder)


def _run_single_mode(config, builder_class) -> int:
    """单区间模式|Single interval mode"""
    output_dir = config.get_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_files = config.get_output_files(config.chrom, config.start, config.end)
    logger_manager = HapTypeLogger(Path(output_files['log']))
    logger = logger_manager.get_logger()

    logger.info("开始单倍型提取|Starting haplotype extraction")
    logger.info(f"区域|Region: {config.chrom}:{config.start}-{config.end}")

    builder = builder_class(logger)
    result = _process_single_region(
        builder, config, logger,
        chrom=config.chrom,
        start=config.start,
        end=config.end,
        name="",
        output_files=output_files,
    )

    if not result:
        logger.error("单倍型提取失败|Haplotype extraction failed")
        return 1

    logger.info("单倍型提取完成|Haplotype extraction completed")
    return 0


def _run_bed_mode(config, builder_class) -> int:
    """BED文件批量模式|BED file batch mode"""
    bed_records = config.parse_bed()
    n_records = len(bed_records)

    if n_records == 0:
        print("BED文件无有效记录|No valid records in BED file")
        return 1

    output_dir = config.get_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    # 共用日志文件
    log_path = output_dir / "hap_type.log"
    logger_manager = HapTypeLogger(log_path)
    logger = logger_manager.get_logger()

    logger.info(f"开始批量单倍型提取|Starting batch haplotype extraction")
    logger.info(f"BED记录数|BED records: {n_records}")
    logger.info(f"VCF文件|VCF file: {config.vcf_file}")

    builder = builder_class(logger)
    success = 0
    skipped = 0
    summary = []

    for i, (chrom, start, end, name) in enumerate(bed_records):
        logger.info(f"[{i + 1}/{n_records}] 区间|Region: {chrom}:{start}-{end} ({name})")

        out_files = config.get_output_files(chrom, start, end, name)
        result = _process_single_region(
            builder, config, logger,
            chrom=chrom, start=start, end=end,
            name=name, output_files=out_files,
        )

        if result is None:
            skipped += 1
        else:
            success += 1
            summary.append(f"  {name}: {result['n_variants']} variants, "
                          f"{result['n_haplotypes']} haplotypes, {result['n_samples']} samples")

    logger.info("=" * 50)
    logger.info(f"批量处理完成|Batch processing completed")
    logger.info(f"成功|Success: {success}/{n_records}")
    if skipped:
        logger.info(f"跳过|Skipped: {skipped}/{n_records}")
    for line in summary:
        logger.info(line)
    logger.info("=" * 50)

    return 0 if skipped < n_records else 1
