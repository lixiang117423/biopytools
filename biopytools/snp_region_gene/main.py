"""SNP区域基因提取主程序|SNP Region Gene Extractor Main Program"""

import argparse
import sys
import os
from .config import SnpRegionConfig
from .utils import SnpRegionLogger, parse_snp_file
from .processor import SnpRegionProcessor


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='SNP区域基因提取工具|SNP Region Gene Extractor Tool\n'
                    '根据SNP位置和指定区域，提取相关基因的CDS和蛋白序列|'
                    'Extract CDS and protein sequences of genes in specified regions around SNPs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--snp',
                       required=True,
                       help='SNP位置文件（格式：Chr01:24770）|SNP position file (format: Chr01:24770)')

    parser.add_argument('-g', '--gff',
                       required=True,
                       help='GFF3注释文件|GFF3 annotation file')

    parser.add_argument('-G', '--genome',
                       required=True,
                       help='基因组FASTA文件|Genome FASTA file')

    # 可选参数|Optional parameters
    parser.add_argument('-l', '--left',
                       type=int,
                       default=0,
                       help='SNP上游距离（bp）|Upstream distance from SNP (bp)')

    parser.add_argument('-r', '--right',
                       type=int,
                       default=0,
                       help='SNP下游距离（bp）|Downstream distance from SNP (bp)')

    parser.add_argument('-p', '--promoter',
                       type=int,
                       default=2000,
                       help='启动子区域距离（bp）|Promoter region distance (bp)')

    parser.add_argument('-o', '--output',
                       default='./snp_region_output',
                       help='输出文件前缀|Output file prefix')

    parser.add_argument('--gffread-path',
                       default='gffread',
                       help='gffread程序路径|gffread program path')

    parser.add_argument('--seqkit-path',
                       default='seqkit',
                       help='seqkit程序路径|seqkit program path')

    parser.add_argument('--keep-temp',
                       action='store_true',
                       help='保留临时文件|Keep temporary files')

    args = parser.parse_args()

    try:
        # 创建配置对象|Create configuration object
        config = SnpRegionConfig(
            snp_file=args.snp,
            gff_file=args.gff,
            genome_file=args.genome,
            left=args.left,
            right=args.right,
            promoter=args.promoter,
            output_prefix=args.output,
            gffread_path=args.gffread_path,
            seqkit_path=args.seqkit_path
        )

        # 验证配置|Validate configuration
        config.validate()

        # 初始化日志|Initialize logger
        log_file = f"{config.output_prefix}.log"
        logger_manager = SnpRegionLogger(log_file=log_file)
        logger = logger_manager.get_logger()

        # 打印配置信息|Print configuration
        logger.info("="*60)
        logger.info("SNP区域基因提取工具|SNP Region Gene Extractor Tool")
        logger.info("="*60)
        logger.info(f"SNP文件|SNP file: {config.snp_file}")
        logger.info(f"GFF文件|GFF file: {config.gff_file}")
        logger.info(f"基因组文件|Genome file: {config.genome_file}")
        logger.info(f"上游距离|Left distance: {config.left} bp")
        logger.info(f"下游距离|Right distance: {config.right} bp")
        logger.info(f"启动子距离|Promoter distance: {config.promoter} bp")
        logger.info(f"输出前缀|Output prefix: {config.output_prefix}")
        logger.info("="*60)

        # 步骤1: 解析SNP文件|Step 1: Parse SNP file
        logger.info("步骤1|Step 1: 解析SNP文件|Parsing SNP file")
        snp_list = parse_snp_file(config.snp_file, logger)
        logger.info(f"找到|Found {len(snp_list)} 个唯一SNP|unique SNPs")

        # 步骤2: 解析GFF3文件|Step 2: Parse GFF3 file
        logger.info("步骤2|Step 2: 解析GFF3文件|Parsing GFF3 file")
        processor = SnpRegionProcessor(config, logger)
        processor.parse_gff3()

        # 步骤3: 提取所有序列|Step 3: Extract all sequences
        logger.info("步骤3|Step 3: 提取所有CDS和蛋白序列|Extracting all CDS and protein sequences")
        if not processor.extract_all_sequences():
            logger.error("序列提取失败|Sequence extraction failed")
            sys.exit(1)

        # 步骤4: 查找每个SNP的相关基因|Step 4: Find related genes for each SNP
        logger.info("步骤4|Step 4: 查找SNP相关基因|Finding related genes for SNPs")
        all_matches = []  # 存储所有匹配结果
        all_mrna_ids = set()  # 收集所有唯一的mRNA ID

        gene_list_lines = []  # gene_list.txt的行内容
        gene_list_lines.append("SNP_Position\tChromosome\tSNP_Pos\tGene_ID\tmRNA_ID\tStrand\tFeature\tDistance")

        for snp_idx, (chrom, pos) in enumerate(snp_list, 1):
            logger.debug(f"处理SNP|Processing SNP {snp_idx}/{len(snp_list)}: {chrom}:{pos}")

            matches = processor.find_genes_for_snp(chrom, pos)

            if matches:
                for match in matches:
                    snp_id = f"{chrom}:{pos}"
                    mrna_id = match['mrna_id']

                    all_matches.append({
                        'snp_id': snp_id,
                        'chrom': chrom,
                        'pos': pos,
                        'gene_id': match['gene_id'],
                        'mrna_id': mrna_id,
                        'strand': match['strand'],
                        'feature': match['feature'],
                        'distance': match['distance']
                    })

                    all_mrna_ids.add(mrna_id)

                    # 添加到gene_list行|Add to gene_list line
                    line = f"{snp_id}\t{chrom}\t{pos}\t{match['gene_id']}\t{mrna_id}\t{match['strand']}\t{match['feature']}\t{match['distance']}"
                    gene_list_lines.append(line)
            else:
                logger.debug(f"SNP {chrom}:{pos} 未找到相关基因|No related genes found")

        logger.info(f"总共找到|Total matches found: {len(all_matches)} 个基因关联|gene associations")
        logger.info(f"唯一mRNA数量|Unique mRNA count: {len(all_mrna_ids)}")

        if not all_mrna_ids:
            logger.warning("没有找到任何相关基因|No related genes found for any SNP")
            logger.info("程序结束|Program finished")
            sys.exit(0)

        # 步骤5: 提取区域序列|Step 5: Extract region sequences
        logger.info("步骤5|Step 5: 提取区域基因序列|Extracting region gene sequences")
        mrna_id_list = list(all_mrna_ids)
        if not processor.extract_region_sequences(mrna_id_list):
            logger.error("区域序列提取失败|Region sequence extraction failed")
            if not args.keep_temp:
                processor.cleanup_temp_files()
            sys.exit(1)

        # 步骤6: 写入gene_list文件|Step 6: Write gene_list file
        logger.info("步骤6|Step 6: 写入基因列表文件|Writing gene list file")
        with open(config.gene_list_output, 'w') as f:
            f.write('\n'.join(gene_list_lines) + '\n')
        logger.info(f"基因列表文件已保存|Gene list file saved: {config.gene_list_output}")

        # 清理临时文件|Clean up temporary files
        if not args.keep_temp:
            logger.info("步骤7|Step 7: 清理临时文件|Cleaning up temporary files")
            processor.cleanup_temp_files()
        else:
            logger.info("步骤7|Step 7: 保留临时文件|Keeping temporary files")

        # 完成|Finished
        logger.info("="*60)
        logger.info("分析完成|Analysis completed successfully!")
        logger.info(f"输出文件|Output files:")
        logger.info(f"  CDS序列|CDS sequences: {config.cds_output}")
        logger.info(f"  蛋白序列|Protein sequences: {config.protein_output}")
        logger.info(f"  基因列表|Gene list: {config.gene_list_output}")
        logger.info(f"  日志文件|Log file: {log_file}")
        logger.info("="*60)

        sys.exit(0)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"程序执行出错|Program execution error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
