"""
处理器模块 | Processor Module
"""

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from .parser import GFFParser
from .extractor import VariantExtractor
from .writer import OutputWriter

class GeneProcessor:
    """基因处理器 | Gene Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.gff_parser = GFFParser(logger)
        self.extractor = VariantExtractor(logger)
        self.writer = OutputWriter(logger)
    
    def process_single_gene(self, gene_id: str):
        """处理单个基因 | Process single gene"""
        self.logger.info(f"\n{'='*100}")
        self.logger.info(f"🧬 开始处理基因 | Processing gene: {gene_id}")
        self.logger.info(f"{'='*100}")
        
        # 1. 从GFF文件获取基因位置 | Get gene location from GFF
        chrom, gene_start, gene_end, strand = self.gff_parser.parse_gene_location(
            self.config.gff_file, gene_id
        )
        
        if chrom is None:
            self.logger.error(f"❌ 未能找到基因 {gene_id}，跳过 | Gene {gene_id} not found, skipping")
            return False
        
        # 2. 计算上下游区域 | Calculate upstream/downstream region
        region_start = max(1, gene_start - self.config.extend_bp)
        region_end = gene_end + self.config.extend_bp
        
        region_info = {
            'gene_id': gene_id,
            'chrom': chrom,
            'strand': strand,
            'gene_start': gene_start,
            'gene_end': gene_end,
            'region_start': region_start,
            'region_end': region_end,
            'extend': self.config.extend_bp
        }
        
        self.logger.info(f"📐 目标区域(±{self.config.extend_bp}bp) | Target region:")
        self.logger.info(f"   {chrom}:{region_start}-{region_end} (链:{strand})")
        self.logger.info(f"   总长度 | Total length: {region_end - region_start + 1} bp")
        
        # 3. 提取变异 | Extract variants
        self.logger.info(f"\n{'='*100}")
        self.logger.info("📊 开始提取变异 | Starting variant extraction")
        self.logger.info(f"{'='*100}")
        
        # 提取外显子区域变异 | Extract exonic variants
        exonic_variants = self.extractor.extract_variants_in_region(
            self.config.exonic_file, chrom, region_start, region_end,
            gene_start, gene_end, strand, is_exonic=True
        )
        
        # 提取所有变异 | Extract all variants
        all_variants = self.extractor.extract_variants_in_region(
            self.config.all_variant_file, chrom, region_start, region_end,
            gene_start, gene_end, strand, is_exonic=False
        )
        
        # 4. 输出文件 | Output files
        self.logger.info(f"\n{'='*100}")
        self.logger.info("💾 正在生成输出文件 | Generating output files")
        self.logger.info(f"{'='*100}")
        
        base_name = gene_id.replace('gene-', '')
        
        if self.config.output_format == 'excel':
            # 输出Excel文件（两个sheet）| Output Excel file (two sheets)
            excel_file = self.config.output_path / f"{base_name}_variants.xlsx"
            self.writer.write_excel_output(exonic_variants, all_variants, excel_file, region_info)
        else:
            # 输出TXT文件（两个独立文件）| Output TXT files (two separate files)
            exonic_txt = self.config.output_path / f"{base_name}_exonic_variants.txt"
            self.writer.write_txt_output(exonic_variants, exonic_txt, region_info, "外显子变异")
            
            all_txt = self.config.output_path / f"{base_name}_all_variants.txt"
            self.writer.write_txt_output(all_variants, all_txt, region_info, "所有变异")
        
        # 5. 输出汇总信息 | Output summary
        self.logger.info(f"\n{'='*100}")
        self.logger.info(f"✅ 基因 {gene_id} 处理完成 | Gene {gene_id} processing completed")
        self.logger.info(f"{'='*100}")
        self.logger.info(f"📊 统计信息 | Statistics:")
        self.logger.info(f"   外显子变异数 | Exonic variants: {len(exonic_variants)}")
        self.logger.info(f"   总变异数 | Total variants: {len(all_variants)}")
        
        return True
    
    def process_all_genes(self):
        """处理所有基因 | Process all genes"""
        total_genes = len(self.config.gene_ids)
        self.logger.info(f"\n🚀 开始批量处理 {total_genes} 个基因 | Starting batch processing of {total_genes} genes")
        self.logger.info(f"🧵 使用线程数 | Using threads: {self.config.threads}")
        
        success_count = 0
        failed_genes = []
        
        if self.config.threads == 1:
            # 单线程处理 | Single-threaded processing
            for gene_id in self.config.gene_ids:
                if self.process_single_gene(gene_id):
                    success_count += 1
                else:
                    failed_genes.append(gene_id)
        else:
            # 多线程处理 | Multi-threaded processing
            with ThreadPoolExecutor(max_workers=min(self.config.threads, total_genes)) as executor:
                future_to_gene = {executor.submit(self.process_single_gene, gene_id): gene_id 
                                 for gene_id in self.config.gene_ids}
                
                for future in as_completed(future_to_gene):
                    gene_id = future_to_gene[future]
                    try:
                        if future.result():
                            success_count += 1
                    except Exception as e:
                        self.logger.error(f"❌ 处理基因 {gene_id} 时发生错误 | Error processing gene {gene_id}: {e}")
                        failed_genes.append(gene_id)
        
        # 最终汇总 | Final summary
        self.logger.info(f"\n{'='*100}")
        self.logger.info("🎉 所有基因处理完成 | All genes processing completed")
        self.logger.info(f"{'='*100}")
        self.logger.info(f"✅ 成功处理 | Successfully processed: {success_count}/{total_genes}")
        if failed_genes:
            self.logger.warning(f"❌ 失败的基因 | Failed genes: {', '.join(failed_genes)}")
        self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
