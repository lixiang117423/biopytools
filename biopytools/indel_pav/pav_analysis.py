"""
PAV分析核心模块 | PAV Analysis Core Module
"""

from typing import List, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
from .data_processing import VCFProcessor

class PAVAnalyzer:
    """PAV分析器 | PAV Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.vcf_processor = VCFProcessor(config, logger)
        self.results = []
    
    def run_pav_analysis(self):
        """运行PAV分析 | Run PAV analysis"""
        self.logger.info("🧬 开始INDEL PAV分析 | Starting INDEL PAV analysis")
        
        # 解析VCF头部
        samples = self.vcf_processor.parse_vcf_header()
        
        # 处理变异
        self.logger.info("🔄 处理INDEL变异 | Processing INDEL variants")
        
        variant_count = 0
        
        # 使用生成器处理大文件
        for variant in self.vcf_processor.process_vcf_variants():
            self.results.append(variant)
            variant_count += 1
            
            if variant_count % 1000 == 0:
                self.logger.info(f"💾 已收集 {variant_count} 个合格INDEL | Collected {variant_count} qualified INDELs")
        
        self.logger.info(f"✅ PAV分析完成，共处理 {len(self.results)} 个INDEL | PAV analysis completed, processed {len(self.results)} INDELs")
        
        # 统计信息
        self.log_statistics()
        
        return samples, self.results
    
    def log_statistics(self):
        """记录统计信息 | Log statistics"""
        if not self.results:
            self.logger.warning("⚠️  没有找到合格的INDEL | No qualified INDELs found")
            return
        
        # 统计INDEL类型
        insertions = sum(1 for r in self.results if r['type'] == 'insertion')
        deletions = sum(1 for r in self.results if r['type'] == 'deletion')
        
        # 统计长度分布
        lengths = [r['length'] for r in self.results]
        avg_length = sum(lengths) / len(lengths)
        max_length = max(lengths)
        min_length = min(lengths)
        
        # 统计染色体分布
        chroms = {}
        for r in self.results:
            chrom = r['chrom']
            chroms[chrom] = chroms.get(chrom, 0) + 1
        
        self.logger.info("📊 INDEL统计信息 | INDEL Statistics:")
        self.logger.info(f"   🔸 插入 | Insertions: {insertions} ({insertions/len(self.results)*100:.1f}%)")
        self.logger.info(f"   🔹 删除 | Deletions: {deletions} ({deletions/len(self.results)*100:.1f}%)")
        self.logger.info(f"   📏 平均长度 | Average length: {avg_length:.1f} bp")
        self.logger.info(f"   📏 长度范围 | Length range: {min_length}-{max_length} bp")
        self.logger.info(f"   🧬 染色体数量 | Chromosome count: {len(chroms)}")
        
        # 显示前5个染色体的INDEL数量
        top_chroms = sorted(chroms.items(), key=lambda x: x[1], reverse=True)[:5]
        self.logger.info(f"   🏆 前5染色体 | Top 5 chromosomes: {', '.join([f'{c}({n})' for c, n in top_chroms])}")

class PAVFilter:
    """PAV结果过滤器 | PAV Result Filter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def filter_results(self, results: List[Dict]) -> List[Dict]:
        """过滤PAV结果 | Filter PAV results"""
        if not results:
            return results
        
        self.logger.info("🔍 应用PAV结果过滤 | Applying PAV result filters")
        
        original_count = len(results)
        filtered_results = []
        
        for result in results:
            # 可以在这里添加更多过滤条件
            # 例如：基于频率、区域等的过滤
            
            # 计算变异频率
            genotypes = result['genotypes']
            alt_count = sum(genotypes)
            frequency = alt_count / len(genotypes) if genotypes else 0
            
            # 这里可以添加频率过滤
            # if frequency < min_frequency or frequency > max_frequency:
            #     continue
            
            filtered_results.append(result)
        
        filtered_count = len(filtered_results)
        removed_count = original_count - filtered_count
        
        if removed_count > 0:
            self.logger.info(f"🚫 过滤移除 {removed_count} 个INDEL | Filtered out {removed_count} INDELs")
        
        self.logger.info(f"✅ 过滤后保留 {filtered_count} 个INDEL | {filtered_count} INDELs retained after filtering")
        
        return filtered_results
