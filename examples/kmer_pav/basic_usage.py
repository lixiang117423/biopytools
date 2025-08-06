#!/usr/bin/env python3
"""
KmerPav基本使用示例 | KmerPav Basic Usage Example
"""

from biopytools.kmer_pav import KmerPavProcessor

def main():
    # 创建处理器 | Create processor
    processor = KmerPavProcessor(
        input_dir="sample_data",
        output_dir="results",
        threads=4
    )
    
    # 运行处理 | Run processing
    processor.run_batch_processing()
    
    print("✅ 处理完成！| Processing completed!")
    print("结果保存在 results/ 目录 | Results saved in results/ directory")

if __name__ == "__main__":
    main()
