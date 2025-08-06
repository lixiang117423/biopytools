#!/usr/bin/env python3
"""
Annovar基本使用示例 | Annovar Basic Usage Example
"""

from biopytools.annovar import AnnovarProcessor

def main():
    # 创建处理器 | Create processor
    processor = AnnovarProcessor(
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
