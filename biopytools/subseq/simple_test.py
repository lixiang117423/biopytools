#!/usr/bin/env python3
"""
简单测试脚本 | Simple test script
"""

import sys
import tempfile
import os
from Bio import SeqIO

# 添加父目录到路径 | Add parent directory to path
sys.path.insert(0, '/share/org/YZWL/yzwl_lixg/software/biopytools/biopytools')

from subseq.core import SequenceExtractor
from subseq.utils import SubseqLogger

class TestConfig:
    def __init__(self):
        self.input_fasta = None
        self.output_fasta = None
        self.id_list_file = None
        self.keep_order = True

def create_test_files():
    """创建测试文件 | Create test files"""
    # 创建FASTA文件 | Create FASTA file
    fasta_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)

    test_sequences = [
        ('seq1', 'ATGCGATCGATCGATCGATCG'),
        ('seq2', 'GCTAGCTAGCTAGCTAGCTAG'),
        ('chr1_seq3', 'TTTTAAAAAAAAAAAAAAA'),
        ('gene1', 'ATGATGATGATGATGATGATG'),
        ('seq4', 'CCCCCCCCCCCCCCCCCCCC'),
    ]

    for seq_id, seq in test_sequences:
        fasta_file.write(f">{seq_id}\n{seq}\n")
    fasta_file.close()

    # 创建ID列表文件 | Create ID list file
    id_list_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    id_list_file.write("seq1\nchr1_seq3\ngene1\nnonexistent_id\n")
    id_list_file.close()

    # 输出文件路径 | Output file path
    output_file = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False).name

    return fasta_file.name, id_list_file.name, output_file

def test_extraction():
    """测试序列提取 | Test sequence extraction"""
    print("🧪 开始序列提取测试 | Starting sequence extraction test")

    # 创建测试文件 | Create test files
    fasta_file, id_list_file, output_file = create_test_files()

    try:
        # 配置 | Configuration
        config = TestConfig()
        config.input_fasta = fasta_file
        config.output_fasta = output_file
        config.id_list_file = id_list_file

        # 设置日志 | Setup logging
        log_dir = os.path.dirname(output_file)
        logger_manager = SubseqLogger(log_dir)
        logger = logger_manager.get_logger()

        # 创建提取器 | Create extractor
        extractor = SequenceExtractor(config, logger)

        # 执行提取 | Perform extraction
        success = extractor.extract_sequences_by_id_list()

        if success and os.path.exists(output_file):
            print("✅ 提取成功 | Extraction successful")

            # 检查输出结果 | Check output results
            extracted_sequences = list(SeqIO.parse(output_file, 'fasta'))
            print(f"📊 提取到 {len(extracted_sequences)} 条序列 | Extracted {len(extracted_sequences)} sequences")

            for seq_record in extracted_sequences:
                print(f"    📄 {seq_record.id}: {len(seq_record.seq)} bp")

            return True
        else:
            print("❌ 提取失败 | Extraction failed")
            return False

    except Exception as e:
        print(f"❌ 测试错误 | Test error: {e}")
        return False
    finally:
        # 清理临时文件 | Clean up temporary files
        for file in [fasta_file, id_list_file, output_file]:
            if os.path.exists(file):
                os.unlink(file)

if __name__ == '__main__':
    print("🚀 序列子集提取模块测试 | Sequence Subsequence Extraction Module Test")
    success = test_extraction()

    if success:
        print("\n🎉 测试通过！| Test passed!")
        sys.exit(0)
    else:
        print("\n❌ 测试失败 | Test failed")
        sys.exit(1)