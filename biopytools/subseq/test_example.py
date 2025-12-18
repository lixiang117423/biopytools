"""
序列子集提取模块测试示例 | Sequence Subsequence Extraction Module Test Example
"""

import os
import tempfile
from Bio import SeqIO
from pathlib import Path

def create_test_data():
    """创建测试数据 | Create test data"""
    # 创建测试FASTA文件 | Create test FASTA file
    test_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)

    test_sequences = [
        ('seq1', 'ATGCGATCGATCGATCGATCG'),
        ('seq2', 'GCTAGCTAGCTAGCTAGCTAG'),
        ('chr1_seq3', 'TTTTAAAAAAAAAAAAAAA'),
        ('chr2_seq4', 'CCCCGGGGCCCCGGGGCCCC'),
        ('gene1', 'ATGATGATGATGATGATGATG'),
        ('gene2', 'GCTGCTGCTGCTGCTGCTGCT'),
        ('long_seq7', 'A' * 2000),  # 长序列
        ('short_seq8', 'ATGC'),    # 短序列
    ]

    for seq_id, seq in test_sequences:
        test_fasta.write(f">{seq_id}\n{seq}\n")

    test_fasta.close()

    # 创建ID列表文件 | Create ID list file
    test_id_list = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    test_id_list.write("seq1\nchr1_seq3\ngene1\nnonexistent_id\n")
    test_id_list.close()

    return test_fasta.name, test_id_list.name

def test_id_list_extraction():
    """测试ID列表提取 | Test ID list extraction"""
    print("🧪 测试ID列表提取功能 | Testing ID list extraction functionality")

    # 创建测试数据 | Create test data
    fasta_file, id_list_file = create_test_data()
    output_file = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False).name

    try:
        # 导入并测试 | Import and test
        from .main import Config, main
        from .utils import SubseqLogger

        # 配置 | Configuration
        config = Config()
        config.input_fasta = fasta_file
        config.output_fasta = output_file
        config.id_list_file = id_list_file
        config.keep_order = True

        # 设置日志 | Setup logging
        logger_manager = SubseqLogger(Path(output_file).parent)
        logger = logger_manager.get_logger()

        # 导入提取器 | Import extractor
        from .core import SequenceExtractor
        extractor = SequenceExtractor(config, logger)

        # 执行提取 | Perform extraction
        success = extractor.extract_sequences_by_id_list()

        if success and os.path.exists(output_file):
            # 检查输出结果 | Check output results
            extracted_sequences = list(SeqIO.parse(output_file, 'fasta'))
            print(f"✅ 成功提取 {len(extracted_sequences)} 条序列 | Successfully extracted {len(extracted_sequences)} sequences")

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

def test_pattern_extraction():
    """测试模式匹配提取 | Test pattern matching extraction"""
    print("\n🧪 测试模式匹配提取功能 | Testing pattern matching extraction functionality")

    # 创建测试数据 | Create test data
    fasta_file, _ = create_test_data()
    output_file = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False).name

    try:
        from .main import Config
        from .utils import SubseqLogger
        from .core import SequenceExtractor

        # 配置 | Configuration
        config = Config()
        config.input_fasta = fasta_file
        config.output_fasta = output_file
        config.pattern = 'chr1_'
        config.pattern_type = 'startswith'
        config.case_sensitive = True

        # 设置日志 | Setup logging
        logger_manager = SubseqLogger(Path(output_file).parent)
        logger = logger_manager.get_logger()

        # 执行提取 | Perform extraction
        extractor = SequenceExtractor(config, logger)
        success = extractor.extract_sequences_by_pattern()

        if success and os.path.exists(output_file):
            # 检查输出结果 | Check output results
            extracted_sequences = list(SeqIO.parse(output_file, 'fasta'))
            print(f"✅ 成功提取 {len(extracted_sequences)} 条匹配序列 | Successfully extracted {len(extracted_sequences)} matching sequences")

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
        for file in [fasta_file, output_file]:
            if os.path.exists(file):
                os.unlink(file)

def test_length_extraction():
    """测试长度筛选提取 | Test length filtering extraction"""
    print("\n🧪 测试长度筛选提取功能 | Testing length filtering extraction functionality")

    # 创建测试数据 | Create test data
    fasta_file, _ = create_test_data()
    output_file = tempfile.NamedTemporaryFile(suffix='.fasta', delete=False).name

    try:
        from .main import Config
        from .utils import SubseqLogger
        from .core import SequenceExtractor

        # 配置 | Configuration
        config = Config()
        config.input_fasta = fasta_file
        config.output_fasta = output_file
        config.min_length = 1000
        config.max_length = 3000

        # 设置日志 | Setup logging
        logger_manager = SubseqLogger(Path(output_file).parent)
        logger = logger_manager.get_logger()

        # 执行提取 | Perform extraction
        extractor = SequenceExtractor(config, logger)
        success = extractor.extract_sequences_by_length()

        if success and os.path.exists(output_file):
            # 检查输出结果 | Check output results
            extracted_sequences = list(SeqIO.parse(output_file, 'fasta'))
            print(f"✅ 成功提取 {len(extracted_sequences)} 条符合长度范围的序列 | Successfully extracted {len(extracted_sequences)} sequences within length range")

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
        for file in [fasta_file, output_file]:
            if os.path.exists(file):
                os.unlink(file)

def run_tests():
    """运行所有测试 | Run all tests"""
    print("🚀 开始运行序列子集提取模块测试 | Starting sequence subsequence extraction module tests")

    tests = [
        ("ID列表提取 | ID list extraction", test_id_list_extraction),
        ("模式匹配提取 | Pattern matching extraction", test_pattern_extraction),
        ("长度筛选提取 | Length filtering extraction", test_length_extraction),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ 测试 '{test_name}' 发生异常 | Test '{test_name}' threw exception: {e}")
            results.append((test_name, False))

    # 总结 | Summary
    print("\n📊 测试结果总结 | Test Results Summary:")
    passed = 0
    total = len(results)

    for test_name, result in results:
        status = "✅ 通过 | PASSED" if result else "❌ 失败 | FAILED"
        print(f"    {test_name}: {status}")
        if result:
            passed += 1

    print(f"\n🎯 总计 | Total: {passed}/{total} 测试通过 | tests passed")

    if passed == total:
        print("🎉 所有测试都通过了！| All tests passed!")
        return True
    else:
        print("⚠️ 部分测试失败，请检查实现 | Some tests failed, please check implementation")
        return False

if __name__ == '__main__':
    run_tests()