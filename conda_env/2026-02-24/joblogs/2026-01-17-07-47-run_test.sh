  # 1. 运行 BWA 测试（会自动分析多重比对）
  ./run_bwa_test.sh

  # 2. 运行 Python 详细分析（需要先运行BWA）
  python3 analyze_multiple_hits.py
