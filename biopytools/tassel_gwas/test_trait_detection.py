#!/usr/bin/env python3
"""
测试表型文件解析逻辑
"""

from pathlib import Path

# Create a test phenotype file similar to your format
test_pheno = Path('./test_trait.txt')
with open(test_pheno, 'w') as f:
    f.write('<Trait>\tfeture_1\tfeture_2\n')
    f.write('DRR054136\t0.3381131038664059\t1.2563186801505348\n')
    f.write('ERR626208\t0.52842002593229\t1.3192560593039047\n')

# Test the trait detection logic
with open(test_pheno, 'r') as f:
    header_line = f.readline().strip()
    columns = header_line.split('\t')
    total_cols = len(columns)

print('检测到总列数 | Detected total columns:', total_cols)
print('表头信息 | Header info:', columns)

if columns[0] != '<Trait>':
    print('⚠️ 第一列不是<Trait>')
else:
    print('✅ 第一列是<Trait>')

# 确定要处理的表型列范围 | Determine trait columns to process
trait_columns = list(range(1, total_cols))  # 从第2列开始
print('开始处理表型列 | Starting to process trait columns:', trait_columns)
print('表型名称 | Trait names:', [columns[i] for i in trait_columns])

# Clean up
test_pheno.unlink()
print('✅ 测试完成，文件已清理')