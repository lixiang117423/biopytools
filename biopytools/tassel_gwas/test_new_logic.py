#!/usr/bin/env python3
"""
测试新的自动识别表型逻辑
"""

from pathlib import Path

# 创建测试表型文件，格式类似你的文件
test_pheno = Path('./test_trait_new.txt')
with open(test_pheno, 'w') as f:
    f.write('<Trait>\tfeture_1\tfeture_2\n')
    f.write('DRR054136\t0.3381131038664059\t1.2563186801505348\n')
    f.write('ERR626208\t0.52842002593229\t1.3192560593039047\n')
    f.write('SRR123456\t0.4123456789012345\t1.3987654321098765\n')

print('测试表型文件格式:')
print('==================')
with open(test_pheno, 'r') as f:
    for i, line in enumerate(f, 1):
        print(f'第{i}行: {line.strip()}')

# 测试表型识别逻辑
print('\n表型识别测试:')
print('============')
with open(test_pheno, 'r') as f:
    header_line = f.readline().strip()
    columns = header_line.split('\t')
    total_cols = len(columns)

print(f'检测到总列数 | Detected total columns: {total_cols}')
print(f'表头信息 | Header info: {columns}')

if columns[0] != '<Trait>':
    print('⚠️ 第一列不是<Trait>')
else:
    print('✅ 第一列是<Trait>')

# 确定要处理的表型列范围 | Determine trait columns to process
trait_columns = list(range(1, total_cols))  # 从第2列开始
print(f'开始处理表型列 | Starting to process trait columns: {trait_columns}')
print(f'表型名称 | Trait names: {[columns[i] for i in trait_columns]}')

# 检查数据行数
data_rows = 0
with open(test_pheno, 'r') as f:
    next(f)  # 跳过表头
    for line in f:
        if line.strip() and not line.startswith('#'):
            data_rows += 1
            cols = len(line.strip().split('\t'))
            if cols != total_cols:
                print(f'❌ 数据行列数({cols})与表头列数({total_cols})不匹配')

print(f'数据行数 | Data rows: {data_rows}')
print(f'表型数量 | Traits: {total_cols - 1}')

# 清理
test_pheno.unlink()
print('\n✅ 测试完成，文件已清理')