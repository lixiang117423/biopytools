# 📥 从CNCB获取测序数据下载链接模块

**批量获取CNCB数据库测序数据下载链接的专业工具 | Professional Tool for Batch Downloading Sequencing Data Links from CNCB Database**

## 📖 功能概述 | Overview

从CNCB获取测序数据下载链接模块是一个专业的生物信息学工具，用于从国家基因组科学数据中心（CNCB）的FTP服务器批量查找和获取测序数据的下载链接。支持SRA、ERA、DRA三大国际数据库，具备智能路径搜索、缓存优化、错误重试等功能，可大大提高数据获取效率。

## ✨ 主要特性 | Key Features

- **🌐 多数据库支持**: 支持SRA、ERA、DRA数据库的Run ID查询
- **⚡ 智能搜索**: 基于缓存机制的路径搜索算法，显著提升查询效率
- **🔄 重试机制**: 内置FTP连接重试和错误恢复机制
- **📊 批量处理**: 支持大规模项目的批量Run ID查询
- **📜 自动脚本**: 自动生成wget下载脚本，支持断点续传
- **📝 详细日志**: 完整的处理过程记录和错误追踪
- **🛡️ 质量控制**: 输入文件验证和输出结果校验
- **🎯 灵活配置**: 丰富的参数配置，适应不同使用场景

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本链接查询
python -m biopytools.get_link_from_CNCB.main projects.txt

# 指定输出文件
python -m biopytools.get_link_from_CNCB.main -i projects.txt -o my_links.txt

# 详细输出模式
python -m biopytools.get_link_from_CNCB.main projects.txt -v
```

### 高级用法 | Advanced Usage

```python
# 直接使用模块
from biopytools.get_link_from_CNCB.main import CNCLinkExtractor

extractor = CNCLinkExtractor(
    input_file="projects.txt",
    output_file="links.txt",
    failed_file="failed.txt",
    verbose=True,
    generate_download_script=True
)

success = extractor.run_extraction()
```

## 📋 命令行参数 | Command Line Parameters

### 必需参数 | Required Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `input_file` | `None` | 📁 输入文件路径（ProjectID和RunID，Tab分隔） |

### 输出配置参数 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `[input]_links.txt` | 📄 输出链接文件路径 |
| `-f, --failed` | `[input]_failed.txt` | ❌ 失败记录文件路径 |
| `--download-script` | `download.sh` | 📜 下载脚本文件名 |

### FTP配置参数 | FTP Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--ftp-host` | `download2.cncb.ac.cn` | 🌐 FTP服务器地址 |
| `--ftp-timeout` | `60` | ⏱️ FTP连接超时时间(秒) |

### 性能配置参数 | Performance Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--retry-attempts` | `3` | 🔄 FTP连接重试次数 |
| `--threads` | `4` | 🧵 并发线程数 |

### 日志配置参数 | Logging Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-v, --verbose` | `False` | 📝 详细输出模式 |
| `--log-file` | `None` | 📊 日志文件路径 |

### 输出选项参数 | Output Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--no-download-script` | `False` | ⏭️ 不生成下载脚本 |
| `--no-executable` | `False` | 🔒 不设置脚本执行权限 |

### 处理步骤说明 | Processing Steps

| 步骤 | 名称 | 描述 |
|------|------|------|
| **1** | 📖 文件验证 | 验证输入文件格式和可读性 |
| **2** | 🌐 FTP连接 | 连接CNCB FTP服务器并测试权限 |
| **3** | 🔍 路径搜索 | 智能搜索Run ID对应的FTP路径 |
| **4** | 📄 文件查找 | 查找可用的数据文件 |
| **5** | 💾 结果保存 | 保存链接和失败记录 |
| **6** | 📜 脚本生成 | 生成自动下载脚本 |
| **7** | 📊 报告生成 | 生成详细处理报告 |

## 📁 输入文件格式 | Input File Format

### 项目Run ID文件 | Project Run ID File

Tab分隔的两列格式文件：

```
ProjectID    RunID
PRJNA123456  SRR12345678
PRJNA123456  SRR12345679
PRJNA789012  ERR12345678
PRJNA789012  DRR12345678
PRJNA345678  SRR34567890
# 注释行会被忽略
```

**文件要求**:
- UTF-8编码格式
- Tab分隔的两列数据
- 第一列为ProjectID，第二列为Run ID
- 支持注释行（以#开头）
- 每个Run ID占据独立一行

**支持的数据前缀**:
- `SRR*`: SRA (Sequence Read Archive)
- `ERR*`: ERA (European Nucleotide Archive)
- `DRR*`: DRA (DNA Data Bank of Japan)

## 📂 输出文件说明 | Output Files Description

### 下载链接文件 | Download Links File

包含所有找到的数据文件下载链接：

```
ftp://download2.cncb.ac.cn/INSDC/SRA/SRR123/SRR12345678/SRR12345678.sra
ftp://download2.cncb.ac.cn/INSDC/SRA/SRR123/SRR12345679/SRR12345679_1.fastq.gz
ftp://download2.cncb.ac.cn/INSDC/SRA/SRR123/SRR12345679/SRR12345679_2.fastq.gz
```

### 失败记录文件 | Failed Records File

记录未能找到文件的Run ID：

```
PRJNA789012	ERR12345679
PRJNA345678	SRR34567891
```

### 自动下载脚本 | Auto Download Script

自动生成的wget下载脚本：

```bash
#!/bin/bash
# Auto-generated download script for CNCB data
# Generated on: 2024-12-17 20:30:45
# Total URLs: 125

wget -c 'ftp://download2.cncb.ac.cn/INSDC/SRA/SRR123/SRR12345678/SRR12345678.sra'
wget -c 'ftp://download2.cncb.ac.cn/INSDC/SRA/SRR123/SRR12345679/SRR12345679_1.fastq.gz'
wget -c 'ftp://download2.cncb.ac.cn/INSDC/SRA/SRR123/SRR12345679/SRR12345679_2.fastq.gz'

echo 'Download script completed'
```

### 详细报告文件 | Detailed Report File

包含完整的处理统计信息：

```
CNCB数据链接提取报告 | CNCB Data Link Extraction Report
============================================================

生成时间 | Generated Time: 2024-12-17 20:30:45

📊 项目统计 | Project Statistics:
------------------------------
总项目数 | Total Projects: 12
总Run ID数 | Total Run IDs: 1,245

✅ 成功统计 | Success Statistics:
------------------------------
成功链接数 | Successful URLs: 1,180
成功率 | Success Rate: 94.8%

❌ 失败统计 | Failure Statistics:
------------------------------
失败ID数 | Failed IDs: 65

失败的Run IDs列表 | List of Failed Run IDs:
PRJNA789012	ERR12345679
PRJNA345678	SRR34567891
...
```

## 🔧 使用示例 | Usage Examples

### 小规模数据获取

```bash
# 获取单个项目的测序数据
python -m biopytools.get_link_from_CNCB.main single_project.txt

# 输出示例：
# 📄 single_project_links.txt - 包含所有下载链接
# 📜 download.sh - 自动下载脚本
# ❌ single_project_failed.txt - 失败记录（如有）
```

### 大规模批量查询

```bash
# 大规模项目批量查询
python -m biopytools.get_link_from_CNCB.main \
    large_project_list.txt \
    -o comprehensive_links.txt \
    -f comprehensive_failed.txt \
    --download-script comprehensive_download.sh \
    --threads 8 \
    --ftp-timeout 120 \
    --retry-attempts 5 \
    -v \
    --log-file comprehensive_process.log

# 使用生成的下载脚本
nohup bash comprehensive_download.sh > download.log 2>&1 &
```

### 高级配置使用

```python
# Python脚本方式使用
from biopytools.get_link_from_CNCB.main import CNCLinkExtractor

# 配置参数
config = {
    'input_file': 'custom_projects.txt',
    'output_file': 'custom_links.txt',
    'failed_file': 'custom_failed.txt',
    'download_script': 'custom_download.sh',
    'ftp_host': 'download2.cncb.ac.cn',
    'ftp_timeout': 90,
    'retry_attempts': 4,
    'max_threads': 6,
    'verbose': True,
    'log_file': 'cncb_process.log',
    'generate_download_script': True,
    'script_executable': True
}

# 创建并运行提取器
extractor = CNCLinkExtractor(**config)
success = extractor.run_extraction()

if success:
    print("✅ 链接提取成功完成！")
    print("📊 统计信息:")
    print(f"  - 成功链接数: {len(extractor.success_urls)}")
    print(f"  - 失败ID数: {len(extractor.failed_ids)}")
else:
    print("❌ 链接提取失败")
```

### 条件查询和过滤

```bash
# 预先过滤Run ID（可选）
awk '$2 ~ /^SRR/' all_projects.txt > sra_only.txt

# 分批处理大量数据
split -l 1000 huge_project_list.txt batch_
for batch in batch_*; do
    python -m biopytools.get_link_from_CNCB.main "$batch" \
        -o "${batch}_links.txt" \
        -f "${batch}_failed.txt" \
        --download-script "${batch}_download.sh"
done

# 合并结果
cat batch_*_links.txt > all_links.txt
cat batch_*_failed.txt > all_failed.txt
```

## ⚡ 性能优化 | Performance Optimization

### 网络配置优化

```bash
# 增加超时时间（网络较慢时）
python -m biopytools.get_link_from_CNCB.main projects.txt --ftp-timeout 180

# 增加重试次数（不稳定网络）
python -m biopytools.get_link_from_CNCB.main projects.txt --retry-attempts 8

# 减少并发（避免过载）
python -m biopytools.get_link_from_CNCB.main projects.txt --threads 2
```

### 数据预筛选

```python
# 预先验证Run ID格式
def validate_run_ids(input_file, output_file):
    valid_lines = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    project_id, run_id = parts
                    if run_id[:3] in ['SRR', 'ERR', 'DRR']:
                        valid_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(valid_lines)

# 使用验证后的文件
validate_run_ids('raw_projects.txt', 'validated_projects.txt')
python -m biopytools.get_link_from_CNCB.main validated_projects.txt
```

### 缓存优化

```python
# 缓存命中统计示例
from biopytools.get_link_from_CNCB.utils import PathCache

# 创建缓存（内部自动使用）
cache = PathCache()
cache.set('SRR12345', '/INSDC/SRA/SRR123')

# 查看缓存统计
stats = cache.get_stats()
print(f"缓存命中率 | Cache hit rate: {stats['hit_rate_percent']}%")
```

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**网络连接问题**
```bash
# 测试FTP连接
ftp download2.cncb.ac.cn

# 检查防火墙设置
telnet download2.cncb.ac.cn 21

# 增加超时时间和重试次数
python -m biopytools.get_link_from_CNCB.main projects.txt \
    --ftp-timeout 180 \
    --retry-attempts 10
```

**找不到文件问题**
```bash
# 检查Run ID格式
grep -E "^(SRR|ERR|DRR)" projects.txt | wc -l

# 验证Run ID是否存在
curl -s "https://www.ncbi.nlm.nih.gov/sra/SRR12345678"

# 使用详细模式查看搜索过程
python -m biopytools.get_link_from_CNCB.main projects.txt -v
```

**权限问题**
```bash
# 检查输出目录权限
ls -la output_directory/

# 手动创建输出目录
mkdir -p output_directory
chmod 755 output_directory
```

### 性能问题 | Performance Issues

**处理速度慢**
```python
# 减少并发线程数（避免服务器限制）
python -m biopytools.get_link_from_CNCB.main projects.txt --threads 2

# 分批处理大量数据
split -l 500 large_list.txt batch_
for batch in batch_*; do
    python -m biopytools.get_link_from_CNCB.main "$batch" &
done
wait
```

**内存使用过高**
```python
# 减少缓存大小（通过配置）
python -m biopytools.get_link_from_CNCB.main projects.txt --threads 1

# 监控内存使用
python -c "
import psutil
import time
while True:
    print(f'Memory: {psutil.virtual_memory().percent}%')
    time.sleep(5)
" &
```

### 调试技巧 | Debugging Tips

```python
# 启用详细日志
python -m biopytools.get_link_from_CNCB.main projects.txt \
    -v \
    --log-file debug.log

# 检查日志文件中的错误
grep "ERROR\|WARNING" debug.log

# 验证输出文件
head -10 projects_links.txt
wc -l projects_links.txt

# 测试生成的下载脚本
bash -n download.sh  # 检查语法
```

## 📊 质量控制建议 | Quality Control Recommendations

### 输入数据验证

```bash
# 1. 格式验证
python -c "
from biopytools.get_link_from_CNCB.utils import InputFileParser
valid, msg = InputFileParser.validate_input_file('projects.txt')
print(f'验证结果 | Validation: {valid}')
print(f'消息 | Message: {msg}')
"

# 2. 内容统计
awk -F'\t' '!/^#/ && NF==2 {print $2}' projects.txt | sort | uniq -c | sort -nr

# 3. 前缀分布
grep -E "^\S+\t(SRR|ERR|DRR)" projects.txt | cut -f2 | cut -c1-3 | sort | uniq -c
```

### 结果验证

```bash
# 1. 链接有效性检查
python -c "
import urllib.request
import ssl

def check_ftp_url(url):
    try:
        context = ssl._create_unverified_context()
        with urllib.request.urlopen(url, timeout=10, context=context) as response:
            return response.status == 200
    except:
        return False

# 检查前10个链接
with open('projects_links.txt') as f:
    urls = f.readlines()[:10]

for url in urls:
    url = url.strip()
    if url:
        valid = check_ftp_url(url)
        print(f'{\"✅\" if valid else \"❌\"} {url}')
"

# 2. 文件大小统计
python -c "
import ftplib
import re

def get_file_size(url):
    match = re.search(r'ftp://[^/]+(.+)', url)
    if match:
        file_path = match.group(1)
        try:
            ftp = ftplib.FTP('download2.cncb.ac.cn')
            ftp.login()
            size = ftp.size(file_path)
            ftp.quit()
            return size
        except:
            return None
    return None

# 统计文件大小
total_size = 0
with open('projects_links.txt') as f:
    for url in f:
        size = get_file_size(url.strip())
        if size:
            total_size += size

print(f'总文件大小 | Total size: {total_size / 1024**3:.2f} GB')
"
```

## 🔗 相关文档 | Related Documentation

- [CNCB官方网站](https://www.cncb.ac.cn/)
- [SRA数据库文档](https://www.ncbi.nlm.nih.gov/sra/docs/)
- [ERA数据库文档](https://www.ebi.ac.uk/ena/browser/)
- [DRA数据库文档](https://ddbj.nig.ac.jp/)
- [FTP协议说明](https://tools.ietf.org/html/rfc959)
- [GNU wget手册](https://www.gnu.org/software/wget/manual/)
- [biopytools其他模块](../README.md)

## 📄 许可证 | License

本模块遵循MIT许可证。详细信息请参见LICENSE文件。

## 🤝 贡献指南 | Contributing

欢迎提交Issue和Pull Request来改进本模块。

## 📞 技术支持 | Support

如有技术问题，请联系：
- 邮箱: yzwl_lixg@outlook.com
- 项目地址: https://github.com/your-org/biopytools

---

**最后更新**: 2024年12月17日
**版本**: 1.0.0
**作者**: biopytools开发团队