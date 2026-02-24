#!/bin/bash

# 1. 定义路径配置
# 原始数据目录
RAW_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/raw/第一批下机数据/fastq"
# 过滤后数据目录
CLEAN_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/clean/第一批下机数据"
# 目标存放目录 (需要拷贝到的地方)
DEST_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/raw/第一批下机数据/fastq2"

# 2. 如果目标目录不存在，则创建
if [ ! -d "$DEST_DIR" ]; then
    echo "目标目录不存在，正在创建: $DEST_DIR"
    mkdir -p "$DEST_DIR"
fi

echo "开始扫描未过滤的文件..."

# 3. 遍历原始目录下的所有 .fq.gz 文件
for raw_file_path in "$RAW_DIR"/*.fq.gz; do
    # 获取文件名 (例如: W727A_1.fq.gz)
    filename=$(basename "$raw_file_path")
    
    # 构造对应的 clean 文件名
    # 逻辑: 将 .fq.gz 替换为 .clean.fq.gz
    # W727A_1.fq.gz -> W727A_1.clean.fq.gz
    clean_filename="${filename%.fq.gz}.clean.fq.gz"
    
    # 定义完整的 clean 文件路径用于检测
    clean_file_path="$CLEAN_DIR/$clean_filename"
    
    # 4. 检测是否存在对应的 clean 文件
    if [ -f "$clean_file_path" ]; then
        # 如果 clean 文件存在，跳过
        # echo "已过滤: $filename -> 跳过"
        :
    else
        # 如果 clean 文件不存在，说明需要拷贝
        echo "[发现未过滤样品] $filename"
        echo "  -> 正在拷贝到 fastq2 文件夹..."
        
        # 执行拷贝命令
        cp "$raw_file_path" "$DEST_DIR/"
        
        # 如果想节省空间，建议把上面那行 cp 改成下面这行软链接:
        # ln -s "$raw_file_path" "$DEST_DIR/$filename"
    fi
done

echo "所有操作完成。"

biopytools fastp \
    -i /share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/raw/第一批下机数据/fastq2 \
    -o /share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/clean/第一批下机数据2 \
    -t 12