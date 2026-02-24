# #!/bin/bash

# # 定义输入和输出目录
# SRA_DIR="/share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/sra"
# OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/raw"

# # 创建输出目录（如果不存在）
# mkdir -p ${OUTPUT_DIR}

# # 遍历所有.sra文件
# for sra_file in ${SRA_DIR}/*.sra
# do
#     # 获取文件名（不含路径和扩展名）
#     base_name=$(basename ${sra_file} .sra)
    
#     echo "处理文件: ${base_name}"
    
#     # 使用fastq-dump转换
#     # --split-3: 自动拆分双端测序文件
#     # --gzip: 压缩输出（可选，如不需要压缩可删除此参数）
#     # --outdir: 指定输出目录
#     fastq-dump --split-3 \
#                --gzip \
#                --outdir ${OUTPUT_DIR} \
#                ${sra_file}
    
#     echo "完成: ${base_name}"
# done

# echo "所有文件处理完成！"


biopytools sra2fastq \
    -i /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/sra \
    -o /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/raw