# =============================================================================
#  09-bio-tools.zsh - 生物信息学工具模块
#  Bioinformatics Tools Module
# =============================================================================

# =============================================================================
# [Aspera SRA/EBI 下载助手函数] Aspera SRA/EBI Download Helper
# =============================================================================
#
# 这个函数简化了从 EBI SRA 数据库使用 Aspera 进行下载的过程。
#
# 用法:
#   sra_download <远程文件路径> [本地存放目录]
#
# 参数:
#   $1 (必需): 远程服务器上的文件路径。
#              例如: /vol1/fastq/SRR305/000/SRR30556600/SRR30556600_2.fastq.gz
#   $2 (可选): 文件要下载到的本地目录。如果省略，默认为当前目录 (./)。
#
# 示例:
#   # 下载文件到当前目录
#   sra_download /vol1/fastq/SRR305/000/SRR30556600/SRR30556600_2.fastq.gz
#
#   # 下载文件到指定目录 ~/data/
#   sra_download /vol1/fastq/SRR305/000/SRR30556600/SRR30556600_2.fastq.gz ~/data/
#
sra_dl() {
  # 检查是否提供了必需的远程路径参数
  if [[ -z "$1" ]]; then
    echo "错误: 请提供远程文件路径。"
    echo "用法: sra_dl <远程文件路径> [本地存放目录]"
    return 1
  fi

  # 定义固定的部分 - 使用统一路径配置
  local aspera_key="${ASPERA_KEY:-~/miniforge3/envs/aspera_v.3.9.6/etc/asperaweb_id_dsa.openssh}"
  local base_ascp_cmd="ascp -T -k 1 -l 300m -P 33001 -i ${aspera_key}"
  local server_info="${ASPERA_SERVER:-era-fasp@fasp.sra.ebi.ac.uk:}"

  # 定义可变的部分
  local remote_file_path="$1"
  # 如果提供了第二个参数（本地目录），则使用它；否则，默认为当前目录 "."
  local local_destination_path="${2:-.}"

  # 拼接并执行最终命令
  echo "==> 准备下载: ${remote_file_path}"
  echo "==> 保存到: ${local_destination_path}"
  
  # 执行命令。使用引号确保带空格的路径也能正常工作。
  eval "${base_ascp_cmd} ${server_info}${remote_file_path} \"${local_destination_path}\""
}

# GTDB-TK 数据库路径已在 00-path-config.zsh 中统一配置
# 默认路径: $HOME/database/gtdbtk/release226
# export GTDBTK_DATA_PATH=$HOME/database/gtdbtk/release226

# 模块加载成功标记
export ZSH_MODULE_BIO_TOOLS_LOADED=1