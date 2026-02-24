# =============================================================================
#  05-utils.zsh - 通用工具函数模块
#  Utilities Functions Module
# =============================================================================

# -----------------------------------------------------------------------------
#  安全删除系统 (Safe Delete System)
# -----------------------------------------------------------------------------
# function rm() {
#   if [ $# -eq 0 ]; then
#     truerm
#     return
#   fi

#   if ! command -v trash-put &> /dev/null; then
#     echo "❌ 错误: 'trash-put' 命令未找到。"
#     echo "   请先安装 'trash-cli' (e.g., pip install trash-cli)。"
#     printf "是否要继续使用系统的 'rm' 命令进行永久删除? (y/N) "
#     read -r response
#     if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
#       echo "⚠️ 警告: 将执行永久删除操作！"
#       truerm "$@"
#     else
#       echo "操作已取消。"
#     fi
#     return
#   fi

#   echo "🤔 您想如何处理以下文件/目录?"
#   for item in "$@"; do
#     echo "   - \"$item\""
#   done

#   printf "请选择: [1]永久删除, [2]移入回收站, [3]取消 > "

#   if [ -n "$ZSH_VERSION" ]; then
#     read -r -k 1 choice
#   else
#     read -r -n 1 choice
#   fi
#   echo

#   case "$choice" in
#     1)
#       echo "🔥 正在永久删除..."
#       truerm -rf "$@"
#       if [ $? -eq 0 ]; then
#         echo "✅ 成功：文件/目录已永久删除。"
#       else
#         echo "❌ 错误：永久删除失败。"
#       fi
#       ;;
#     2)
#       echo "♻️ 正在移动到回收站..."
#       trash-put "$@"
#       if [ $? -eq 0 ]; then
#         echo "✅ 成功：文件/目录已移动到回收站。"
#         echo "   - 查看回收站: trash-list"
#         echo "   - 恢复文件:   trash-restore"
#         echo "   - 清空回收站: trash-empty"
#       else
#         echo "❌ 错误：移动文件到回收站失败。"
#       fi
#       ;;
#     3 | *)
#       echo "🚫 操作已取消。"
#       ;;
#   esac
# }

function rm() {
  if [ $# -eq 0 ]; then
    truerm
    return
  fi

  if ! command -v trash-put &> /dev/null; then
    echo "❌ 错误: 'trash-put' 命令未找到。"
    echo "   请先安装 'trash-cli' (e.g., pip install trash-cli)。"
    printf "是否要继续使用系统的 'rm' 命令进行永久删除? (y/N) "
    read -r response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
      echo "⚠️ 警告: 将执行永久删除操作！"
      truerm "$@"
    else
      echo "操作已取消。"
    fi
    return
  fi

  echo "🤔 您想如何处理以下文件/目录?"
  for item in "$@"; do
    echo "   - \"$item\""
  done

  printf "请选择: [1]永久删除, [2]移入回收站, [3]取消 > "
  
  # 修改这里：移除单字符读取参数，允许用户输入完整内容并按回车
  read -r choice

  case "$choice" in
    1)
      echo "🔥 正在永久删除..."
      truerm -rf "$@"
      if [ $? -eq 0 ]; then
        echo "✅ 成功：文件/目录已永久删除。"
      else
        echo "❌ 错误：永久删除失败。"
      fi
      ;;
    2)
      echo "♻️  正在移动到回收站..."
      trash-put "$@"
      if [ $? -eq 0 ]; then
        echo "✅ 成功：文件/目录已移动到回收站。"
        echo "   - 查看回收站: trash-list"
        echo "   - 恢复文件:   trash-restore"
        echo "   - 清空回收站: trash-empty"
      else
        echo "❌ 错误：移动文件到回收站失败。"
      fi
      ;;
    3 | *)
      echo "🚫 操作已取消。"
      ;;
  esac
}

# -----------------------------------------------------------------------------
#  实用工具函数 (Utility Functions)
# -----------------------------------------------------------------------------

# 不区分大小写的grep
grepi() {
  grep -i "$@"
}

# 创建目录并立即进入
mkcd() {
  mkdir -p "$1" && cd "$1"
}

# 简化mamba环境激活
mma() {
  if [ -z "$1" ]; then
    echo "用法: mma <环境名称>"
    mamba env list
    return 1
  fi
  conda activate "$1"
}

# =============================================================================
# 函数：mmca (Mamba Create and Activate)
# 功能：创建一个新的 mamba 环境并立即激活它。
# 用法：
#   mmca <env_name>
#   mmca <env_name> python=3.10 numpy pandas
# =============================================================================
function mmca() {
  # 检查是否提供了环境名称
  if [ -z "$1" ]; then
    echo "用法错误: 请提供一个环境名称。"
    echo "例如: mmca my_env"
    return 1 # 返回错误码
  fi

  # 将第一个参数（环境名称）保存到变量中
  local env_name=$1

  echo "--> 正在创建 Mamba 环境: $env_name..."
  # 使用 "$@" 将所有参数（环境名和包名）传递给 mamba create
  # 使用 -y 自动确认安装
  # 使用 && 确保只有在创建成功后才执行激活命令
  if mamba create -n "$@" -y && mamba activate "$env_name"; then
    echo "--> 环境 '$env_name' 创建并激活成功！"
  else
    echo "--> 操作失败。请检查上面的错误信息。"
    return 1 # 返回错误码
  fi
}

# 智能解压函数 v2.1 - 修复了 shell 兼容性和文件格式识别问题
function x() {
    # 默认设置和参数解析
    local keep_original=true
    local target_dir=""
    local files=()
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -k|--keep)
                keep_original=true
                shift
                ;;
            -r|--remove|--delete)
                keep_original=false
                shift
                ;;
            -h|--help)
                cat << 'EOF'
用法: x <压缩文件> [目标目录] [选项]

解压各种格式的压缩文件，默认保留原文件

参数:
  <压缩文件>     要解压的文件路径
  [目标目录]     解压到指定目录（可选）

选项:
  -k, --keep     保留原压缩文件（默认行为）
  -r, --remove   解压后删除原压缩文件
  -h, --help     显示此帮助信息

支持的格式:
  tar系列: tar.gz, tgz, tar.bz2, tbz2, tar.xz, txz, tar.Z, tar.lzma, tar.lz4, tar.zst, tar
  常用格式: zip, rar, 7z
  单文件: gz, bz2, xz, Z, lzma, lz4, zst
  生物信息: fastq.gz, fasta.gz, gff3.gz
  系统包: deb, rpm
EOF
                return 0
                ;;
            -*)
                echo "❌ 未知选项: $1" >&2
                echo "使用 x --help 查看帮助信息" >&2
                return 1
                ;;
            *)
                if [ -f "$1" ]; then
                    files+=("$1")
                elif [ -z "$target_dir" ] && [ -d "$1" ]; then
                    target_dir="$1"
                else
                    files+=("$1")
                fi
                shift
                ;;
        esac
    done
    
    if [ ${#files[@]} -eq 0 ]; then
        echo "❌ 错误: 请提供要解压的文件" >&2
        echo "使用 x --help 查看帮助信息" >&2
        return 1
    fi
    
    local success_count=0
    local total_files=${#files[@]}

    for file in "${files[@]}"; do
        echo ""
        echo "📦 处理文件 ($((success_count + 1))/$total_files): $(basename "$file")"
        echo "===================================================================================="

        # 手动构建绝对路径
        local file_abs
        if [[ "$file" = /* ]]; then
            file_abs="$file"
        else
            file_abs="$PWD/$file"
        fi

        # 检查文件是否存在
        if [ ! -f "$file_abs" ]; then
            echo "❌ 错误: 文件不存在 - '$file'" >&2
            continue
        fi

        local original_dir
        original_dir=$(pwd)

        # 获取文件所在目录
        local file_dir
        file_dir=$(dirname "$file_abs")

        local work_dir="$original_dir"
        if [ -n "$target_dir" ]; then
            # 如果用户指定了目标目录，使用指定的目录
            mkdir -p "$target_dir"
            if [[ "$target_dir" = /* ]]; then
                work_dir="$target_dir"
            else
                work_dir="$PWD/$target_dir"
            fi
        else
            # 否则使用原文件所在的目录
            work_dir="$file_dir"
        fi
        
        echo "📂 解压位置: $work_dir"
        
        cd "$work_dir" || {
            echo "❌ 错误: 无法切换到目录 '$work_dir'" >&2
            continue
        }
        
        local extract_success=false
        # FIX: 使用更兼容的方式转换为小写
        local file_lower
        file_lower=$(basename "$file_abs" | tr '[:upper:]' '[:lower:]')

        # 根据文件扩展名进行匹配，优先匹配长扩展名
        case "$file_lower" in
            # --- 生物信息学常用格式 ---
            *.fastq.gz|*.fq.gz)
                echo "🧬 检测到 FASTQ 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.gz}")
                gzip -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.fasta.gz|*.fa.gz|*.fas.gz)
                echo "🧬 检测到 FASTA 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.gz}")
                gzip -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.gff3.gz|*.gff.gz|*.gtf.gz)
                echo "🧬 检测到 GFF/GTF 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.gz}")
                gzip -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            # --- tar 系列格式 ---
            *.tar.gz|*.tgz)
                echo "📦 检测到 tar.gz 压缩包..."
                tar -xzf "$file_abs" && extract_success=true
                ;;
            *.tar.bz2|*.tbz2)
                echo "📦 检测到 tar.bz2 压缩包..."
                tar -xjf "$file_abs" && extract_success=true
                ;;
            *.tar.xz|*.txz)
                echo "📦 检测到 tar.xz 压缩包..."
                tar -xJf "$file_abs" && extract_success=true
                ;;
            *.tar.zst)
                echo "📦 检测到 tar.zst 压缩包..."
                zstd -dc "$file_abs" | tar -xf - && extract_success=true
                ;;
            *.tar.z)
                echo "📦 检测到 tar.Z 压缩包..."
                tar -xZf "$file_abs" && extract_success=true
                ;;
            *.tar.lzma)
                echo "📦 检测到 tar.lzma 压缩包..."
                lzma -dc "$file_abs" | tar -xf - && extract_success=true
                ;;
            *.tar.lz4)
                echo "📦 检测到 tar.lz4 压缩包..."
                lz4 -dc "$file_abs" | tar -xf - && extract_success=true
                ;;
            *.tar)
                echo "📦 检测到 tar 文档..."
                tar -xf "$file_abs" && extract_success=true
                ;;
            # --- 单文件压缩格式 ---
            *.gz)
                echo "🗜️  检测到 gzip 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.gz}")
                gzip -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.bz2)
                echo "🗜️  检测到 bzip2 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.bz2}")
                bzip2 -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.xz)
                echo "🗜️  检测到 xz 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.xz}")
                xz -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.z)
                echo "🗜️  检测到 compress 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.Z}")
                uncompress -c "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.lzma)
                echo "🗜️  检测到 lzma 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.lzma}")
                lzma -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.lz4)
                echo "🗜️  检测到 lz4 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.lz4}")
                lz4 -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            *.zst)
                echo "🗜️  检测到 zstd 压缩文件..."
                local output_file
                output_file=$(basename "${file_abs%.zst}")
                zstd -dc "$file_abs" > "$output_file" && extract_success=true
                ;;
            # --- 常用压缩包格式 ---
            *.zip)
                echo "📂 检测到 ZIP 压缩包..."
                unzip -o "$file_abs" && extract_success=true
                ;;
            *.rar)
                echo "📂 检测到 RAR 压缩包..."
                unrar x -o+ "$file_abs" && extract_success=true
                ;;
            *.7z)
                echo "📂 检测到 7z 压缩包..."
                7z x -o"$work_dir" "$file_abs" && extract_success=true
                ;;
            # --- 系统包格式 ---
            *.deb)
                echo "📦 检测到 DEB 软件包..."
                ar x "$file_abs" && extract_success=true
                ;;
            *.rpm)
                echo "📦 检测到 RPM 软件包..."
                rpm2cpio "$file_abs" | cpio -idmv && extract_success=true
                ;;
            *)
                echo "❌ 错误: 不支持的文件格式 '$(basename "$file_abs")'" >&2
                echo "💡 支持的格式请使用 'x --help' 查看"
                cd "$original_dir"
                continue
                ;;
        esac
        
        if [ "$extract_success" = true ]; then
            echo "✅ 解压成功: $(basename "$file_abs")"
            if [ "$keep_original" = false ]; then
                echo "🗑️  删除原文件: $(basename "$file_abs")"
                command rm -f "$file_abs"
            else
                echo "💾 保留原文件: $(basename "$file_abs")"
            fi
            ((success_count++))
        else
            echo "❌ 解压失败: $(basename "$file_abs")" >&2
        fi
        
        cd "$original_dir"
    done
    
    echo ""
    echo "🎉 处理完成！"
    echo "📊 成功解压: $success_count/$total_files 个文件"
    
    return $((total_files - success_count))
}

# 历史命令检索 - 搜索 shell 内存历史
function hgg() {
    local pattern="$1"
    local count="${2:-10}"

    # 只使用内存历史，但增加搜索范围
    history 1 | grep "$pattern" | tail -n "$count" | while read -r line; do
        # 清理格式，只显示命令部分
        echo "$line" | sed -E 's/^ *[0-9]+\*? *//'
    done
}

# 历史命令检索 - 搜索 ~/history_commands.txt
function hg() {
    if [ -z "$1" ]; then
        echo "用法: hg <搜索关键词> [显示条数]"
        echo "示例: hg ls      # 搜索包含 ls 的命令"
        echo "      hg git 20  # 搜索包含 git 的命令，显示 20 条"
        return 1
    fi

    local pattern="$1"
    local count="${2:-20}"

    if [ ! -f "$LOG_COMMANDS_FILE" ]; then
        echo "日志文件不存在: $LOG_COMMANDS_FILE"
        return 1
    fi

    # 搜索并显示后 N 条匹配记录（带颜色），最新的在最下面
    grep -i "$pattern" "$LOG_COMMANDS_FILE" | tail -n "$count" | _colorize_log
}

# hg 函数的补全功能
# 支持补全当前目录文件名
_hg() {
    local curcontext="$curcontext" state line
    typeset -A opt_args

    _arguments -C \
        '1: :->first_arg' \
        '2: :->second_arg'

    case $state in
        first_arg)
            # 第一个参数：既可以是文件名，也可以是搜索关键词
            # 提供当前目录文件补全
            _files
            ;;
        second_arg)
            # 第二个参数：数字（显示条数）
            _numbers
            ;;
    esac
}

# 注册补全函数
compdef _hg hg

# git命令
function gacp() {
    git add "$1"
    git commit -m "$2"
    git push origin main
}

# grep ">"
function fh() {
  # 如果第一个参数是 -h 或 --help，显示帮助信息
  if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    # 使用 'heredoc' 格式化输出帮助文本
    cat <<EOF
Usage: fahead [-c | --count] [file1.fa file2.fa ...]

Description:
  Quickly finds and displays FASTA headers (lines starting with '>')
  from specified files or from standard input if no files are given.

Options:
  -c, --count   Instead of printing the headers, print the count of headers.
  -h, --help    Show this help message.

Examples:
  # Show headers from a single file
  fahead my_genome.fa

  # Count headers in a file
  fahead -c my_genome.fa

  # Read from a pipe and count headers
  cat *.fa | fahead -c
EOF
    return 0 # 成功退出函数
  fi

  # 如果第一个参数是 -c 或 --count，则进行计数
  if [[ "$1" == "-c" || "$1" == "--count" ]]; then
    shift # 移除 -c 参数，剩下的 $@ 就是文件名了
    # 对剩下的文件执行 grep 并通过管道传给 wc -l 来计数
    grep ">" "$@" | wc -l
  else
    # 默认行为：直接 grep
    grep ">" "$@"
  fi
}

# fh 函数的补全功能
_fh() {
  local curcontext="$curcontext" state line
  typeset -A opt_args

  _arguments -C \
    '1: :->first_arg' \
    '2: :->rest_args'

  case $state in
    first_arg)
      # 第一个参数可以是选项或文件名
      _alternative \
        'options:option:(-c --count -h --help)' \
        'files:fasta file:_files -g "*.fa(-.)"' \
        'files:fasta file:_files -g "*.fasta(-.)"' \
        'files:fasta file:_files -g "*.fna(-.)"' \
        'files:fasta file:_files -g "*.ffn(-.)"' \
        'files:fasta file:_files -g "*.faa(-.)"' \
        'files:fasta file:_files -g "*.frn(-.)"'
      ;;
    rest_args)
      # 后续参数只补全文件名
      _files -g "*.fa(-.)"
      _files -g "*.fasta(-.)"
      _files -g "*.fna(-.)"
      _files -g "*.ffn(-.)"
      _files -g "*.faa(-.)"
      _files -g "*.frn(-.)"
      ;;
  esac
}

# 注册补全函数
compdef _fh fh

function lnbin() {
    # 检查是否提供了至少一个参数
    if [[ $# -lt 1 ]]; then
        echo "${RED}错误: 请提供一个要链接的程序名。${RESET}"
        echo "用法: lnbin <程序名> [可选的目标目录，默认为 ~/bin]"
        return 1
    fi

    local program_name="$1"
    # 使用参数扩展来设置默认的目标目录为 ~/bin
    # 如果提供了第二个参数，则使用它；否则，使用 $HOME/bin
    local dest_dir="${2:-$HOME/.local/bin}"
    
    # --- 步骤 1: 查找程序的完整路径 ---
    echo "${BLUE}--> 正在查找 '${program_name}'...${RESET}"
    local source_path
    source_path=$(command -v "$program_name")

    # 检查 'command -v' 是否成功找到路径
    if [[ -z "$source_path" ]]; then
        echo "${RED}错误: 在你的 PATH 中找不到程序 '${program_name}'。${RESET}"
        return 1
    fi
    echo "${GREEN}✔ 找到程序: ${source_path}${RESET}"

    # --- 步骤 2: 创建符号链接 ---
    # 确保目标目录存在，-p 选项表示如果不存在则创建，且不会因已存在而报错
    mkdir -p "$dest_dir"

    local dest_path="$dest_dir/$program_name"
    echo "${BLUE}--> 准备在 ${dest_path} 创建链接...${RESET}"

    # 检查目标位置是否已经存在文件或链接
    if [[ -e "$dest_path" ]]; then
        echo "${YELLOW}警告: '${dest_path}' 已经存在。跳过创建。${RESET}"
        return 0
    fi

    # 执行链接命令并检查结果
    if ln -s "$source_path" "$dest_path"; then
        echo "${GREEN}✅ 成功创建符号链接！${RESET}"
        echo "   链接: ${dest_path}"
        echo "   指向: ${source_path}"
    else
        echo "${RED}错误: 创建符号链接失败。${RESET}"
        echo "请检查你是否有 '${dest_dir}' 目录的写入权限。"
        echo "如果目标是系统目录 (如 /usr/local/bin)，你可能需要使用 'sudo'。"
        return 1
    fi
}

# 从 off-line 服务器同步文件的 rsync 快捷函数
# 用法: rsync_offline <远程文件路径> [本地目标路径]
function copy() {
  # 检查是否提供了至少一个参数（远程文件路径）
  if [[ -z "$1" ]]; then
    echo "错误: 请提供远程文件的路径。"
    echo "用法: rsync_offline <远程文件路径> [本地目标路径]"
    echo "示例: rsync_offline /share/org/YZWL/.../file.txt ./"
    return 1
  fi

  # 第一个参数是远程路径
  local remote_path="$1"
  # 第二个参数是本地目标路径，如果未提供，则默认为当前目录 "."
  local local_dest="${2:-.}"

  echo "==> 从 off-line:${remote_path} 同步到 ${local_dest}"
  # 执行 rsync 命令，并添加 --progress 以显示进度
  rsync -avz --partial --progress "off-line:${remote_path}" "${local_dest}"
}

# 函数: xg (xz | grep)
# 用法: xg <要搜索的关键词> [grep的其他参数]
# 示例: xg psoja
#       xg TTTG -i
xg() {
  # 检查是否提供了搜索关键词
  if [ -z "$1" ]; then
    echo "错误: 请提供一个搜索关键词。"
    echo "用法: xg <要搜索的关键词>"
    return 1
  fi
  
  # 第一个参数作为搜索模式
  local pattern="$1"
  
  # "$@" 会将所有其他参数 (如 -i, --color) 传递给 grep
  shift
  
  # 执行核心命令
  xz | grep --color=auto "$pattern" "$@"
}

cpp() {
  # 检查是否提供了文件名作为参数
  if [ -z "$1" ]; then
    echo "用法: cprp <文件或目录名>"
    return 1
  fi
  
  # 检查文件或目录是否存在
  if [ ! -e "$1" ]; then
    echo "错误: '$1' 不存在。"
    return 1
  fi
  
  # 获取绝对路径
  abs_path=$(realpath "$1")
  
  # 尝试多种剪贴板工具
  if command -v xclip >/dev/null 2>&1; then
    # Linux with X11
    echo -n "$abs_path" | xclip -selection clipboard
    echo "✅ 已复制到剪贴板 (xclip): $abs_path"
  elif command -v xsel >/dev/null 2>&1; then
    # Linux with X11 (alternative)
    echo -n "$abs_path" | xsel --clipboard
    echo "✅ 已复制到剪贴板 (xsel): $abs_path"
  elif command -v pbcopy >/dev/null 2>&1; then
    # macOS
    echo -n "$abs_path" | pbcopy
    echo "✅ 已复制到剪贴板 (pbcopy): $abs_path"
  elif command -v clip.exe >/dev/null 2>&1; then
    # WSL (Windows Subsystem for Linux)
    echo -n "$abs_path" | clip.exe
    echo "✅ 已复制到剪贴板 (clip.exe): $abs_path"
  else
    # 如果没有剪贴板工具，至少输出路径
    echo "⚠️  警告: 未找到剪贴板工具 (xclip/xsel/pbcopy/clip.exe)"
    echo "📋 路径已输出到终端: $abs_path"
    # 可选：将路径保存到临时文件
    echo "$abs_path" > /tmp/cprp_last_path.txt
    echo "💾 路径已保存到: /tmp/cprp_last_path.txt"
  fi
}

myzip() {
    # 1. 获取用户输入的文件名，并去除末尾可能存在的 "/" (比如文件夹补全时会有/)
    local target="${1%/}"
    
    # 2. 检查输入是否存在
    if [ ! -e "$1" ]; then
        echo "错误: 文件或目录 '$1' 不存在。"
        return 1
    fi

    # 3. 执行压缩命令
    # 格式: tar -czvf 文件名.tar.gz 文件名
    echo "正在压缩: ${target} -> ${target}.tar.gz ..."
    tar -czvf "${target}.tar.gz" "$1"
    
    echo "完成。"
}

# singularity拉取镜像
pull() {
    if [ -z "$1" ]; then
        echo "错误: 请提供镜像名称 (例如: pull aryeelab/hicpro)"
        return 1
    fi
    
    echo "正在通过 docker.1ms.run 镜像源拉取: $1 ..."
    singularity pull "docker://docker.1ms.run/$1"
}

# 复制绝对路径到剪切板
kk() {
    local path=$(realpath "$@")
    echo "$path"
    printf "\033]52;c;$(printf "%s" "$path" | ${MINIFORGE3_DIR:-$HOME/miniforge3}/bin/python3 -c "import sys, base64; print(base64.b64encode(sys.stdin.buffer.read()).decode(), end='')")\a"
}

# 增强 pwd 命令：显示当前路径并自动复制到剪贴板
pwd() {
    local path=$(builtin pwd)
    echo "$path"
    printf "\033]52;c;$(printf "%s" "$path" | ${MINIFORGE3_DIR:-$HOME/miniforge3}/bin/python3 -c "import sys, base64; print(base64.b64encode(sys.stdin.buffer.read()).decode(), end='')")\a"
}

# -----------------------------------------------------------------------------
#  命令历史记录 (Command History Logger)
# -----------------------------------------------------------------------------
# 自动记录命令执行的时间、目录和命令内容到 ~/history_commands.txt

LOG_COMMANDS_FILE="$HOME/history_commands.txt"

# ANSI 颜色代码
COLOR_TIME="\033[32m"      # 绿色 - 时间
COLOR_PATH="\033[34m"      # 蓝色 - 路径
COLOR_CMD="\033[33m"       # 黄色 - 命令
COLOR_SEP="\033[90m"       # 灰色 - 分隔符
COLOR_RESET="\033[0m"      # 重置颜色

# 确保日志文件存在
touch "$LOG_COMMANDS_FILE" 2>/dev/null || true

# zsh 使用 preexec 钩子
if [ -n "$ZSH_VERSION" ]; then
    # 记录命令执行的函数（追加模式，更高效）
    log_command_zsh() {
        local cmd="$1"
        # 排除空命令和一些不需要记录的命令
        if [ -n "$cmd" ]; then
            local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
            # 路径简化：用 ~ 替代 $HOME
            local display_path="${PWD/#$HOME/~}"
            # 纯文本格式存储 - 追加到文件开头（使用临时文件）
            local new_line="[${timestamp}] [${display_path}] → ${cmd}"

            # 更高效的实现：直接追加到文件末尾
            # 如需新记录在前，可定期使用 sort/tac 命令处理
            echo "$new_line" >> "$LOG_COMMANDS_FILE"
        fi
    }

    # 添加 preexec 钩子
    autoload -U add-zsh-hook
    add-zsh-hook preexec log_command_zsh
fi

# 给日志添加颜色（内部函数）
# 纯文本格式: [时间] [路径] → 命令
_colorize_log() {
    perl -pe "s/\[([^]]+)\] \[([^]]+)\] → (.*)/${COLOR_SEP}[${COLOR_TIME}\1${COLOR_SEP}] [${COLOR_PATH}\2${COLOR_SEP}] ${COLOR_SEP}→${COLOR_RESET} ${COLOR_CMD}\3${COLOR_RESET}/"
}

# 查看历史命令的快捷函数（带颜色高亮）
hlog() {
    if [ -f "$LOG_COMMANDS_FILE" ]; then
        # 显示前 N 条记录（默认 20 条，最新记录在前）
        local count="${1:-20}"
        head -n "$count" "$LOG_COMMANDS_FILE" | _colorize_log
    else
        echo "日志文件不存在: $LOG_COMMANDS_FILE"
    fi
}

# 搜索历史命令
hloggrep() {
    if [ -z "$1" ]; then
        echo "用法: hloggrep <搜索关键词>"
        return 1
    fi
    if [ -f "$LOG_COMMANDS_FILE" ]; then
        grep -i "$1" "$LOG_COMMANDS_FILE" | _colorize_log
    else
        echo "日志文件不存在: $LOG_COMMANDS_FILE"
    fi
}

# -----------------------------------------------------------------------------
#  模块加载成功标记
# -----------------------------------------------------------------------------
export ZSH_MODULE_UTILS_LOADED=1
