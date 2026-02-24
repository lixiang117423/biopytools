#!/bin/bash
# 完整的sub函数 - 带历史记录功能和脚本备份

# =============================================================================
#  内部通用函数：统一处理作业提交逻辑
# =============================================================================
_sub_common() {
  local DEFAULT_QUEUE="$1"
  local DEFAULT_CPUS="$2"
  local DEFAULT_MEM="$3"
  shift 3  # 移除前三个参数，剩下的传给逻辑

  local JOB_NAME="" QUEUE="$DEFAULT_QUEUE" CPUS="$DEFAULT_CPUS" MEM="$DEFAULT_MEM" HOSTS=1 GPUS=0
  local OUT_FILE="" ERR_FILE="" SCRIPT_TO_RUN=""
  
  # 🏠 全局日志配置
  local DEFAULT_LOG_DIR="${HOME}/joblogs"
  
  # 📝 保存原始参数用于历史记录
  local ORIGINAL_ARGS=("$@")
  
  local usage="
Usage: sub [job_name] [OPTIONS] /path/to/your/script.sh
   or: sub [OPTIONS] -j <job_name> /path/to/your/script.sh

  'csub' 命令的封装,简化作业提交。🚀

必填项: 📌
  job_name            作业名称 (位置参数,推荐用法)
  /path/to/script.sh  要运行的脚本路径 (最后一个参数)

选项: ⚙️
  -j <name>       作业名称 (选项参数,向后兼容)
  -q <queue>      队列名称 (默认: ${QUEUE})
  -n <number>     CPU 核心数 (默认: ${CPUS})
  -m <gigabytes>  内存大小 (GB) (默认: ${MEM})
  -g <number>     GPU 数量 (默认: ${GPUS})
  -h <number>     主机/节点数 (默认: ${HOSTS})
  -o <path>       标准输出文件 (默认: ${DEFAULT_LOG_DIR}/YYYY-MM-DD-JOB_NAME.out)
  -e <path>       标准错误文件 (默认: ${DEFAULT_LOG_DIR}/YYYY-MM-DD-JOB_NAME.err)
  -l <dir>        指定日志目录 (默认: ${DEFAULT_LOG_DIR})
  --local-log     使用当前目录的logs子目录
  --merge-log     合并输出到单个.log文件 (YYYY-MM-DD-JOB_NAME.log)
  --no-backup     不备份脚本文件
  --help          显示此帮助信息 ❓

日志管理: 📋
  默认日志位置: ${DEFAULT_LOG_DIR}/
  脚本备份位置: ${DEFAULT_LOG_DIR}/YYYY-MM-DD-HH-MM-script_name.sh
  查看日志: ls ${DEFAULT_LOG_DIR}/
  清理旧日志: rm ${DEFAULT_LOG_DIR}/2023-*
  查看提交历史: sub_logs history
  
  支持两种日志格式:
  - 分离模式: YYYY-MM-DD-jobname.out/.err (默认)
  - 合并模式: YYYY-MM-DD-jobname.log (--merge-log)

示例: 💡
  sub my_job script.sh                    # 日志: ~/joblogs/2023-10-27-my_job.out/err
                                          # 备份: ~/joblogs/2023-10-27-15-30-script.sh
  sub my_job --merge-log script.sh        # 日志: ~/joblogs/2023-10-27-my_job.log
  sub my_job --local-log script.sh        # 日志: ./logs/2023-10-27-my_job.out/err
  sub my_job --no-backup script.sh        # 不备份脚本文件
  sub my_job -l /data/logs script.sh      # 日志: /data/logs/2023-10-27-my_job.out/err
  sub my_job -o custom.out script.sh      # 自定义输出文件
"

  if [[ $# -eq 0 ]] || [[ "$1" == "--help" ]]; then
    echo "$usage"
    return 0
  fi

  local LOG_DIR="$DEFAULT_LOG_DIR"
  local USE_LOCAL_LOG=false
  local MERGE_LOG=false
  local NO_BACKUP=false

  # 解析非选项参数
  local args=()
  while [[ $# -gt 0 ]]; do
    case $1 in
      --local-log)
        USE_LOCAL_LOG=true
        shift
        ;;
      --merge-log)
        MERGE_LOG=true
        shift
        ;;
      --no-backup)
        NO_BACKUP=true
        shift
        ;;
      --help)
        echo "$usage"
        return 0
        ;;
      -*)
        args+=("$1")
        if [[ "$1" =~ ^-[jqnmhoelg]$ ]]; then
          shift
          args+=("$1")
        fi
        shift
        ;;
      *)
        args+=("$1")
        shift
        ;;
    esac
  done

  set -- "${args[@]}"

  # 处理位置参数中的作业名称
  if [[ $# -gt 0 && "$1" != -* ]]; then
    JOB_NAME="$1"
    shift
  fi

  # 处理选项
  while getopts ":j:q:n:m:h:o:e:g:l:" opt; do
    case ${opt} in
      j)
        if [[ -n "$JOB_NAME" ]]; then
          echo "⚠️ 警告: 作业名称已通过位置参数设置为 '$JOB_NAME'，忽略 -j 选项" >&2
        else
          JOB_NAME="$OPTARG"
        fi
        ;;
      q) QUEUE="$OPTARG" ;;
      n) CPUS="$OPTARG" ;;
      m) MEM="$OPTARG" ;;
      h) HOSTS="$OPTARG" ;;
      o) OUT_FILE="$OPTARG" ;;
      e) ERR_FILE="$OPTARG" ;;
      g) GPUS="$OPTARG" ;;
      l) LOG_DIR="$OPTARG" ;;
      \?) echo "❌ 无效选项: -$OPTARG" >&2; echo "$usage"; return 1 ;;
      :) echo "❌ 选项 -$OPTARG 需要一个参数" >&2; echo "$usage"; return 1 ;;
    esac
  done
  shift $((OPTIND -1))

  SCRIPT_TO_RUN="$1"
  
  # 验证必填参数
  if [[ -z "$JOB_NAME" ]]; then
    echo "❌ 错误: 必须指定作业名称。" >&2
    echo "  方式1: sub <job_name> <script>" >&2
    echo "  方式2: sub -j <job_name> <script>" >&2
    echo "$usage"
    return 1
  fi
  if [[ -z "$SCRIPT_TO_RUN" ]]; then
    echo "❌ 错误: 必须在最后指定要运行的脚本。" >&2
    echo "$usage"
    return 1
  fi

  # 验证脚本文件是否存在
  if [[ ! -f "$SCRIPT_TO_RUN" ]]; then
    echo "❌ 错误: 脚本文件不存在: $SCRIPT_TO_RUN" >&2
    return 1
  fi

  # 📅 获取当前日期和时间
  local TODAY
  # TODAY=$(date +'%Y-%m-%d')
  TODAY=$(date +'%Y-%m-%d-%H-%M')
  local SUBMIT_TIME
  SUBMIT_TIME=$(date +'%Y-%m-%d %H:%M:%S')
  local WORK_DIR
  WORK_DIR=$(pwd)

  # 🗂️ 决定日志目录
  if [[ "$USE_LOCAL_LOG" == true ]]; then
    LOG_DIR="./logs"
  fi

  # 📂 设置默认日志文件路径
  if [[ "$MERGE_LOG" == true ]]; then
    # 合并模式：使用单个.log文件
    if [[ -z "$OUT_FILE" ]]; then OUT_FILE="${LOG_DIR}/${TODAY}-${JOB_NAME}.log"; fi
    ERR_FILE="$OUT_FILE"  # 错误也输出到同一个文件
  else
    # 分离模式：使用.out/.err文件
    if [[ -z "$OUT_FILE" ]]; then OUT_FILE="${LOG_DIR}/${TODAY}-${JOB_NAME}.out"; fi
    if [[ -z "$ERR_FILE" ]]; then ERR_FILE="${LOG_DIR}/${TODAY}-${JOB_NAME}.err"; fi
  fi
  
  # 创建日志目录
  if ! mkdir -p "$LOG_DIR" 2>/dev/null; then
    echo "❌ 错误: 无法创建日志目录 $LOG_DIR" >&2
    echo "💡 建议: 使用 --local-log 或 -l 指定其他目录" >&2
    return 1
  fi

  # 🗃️ 准备脚本备份路径
  local SCRIPT_BACKUP=""
  if [[ "$NO_BACKUP" != true ]]; then
    local script_basename=$(basename "$SCRIPT_TO_RUN")
    SCRIPT_BACKUP="${LOG_DIR}/${TODAY}-${script_basename}"
  fi

  local R_STRING="rusage[mem=${MEM}G]"
  if [[ ${GPUS} -gt 0 ]]; then
      R_STRING+=",ngpus_excl_p=${GPUS}"
  fi
  R_STRING+=" span[hosts=${HOSTS}]"

  # 收集所有输出到变量
  local full_output=""

  full_output+="---> 🚀 提交作业,参数如下:
"
  full_output+="     📝 作业名称: ${JOB_NAME}
     🔄 队列    : ${QUEUE}
     💻 CPU核心 : ${CPUS}
     🧠 内存    : ${MEM}G
     🎮 GPU数量 : ${GPUS}
     🖥️  主机数 : ${HOSTS}
     📁 日志目录: ${LOG_DIR}/
"

  if [[ "$MERGE_LOG" == true ]]; then
    full_output+="    📄 合并日志: ${OUT_FILE}
"
  else
    full_output+="     📤 输出日志: ${OUT_FILE}
     📥 错误日志: ${ERR_FILE}
"
  fi

  if [[ "$NO_BACKUP" != true ]]; then
    full_output+="     📋 脚本备份: ${SCRIPT_BACKUP}
"
  fi

  full_output+="     📜 执行脚本: ${SCRIPT_TO_RUN}
"
  full_output+="--------------------------------------------------
"

  # 🚀 执行作业提交
  local csub_output
  csub_output=$(csub -J "${JOB_NAME}" -q "${QUEUE}" -n "${CPUS}" \
       -R "${R_STRING}" \
       -o "${OUT_FILE}" -e "${ERR_FILE}" "bash ${SCRIPT_TO_RUN}" 2>&1)
  local csub_exit_code=$?

  # 📋 备份脚本文件（仅在提交成功时）
  if [[ $csub_exit_code -eq 0 && "$NO_BACKUP" != true ]]; then
    if cp "$SCRIPT_TO_RUN" "$SCRIPT_BACKUP" 2>/dev/null; then
      full_output+="✅ 脚本已备份到: $SCRIPT_BACKUP
"
    else
      full_output+="⚠️ 警告: 无法备份脚本文件到 $SCRIPT_BACKUP
"
    fi
  elif [[ "$NO_BACKUP" == true ]]; then
    SCRIPT_BACKUP=""  # 清空备份路径，用于历史记录
  fi

  # 📝 记录提交历史（获取输出消息）
  local history_msg
  history_msg=$(record_submission_history "$SUBMIT_TIME" "$WORK_DIR" "$csub_exit_code" "$csub_output" \
    "$JOB_NAME" "$QUEUE" "$CPUS" "$MEM" "$GPUS" "$HOSTS" "$LOG_DIR" "$OUT_FILE" "$ERR_FILE" "$SCRIPT_TO_RUN" "$SCRIPT_BACKUP" "$MERGE_LOG" "${ORIGINAL_ARGS[@]}")
  full_output+="$history_msg
"

  # 添加csub输出
  if [[ -n "$csub_output" ]]; then
    full_output+="$csub_output
"
  fi

  # 统一输出所有信息
  echo "$full_output"

  # 复制所有输出到剪贴板
  printf "\033]52;c;$(printf "%s" "$full_output" | ${MINIFORGE3_DIR:-$HOME/miniforge3}/bin/python3 -c "import sys, base64; print(base64.b64encode(sys.stdin.buffer.read()).decode(), end='')")\a"

  return $csub_exit_code
}

# =============================================================================
#  公共接口函数 - 调用内部通用函数
# =============================================================================

# 默认 c02 队列 (CPUS=64, MEM=300G)
function sub() {
  _sub_common "c02" 64 300 "$@"
}

# c01 队列 (CPUS=88, MEM=500G)
function sub2c01() {
  _sub_common "c01" 88 500 "$@"
}

# c02 队列 (CPUS=64, MEM=300G) - 与 sub 相同
function sub2c02() {
  _sub_common "c02" 64 300 "$@"
}

# 📝 记录提交历史函数
function record_submission_history() {
  local submit_time="$1"
  local work_dir="$2"
  local exit_code="$3"
  local csub_output="$4"
  local job_name="$5"
  local queue="$6"
  local cpus="$7"
  local mem="$8"
  local gpus="$9"
  local hosts="${10}"
  local log_dir="${11}"
  local out_file="${12}"
  local err_file="${13}"
  local script_to_run="${14}"
  local script_backup="${15}"
  local merge_log="${16}"

  # 🔄 转换为绝对路径
  if [[ "$script_to_run" != /* ]]; then
    script_to_run="$work_dir/$script_to_run"
  fi

  shift 16
  local original_args=("$@")
  
  local history_file="${HOME}/submitted_jobs.txt"
  
  # 确定提交状态
  local submit_status="失败"
  if [[ $exit_code -eq 0 ]]; then
    submit_status="成功"
  fi
  
  # 尝试提取作业ID
  local job_id="未知"
  if [[ -n "$csub_output" ]]; then
    # 尝试从csub输出中提取作业ID (格式可能需要根据实际情况调整)
    job_id=$(echo "$csub_output" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+' | head -1)
    if [[ -z "$job_id" ]]; then
      job_id=$(echo "$csub_output" | grep -oE '[0-9]{6,}' | head -1)
    fi
    if [[ -z "$job_id" ]]; then
      job_id="未提取到"
    fi
  fi
  
  # 重构完整命令
  local full_command="sub"
  for arg in "${original_args[@]}"; do
    if [[ "$arg" =~ [[:space:]] ]]; then
      full_command+=" \"$arg\""
    else
      full_command+=" $arg"
    fi
  done
  
  # 创建历史文件目录（如果不存在）
  local history_dir=$(dirname "$history_file")
  mkdir -p "$history_dir" 2>/dev/null
  
  # 清理变量中的换行符和多余空白
  script_to_run=$(echo "$script_to_run" | tr -d '\n\r' | xargs)
  out_file=$(echo "$out_file" | tr -d '\n\r' | xargs)
  err_file=$(echo "$err_file" | tr -d '\n\r' | xargs)
  [[ -n "$script_backup" ]] && script_backup=$(echo "$script_backup" | tr -d '\n\r' | xargs)

  # 写入历史记录
  {
    echo "=============================================================================================="
    echo "⏰ 提交时间: $submit_time"
    echo "📂 工作目录: $work_dir"
    echo "🚀 提交作业参数如下:"
    echo "    📝 作业名称: ${job_name}"
    echo "    🔄 队列    : ${queue}"
    echo "    💻 CPU核心 : ${cpus}"
    echo "    🧠 内存    : ${mem}G"
    echo "    🎮 GPU数量 : ${gpus}"
    echo "    🖥️  主机数  : ${hosts}"
    echo "    📁 日志目录: ${log_dir}/"
    if [[ "$merge_log" == true ]]; then
      echo "    📄 合并日志: ${out_file}"
    else
      echo "    📤 输出日志: ${out_file}"
      echo "    📥 错误日志: ${err_file}"
    fi
    if [[ -n "$script_backup" ]]; then
      echo "    📋 脚本备份: ${script_backup}"
    fi
    echo "    📜 执行脚本: ${script_to_run}"
    echo "💻 完整命令: $full_command"
    echo "✅ 提交状态: $submit_status"
    echo "🆔 作业ID: $job_id"
    # echo ""
  } >> "$history_file"
  
  # 验证写入是否成功
  if [[ -f "$history_file" ]]; then
    echo "📝 已记录到提交历史: $history_file"
  else
    echo "⚠️ 警告: 无法写入提交历史文件: $history_file"
  fi
}

# 🛠️ 辅助函数：日志管理
function sub_logs() {
  local DEFAULT_LOG_DIR="${HOME}/joblogs"
  local action="$1"
  
  case "$action" in
    "list"|"ls")
      echo "📋 日志目录内容:"
      ls -la "$DEFAULT_LOG_DIR" 2>/dev/null || echo "日志目录不存在: $DEFAULT_LOG_DIR"
      ;;
    "today")
      local today=$(date +'%Y-%m-%d')
      echo "📅 今天的日志 ($today):"
      ls -la "$DEFAULT_LOG_DIR"/*"$today"*.{out,err,log,sh} 2>/dev/null || echo "今天还没有日志"
      ;;
    "clean")
      read -p "🗑️ 确定要清理7天前的日志吗? [y/N]: " -n 1 -r
      echo
      if [[ $REPLY =~ ^[Yy]$ ]]; then
        find "$DEFAULT_LOG_DIR" -name "20*-*-*-*.out" -mtime +7 -delete 2>/dev/null
        find "$DEFAULT_LOG_DIR" -name "20*-*-*-*.err" -mtime +7 -delete 2>/dev/null
        find "$DEFAULT_LOG_DIR" -name "20*-*-*-*.log" -mtime +7 -delete 2>/dev/null
        find "$DEFAULT_LOG_DIR" -name "20*-*-*-*.sh" -mtime +7 -delete 2>/dev/null
        echo "✅ 已清理7天前的日志和脚本备份"
      fi
      ;;
    "tail")
      local job_name="$2"
      if [[ -n "$job_name" ]]; then
        local today=$(date +'%Y-%m-%d')
        # 优先查找.log文件，然后是.out文件
        local log_file="$DEFAULT_LOG_DIR/$today-$job_name.log"
        if [[ ! -f "$log_file" ]]; then
          log_file="$DEFAULT_LOG_DIR/$today-$job_name.out"
        fi
        
        if [[ -f "$log_file" ]]; then
          echo "📖 正在查看: $log_file"
          tail -f "$log_file"
        else
          echo "❌ 日志文件不存在:"
          echo "   - $DEFAULT_LOG_DIR/$today-$job_name.log"
          echo "   - $DEFAULT_LOG_DIR/$today-$job_name.out"
        fi
      else
        echo "用法: sub_logs tail <job_name>"
      fi
      ;;
    "script")
      local job_name="$2"
      if [[ -n "$job_name" ]]; then
        local today=$(date +'%Y-%m-%d')
        echo "🔍 搜索脚本备份 (模糊匹配 '$job_name'):"
        local found_scripts=()
        while IFS= read -r -d '' script_file; do
          found_scripts+=("$script_file")
        done < <(find "$DEFAULT_LOG_DIR" -name "*$job_name*.sh" -print0 2>/dev/null)
        
        if [[ ${#found_scripts[@]} -eq 0 ]]; then
          echo "❌ 未找到包含 '$job_name' 的脚本备份"
        elif [[ ${#found_scripts[@]} -eq 1 ]]; then
          echo "📋 找到脚本: ${found_scripts[0]}"
          echo "📖 脚本内容:"
          echo "----------------------------------------"
          cat "${found_scripts[0]}"
        else
          echo "📋 找到多个脚本备份:"
          for i in "${!found_scripts[@]}"; do
            echo "  $((i+1)). $(basename "${found_scripts[$i]}")"
          done
          echo ""
          read -p "请选择要查看的脚本 [1-${#found_scripts[@]}]: " choice
          if [[ "$choice" =~ ^[0-9]+$ ]] && [[ $choice -ge 1 ]] && [[ $choice -le ${#found_scripts[@]} ]]; then
            local selected_script="${found_scripts[$((choice-1))]}"
            echo "📖 脚本内容: $selected_script"
            echo "----------------------------------------"
            cat "$selected_script"
          else
            echo "❌ 无效选择"
          fi
        fi
      else
        echo "用法: sub_logs script <job_name_pattern>"
      fi
      ;;
    "history")
      local count="$2"
      local history_file="${HOME}/submitted_jobs.txt"
      
      if [[ ! -f "$history_file" ]]; then
        echo "📝 还没有提交历史记录"
        return 0
      fi
      
      echo "📋 作业提交历史:"
      echo ""
      
      if [[ "$count" == "all" ]]; then
        cat "$history_file"
      elif [[ "$count" =~ ^[0-9]+$ ]]; then
        # 显示最近N条记录
        local total_blocks=$(grep -c "^==============" "$history_file")
        if [[ $total_blocks -le $count ]]; then
          cat "$history_file"
        else
          # 计算要跳过的块数
          local skip_blocks=$((total_blocks - count))
          awk -v skip="$skip_blocks" '
            /^==============/ { block_count++ }
            block_count > skip { print }
          ' "$history_file"
        fi
      else
        # 默认显示最近5条记录
        local total_blocks=$(grep -c "^==============" "$history_file")
        if [[ $total_blocks -le 5 ]]; then
          cat "$history_file"
        else
          local skip_blocks=$((total_blocks - 5))
          awk -v skip="$skip_blocks" '
            /^==============/ { block_count++ }
            block_count > skip { print }
          ' "$history_file"
        fi
      fi
      ;;
    "search")
      local keyword="$2"
      local history_file="${HOME}/submitted_jobs.txt"
      
      if [[ ! -f "$history_file" ]]; then
        echo "📝 还没有提交历史记录"
        return 0
      fi
      
      if [[ -z "$keyword" ]]; then
        echo "用法: sub_logs search <关键词>"
        return 1
      fi
      
      echo "🔍 搜索包含 '$keyword' 的提交记录:"
      echo ""
      
      # 使用awk搜索包含关键词的记录块
      awk -v keyword="$keyword" '
        /^==============/ { 
          if (found) print block
          block = $0 "\n"
          found = 0
          next
        }
        { 
          block = block $0 "\n"
          if (tolower($0) ~ tolower(keyword)) found = 1
        }
        END { if (found) print block }
      ' "$history_file"
      ;;
    *)
      echo "🔧 日志管理工具"
      echo "用法: sub_logs <action>"
      echo ""
      echo "动作:"
      echo "  list/ls     - 列出所有日志文件"
      echo "  today       - 显示今天的日志和脚本备份"
      echo "  clean       - 清理7天前的日志和脚本备份"
      echo "  tail <job>  - 实时查看指定作业的日志"
      echo "  script <job> - 查看指定作业的脚本备份"
      echo "  history [N] - 显示最近N条提交历史 (默认5条)"
      echo "  history all - 显示所有提交历史"
      echo "  search <词> - 搜索提交历史"
      echo ""
      echo "示例:"
      echo "  sub_logs today"
      echo "  sub_logs tail my_job"
      echo "  sub_logs script my_job"
      echo "  sub_logs history 10"
      echo "  sub_logs search interproscan"
      echo ""
      echo "💡 快速查看错误日志:"
      echo "  tail -f ~/joblogs/2023-10-27-my_job.err"
      echo ""
      echo "📝 提交历史文件位置: ~/submitted_jobs.txt"
      echo "📋 脚本备份位置: ~/joblogs/YYYY-MM-DD-HH-MM-script.sh"
      ;;
  esac
}

# 🚀 快速别名
alias sublogs='sub_logs'
alias subtoday='sub_logs today'
alias subhistory='sub_logs history'

#
# 卡片式/树状 cj 命令格式化函数 (终极定制版 V2)
#
# 功能:
# - 使用任务名称作为主标题，任务ID作为详细信息，更符合直觉。
# - 主标题颜色已更换为亮绿色。
# - 增加了丰富的 Emoji 和颜色，让界面生动、信息直观。
#
function cj_formatted() {
  local now
  now=$(date +%s)
  
  local cj_output
  cj_output=$(command cjobs "$@")
  
  if [ "$(echo "$cj_output" | sed '/^$/d' | wc -l)" -le 1 ]; then
      echo "✅ 没有正在运行或排队的任务。"
      return
  fi

  echo "$cj_output" | gawk -v now="$now" '
    BEGIN {
        # --- ANSI 颜色代码 (这是最终版) ---
        C_GREEN_B="\033[1;32m"  # 亮绿色 (粗体) - 用于主标题
        C_CYAN_B="\033[1;36m"   # 亮青色 (粗体) - 用于任务状态
        C_RESET="\033[0m"       # 重置所有颜色

        # 月份映射
        split("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec", m, " ");
        for (i=1; i<=12; i++) M[m[i]] = i;
    }

    # 跳过原始表头，处理数据行
    NR > 1 {
        # 计算运行时间
        month_str=$8; day=$9; split($10, t, ":");
        year=strftime("%Y", now);
        datespec=sprintf("%d %d %d %d %d 0", year, M[month_str], day, t[1], t[2]);
        submit_ts=mktime(datespec);
        if (submit_ts > now) {
            datespec=sprintf("%d %d %d %d %d 0", year-1, M[month_str], day, t[1], t[2]);
            submit_ts=mktime(datespec);
        }
        duration = now - submit_ts;
        if (duration < 0) duration = 0;
        
        run_time_str = (duration < 86400) ? int(duration/3600)"h" : int(duration/86400)"d";

        # 根据任务状态选择 Emoji
        status_emoji = "⚙️"; # 默认符号
        if ($3 == "RUN")      { status_emoji = "🏃"; } 
        else if ($3 == "PEND")  { status_emoji = "⏳"; }
        else if ($3 == "DONE")  { status_emoji = "✅"; }
        else if ($3 == "EXIT")  { status_emoji = "❌"; }

        # --- 打印格式化的卡片 (这是最终版布局) ---
        # 主标题行: [任务名称] 状态 (运行时长)
        printf("%s[%s]%s %s %s%s (%s)%s\n", C_GREEN_B, $7, C_RESET, status_emoji, C_CYAN_B, $3, run_time_str, C_RESET);
        
        # 详细信息行: 任务ID, 用户, 队列, 主机路由
        printf("  ├─ 🆔  任务ID:   %s\n", $1);
        # printf("  ├─ 👤  用户:     %s\n", $2);
        # printf("  ├─ 👤  用户:     Xiang LI\n");
        printf("  ├─ 💻  队列:     %s\n", $4);
        # printf("  └─ 🔗  主机路由: %s -> %s\n", $5, $6);
        
        # 在每个卡片后打印一个空行作为分隔
        # print "";
    }
  '
}

function cjj() {
    # 1. 定义日志文件路径 - 使用统一配置
    local log_file="${SUBMITTED_JOBS_FILE:-$HOME/submitted_jobs.txt}"
    
    # 2. 打印表头
    if [ ! -f "$log_file" ]; then
        echo "❌ 错误: 找不到日志文件 $log_file"
        return 1
    fi
    
    echo "================================================================================================================"
    echo "🚀 当前运行任务详情汇总"
    echo "📂 数据来源: $log_file"
    echo "================================================================================================================"

    # 3. 核心修改：使用管道 + while read 逐行读取 ID，防止变量合并
    #    cjobs | 跳过第一行 | 打印第一列 | 过滤非数字字符(防止颜色代码干扰)
    cjobs | awk 'NR>1 {print $1}' | sed 's/[^0-9]*//g' | while read -r job_id; do
        
        # 防止空行
        if [ -z "$job_id" ]; then continue; fi

        # 4. 针对每个 job_id 去文件中倒序查找
        tac "$log_file" | awk -v target="$job_id" '
            BEGIN { 
                # 以长分隔线为记录分隔符
                RS="==============================================================================================" 
            }
            
            # 匹配规则：使用正则表达式匹配 ID，允许冒号前后有空格
            # $0 是当前的一个完整记录块
            $0 ~ "🆔 作业ID\\s*:\\s*" target {
                
                # 提取各个字段 (正则优化：允许冒号后有空格)
                match($0, /📝 作业名称\s*:\s*([^\n]+)/, name)
                match($0, /✅ 提交状态\s*:\s*([^\n]+)/, status)
                match($0, /📂 工作目录\s*:\s*([^\n]+)/, dir)
                match($0, /🔄 队列\s*:\s*([^\n]+)/, queue)
                match($0, /📜 执行脚本\s*:\s*([^\n]+)/, script)
                match($0, /💻 完整命令\s*:\s*([^\n]+)/, cmd)
                match($0, /📤 输出日志\s*:\s*([^\n]+)/, out_log)
                match($0, /📥 错误日志\s*:\s*([^\n]+)/, err_log)
                
                # 打印详情
                # print "🆔 作业ID  : " target
                print "📝 作业名称: " (name[1] ? name[1] : "N/A")
                # print "🔄 队列名称: " (queue[1] ? queue[1] : "N/A")
                print "📂 工作目录: " (dir[1] ? dir[1] : "N/A")
                print "📜 执行脚本: " (script[1] ? script[1] : "N/A")
                print "📤 输出日志: " (out_log[1] ? out_log[1] : "N/A")
                print "📥 错误日志: " (err_log[1] ? err_log[1] : "N/A")
                
                found=1
                exit # 找到就退出awk
            }
            
            END {
                # 如果没找到
                if (!found) {
                    print "🆔 作业ID  : " target
                    print "⚠️  (未在日志文件中找到该任务的提交记录)"
                }
            }
        '
        echo "----------------------------------------------------------------------------------------------------------------"
    done
}

# 批量提交脚本
# ==================== 批量作业提交工具 ====================
# 用法: batch_sub -i <cmd_list> [-j <prefix>] [-n <cpus>] [-m <mem>]
# 将此函数添加到 ~/.zshrc 中，然后执行 source ~/.zshrc
# =========================================================

batch_sub() {
    # 默认配置
    local CPU_NUM=64
    local MEM_GB=300
    local INPUT_FILE=""
    local JOB_PREFIX=""
    local SLEEP_TIME=5
    local DRY_RUN=false
    
    # 帮助信息
    local usage="
📋 批量作业提交工具 (Batch Job Submitter)

用法: batch_sub -i <cmd_list> [-j <prefix>] [-n <cpus>] [-m <mem>] [-s <sleep>] [-d]

必选参数:
  -i <file>    包含命令的列表文件 (每行一条命令)

可选参数:
  -j <name>    作业名前缀 (推荐使用，便于作业管理)
                 不指定时自动从命令中提取名称
  -n <number>  CPU 核心数 (默认: 64)
  -m <gb>      内存大小 GB (默认: 300)
  -s <sec>     提交间隔秒数 (默认: 5)
  -d           试运行模式 (仅生成脚本，不实际提交)
  -h           显示此帮助信息

命名规则 (优先级从高到低):
  1. 用户指定前缀 (-j MyTask)     → MyTask_1, MyTask_2...
  2. 自动提取染色体 (Chr01...)     → Job_Chr01, Job_Chr02...
  3. 提取脚本文件名 (run.sh)       → run_1, run_2...
  4. 纯行号回退                    → Job_Line1, Job_Line2...

示例:
  # 1. 基础用法 (自动命名)
  batch_sub -i commands.txt

  # 2. 指定作业名和资源
  batch_sub -i tasks.txt -j Mapping -n 40 -m 200

  # 3. 试运行 (检查脚本生成，不提交)
  batch_sub -i commands.txt -d

  # 4. 慢速提交 (避免过载)
  batch_sub -i commands.txt -s 1.0
"

    # 参数解析
    while getopts ":i:n:m:j:s:dh" opt; do
        case $opt in
            i) INPUT_FILE="$OPTARG" ;;
            n) CPU_NUM="$OPTARG" ;;
            m) MEM_GB="$OPTARG" ;;
            j) JOB_PREFIX="$OPTARG" ;;
            s) SLEEP_TIME="$OPTARG" ;;
            d) DRY_RUN=true ;;
            h) echo "$usage"; return 0 ;;
            \?) echo "❌ 无效选项: -$OPTARG"; echo "$usage"; return 1 ;;
            :) echo "❌ 选项 -$OPTARG 需要参数"; echo "$usage"; return 1 ;;
        esac
    done

    # 输入验证
    if [[ -z "$INPUT_FILE" ]]; then
        echo "❌ 错误: 缺少必需参数 -i (命令列表文件)"
        echo "$usage"
        return 1
    fi

    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "❌ 错误: 文件不存在: '$INPUT_FILE'"
        return 1
    fi

    if [[ ! -r "$INPUT_FILE" ]]; then
        echo "❌ 错误: 文件无读取权限: '$INPUT_FILE'"
        return 1
    fi

    # 参数范围检查
    if ! [[ "$CPU_NUM" =~ ^[0-9]+$ ]] || (( CPU_NUM <= 0 )); then
        echo "❌ 错误: CPU 数量必须是正整数 (当前: $CPU_NUM)"
        return 1
    fi

    if ! [[ "$MEM_GB" =~ ^[0-9]+$ ]] || (( MEM_GB <= 0 )); then
        echo "❌ 错误: 内存大小必须是正整数 (当前: $MEM_GB GB)"
        return 1
    fi

    # 环境检查
    if ! command -v sub &> /dev/null && [[ "$DRY_RUN" == false ]]; then
        echo "⚠️  警告: 未找到 'sub' 命令，将自动启用试运行模式 (-d)"
        DRY_RUN=true
    fi

    # 初始化
    local TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    local TEMP_DIR="./batch_jobs_${TIMESTAMP}"
    
    if ! mkdir -p "$TEMP_DIR"; then
        echo "❌ 错误: 无法创建目录 '$TEMP_DIR'"
        return 1
    fi

    echo "=========================================="
    echo "📂 工作目录: $TEMP_DIR"
    echo "📄 输入文件: $INPUT_FILE"
    echo "⚙️  资源配置: CPU=$CPU_NUM 核, 内存=${MEM_GB}GB"
    [[ -n "$JOB_PREFIX" ]] && echo "🏷️  作业前缀: $JOB_PREFIX"
    [[ "$DRY_RUN" == true ]] && echo "🔍 模式: 试运行 (不实际提交)"
    echo "=========================================="

    local LINE_NUM=0
    local SUBMIT_COUNT=0
    local SKIP_COUNT=0
    local FAIL_COUNT=0
    local -a FAILED_JOBS=()

    # 主处理循环
    while IFS= read -r RAW_LINE || [[ -n "$RAW_LINE" ]]; do
        ((LINE_NUM++))
        
        # 去除首尾空格
        local CMD=$(echo "$RAW_LINE" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')

        # 跳过空行和注释
        if [[ -z "$CMD" ]] || [[ "$CMD" == \#* ]]; then
            ((SKIP_COUNT++))
            continue
        fi

        # 作业命名逻辑
        local JOB_NAME=""

        # 策略1: 用户指定前缀
        if [[ -n "$JOB_PREFIX" ]]; then
            JOB_NAME="${JOB_PREFIX}_$((SUBMIT_COUNT + 1))"

        # 策略2: 自动匹配染色体 (Chr01, Chr1, chr01 等)
        elif [[ "$CMD" =~ [Cc]hr0*([0-9]+) ]]; then
            local CHR_NUM="${BASH_REMATCH[1]}"
            JOB_NAME="Job_Chr$(printf "%02d" $CHR_NUM)"

        # 策略3: 提取脚本文件名
        else
            local SCRIPT_FILE=$(echo "$CMD" | grep -oE '[^[:space:]]+\.(sh|py|pl|R|jar|rb|js)' | head -1)
            
            if [[ -n "$SCRIPT_FILE" ]]; then
                local BASE=$(basename "$SCRIPT_FILE")
                local NAME_NO_EXT="${BASE%.*}"
                # 清理文件名中的特殊字符
                NAME_NO_EXT=$(echo "$NAME_NO_EXT" | tr -c '[:alnum:]_-' '_')
                JOB_NAME="${NAME_NO_EXT}_$((SUBMIT_COUNT + 1))"
            else
                # 策略4: 回退到行号
                JOB_NAME="Job_Line${LINE_NUM}"
            fi
        fi

        # 生成作业脚本
        local SCRIPT_NAME="${TEMP_DIR}/${JOB_NAME}.sh"
        
        # 如果文件已存在，添加后缀避免覆盖
        local SUFFIX=1
        while [[ -f "$SCRIPT_NAME" ]]; do
            SCRIPT_NAME="${TEMP_DIR}/${JOB_NAME}_v${SUFFIX}.sh"
            ((SUFFIX++))
        done

        cat > "$SCRIPT_NAME" <<EOF
#!/bin/bash
#================================================
# 自动生成的作业脚本
# 原始命令文件: $INPUT_FILE
# 原始行号: $LINE_NUM
# 生成时间: $(date '+%Y-%m-%d %H:%M:%S')
#================================================

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任何一个失败都返回失败

echo "🚀 作业开始: $JOB_NAME"
echo "⏰ 开始时间: \$(date '+%Y-%m-%d %H:%M:%S')"
echo "📍 执行命令: $CMD"
echo "=========================================="

# 执行原始命令
$CMD

EXIT_CODE=\$?

echo "=========================================="
echo "⏰ 结束时间: \$(date '+%Y-%m-%d %H:%M:%S')"
if [[ \$EXIT_CODE -eq 0 ]]; then
    echo "✅ 作业成功: $JOB_NAME"
else
    echo "❌ 作业失败: $JOB_NAME (退出码: \$EXIT_CODE)"
fi

exit \$EXIT_CODE
EOF

        chmod +x "$SCRIPT_NAME"

        # 提交作业
        echo "[$((SUBMIT_COUNT + 1))] 📤 提交: $JOB_NAME (来自第 $LINE_NUM 行)"
        
        if [[ "$DRY_RUN" == true ]]; then
            echo "   └─ [试运行] 脚本已生成: $(basename "$SCRIPT_NAME")"
        else
            if sub "$JOB_NAME" -n "$CPU_NUM" -m "$MEM_GB" "$SCRIPT_NAME" 2>&1; then
                echo "   └─ ✅ 提交成功"
            else
                echo "   └─ ❌ 提交失败"
                ((FAIL_COUNT++))
                FAILED_JOBS+=("$JOB_NAME (行 $LINE_NUM)")
            fi
            
            sleep "$SLEEP_TIME"
        fi

        ((SUBMIT_COUNT++))

    done < "$INPUT_FILE"

    # 总结报告
    echo ""
    echo "=========================================="
    echo "📊 批量提交完成"
    echo "=========================================="
    echo "📄 总行数: $LINE_NUM"
    echo "✅ 已提交: $SUBMIT_COUNT 个作业"
    echo "⏭️  已跳过: $SKIP_COUNT 行 (空行/注释)"
    
    if (( FAIL_COUNT > 0 )); then
        echo "❌ 失败数: $FAIL_COUNT"
        echo ""
        echo "失败作业列表:"
        for job in "${FAILED_JOBS[@]}"; do
            echo "  - $job"
        done
        echo ""
        echo "⚠️  建议: 检查失败作业的脚本文件在 $TEMP_DIR"
    fi
    
    echo "📂 脚本目录: $TEMP_DIR"
    
    if [[ "$DRY_RUN" == true ]]; then
        echo ""
        echo "💡 提示: 这是试运行模式，未实际提交作业"
        echo "   去掉 -d 参数可真实提交"
    fi
    
    echo "=========================================="

    # 返回状态码
    if (( FAIL_COUNT > 0 )); then
        return 1
    else
        return 0
    fi
}

# 模块加载成功标记
export ZSH_MODULE_CLUSTER_TOOLS_LOADED=1