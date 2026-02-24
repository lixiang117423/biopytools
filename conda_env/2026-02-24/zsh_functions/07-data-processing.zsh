# =============================================================================
#  07-data-processing.zsh - 数据处理函数模块
#  Data Processing Functions Module
# =============================================================================

# =============================================================================
# [基础数据筛选功能] Basic Data Filtering Functions
# =============================================================================

# 过滤文件中的特定列 - 支持多种比较操作符
myfilter() {
    if [ $# -lt 4 ]; then
        echo "用法: myfilter <列号> <操作符> <值> <文件>"
        echo "操作符支持:"
        echo "  eq  - 等于 (==)"
        echo "  ne  - 不等于 (!=)"
        echo "  gt  - 大于 (>)"
        echo "  lt  - 小于 (<)"
        echo "  ge  - 大于等于 (>=)"
        echo "  le  - 小于等于 (<=)"
        echo "  match    - 包含模式 (~)"
        echo "  nomatch  - 不包含模式 (!~)"
        echo ""
        echo "例如: myfilter 3 gt 100 data.txt"
        echo "     myfilter 2 eq hello data.txt"
        echo "     myfilter 1 match error data.txt"
        return 1
    fi
    
    local column="$1"
    local operator="$2"
    local value="$3"
    local file="$4"
    
    # 将文字操作符转换为符号
    local awk_op
    case "$operator" in
        "eq"|"="|"==")
            awk_op="=="
            ;;
        "ne"|"!=")
            awk_op="!="
            ;;
        "gt"|">")
            awk_op=">"
            ;;
        "lt"|"<")
            awk_op="<"
            ;;
        "ge"|">=")
            awk_op=">="
            ;;
        "le"|"<=")
            awk_op="<="
            ;;
        "match"|"~")
            awk_op="~"
            ;;
        "nomatch"|"!~")
            awk_op="!~"
            ;;
        *)
            echo "错误: 不支持的操作符 '$operator'"
            echo "支持的操作符: eq ne gt lt ge le match nomatch"
            return 1
            ;;
    esac
    
    # 检查文件是否存在
    if [ ! -f "$file" ]; then
        echo "错误: 文件 '$file' 不存在"
        return 1
    fi
    
    # 根据操作符类型处理
    case "$awk_op" in
        "~"|"!~")
            # 正则匹配
            awk -v col="$column" -v op="$awk_op" -v val="$value" '
                {
                    if (op == "~") {
                        if ($col ~ val) print
                    } else if (op == "!~") {
                        if ($col !~ val) print
                    }
                }
            ' "$file"
            ;;
        *)
            # 数值或字符串比较
            awk -v col="$column" -v op="$awk_op" -v val="$value" '
                {
                    if (op == "==") {
                        if ($col == val) print
                    } else if (op == "!=") {
                        if ($col != val) print
                    } else if (op == ">") {
                        if ($col > val) print
                    } else if (op == "<") {
                        if ($col < val) print
                    } else if (op == ">=") {
                        if ($col >= val) print
                    } else if (op == "<=") {
                        if ($col <= val) print
                    }
                }
            ' "$file"
            ;;
    esac
}

# myfilter的便捷别名
myfgt() { myfilter "$1" gt "$2" "$3"; }      # 大于
myflt() { myfilter "$1" lt "$2" "$3"; }      # 小于
myfeq() { myfilter "$1" eq "$2" "$3"; }      # 等于
myfne() { myfilter "$1" ne "$2" "$3"; }      # 不等于
myfge() { myfilter "$1" ge "$2" "$3"; }      # 大于等于
myfle() { myfilter "$1" le "$2" "$3"; }      # 小于等于
myfmatch() { myfilter "$1" match "$2" "$3"; }   # 包含模式
myfnmatch() { myfilter "$1" nomatch "$2" "$3"; } # 不包含模式

# 按列名过滤非零行
filter_col_nonzero() {
  if [ -z "$1" ]; then
    echo "Usage: filter_col_nonzero <ColumnName> [filename]" >&2
    echo "Example: filter_col_nonzero FPKM data.txt" >&2
    echo "If no filename is provided, it reads from standard input (stdin)." >&2
    return 1
  fi

  local target_col_name="$1"
  local input_file="$2"

  awk -F'\t' -v target_col="$target_col_name" '
    NR==1 {
      for(i=1; i<=NF; i++) {
        if($i == target_col) {
          col_num = i
        }
      }
      if (col_num == 0) {
        print "Error: Column '\''" target_col "'\'' not found." > "/dev/stderr"
        exit 1
      }
      print
    }
    NR > 1 && col_num > 0 && $col_num != 0
  ' "$input_file"
}

# =============================================================================
# [Bash版dplyr数据处理函数库] Bash dplyr-style Data Processing Functions
# =============================================================================

# SELECT - 选择指定列
dselect() {
    if [ $# -eq 0 ]; then
        echo "用法: dselect <列号1,列号2,...> [文件]"
        echo "例如: dselect 1,3,5 data.txt"
        echo "     cat data.txt | dselect 2,4"
        return 1
    fi
    
    local columns="$1"
    local file="$2"
    
    # 将逗号分隔的列号转换为awk格式
    local awk_cols=$(echo "$columns" | sed 's/,/ "," /g' | sed 's/^/\$/; s/$/ ","/' | sed 's/, *$//')
    
    if [ -n "$file" ]; then
        awk -v OFS="\t" "{print $awk_cols}" "$file"
    else
        awk -v OFS="\t" "{print $awk_cols}"
    fi
}

# FILTER - 筛选行（dplyr风格接口）
dfilter() {
    if [ $# -lt 3 ]; then
        echo "用法: dfilter <列号> <操作符> <值> [文件]"
        echo "操作符: eq ne gt lt ge le match nomatch"
        echo "例如: dfilter 3 gt 100 data.txt"
        echo "     cat data.txt | dfilter 2 eq hello"
        return 1
    fi
    
    local column="$1"
    local operator="$2"
    local value="$3"
    local file="$4"
    
    if [ -n "$file" ]; then
        myfilter "$column" "$operator" "$value" "$file"
    else
        # 从stdin读取
        myfilter "$column" "$operator" "$value" /dev/stdin
    fi
}

# ARRANGE - 排序
darrange() {
    if [ $# -eq 0 ]; then
        echo "用法: darrange <列号> [asc|desc] [文件]"
        echo "例如: darrange 3 desc data.txt"
        echo "     cat data.txt | darrange 2 asc"
        return 1
    fi
    
    local column="$1"
    local order="${2:-asc}"
    local file="$3"
    
    local sort_opts=""
    if [ "$order" = "desc" ]; then
        sort_opts="-r"
    fi
    
    # 检查是否为数值排序
    if [ -n "$file" ]; then
        # 检测列是否为数值
        local is_numeric=$(awk -v col="$column" 'NR==2 {print ($col ~ /^[0-9.+-]+$/)}' "$file")
        if [ "$is_numeric" = "1" ]; then
            sort -k"$column" -n $sort_opts "$file"
        else
            sort -k"$column" $sort_opts "$file"
        fi
    else
        # 从stdin读取时，假设第一行数据来判断
        local temp_file=$(mktemp)
        cat > "$temp_file"
        local is_numeric=$(awk -v col="$column" 'NR==2 {print ($col ~ /^[0-9.+-]+$/)}' "$temp_file")
        if [ "$is_numeric" = "1" ]; then
            sort -k"$column" -n $sort_opts "$temp_file"
        else
            sort -k"$column" $sort_opts "$temp_file"
        fi
        rm -f "$temp_file"
    fi
}

# MUTATE - 创建新列
dmutate() {
    if [ $# -lt 2 ]; then
        echo "用法: dmutate '<新列表达式>' [文件]"
        echo "例如: dmutate '\$4 = \$2 + \$3' data.txt  # 第4列 = 第2列 + 第3列"
        echo "     dmutate '\$5 = \$1 * 100' data.txt   # 第5列 = 第1列 * 100"
        echo "     cat data.txt | dmutate '\$3 = \$1 / \$2'"
        return 1
    fi
    
    local expression="$1"
    local file="$2"
    
    if [ -n "$file" ]; then
        awk "{$expression; print}" "$file"
    else
        awk "{$expression; print}"
    fi
}

# SUMMARISE - 汇总统计
dsummarise() {
    if [ $# -eq 0 ]; then
        echo "用法: dsummarise <统计函数> <列号> [文件]"
        echo "统计函数: sum mean count min max"
        echo "例如: dsummarise sum 3 data.txt"
        echo "     cat data.txt | dsummarise mean 2"
        return 1
    fi
    
    local func="$1"
    local column="$2"
    local file="$3"
    
    local awk_script=""
    case "$func" in
        "sum")
            awk_script="BEGIN{sum=0} {sum+=\$$column} END{print sum}"
            ;;
        "mean"|"avg")
            awk_script="BEGIN{sum=0; count=0} {sum+=\$$column; count++} END{if(count>0) print sum/count; else print 0}"
            ;;
        "count")
            awk_script="END{print NR}"
            ;;
        "min")
            awk_script="NR==1{min=\$$column} {if(\$$column<min) min=\$$column} END{print min}"
            ;;
        "max")
            awk_script="NR==1{max=\$$column} {if(\$$column>max) max=\$$column} END{print max}"
            ;;
        *)
            echo "错误: 不支持的统计函数 '$func'"
            return 1
            ;;
    esac
    
    if [ -n "$file" ]; then
        awk "$awk_script" "$file"
    else
        awk "$awk_script"
    fi
}

# DISTINCT - 去重
ddistinct() {
    local column="$1"
    local file="$2"
    
    if [ -z "$column" ]; then
        # 整行去重
        if [ -n "$file" ]; then
            sort "$file" | uniq
        else
            sort | uniq
        fi
    else
        # 基于特定列去重
        if [ -n "$file" ]; then
            awk -v col="$column" '!seen[$col]++' "$file"
        else
            awk -v col="$column" '!seen[$col]++'
        fi
    fi
}

# COUNT - 计数分组
dcount() {
    if [ $# -eq 0 ]; then
        echo "用法: dcount <列号> [文件]"
        echo "例如: dcount 2 data.txt  # 统计第2列各值的出现次数"
        echo "     cat data.txt | dcount 1"
        return 1
    fi
    
    local column="$1"
    local file="$2"
    
    if [ -n "$file" ]; then
        awk -v col="$column" '{count[$col]++} END{for(i in count) print i "\t" count[i]}' "$file" | sort
    else
        awk -v col="$column" '{count[$col]++} END{for(i in count) print i "\t" count[i]}' | sort
    fi
}

# SLICE - 按行号选择
dslice() {
    if [ $# -eq 0 ]; then
        echo "用法: dslice <起始行>:<结束行> [文件]"
        echo "     dslice <行号> [文件]"
        echo "例如: dslice 5:10 data.txt  # 选择第5到10行"
        echo "     dslice 3 data.txt      # 选择第3行"
        echo "     cat data.txt | dslice 1:5"
        return 1
    fi
    
    local range="$1"
    local file="$2"
    
    if [[ "$range" == *":"* ]]; then
        # 范围选择
        local start=$(echo "$range" | cut -d: -f1)
        local end=$(echo "$range" | cut -d: -f2)
        if [ -n "$file" ]; then
            sed -n "${start},${end}p" "$file"
        else
            sed -n "${start},${end}p"
        fi
    else
        # 单行选择
        if [ -n "$file" ]; then
            sed -n "${range}p" "$file"
        else
            sed -n "${range}p"
        fi
    fi
}

# HEAD/TAIL - 查看头部/尾部
dhead() {
    local n="${1:-10}"
    local file="$2"
    
    if [ -n "$file" ]; then
        head -n "$n" "$file"
    else
        head -n "$n"
    fi
}

dtail() {
    local n="${1:-10}"
    local file="$2"
    
    if [ -n "$file" ]; then
        tail -n "$n" "$file"
    else
        tail -n "$n"
    fi
}

# =============================================================================
# [实用工具函数] Utility Functions
# =============================================================================

# 数据预览
dglimpse() {
    local file="$1"
    if [ -z "$file" ]; then
        echo "用法: dglimpse <文件>"
        return 1
    fi
    
    echo "=== 文件信息 ==="
    echo "行数: $(wc -l < "$file")"
    echo "列数: $(awk '{print NF; exit}' "$file")"
    echo ""
    echo "=== 前5行 ==="
    head -5 "$file"
    echo ""
    echo "=== 各列统计 ==="
    local ncols=$(awk '{print NF; exit}' "$file")
    for i in $(seq 1 $ncols); do
        echo "列 $i:"
        echo "  唯一值数: $(awk -v col=$i '{print $col}' "$file" | sort | uniq | wc -l)"
        echo "  示例值: $(awk -v col=$i 'NR<=3 {printf "%s ", $col}' "$file")"
        echo ""
    done
}

# 快速统计
dquick_stats() {
    local column="$1"
    local file="$2"
    
    if [ -z "$column" ]; then
        echo "用法: dquick_stats <列号> [文件]"
        return 1
    fi
    
    echo "=== 第${column}列统计 ==="
    echo "总数: $(dsummarise count "$column" "$file")"
    echo "最小值: $(dsummarise min "$column" "$file")"
    echo "最大值: $(dsummarise max "$column" "$file")"
    echo "平均值: $(dsummarise mean "$column" "$file")"
    echo "唯一值: $(if [ -n "$file" ]; then awk -v col="$column" '{print $col}' "$file"; else awk -v col="$column" '{print $col}'; fi | sort | uniq | wc -l)"
}

# 帮助信息
dhelp() {
    echo "=== Shell工具函数库 ==="
    echo ""
    echo "📊 数据处理 (dplyr风格):"
    echo "  dselect   - 选择列"
    echo "  dfilter   - 筛选行"
    echo "  ddistinct - 去重"
    echo "  dslice    - 按行号选择"
    echo "  darrange  - 排序"
    echo "  dmutate   - 创建新列"
    echo "  dsummarise- 汇总统计"
    echo "  dcount    - 计数分组"
    echo "  dhead     - 查看头部"
    echo "  dtail     - 查看尾部"
    echo "  dglimpse  - 数据预览"
    echo "  dquick_stats - 快速统计"
    echo ""
    echo "📋 传统数据筛选:"
    echo "  myfilter  - 基础筛选（支持eq,ne,gt,lt,ge,le,match,nomatch）"
    echo "  filter_col_nonzero - 按列名过滤非零行"
    echo ""
    echo "使用 <函数名> -h 或 <函数名> 可查看具体用法"
    echo "数据处理函数都支持管道操作!"
}

# 模块加载成功标记
export ZSH_MODULE_DATA_PROCESSING_LOADED=1