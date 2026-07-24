# 超算日志系统分离规范|Job Scheduler Log Separation

> 📖 本文件为 CLAUDE.md 的按需参考文档,**配置超算(作业调度)日志时读取**。核心日志格式见 CLAUDE.md §二。
> On-demand reference for CLAUDE.md; read when configuring logging for the HPC job scheduler.

---

在超算系统（如使用 `sub` 函数提交任务）运行时，必须正确配置 stdout 和 stderr handler，以实现日志分离：

> **原则|Principle:** INFO → stdout → .out 文件，WARNING/ERROR → stderr → .err 文件

## 1. 标准配置|Standard Configuration

```python
import logging
import sys

logger = logging.getLogger("ModuleName")
logger.setLevel(logging.DEBUG)
logger.handlers.clear()
logger.propagate = False  # 避免重复|Avoid duplicates

formatter = logging.Formatter(
    '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# stdout handler - INFO级别|stdout handler - INFO level
# → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.INFO)
stdout_handler.setFormatter(formatter)
logger.addHandler(stdout_handler)

# stderr handler - WARNING及以上|stderr handler - WARNING and above
# → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
stderr_handler = logging.StreamHandler(sys.stderr)
stderr_handler.setLevel(logging.WARNING)
stderr_handler.setFormatter(formatter)
logger.addHandler(stderr_handler)

# 文件handler - 所有级别|File handler - all levels
# → 本地完整日志|→ Local complete log
file_handler = logging.FileHandler('99_logs/pipeline.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
```

## 2. 日志流向|Log Flow

```
Python Logger (ModuleName)
    │
    ├─ stdout handler (INFO)
    │   └─→ 系统stdout → 超算捕获 → YYYY-MM-DD-HH-MM-jobname.out
    │       内容：正常流程日志|Content: Normal operation logs
    │
    ├─ stderr handler (WARNING+)
    │   └─→ 系统stderr → 超算捕获 → YYYY-MM-DD-HH-MM-jobname.err
    │       内容：错误和警告|Content: Errors and warnings
    │
    └─ file handler (DEBUG+)
        └─→ 本地文件 → 99_logs/pipeline.log
            内容：完整调试日志|Content: Complete debug logs
```

## 3. 常见错误|Common Mistakes

❌ **错误1：所有日志都输出到stderr**
```python
# 导致 .out 文件为空，.err 文件包含所有信息
handler = logging.StreamHandler(sys.stderr)  # ❌
```

❌ **错误2：stdout 和 stderr 都输出INFO**
```python
# 导致 .out 和 .err 内容重复
stdout_handler.setLevel(logging.INFO)   # ✅
stderr_handler.setLevel(logging.INFO)    # ❌ 应该是 WARNING
```

❌ **错误3：忘记 propagate=False**
```python
# 导致日志重复输出
logger.propagate = True  # ❌ 会传播到root logger
logger.propagate = False  # ✅ 避免重复
```
