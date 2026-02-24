# =============================================================================
#  统一路径配置文件
#  Unified Path Configuration
# =============================================================================
# 此文件定义所有常用的路径，便于统一管理和修改
# 所有路径均基于 $HOME，提升可移植性
# =============================================================================

# =============================================================================
#  用户基础路径 (User Base Paths)
# =============================================================================
export USER_HOME="${HOME}"
export USER_BASE="${USER_BASE:-$HOME}"

# =============================================================================
#  软件工具路径 (Software & Tools Paths)
# =============================================================================
export MINIFORGE3_DIR="${MINIFORGE3_DIR:-${USER_BASE}/miniforge3}"
export SOFTWARE_DIR="${SOFTWARE_DIR:-${USER_BASE}/software}"
export BIOPYTOOLS_DIR="${BIOPYTOOLS_DIR:-${USER_BASE}/software/biopytools}"
export SCRIPTS_DIR="${SCRIPTS_DIR:-${SOFTWARE_DIR}/scripts}"

# =============================================================================
#  数据库路径 (Database Paths)
# =============================================================================
export DATABASE_DIR="${DATABASE_DIR:-${USER_BASE}/database}"
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-${DATABASE_DIR}/gtdbtk/release226}"

# =============================================================================
#  日志与备份路径 (Log & Backup Paths)
# =============================================================================
export JOB_LOGS_DIR="${JOB_LOGS_DIR:-${USER_HOME}/joblogs}"
export SUBMITTED_JOBS_FILE="${SUBMITTED_JOBS_FILE:-${USER_BASE}/submitted_jobs.txt}"
export BACKUP_BASE_DIR="${BACKUP_BASE_DIR:-${USER_HOME}/conda_env_backups}"

# =============================================================================
#  ZSH 配置路径 (ZSH Configuration Paths)
# =============================================================================
export ZSH_FUNCTIONS_DIR="${ZSH_FUNCTIONS_DIR:-${USER_HOME}/zsh/functions}"
export ZSH_CUSTOM_PLUGINS_DIR="${ZSH_CUSTOM_PLUGINS_DIR:-${USER_HOME}/.oh-my-zsh/custom/plugins}"

# =============================================================================
#  本地二进制路径 (Local Binary Paths)
# =============================================================================
export LOCAL_BIN_DIR="${LOCAL_BIN_DIR:-${USER_HOME}/.local/bin}"

# =============================================================================
#  Conda/Mamba 路径 (Conda/Mamba Paths)
# =============================================================================
export CONDA_EXE="${CONDA_EXE:-${MINIFORGE3_DIR}/bin/conda}"
export MAMBA_EXE="${MAMBA_EXE:-${MINIFORGE3_DIR}/bin/mamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${MINIFORGE3_DIR}}"

# =============================================================================
#  Aspera 配置 (Aspera Configuration for SRA Download)
# =============================================================================
export ASPERA_KEY="${ASPERA_KEY:-${MINIFORGE3_DIR}/envs/aspera_v.3.9.6/etc/asperaweb_id_dsa.openssh}"
export ASPERA_SERVER="${ASPERA_SERVER:-era-fasp@fasp.sra.ebi.ac.uk:}"

# =============================================================================
#  路径配置已加载标记 (Path Configuration Loaded Marker)
# =============================================================================
export ZSH_PATH_CONFIG_LOADED=1
